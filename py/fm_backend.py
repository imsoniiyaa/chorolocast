#This backend code helps to run the choronos and timesfm model(foundation models) and send it to R for verifying
import numpy as np

LEVELS = [round(0.1 * i, 1) for i in range(1, 10)]

#First method - choronos model
#Using amazon prediction model since it is smart and can predict the future values based on the past values
class ChronosBolt:
    name = "chronos"

    def __init__(self, model_id="amazon/chronos-bolt-base", device="cpu"): 
        import torch
        from chronos import BaseChronosPipeline

        self.torch = torch
        self.pipe = BaseChronosPipeline.from_pretrained(
            model_id, device_map=device, torch_dtype=torch.float32
        )

    def predict(self, ctxs, horizon):
        outs = []
        for i in range(0, len(ctxs), 32):
            batch = [self.torch.tensor(np.asarray(c, dtype=np.float32)) for c in ctxs[i : i + 32]]
            q, _ = self.pipe.predict_quantiles(
                batch, prediction_length=int(horizon), quantile_levels=LEVELS
            )
            outs.append(q.float().numpy())
        return np.concatenate(outs)

#Second method - timesfm model
#TimesFM is a foundation model for time series forecasting developed by Google Research.
#It is designed to handle a wide range of time series forecasting tasks, including long-term forecasting,
#multivariate forecasting, and probabilistic forecasting. 
# TimesFM uses a transformer-based architecture to capture complex temporal dependencies in time series data.
class TimesFM25:
    name = "timesfm"

    def __init__(self, ctx_len=512):
        import timesfm
        import torch

        torch.set_float32_matmul_precision("high")
        self.m = timesfm.TimesFM_2p5_200M_torch.from_pretrained(
            "google/timesfm-2.5-200m-pytorch", torch_compile=False  # avoid torch.compile issues on older torch / CPU
        )
        self.m.compile(
            timesfm.ForecastConfig(
                max_context=int(np.ceil(ctx_len / 32) * 32),
                max_horizon=128,
                normalize_inputs=True,
                use_continuous_quantile_head=True,
                force_flip_invariance=True,
                infer_is_positive=True,
                fix_quantile_crossing=True,
            )
        )

    def predict(self, ctxs, horizon):
        outs = []
        for i in range(0, len(ctxs), 32):
            _, q = self.m.forecast(
                horizon=int(horizon), inputs=[np.asarray(c, dtype=np.float32) for c in ctxs[i : i + 32]]
            )
            outs.append(np.asarray(q)[:, : int(horizon), 1:10])  # index 0 = mean, 1..9 = P10..P90
        return np.concatenate(outs)

#If R requests a model choronos or timesfm, this function will load the model and return it to R
#If the model is not recognized, it will raise a ValueError
def load_backend(name, ctx_len=512):
    if name == "chronos":
        return ChronosBolt()
    if name == "timesfm":
        return TimesFM25(int(ctx_len))
    raise ValueError(f"unknown model: {name}")

#This function will take the model, context and horizon and return the prediction to R
def predict(backend, ctxs, horizon):
    ctxs = [np.asarray(c, dtype=np.float32) for c in ctxs]
    return np.asarray(backend.predict(ctxs, int(horizon)), dtype=np.float64)

#This function will take the model, context and horizon and return the prediction to R in a flat list
def predict_flat(backend, ctxs, horizon):
    return predict(backend, ctxs, horizon).ravel().tolist()