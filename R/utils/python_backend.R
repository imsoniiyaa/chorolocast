get_py_dir <- function() {
    cands <- c(
        file.path(getwd(), "py"),
        file.path(getwd(), "..", "py")
    )

    cands <- cands[
        file.exists(file.path(cands, "fm_backend.py"))
    ]

    if (length(cands)) {
        normalizePath(cands[1])
    } else {
        NA_character_
    }
}
backend_python <- function(model, ctx_len) {
    PY_DIR <- get_py_dir()
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("install.packages('reticulate')")
    }

    if (is.na(PY_DIR)) {
        stop("cannot find py/fm_backend.py; set CHLOROCAST_PY_DIR")
    }

    m <- reticulate::import_from_path(
        "fm_backend",
        path = PY_DIR
    )

    b <- m$load_backend(
        model,
        as.integer(ctx_len)
    )

    list(
        name = model,
        predict = function(ctxs, horizon) {
            flat <- m$predict_flat(
                b,
                lapply(ctxs, function(cx) as.numeric(if (is.list(cx)) cx$y else cx)),
                as.integer(horizon)
            )

            aperm(
                array(
                    as.numeric(unlist(flat)),
                    dim = c(9, horizon, length(ctxs))
                ),
                c(3, 2, 1)
            )
        }
    )
}
