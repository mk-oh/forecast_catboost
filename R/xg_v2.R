

######################################################
#### Function for Xgboost_ts



################### Forecast #####################################################
############################################################

############################################################
xgboost.forecast <- function (y,
                              #- priod
                              h = 42,
                              #- data lag set
                              maxlag = max(42, 4 * frequency(y)),
                              #- External Value Use (require : Matrix type)
                              xreg = NULL,
                              pred_xreg = NULL,
                              #- Model Parameter - #
                              nrounds = 500,
                              params = NULL,
                              nrounds_method = c("cv",# "v",
                                                 "none"),
                              nfold = ifelse(length(y) > 720, 10, 5),

                              #- Data. Preprocess
                              # 1. Boxcox Trans
                              boxcox = TRUE,
                              lambda = BoxCox.lambda(y),
                              # 2. Outlier
                              ts_clean = TRUE,
                              ts_lambda = 1,
                              # 3. Data lag using Type
                              season_type = c("decompose", 'dummy', "none"),
                              verbose = TRUE, ...)  {

  nrounds_method = match.arg(nrounds_method)
  season_type = match.arg(season_type)


  ###### Check Data format

  if(is.null(params)) {
    params = list(depth = 6,
                  boosting_type = 'gbdt',
                  eta = 0.15,
                  gamma = 1,
                  lambda = 0.9,
                  eval_metric  = 'mae')
  }
  #########

  if (!"ts" %in% class(y)) {
    stop("y class check : y need to ts class")
  }

  if (length(y))

    if (!is.null(xreg)) {
      xreg <- as.matrix(xreg)
      if (!is.numeric(xreg) | !is.matrix(xreg)) {
        stop("xreg check : class is numeric vector or a matrix")
      }
    }


  #### Outlier preprocess - (require : forecast package)
  if(ts_clean == TRUE) {
    y <- tsclean(y, ts_lambda)
  }

  temp_y <- y
  origxreg <- xreg

  f <- stats::frequency(y)
  orig_n <- length(y)

  #### Check a length of y
  if (orig_n < 46) {
    stop("must be data length > 45")
  }


  #### Check Maxlag period / trans.

  if (maxlag > (orig_n - f )) {
    warning(paste("y is too short ", maxlag, "-> ",
                  orig_n - f - round(f/2)))
    maxlag <- orig_n - f - round(f/2)
  }

  #### Log Translate  ---->  Crush...

  # if (log == TRUE) {
  #   y <- log1p(y)
  # }

  #### BoxCox Translate  (require : forecast package )
  if (boxcox == TRUE) {
    lambda <- BoxCox.lambda(y)
  } else {
    lambda <- 1
  }

  target_y <- Boxcox(y, lambda = lambda)


  if (maxlag != round(maxlag)) {
    maxlag <- round(maxlag)
    if (verbose) {
      message(paste("maxlag is need to int, change to ", maxlag))
    }
  }


  if (season_type == "decompose") {
    decomp <- decompose(target_y, type = "multiplicative")
    target_y <- seasadj(decomp)

  }


  n <- orig_n - maxlag

  # Target split -> y range : [1 ~ maxlag]
  y2 <- ts(target_y[-(1:(maxlag))],
           start = time(target_y)[maxlag + 1],
           frequency = f)

  if (season_type == "dummy" & f > 1) {
    ncolx <- maxlag + f - 1
  }

  if (season_type == "decompose" | (season_type == "none" | f == 1)) {
    ncolx <- maxlag
  }



  # if (!"ts" %in% class(y)) {
  #    stop("class check : y is time-seires")
  # }

  # if (!is.null(xreg)) {
  #     if (class(xreg) == "ts" | "data.frame" %in% class(xreg)) {
  #         xreg <- as.matrix(xreg)
  #     }
  #     if (!is.numeric(xreg) | !is.matrix(xreg)) {
  #         stop("xreg check : class is numeric vector or a matrix")
  #     }
  # }

  # # if(ts_clean == TRUE) {
  # #     y <- tsclean(y)
  # # }

  # # orign <- length(y)
  # # target_y <- Boxcox(y, lambda = lambda)
  # # K = max(1, min(round(f/4 - 1), 10)


  # # if (length(y) < 46) {
  # #     stop("must be data length > 45")
  # # }


  # # if (maxlag > (orign - f - round(f/4))) {
  # #     warning(paste("y is too short ", maxlag, "--> ",
  # #         orign - f - round(f/4)))
  # #     maxlag <- orign - f - round(f/4)
  # # }

  # # if (maxlag != round(maxlag)) {
  # #     maxlag <- round(maxlag)
  # #     if (verbose) {
  # #         message(paste("maxlag is need to int, change to ", maxlag))
  # #     }
  # # }


  # if (season_type == "decompose") {
  #     decomp <- decompose(target_y, type = "multiplicative")
  #     target_y <- seasadj(decomp)
  # }

  # if (season_type == "dummy" & f > 1) {
  #     ncolx <- maxlag + f - 1
  # }

  # xreg_temp <- xreg
  # n <- orign - maxlag

  # y2 <- ts(target_y[-(1:(maxlag))],
  #         start = time(target_y)[maxlag + 1],
  #         frequency = f)

  # if (nrounds_method == "cv" & n < 90) {
  #     warning("Data is too short. need to > 90")
  #     nrounds_method <- "v"
  # }

  # if (season_type == "decompose") {
  #     ncolx <- maxlag
  # }

  # if (season_type == "none" | f == 1) {
  #     ncolx <- maxlag
  # }



  x <- matrix(0, nrow = n, ncol = ncolx)
  x[, 1:maxlag] <- lag_y(target_y, maxlag, keeporig = FALSE)


  if (f == 1 || season_type == "decompose" || season_type == "none") {
    colnames(x) <- c(paste0("lag", 1:maxlag))
  }


  if (f > 1 & season_type == "dummy") {
    tmp <- data.frame(y = 1,
                      x = as.character(rep_len(1:f, n))
    )
    seasons <- stats::model.matrix(y ~ x, data = tmp)[, -1]
    x[, maxlag + 1:(f - 1)] <- seasons
    colnames(x) <- c(paste0("lag", 1:maxlag),
                     paste0("season", 2:f))
  }


  if (!is.null(xreg)) {
    xreg <- lag_xreg(xreg, maxlag = maxlag)
    x <- cbind(x, xreg[, , drop = FALSE])
  }


  if (nrounds_method == "cv") {
    if (verbose) {
      message("Start Cross-Validation")
    }
    cv <- xgb.cv(params= params,data = x, label = y2, nrounds = nrounds,
                 nfold = nfold, early_stopping_rounds = 10, maximize = FALSE,
                 verbose = verbose, ...)
    nrounds_use <- cv$best_iteration
  } else {
    if (nrounds_method == "v") {
      nrounds_use <- validate_xgbar(y, xreg = xreg, early_stopping_rounds = 10, params = params)$best_nrounds
    } else {
      nrounds_use <- nrounds
    }
  }


  if (verbose) {
    message("Fitting model")
  }



  model <- xgboost(data = x, label = y2, nrounds = nrounds,
                   params = params, early_stopping_rounds = 20,
                   verbose = verbose)

  fitted <- ts(c(rep(NA, maxlag),
                 predict(model, newdata = x)),
               frequency = f,
               start = min(time(target_y))
  )

  diffs <- 0

  if (season_type == "decompose") {
    fitted <- fitted * decomp$seasonal
  }


  fitted <- InvBoxcox(fitted, lambda = lambda)
  method <- paste0("xg(", maxlag, ", ", diffs, ", ")

  if (f == 1 | season_type == "none") {
    method <- paste0(method, "'non-seasonal')")
  } else {
    method <- paste0(method, "'", season_type, "')")
  }


  output <- list(y = temp_y, y2 = y2, x = x, model = model,
                 fitted = fitted, maxlag = maxlag, season_type = season_type,
                 diffs = diffs, lambda = lambda, method = method)
  if (season_type == "decompose") {
    output$decomp <- decomp
  }


  if (!is.null(xreg)) {
    output$xreg_temp = xreg_temp
    output$ncolxreg <- ncol(xreg_temp)
  }


  pred <- forecast_ap(output, h = h, xreg = pred_xreg, lambda = lambda)

  return(list(model, pred)
  )
}
