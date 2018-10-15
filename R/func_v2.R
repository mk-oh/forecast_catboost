#############################
lag_y <- function (y, maxlag, keeporig = TRUE) {
  if (!is.vector(y) & !is.ts(y)) {
    stop("x must be a vector or time series")
  }

  x <- as.vector(y)
  n <- length(y)
  z <- matrix(0,
              nrow = (n - maxlag),
              ncol = maxlag + 1)

  for (i in 1:ncol(z)) {
    z[, i] <- x[(maxlag + 2 - i):(n + 1 - i)]
  }

  varname <- "x"
  colnames(z) <- c(varname, paste0(varname, "_lag", 1:maxlag))
  if (!keeporig) {  z <- z[, -1] }

  return(z)
}

########################################################
lag_xreg <- function (x, maxlag) {
  if (!is.matrix(x)) {
    stop("X needs to be a matrix")
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("Var", 1:ncol(x))
  }

  n <- nrow(x)
  m <- matrix(0, nrow = (n - maxlag), ncol = (maxlag + 1) *
                ncol(x))

  for (i in 1:ncol(x)) {
    m[, 1:(maxlag + 1) + (i - 1) * (maxlag + 1)] <- lag_y(x[, i], maxlag = maxlag)
  }
  thenames <- character()
  for (i in 1:ncol(x)) {
    thenames <- c(thenames, paste0(colnames(x)[i], "_lag", 0:maxlag))
  }
  colnames(m) <- thenames
  return(m)
}

#---- Box-Cox Translataion

Boxcox <- function (y, lambda = 1) {
  if (lambda != 0) {
    trans_y <- sign(y) * (((abs(y) + 1)^lambda - 1) / lambda)
  } else {
    trans_y = sign(y) * (log(abs(y) + 1))
  }
  return(trans_y)
}

InvBoxcox <- function (y, lambda = 1) {
  if (lambda != 0) {
    y <- ((abs(y) * lambda + 1)^(1/lambda) - 1) * sign(y)
  } else {
    y <- (exp(abs(y)) - 1) * sign(y)
  }
  return(y)
}

#---- Predict /w Roll-up

rollup_pred <- function(x, y,
                        model, xregpred,
                        object= object, i, f = 7) {
  newrow <- c(y[length(y)], x[nrow(x), 1:(object$maxlag - 1)])
  if (object$maxlag == 1) {
    newrow = newrow[-1]
  }
  if (f > 1 & object$season_type == "dummy") {
    newrow <- c(newrow, x[(nrow(x) + 1 - f),
                          (object$maxlag + 1):(object$maxlag + f - 1)])
  }
  if (!is.null(xregpred)) {
    newrow <- c(newrow, xregpred)
  }
  newrow <- matrix(newrow, nrow = 1)
  colnames(newrow) <- colnames(x)

  pred <- predict(model, newdata = newrow)
  return(list(x = rbind(x, newrow), y = c(y, pred)))
}


#############

validate_xgbar <- function (y, xreg = NULL, nrounds = 50, params = NULL,...) {
  n <- length(y)
  split_n <- round(0.8 * n)

  trainy <- ts(y[1:split_n],
               start = start(y),
               frequency = frequency(y))
  testy <- y[(split_n + 1):n]

  h <- length(testy)

  if (!is.null(xreg)) {
    trainxreg <- xreg[1:split_n, ]
    testxreg <- xreg[(split_n + 1):n, ]
  }




  xg_valid <- function(nrounds) {
    if (!is.null(xreg)) {
      train_md <- xgboost.ts(trainy,
                             xreg = xreg,
                             nrounds_method = "manual",
                             params = params,
                             nrounds = nrounds)
    } else {
      train_md <- xgboost.ts(trainy,
                             nrounds_method = "manual",
                             parms = parms,
                             nrounds = nrounds)
    }

    fc <- forecast(train_md, h = h)
    result <- accuracy(fc, testy)[2, 6]
    return(result)
  }

  mases <- sapply(as.list(1:nrounds), xg_valid)
  best_nrounds <- min(which(mases == min(mases)))
  output <- list(best_nrounds = best_nrounds, best_mase = min(mases))
  return(output)
}




######################################################
#### Function for Xgboost_ts
################### Forecast #####################################################

forecast_ap <- function (object, h = 42, xreg = NULL, lambda = object$lambda,
                         boxcox = TRUE, verbose = TRUE,...
) {

  #----. check a XREG data for predict
  if (!is.null(xreg)) {
    if (is.null(object$ncolxreg)) {
      stop("model isn't used to  with xreg.")
    }
    if (!is.numeric(xreg) | !is.matrix(xreg)) {
      stop("check - xreg")
    }
    if (class(xreg) == "ts" | "data.frame" %in% class(xreg)) {
      xreg <- as.matrix(xreg)
    }

    if (ncol(xreg) != object$ncolxreg) {
      stop("xreg != the original xreg")
    }
    h <- nrow(xreg)
    xreg2 <- lag_xreg(rbind( object$origxreg,xreg), maxlag = object$maxlag)
    # xreg2 <- lag_xreg(rbind(xreg, object$origxreg), maxlag = object$maxlag)
    nn <- nrow(xreg2)
    xreg3 <- xreg2[(nn - h + 1):nn, ]
  }

  #----. non-xreg is Pass
  if (is.null(xreg)) {
    xreg3 <- NULL

  }

  f <- frequency(object$y)
  season_type <- object$season_type

  htime <- time(ts(rep(0, h),
                   frequency = f,
                   start = max(time(object$y)) + 1/f))


  x <- object$x
  y <- object$y2

  for (i in 1:h) {
    tmp <-  rollup_pred(x, y, model = object$model,
                        xregpred = xreg3[i, ], object= object, i = i)
    x <- tmp$x
    y <- tmp$y
  }

  y <- ts(y[-(1:length(object$y2))],
          frequency = f,
          start = max(time(object$y)) + 1/f)

  if (object$diffs > 0) {
    for (i in 1:object$diffs) {
      y <- ts(cumsum(y), start = start(y), frequency = f)
    }
    y <- y + Mod(object$y[length(object$y)], lambda = lambda)
  }

  if (season_type == "decompose") {
    multipliers <- utils::tail(object$decomp$seasonal, f)
    if (h < f) {
      multipliers <- multipliers[1:h]
    }
    y <- y * as.vector(multipliers)
  }

  y <- InvBoxcox(y, lambda = lambda)
  output <- round(y)
  return(output)
}
