##############################
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

forward <- function(x, y, model, xregpred, object= object, i, f = 7) {
        newrow <- c(y[length(y)], x[nrow(x), 1:(object$maxlag - 1)])
        if (object$maxlag == 1) {
           newrow = newrow[-1]
        }
        if (f > 1 & object$seas_method == "dummy") {
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

validate_xgbar <- function (y, xreg = NULL, nrounds = 50, ...) {
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
                nrounds = nrounds)
        } else {
            train_md <- xgboost.ts(trainy, 
                nrounds_method = "manual", 
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

forecast_ap <- function (object, h = 42, xreg = NULL, lambda = object$lambda) {
    if (!is.null(xreg)) {
        if (is.null(object$ncolxreg)) {
            stop("model isn't used to  with xreg.")
        }

        if (class(xreg) == "ts" | "data.frame" %in% class(xreg)) {
            xreg <- as.matrix(xreg)
        }

        if (!is.numeric(xreg) | !is.matrix(xreg)) {
            stop("check - xreg")
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

    if (is.null(xreg)) {
        xreg3 <- NULL

    }
   
    f <- frequency(object$y)
    seas_method <- object$seas_method
   
    htime <- time(ts(rep(0, h), 
        frequency = f,
        start = max(time(object$y)) + 1/f))
   
 
    x <- object$x
    y <- object$y2

    for (i in 1:h) {
        tmp <-  forward(x, y, model = object$model,  
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

    if (seas_method == "decompose") {
        multipliers <- utils::tail(object$decomp$seasonal, f)
        if (h < f) {
            multipliers <- multipliers[1:h]
        }
        y <- y * as.vector(multipliers)
    }

    y <- InvMod(y, lambda = lambda)
    output <- round(y)
    return(output)
}

############################################################
xgboost.forecast <- function (y, h = 42,xreg = NULL, pred_xreg = NULL,
    maxlag = max(42, 4 * frequency(y)), nrounds = 100, nrounds_method = c("manual","cv", "v" ), 
    nfold = ifelse(length(y) > 365, 10, 5), lambda = BoxCox.lambda(y), verbose = FALSE, 
    ts_clean = FALSE, seas_method = c("decompose", 'dummy', "none"),
    K = max(1, min(round(f/4 - 1), 10)), ...)  {

    nrounds_method = match.arg(nrounds_method)
    seas_method = match.arg(seas_method)

###### Check Data format
    if (!"ts" %in% class(y)) {
       stop("class check : y is time-seires")
    }

    if (!is.null(xreg)) {
        if (class(xreg) == "ts" | "data.frame" %in% class(xreg)) {
            xreg <- as.matrix(xreg)
        }
        if (!is.numeric(xreg) | !is.matrix(xreg)) {
            stop("xreg check : class is numeric vector or a matrix")
        }
    }

    if(ts_clean == TRUE) {
        y <- tsclean(y)
    }

    f <- stats::frequency(y)
    temp_y <- y
    orign <- length(y)
    target_y <- Mod(y, lambda = lambda)

    if (length(y) < 46) {
        stop("must be data length > 45")
    }


    if (maxlag > (orign - f - round(f/4))) {
        warning(paste("y is too short ", maxlag, "--> ", 
            orign - f - round(f/4)))
        maxlag <- orign - f - round(f/4)
    }

    if (maxlag != round(maxlag)) {
        maxlag <- round(maxlag)
        if (verbose) {
            message(paste("maxlag is need to int, change to ", maxlag))
        }
    }


    if (seas_method == "decompose") {
        decomp <- decompose(target_y, type = "multiplicative")
        target_y <- seasadj(decomp)
    }

    if (seas_method == "dummy" & f > 1) {
        ncolx <- maxlag + f - 1
    }

    xreg_temp <- xreg
    n <- orign - maxlag

    y2 <- ts(target_y[-(1:(maxlag))], 
            start = time(target_y)[maxlag + 1],
            frequency = f)

    if (nrounds_method == "cv" & n < 90) {
        warning("Data is too short. need to > 90")
        nrounds_method <- "v"
    }

    if (seas_method == "decompose") {
        ncolx <- maxlag
    }

    if (seas_method == "none" | f == 1) {
        ncolx <- maxlag
    }


    x <- matrix(0, nrow = n, ncol = ncolx)
    x[, 1:maxlag] <- lag_y(target_y, maxlag, keeporig = FALSE)
    

    if (f == 1 || seas_method == "decompose" || seas_method == "none") {
        colnames(x) <- c(paste0("lag", 1:maxlag))
    }


    if (f > 1 & seas_method == "dummy") {
        tmp <- data.frame(y = 1, 
            x = as.character(rep_len(1:f, n))
            )
        seasons <- model.matrix(y ~ x, data = tmp)[, -1]
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
        cv <- xgb.cv(params= NULL,data = x, label = y2, nrounds = nrounds, 
            nfold = nfold, early_stopping_rounds = 5, maximize = FALSE, 
            verbose = verbose, ...)
        nrounds_use <- cv$best_iteration
    } else {
        if (nrounds_method == "v") {
            nrounds_use <- validate_xgbar(y, xreg = xreg, ...)$best_nrounds
        } else {
            nrounds_use <- nrounds
        }
    }

    if (verbose) {
        message("Fitting model")
    }

    model <- xgboost(data = x, label = y2, nrounds = nrounds_use, verbose = verbose)
    fitted <- ts(c(rep(NA, maxlag),
                 predict(model, newdata = x)), 
                 frequency = f,
                 start = min(time(target_y))
                 )


    diffs <- 0

    if (seas_method == "decompose") {
        fitted <- fitted * decomp$seasonal
    }


    fitted <- InvMod(fitted, lambda = lambda)
    method <- paste0("xg(", maxlag, ", ", diffs, ", ")
    
    if (f == 1 | seas_method == "none") {
        method <- paste0(method, "'non-seasonal')")
    } else {
        method <- paste0(method, "'", seas_method, "')")
    }

    output <- list(y = temp_y, y2 = y2, x = x, model = model, 
        fitted = fitted, maxlag = maxlag, seas_method = seas_method, 
        diffs = diffs, lambda = lambda, method = method)
    if (seas_method == "decompose") {
        output$decomp <- decomp
    }
   
    if (!is.null(xreg)) {
        output$xreg_temp = xreg_temp
        output$ncolxreg <- ncol(xreg_temp)
    }
    class(output) <- "xg_ts"

    pred <- forecast_ap(output, h = h, xreg = pred_xreg, lambda = lambda)

    return(list(model, 
        pred))
}