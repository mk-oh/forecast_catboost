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
