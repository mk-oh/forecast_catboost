\name{xgboost.forecast}
\alias{xgboost.forecast}
\title{Analysis with xgboost}
\usage{
cat.forecast(y, h = 42, maxlag = max(42, 4 * frequency(y)), xreg = NULL, 
    pred_xreg = NULL, nrounds = 500, params = NULL, nrounds_method = c("cv", 
        "none"), nfold = ifelse(length(y) > 720, 10, 5), boxcox = TRUE, 
    lambda = BoxCox.lambda(y), ts_clean = TRUE, ts_lambda = 1, 
    season_type = c("decompose", "dummy", "none"), verbose = TRUE)
}


\description{
 'Timeseries analysis using XGboost'.
}
\examples{
catforecast:::xgboost.forecast(AirPassengers)
}
