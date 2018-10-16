\name{cat.forecast}
\alias{cat.forecast}
\title{Analysis with catboost}
\usage{
cat.forecast(y, h = 42, xreg = NULL, pred_xreg = NULL, params = NULL, 
    maxlag = max(42, 4 * frequency(y)), valid = TRUE, boxcox = TRUE, 
    ts_clean = TRUE, ts_lambda = 1, season_type = c("dummy", 
        "decompose", "none"), verbose = TRUE, ...)
}


\description{
'Timeseries analysis using Catboost'.
}
\examples{
catforecast:::cat.forecast(AirPassengers)
}
