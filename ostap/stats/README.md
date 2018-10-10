# Core 

* [ostap.stats](README.md)

Utilities to make certain statistical calculations:
  -  Utility to combine several measurements, including correlated and uncorrelated uncertainties
  -  counters :
      - `SE`  : `Ostap::StatEntity`, regular counter, that keeps total sum, number of entries, min/max-values, mean and rms 
      - `WSE` : `Ostap::WStatEntity`, counter with the *weight*, also keeps the total sum, number of entries, min/max values, mean and rms, and also statistics of *weights* 
      - `NSE` : `Ostap::NStatEntity`, *running counter*, useful to keep the statistics of last N-entries 
  - calculation of various *statistic* for generic functions/distributions. All utilities exist in two variants: classes and standalone  functions, *e.g.* `Mode` and `mode`, `Mean` and `mean`, `Median` and `median`, etc...
      - moments 
      - central moments 
      - mean 
      - variance and rms 
      - skewness and kurtosis 
      - median 
      - mode 
      - *width/FWHM*
      - quantiles
      - symmetric and asymmetric *confidence intervals* 
  - calculation of statistics for variables/expressions from `TTree`, `RooAbsData` or `ROOT::DataFrame`
      - `data_get_moment`     : calculate the moments 
      - `data_moment`         : calculate the moments with uncertainty
      - `data_central_moment` : calculate the central moments with uncertainty
      - `data_mean`           : calculate the mean with uncertainty
      - `data_variance`       : calculate the variance with uncertainty
      - `data_dispersion`     : calculate the dispersion with uncertainty
      - `data_rms`            : calculate the RMS with uncertainty
      - `data_skewness`       : calculate the skewness with uncertainty
      - `data_kurtosis`       : calculate the (excess) kurtosis with uncertainty
      - `data_quantile`       : calculate the quantile 
      - `data_median`         : calcualte the median
      - `data_quantiles`      : calculate the quantiles  
      - `data_interval`       : calculate the interval 
      - `data_terciles`       : calcualte two terciles 
      - `data_quartiles`      : calculate three quartiles 
      - `data_quintiles`      : calculate four  quintiles 
      - `data_deciles`        : calculate nine  deciles
  - *U-statistics*, useful for *goodness-of-fit* tests
  - 2D-decorrelation transformation

