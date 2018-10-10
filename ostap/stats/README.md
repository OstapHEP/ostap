# Statistics 

* [ostap.stats](README.md)

Utilities to make certain statistical calculations:
  - [combine.py](combine.py): utility to combine several measurements, including correlated and uncorrelated uncertainties
      - [P.Avery *Combining measurements with correlated errors*, CBX 95 55](http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz)
      - [P.Avery *Combining measurements with correlated errors*, CBX 95 55](http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors)
  - [counters.py](counters.py) :
      - `SE`  : `Ostap::StatEntity`, regular counter,  keeps 
           - total sum
           - number of entries
           - min/max-values
           - mean 
           - rms
           - (efficiency) 
      - `WSE` : `Ostap::WStatEntity`, counter with the *weight*, also keeps 
           - the total sum
           - number of entries
           - min/max values
           - weighted mean 
           - weighted rms
           - full statistics of *weights* 
      - `NSE` : `Ostap::NStatEntity`, *running counter*, useful to keep the statistics for the last N-entries 
  - [moments.py](moments.py): calculation of various *statistic* for generic functions/distributions. All utilities exist in two variants: classes and standalone  functions, *e.g.* `Mode` and `mode`, `Mean` and `mean`, `Median` and `median`, etc...
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
  - [statvars.py](statvars.py): calculation of statistics for variables/expressions from `TTree`, `RooAbsData` or `ROOT::DataFrame`
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
  - [ustat.py](ustat.py): *U-statistics*, useful for *goodness-of-fit* tests
      - [M.Williams *How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics*](https://doi.org/10.1088/1748-0221/5/09/P09004)
      - [M.Williams *How good are your fits? Unbinned multivariate goodness-of-fit tests in high energy physics*](http://arxiv.org/abs/arXiv:1003.1768)
  - [corr2d.py](corr2d.py) : 2D-decorrelation transformation for pair of variables from `TTree` or `RooAbsData`

