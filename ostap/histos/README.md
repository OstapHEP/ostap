# Histos

* [ostap.histos](README.md)

Collection of utilities and decorators for easy manipulations with the histograms & graphs
- bin-access, looping, ranges, etc...
- operators for historgams, in particualr for any historgams `h`, `h1` and `h2`, 
any *function*/*callable* `f`  and any *scalar constant* `c` there are valid actions:
```
h1 + h2 , h1 - h2, h1 * h2 , h1 / h2 
h  + f  , h  - f , h  * f  , h  / f 
f  + h  , h  - f , h  * f  , f  / h
h  + c  , h  - c , h  * c  , h  / c 
c  + h  , h  - c , h  * c  , f  / c
```
- functions: ` h1 = sin ( h ) `
- histogram as functions (local interpolation): ` v = h ( 0.1 ) `
- specific methods to get variosu efficienciy estimates, including *binomial efficiency* and various types of binomial efficiency intervals
- histogram parameterizations 
     - as *Legendre* sum
     - as *Chebyshev* sum
     - as *Bernstein* sum with various variant
         - regular *Bernstein* sum 
         - *Bernstein* sum, forced to be positive  
         - *Bernstein* sum, forced to be positive and monotonic 
         - *Bernstein* sum, forced to be positive, monotonic with the fixed sign of the second derivative 
         - *Bernstein* sum, forced to be positive   with the fixed sign of the second derivative 
     - as *b-spline*
         - regular *b-spline* 
         - positive *b-spline* 
         - positive monotonic *b-spline* 
         - positive monotonic *b-spline* with the fixed sign of the second derivative  
         - positive *b-spline* with the fixed sign of the second derivative  
     - approximation using `RooAbsPdf`    
- many methods for historgam comparison,  inclusing bin-by-bin and shape comparisons 
- summation and integration methods 
- transformation to various type of graphs 
- ... 
       