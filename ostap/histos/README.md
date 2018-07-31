# Histos

* [ostap.histos](README.md)

Collection of utilities and decorators for easy manipulations with the histograms.
 - operators for historgams, in particualr for any historgams `h`, `h1` and `h2`, 
any *function*/*callable* `f`  and any *scalar constant* `c` there are valid actions:
```
h1 + h2 , h1 - h2, h1 * h2 , h1 / h2 
h  + f  , h  - f , h  * f  , h  / f 
f  + h  , h  - f , h  * f  , f  / h
h  + c  , h  - c , h  * c  , h  / c 
c  + h  , h  - c , h  * c  , f  / c
```
