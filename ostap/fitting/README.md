# Fitting

* [ostap.fitting](README.md)

Collection of various utilities that simplify  the communications with [`RooFit`](https://root.cern.ch/roofit)
 - decorations for many  native `RooFit` classes, namely
   - variables   : `RooAbsReal`, `RooRealVar`, ... 
   - collections : `RooArgSet`, `RooArgList` , ...
   - datasets    : `RooAbsData`, `RooDataSet`, `RooDataHist`
   - minuit 
   - auxillary classes : `RooFitResults`, ... 
 - utilities for easy creation of fitting models and easy fitting and plotting  
   - base classes : `PDF`, `MASS`, ... 
   - compound ( signal(s) + background(s) ) models:
       - `Fit1D`
       - `Fit2D`, `Fit2DSym`
       - `Fit3D`, `Fit3DSym`, `Fit3DMix`
   - useful wrappers for bare `RooAbsPdf` objects: 
       - `Generic1D_pdf` 
       - `Generic2D_pdf` 
       - `Generic3D_pdf` 
   - other: `PyPDF`, *resolution*, *convolution*, *efficiencies* , etc..
- zillions of specialized pdfs:
   - *signal/peak-like*  
   - generic smooth background 
   - physics-inspired smooth background models 
   - many widely knows distibutions 


 - [rootfit.py](roofit.py): *head* module for varioud decorations of `RooFit`-objects 
     - [variables.py](variables.py) -  collections of decorations for `RooAbsReal`, `RooRealVar` and related clases
        * the classes get many new methods: 
```python
var = ...
var.fix(1.0)  ## fix it at value=1.0
var.fix()     ## fix it at the curreent value 
var.release()

ve = var.ve    () ## convert to ValueWithError
ve = var.asVe  () ## ditto  
ve = var.as_VE () ## ditto 

if 2.5 in var3 :  ## check if  2.5 is within min/max range 
      ... 
mn , mx = var. minmax() ## get min/max range, is applicable 
mn , mx = var.xminmax() ## ditto

h  = var.histo( bins = 50 ) ## create the correspondiong histogram```


         * trivial math-operations  (`ValueWithError` as return value)

```python
var + 1.0 
var - 2.1 
var * 3.2
var / 4.6
5.1 + var 
6.2 - var 
7.3 * var 
8.4 / var 
var ** 9.5 
1.2 ** var```
         
          * also trivial properties are defined:
             *  `value`  (returns `float` or `ValueWithError`, depending on the type) 
             *  `error`  (returns `float`)
```python
var.value 
var.error
var.value = 10```

        * `SETVAR` - useful context manager to preserve the value of `RooRealVar`

     
