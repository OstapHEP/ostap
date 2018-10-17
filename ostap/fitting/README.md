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
        -   the classes get many new methods
             - fixing and releasing parameters 
```
var = ..
var.fix(1.0)  ## fix it at value=1.0
var.fix()     ## fix it at the curreent value 
var.release()
```



