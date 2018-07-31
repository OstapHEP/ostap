# Fitting

* [ostap.fitting](README.md)

Collection of various utilities that simplify  the communications with [`RooFit`](https://root.cern.ch/roofit)
 - decorations for many  native `RooFit` classes, namely
   - variables   : `RooAbsReal`, `RooRealVar`, ... 
   - collections : `RooArgSet`, `RooArgList` , ...
   - datasets    : `RooAbsData`, `RooDataSet`, `RooDataHist`
   - minuit 
   - auxilalry classes : `RooFitResults`, ... 
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
   - other: `PyPDF`, *resolution*. *efficiencies* , etc..
- zillions of specialized pdfs:
   - *signal/peak-like*  
   - generic smooth background 
   - physics-inspired smooth background models 
   - many widely knows distibutions 


