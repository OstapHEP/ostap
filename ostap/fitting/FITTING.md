# Fitting

Of course, the actual fits are of performed by [`RooFit`](https://root.cern.ch/roofit).
Ostap provides just user-friendly wrappers for fitting.
Almost all Ostap fitting utilities are implemented as subclasses of class `ostap.basic.PDF`.

## 1D-fits 

The 1D-fits are rather trivial 
```python
dataset = ...
pdf     = ... ## Ostap-based pdf, that inherits from class ostap.basic.PDF
## perform the fit
pdf.fitTo ( dataset )
```
A lot of optional arguments are defined
 - `draw`    : draw the fit results?
 - `nbins`   : in case of drawing, the binning scheme to be used 
 - `silent`  : silent fit?
 - `refit`   : refit in case of failure or fit problems?
 - `timer`   : measure the CPU performace
 - `args`    : tuple of `ROOT.RooCmdArg` arguments for `ROOT.RooAbsPdf.fitTo`
 - ...  
A lot of control flags for `RooFit` can be provided via additional case-insensitive keyword
arguments:
 - `verbose`   : converted to `ROOT.RooFit.Verbose ( ... )`
 - `strategy` or `minuit_strategy`  : converted to `ROOT.RooFit.Strategy ( ... )`
 - `printlevel` or `minuit_print`  : converted to `ROOT.RooFit.PrintLevel ( ... )`
 - `printallerrors` or `print_errors`  : converted to `ROOT.RooFit.PrintEvalErrors ( ... )`
 - `warnings` : converted to `ROOT.RooFit.Warnings ( ... )`
 - `weighted` or `sumw2` or `sumw2error`: converted to `ROOT.RooFit.SumW2Error ( ... )`
 - `extended` : converted to `ROOT.RooFit.Extended( ... )`
 - `ncpu` : converted to `ROOT.RooFit.NumCPU( ... )`
 - `range` : converted to `ROOT.RooFit.Range( ... )`
 - `hesse` : converted to `ROOT.RooFit.Hesse( ... )`
 - `initialhesse` : converted to `ROOT.RooFit.InitialHesse( ... )`
 - `minimizer` : converted to `ROOT.RooFit.Minimizer( ... )`
 - `minos` : converted to `ROOT.RooFit.Minos ( ... )`
 - `clone` or `clonedata` : converted to `ROOT.RooFit.CloneData ( ... )`
 - `constrain` : converted to `ROOT.RooFit.ExternalConstraints ( ... )`


The return value is a tuple of fit results `ROOT.RooFitResult` and the `ROOT.RooPlot`  (in case drawing is activated, `None` otherwise):
```python
result , frame = pdf.fitTo ( dataset , draw = True , nbins = 100 )
## draw again:
frame.draw() 
```
Alternatively the fit resutl can be visualized later:
```python
result , _ = pdf.fitTo ( dataset )
frame = pdf.draw ( dataset , nbins = 100 , ... ) 
```

For better control of visualization see [here](../plotting/FITDRAW.md).



