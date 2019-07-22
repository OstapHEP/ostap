# Simultaneous fit 

* [ostap.fitting.simfit](SIMFIT.md)

Collection of utilities that simplify Simultaneous fit using [`RooFit`](https://root.cern.ch/roofit).

## `SimFit`

The helper PDF-like object that builds `RooSimultaneous` from set of Ostap's based pdfs 
```python

sample  = ROOT.RooCategory('Sample','sample','A','B','C')
pdf_A   = ...  ## the regular Ostap's PDF, can be 1D,2D or 3D 
pdf_B   = ...  ## the regular Ostap's PDF, can be 1D,2D or 3D 
pdf_C   = ...  ## the regular Ostap's PDF, can be 1D,2D or 3D 

sim_pdf = SimFit ( sample , { 'A' : pdf_A , 'B' : pdf_B , 'C' : pdf_C } , name = 'Sim' )
                        
```
### makeing fit 

Making the fit is trivial:
```python

dataset = ... ## data set for simultaneous fit

r , f   = sim_pdf.fitTo ( dataset ) 
```
All the argumens of `SimFit.fitTo` methdods are the same as 
for `PDF.fitTo`, e.g.
```python
r , f   = sim_pdf.fitTo ( dataset , draw = True , silent = True , sumw2 = True , ....) 
```
 
### drawing the results 

Drawing the results is delegated to the individual fit components (`pdf_A`, `pdf_B`, `pdf_C` as in this example). The methdod `SimFit.draw` require the 
first argument to code the  sampel category and the drawing variable, e.g, assuming
`pdf_A` is 1D-pdf: 
```python
frame = sim_fit.draw ( 'A'  , dataset , nbins = 100 )  ## assumong pdf_A is 1D-pdf
```
For 2D and 3D-pdfs, one needs to indicate also the drawoing variable. 
e.g. assume `pdf_B` is 2D-pdf, `PDF2`,
 with  1st variable named `mass1` and the second variable 
named `mass2`: 
```python
frame = sim_fit.draw ( 'B/1'          , dataset , nbins = 100 ) 
frame = sim_fit.draw ( 'B/x'          , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( 'B/X'          , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( 'B/mass1'      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B', 1 )      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B','1')      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B','x')      , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( ('B','X')      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B','mass1')  , dataset , nbins = 100 ) ## ditto 
```
and the same for the second variable 
```python
frame = sim_fit.draw ( 'B/2'          , dataset , nbins = 100 ) 
frame = sim_fit.draw ( 'B/y'          , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( 'B/Y'          , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( 'B/mass2'      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B', 2 )      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B','2')      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B','y')      , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( ('B','y')      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('B','mass2')  , dataset , nbins = 100 ) ## ditto 
```
In similar way  for 3D-pdfs, one can also plot the third variable, 
e.g. assume `pdf_C` is 3D-pdf, `PDF3`:
```python
frame = sim_fit.draw ( 'C/3'          , dataset , nbins = 100 ) 
frame = sim_fit.draw ( 'C/z'          , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( 'C/Z'          , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( 'C/mass3'      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('C', 3 )      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('C','3')      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('C','z')      , dataset , nbins = 100 ) ## ditto 
frame = sim_fit.draw ( ('C','Z')      , dataset , nbins = 100 ) ## ditto
frame = sim_fit.draw ( ('C','mass3')  , dataset , nbins = 100 ) ## ditto 
```


## Prepare data set for simultaneous fit :

### `combine_data` 

Helper function to create "combined" data set 
```python

dataset_A = ... # regular RooDataSet
dataset_B = ... # regular RooDataSet
dataset_C = ... # regular RooDataSet

## variables to be combined into final dataset
vars     = ROOT.RooArgSet ( var1 , var2 , .. , varN )
dataset  = combined_data ( sample , vars , { 'A' : dataset_A , 'B' : dataset_B , 'C' : dataset_C } ) 

```


This function can not combine the weighted datasets!  If weighted datastes are needed, combined non-weighed samples,  but provide the weight variable:

```python

## variables to be combined into final dataset
vars     = ROOT.RooArgSet ( var1 , var2 , SS_sw )
dataset  = combined_data ( sample , vars , 
   { 'A' : dataset_A , 'B' : dataset_B , 'C' : dataset_C } , 
   args = ( ROOT.RooFit.WeighrtVar ( SS_sw.name ) , ) )
```


### `combine_hdata` 

Helper function to create "combined" data set from the the 1D-historgams 
```python

h_A     = ...
h_B     = ...
h_C     = ...
vars    = ROOT.RooArgSet ( mass ) 
dataset = combine_hdata ( sample , vars , {'A' : h_A , 'B' : h_B , 'C' : h_C } )
```

