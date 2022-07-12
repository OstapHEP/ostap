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
   - base classes : `PDF`, `PEAK`, ... 
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


## [rootfit.py](roofit.py)

The main module for various decorations of `RooFit`-objects 

### [variables.py](variables.py) 
Collections of decorations for `RooAbsReal`, `RooRealVar` and related clases:
      - trivial math-operations  (`ValueWithError` as return value) 
      - other useful methods and properties 
      - useful context manager `SETVAR` to preserve the value of `RooRealVar`
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
1.2 ** var

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

h = var.histo( bins = 50 ) ## create the correspondiong histogram

var.value 
var.error
var.value = 10

print float(var)
with SETVAR ( var ) : 
    var.value = 10 
    print float(var)
print float(var) 
```

### [dataset.py](dataset.py) 
Collection of decorations for `RooAbsData`, `RooDataSet` and related clases:
      * iterators, getters, slices, trivial operators
      * drawing and projections 
      * statistic for the variables and expressions 
 
```python
data    = ...
entry   = data[14]            ## get certain tenr
entries = data[:500]          ## get the first 500  events 
entries = data[1:-1:2]        ## get every second event
entries = data.sample ( 100 ) ## select  (randomly) 100 events from dataset
entries = data.shuffle ()     ## get the shuffled data 

if 'pt' in data : ...         ## check existcen of variable in dataset      

data2 = ... # another dataset (of the same structure) 

d = data + data2   ## merge  two data sets of the same structure  
data += data2      ## append dataste of the same structure 

data3 = ... # another dataset (of the same lenfth) 

d = data * data3   ## merge two dataset of the same length 

d1 = data * 0.1    ## get (random) 10% the dataset 
d2 = data  / 15    ## get (random) 1/15 of the dataset 
d3 = data  % 20    ## get (random) 1/20 of the dataset  

s1     = data.statVar('pt*pt/mass', 'ez>1000') ## get the statistics for the expression
s2     = data.statCov('pt', 'ez/mass')         ## get the correlations 
s3     = data.sumVar ('pt*weight/ez')          ## the sum for the expression
mn, mx = data.vminmax('pt', 'ez>100')          ## get min/max for the expression
 
h1 = ... ## book the histo
h  = data.project ( h1 , 'pt/ez' , 'mass>100' ) ## make a projection to the histogram 

data.draw('pt/ez','mass>100') ## draw the expression  
 
dsw = data.makeWeighted ( '1/eff' ) ## create the weighted dataset   

s1     = data.statVar('pt*pt/mass', 'ez>1000') ## get the statistics for the expression
s2     = data.statCov('pt', 'ez/mass')         ## get the correlations 
s3     = data.sumVar ('pt*weight/ez')          ## the sum for the expression
mn, mx = data.vminmax('pt', 'ez>100')          ## get min/max for the expression

print data.moment         ( 'mass*mass' , 3   ) ## 3rd moment
print data.central_moment ( 'mass*mass' , 3   ) ## 3rd central moment
print data.skewness       ( 'mass*mass'       ) ## skewness 
print data.kurtosis       ( 'mass*mass'       ) ## kurtosis 

print data.mean           ( 'mass*mass'       ) ## mean
print data.variance       ( 'mass*mass'       ) ## variance 
print data.dispersion     ( 'mass*mass'       ) ## variance 
print data.rms            ( 'mass*mass'       ) ## rms 
print data.quantile       ( 0.3 , 'mass*mass' ) ## 30% quantile 
print data.quantiles      ( ( 0.1 , 0.9 ) , 'mass*mass' ) ## 10 and 90% quantiles 
print data.interval       (  0.05 , 0.95  , 'mass*mass' ) ## 90% interval 
print data.median         (  'mass*mass'      ) ## median
print data.terciles       (  'mass*mass'      ) ## terciles 
print data.quartiles      (  'mass*mass'      ) ## quartiles
print data.quintiles      (  'mass*mass'      ) ## quintiles
print data.deciles        (  'mass*mass'      ) ## deciles

print data.branches() , data.leaves() ## list variables 
print data                            ## print dataset as table 
print data.table  ()                  ## print dataset (another format) 
print data.pprint ()                  ## print dataset (another format)
```

### [roocollections.py](roocollections.py) 
Collection of decorations for `RooArgSet`, `RooArgList` and related clasese
    * iterators, getters, printing, some operators, ...
```python
argset = ...
for a in argset : print a 

argset['pt'] ## get ittem 
argset.pt    ## get attribute/ditto  

if  'pt' in argset : ...  ## check the presence 

arglist = ...
for a in arglist : print a 

arglist[1] ## get the item 
```

### [roofitresult.py](roofitresult.py) 
Collection of decorations for `RooFitResult` clasese
    * printing, access to fit parameters, operations with fit parameters
```python
result = ... ## RooFitResult instance 

print result.params()                      ## get parameters 
print result.params( float_only = False )  ## get parameters 

results.params()['Signal'][0]  ## parameter value (+error) as ValueWithError
results.params()['Signal'][1]  ## parameter as RooAbsReal/RooRealVar

results.param('Signal')[0]  ## parameter value (+error) as ValueWithError
results.param('Signal')[1]  ## parameter as RooAbsReal/RooRealVar

print result.Signal  ## get parameter as attribute 

for p in reuslts : print p                       ## loop over parameters 

for name,v in results.iteritems() : print name,v ## loop over parameters

c = result.corr ('Signal','Mean') ## correlation coefficient betweem two parameters 

cov = result.cov_matrix('Signal', 'Mean', 'Sigma') ## get correlarion sub-matrix 
cov = result.cov('Mean', 'Sigma') ## get covariance for two elements 


v = result.sum      ('Signal','Background') ## calculate sum for two parameters 
v = result.multiply ('Signal','Background') ## calculate 
v = result.divide   ('Signal','Background') ## calculate 
v = result.subtract ('Signal','Background') ## calculate 
v = result.fraction ('Signal','Background') ## calculate Signal/(Signal+Background)
```

## [simfit.py](simfit.py)

Collection of utilities that simplify Simultaneous fit using [`RooFit`](https://root.cern.ch/roofit).
See [here](SIMFIT.md) for more details. 