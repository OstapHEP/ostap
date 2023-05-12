# Ho to apply (Gaussian) constraints in the fit?

## The most simple case: add Gaussian constraint for a single variable

Assuming one needs to fit a data with a model consisting of a sum of "signal" 
Gaussian component and "flat" background:
```
import ostap.fitting.models as Models 

xvar  = ROOT.RooRealVar ( 'x' , 'x-variable' , 0 , 10 ) ## observable

gauss = Models.Gauss_pdf ( 'Gauss' , 
                           xvar  = xvar , 
                           sigma = ( 1 , 0.1 , 2 ) ,  ## ( value, min , max )
                           mean  = ( 5 , 4   , 6 ) )  ## ( value, min , max )
 
model = Models.Fit1D ( signal = gauss , background = 'flat' ) 
```
And assume one knows the value for Gaussian `sigma`  parameter from e.g. simulation to be $1.1\pm0.2$.

### Create the external Gaussian constraint

```
from ostap.core.pyrouts import VE ## Type  that represent value with error  
constraint = gauss.soft_constraint ( gauss.sigma  , VE ( 1.1 , 0.2**2 ) )
```

### Use the Gaussian constraint in the fit:
```
dataset = ...
result , frame = model.fitTo ( dataset , constraints = [ constraint ] )
```

## Use a multivariate Gaussianb constraints

In case one knows several parameters with their covarinace matrix, one can create multivarinate 
Gaussian constraint. Multivarinat constarin can be defiend using the parameyer values and their covariant matrix:
```
vars       = gauss.mean , gauss.sigma
values     = 5.1 , 1.1 
covariance = ... 
constraint = gauss.soft_multivar_constraint ( vars ,  ( values , covariance ) )
```
or just from `RooFitResult` object, e.g form previoss for referenc/conrl dataset:
```

## fit control dataset 
dataset_control = ...
result_control , _ = model.fitTo ( dataset_control ) 

## create constraint 
constraint = gauss.soft_multivar_constraint ( vars ,  result_control )

## use constrait:

dataset = ...
result , frame = model.fitTo ( dataset , constraints = [ constraint ] )

``` 




