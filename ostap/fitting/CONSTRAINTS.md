# Ho to apply (Gaussian) constraints in the fit?

## The most simple case: Gaussian constraint for single variable

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

## Multivariate Gaussian constraint

In case one knows several parameters with their covarinace matrix, one can create multivarinate 
Gaussian constraint. Multivarinat constarin can be defiend using the parameyer values and their covariant matrix:
```
vars       = gauss.mean , gauss.sigma
values     = 5.1 , 1.1 
covariance = ... 
constraint = gauss.soft_multivar_constraint ( vars ,  ( values , covariance ) )
```
or just from the `RooFitResult` object, e.g from previos fit to the reference/control dataset:
```


dataset_control = ...                                 ## control dataset
result_control , _ = model.fitTo ( dataset_control )  ## fit to control dataset 

constraint = gauss.soft_multivar_constraint ( vars ,  result_control ) ## create constraint 

dataset = ...
result , frame = model.fitTo ( dataset , constraints = [ constraint ] ) ## use contraint!

``` 


## practical example: Determine the resolution scale factor

Assume that resolution on variable `x` can be different slightly
 different for nominal and reference/control dataset. In this case still
can use the constrains but add additional scale fatcor for the resolution
```
scale = ROOT.RooRealVar ( 'scale' , 'scale-factor for resolution' , 1.0 , 0.5 , 2.0 ) ## resolution scale factor

corrected_sigma = gauss.vars_multiply ( gauss_sigma , scale , name = 'sigma_corrected' ) ## corrected_sigma = sigma * scale 

corrected_gauss = M.Gauss_pdf ( 'CorrectedGauss'          , 
                                 xvar = x                 , 
                                 sigma  = corrected_sigma , ## attention! 
                                 mean   = gauss.mean      ) ## attention! 
 
corrected_model = M.Fit1D ( signal = corrected_gauss , background = 'flat' , suffix = 'corrected' ) 
 
``` 
And then one can fit to the nominal dataset with constrains 
from the reference/control dataset to get 
value of scale factor:
```
result1 , _ = corrected_model.fitTo ( dataset , constraints = [ constraint ] )
```
In case one know the scale factor, e.g $1.0\pm0.1$, 
this knowledge can be explointed via the external constraint
```
scale_constraint = gauss.soft_constraint ( scale , VE ( 1.0 , 0.1**2 ) ) ## create consraint for scale parameter 

result2 , _ = corrected_model.fitTo ( dataset , constraints = [ constraint, scale_constraint  ] ) ## use two constraints 

```




