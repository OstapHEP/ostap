# Fill `RooDataSet` from `ROOT.Tree`

Almost all [`RooFit`](https://root.cern.ch/roofit) actions deal with `RooDataSet`. Typically one fils `RooDataSet` objects form `ROOT.TTree`. There are several  ways to do it 


1.  Use the corresponsing constructor of `RooDataSet`: native way, very efficient, but limited to the trvivial scalar columns and rather limited filtering functionality   

2.  Via explicit loop over `ROOT.TTree` entries and using `RooDataSet.add` method: very powerful and flexible method, can deals with arboyrary complex columns, expressions,  functions, etc...  but very CPU  inefficient for large data. Large speedup can be obtaine dvia efficicient filetring using `withCuts`-iterator

3.  Using `Selectors` from [`selectors.py`](./selectors.py) module: set of very flexible and powerful technique. For simple cases is falls to (1), but in additon, it allowes to perform very  efficient filtering, exloit the powerful expressions and functions, inclusing `ROOT.RDataFrame` filtering. With the proper usage it could be very CPU efficient.  In addition it allows powerful paralellization, exploiting multicore and multiprocessor architechture

## Simple constructor of `RooDataSet`  and `make_dataset` method 

```python
tree      = ...                    ## the tree
varset    = ROOT.RooArgSet ( ... ) ## the booked variables 
cuts      = 
dataset   = RooDataSet     ( 'ds1' , 'my data' , tree , varset , cuts ) 
```
This is very CPU efficient method, but it has certain limitations: 
 - it deals only with the simple scalar variables in `ROOT.TTree`
 - `cuts` operates only on the variables that are supposed to be in dataset. One cannot use requiremenet on other variables. All variables used in `cuts` must be seelcted for the dataset. It can make the dataset heavier than needed, and it  results in no-optimal CPU performace for such case. 

There is a bit more useful, flexible and powerful variant of it, using `make_dataset` method 

```python
tree = ...
ds = tree.make_dataset ( variables = [ 'pt', 'p' , 'eta' ] , 
                         selection = "abs(ID)==121"        )
``` 
In this approach `selection` can include any  other (scalar) variables in `ROOT.TTree`. 

Also in this approach one can insert another (derived) variables:
```python
tree = ...
ds = tree.make_dataset ( variables = [ ( 'pt2' , 'pT squared' , 0 , 10 , 'pt*pt' ) , 
                         'p' , 'eta' ] , 
                         selection = "abs(ID)==121"        )
``` 
In this example the calculated/derived variable `pt2` will be added into `RooDataSet`.
Any valid `RooFormulaVar` expressions can be used in this approach to define new variables. 
```python
tree = ...
from ostap.fitting  .selectors import Variable 
ds = tree.make_dataset ( variables = [ Variable( var = 'pt2' , 
           description = 'pT squared' , vmin  = 0 , vmax = 10 , acessor = 'pt*pt' ) , 
                         'p' , 'eta' ] , 
                         selection = "abs(ID)==121"        )
``` 

## Explicit loop over `ROOT.TTree` entries 

```python
m1     = ROOT.RooRealVar('m1','mass of 1st particle',3,3.5)
m2     = ROOT.RooRealVar('m2','mass of 2nd particle',4.5,6)
varset = ROOT.RooArgSet ( m1 , m2 )
ds     = ROOT.RooDataSet ( 'ds' , 'data' , varset )  

tree   = ...
for entry in tree :

    ... put some cuts here... 

    ... get/calculate 'm1'
    m1_value = ...
    if not m1_value in m1 : continue 
    m1.setVal ( m1_value ) 

    ... get/calculate 'm2'
    m2_value = ...
    if not m2_value in m2 : continue 
    m2.setVal ( m2_value ) 

    ds.add ( varset )  
```
Definitley this is very powerful method, and arbitrary complicated columns,
functions and expression can be used here..
The technique can be speedup by embedding `cuts` direclty into the loop:

```python

tree   = ...
for entry in tree.withCut ( '( pt>10 ) && (y<5) && abs(ID)== 121' ) :

    ... get/calculate 'm1'
    m1_value = ...
    if not m1_value in m1 : continue 
    m1.setVal ( m1_value ) 

    ... get/calculate 'm2'
    m2_value = ...
    if not m2_value in m2 : continue 
    m2.setVal ( m2_value ) 

    ds.add ( varset )  
```
In case of  harsh cuts it can be very CPU efficient, since cuts are applied 
on `ROOT.TTree` level before python enters the game. 


## Using selectors
This method could be CPU efficient, combining with the great flexibility and power. 

```python
## prepare the list of variables 
variables = [  
    'mfit_Lbs_cc'                   ,
    'pt23'                          ,
    Variable ( 'tmva_3_Lb_BDTG_response' , '' , -0.5 , 0.5 ) ,
    ( 'm12cc', '', 0,      7 )      ,
    'm13cc'                         ,
    'm23cc'                         ,
    ( 'ptk'  , '' , 0 , 4   , 'pt_kaon_lb[0]'  )  ,
    ( 'etak' , '' , 2 , 4.9 , 'eta_kaon_lb[0]' )  ,  
    Variable( 'ququ' , '' , 0 , 10 , acessor = lambda s : s.ququ+s.pt/s.eta )  ,  
    ]

selector = SelectorWithVars ( 
   variables = variables ,
   selection = "...  the cuts are here ... " ,      
   cuts      = lambda s : 1000 < s.run_number < 1000000 ) ## global python cuts 

tree = ...
tree.process ( selector ) 

dataset = selector.data 
```
As one sees, it is a bit verbose technique, 
but one can define arbitrarily complicated columns, functions and expressions.
Note that variables can be renames, and oen can get scalar content from 
array-like variable. 
The  variable is called *trivial* if 
 - it directly refers to the scalar item in `ROOT.TTree`
 - `accessor` function for `Variable` class is a valid `TTreeFormula/RooFormulaVar` expression

```python
## the list of trivial variables 
variables = [  
    'mfit_Lbs_cc'                   ,
    'pt23'                          ,
    Variable ( 'tmva_3_Lb_BDTG_response' , '' , -0.5 , 0.5 ) ,
    ( 'm12cc', '', 0,      7 ),
    Variable ( 'aaa' , '' , -0.5 , 0.5 , 'a*b/c+sin(d) )' ) 
    ]
```

If all variables are  *trivial* and not global python cuts are defined,  very agressive optimization is applied, that, in particular, can include the intermediate 
operations with `ROOT.RDataFrame.Define`, `ROOT.RDataFrame.Filter` and `ROOT.RDataFrame.Snapshot`.

## Using parallel selectors



## Special cases 

There are few special cases, where one needs to add non-trivial variable into `RooDataSet`

1. Adding result of `TMVA` or `TMVA/Chopping` decision 
2. Adding result of `reweighting` procedure 

Both cases are non-trivial and CPU expensive. 
For them two-step procesure is suggested 
1. Fill `RooDataSet` in regular way without thie variables (but keeping in `RooDataSet` the varibales needed to calculate `TMVA`, `TMVA/Chopping` or `reweighting` results
2. Apply dedicated method to add new variables to `RooDataSet`

###  `TMVA` or `TMVA/Chopping` 

Here the dedicated method is `addTMVAResponse`

```python
dataset       = ...
weights_files = ... ## tar-file from `TMVA`
addTMVAReponse ( dataset , 
                 input_variables       , ## list of input variables 
                 weights_files         , ## input files with TMVA weights 
                 prefix  = 'TMVA_'     , 
                 suffix  = '_response' ,
                 verbose = True        )
```
 or `addChoppingResponse` for `TMVA/Chopping`
```python
dataset       = ...
weights_files = ... ## tar-file from `TMVA`
addChoppingReponse ( dataset , 
                     chopper = ' ... '     , ## chopper categrory/fnuction 
                     N       = 11          , ## number of chopping categories 
                     input_variables       , ## list of input variables 
                     weights_files         , ## input files with TMVA weights 
                     prefix  = 'TMVA_'     , 
                     suffix  = '_response' ,
                     verbose = True        )
```
For both methods the parallelization is applicable 


###  `reweighting`

```python
weighter = ...
dataset  = ...
dataset.add_reweigthing ( weighter , name = 'weight')
```
 
