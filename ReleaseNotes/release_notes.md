## New features 

  1. Add `FlattePS_pdf` - similar to `BWPS_pdf`
  1. Make few important steps towards ROOT 6.23/01 ("New PYROOT"). Full  adaptation is not yet achieved, there   are some pending problems with effective inheritance from C++ classes (namely `TSelector`, `PyPDF`, `PyVAR`, ...). There are also some puzzling  crashes... 
  1. rename  tests, make test selection more transparent and easy to navigate back 
  1. rename `ostap/fitting/selectors.py` to `ostap/fitting/pyselectors.py` to avoid the name clash for python3 
  1. add helper script `pplaunch` to launch remote pp-servers via ssh tunnels
  1. update `PyVar`, `PyVar2`, `PyPdf`,`PyPdf2`
  1. (almost) complete update for new PyROOT 
  1. fix `test_fitting_minuit_weighted` - thanks to Dima Golubkov
  1. make more  coherent treatment of ROOT issues 
  1. fix `minuit` for new PyROOT  (signature of `FCN` is different!)
  1. minor update for `minuit` : from now allow access  by parameter name:
```
minuit        = ...
minuit['p2']  = 10 
minuit.minos   ('p1','p2',...)
minuit.release ('p3')
```   
  1. update `pptunnel` + `pplaunch` with better and more informative output 
  1. add `Ostap::Math::BW3L`, `Ostap::Models::BW3L`,  `BW3L_pdf` and extend  test `ostap/fitting/tests/tests_fitting_breiwigner.py` - resurrected version of the Breit-Wigner profile from 3-body decays 
  1. add check for `more_itertools`, provide   replacement for `chunked` when `more_itertools`  is not available 
  1. update `Ostap::Math::NSphere` and `Ostap::Math::Positive` such that for null-parameters the reusltin poisitve function is a constant. The trick is based on properties of Chebyshev polynomials of 1st,. 2nd, 3rd and 4th kind.   
  1. `ostap.logger.table` add parameter `alignment` that specifies the column alignment.
  1. `ostap.parallel.task` add parameter `batch`, that allows to execute the tasks in `batch` mode 
  1. re-add  generic Breit-wigner channel `Ostap::Math::ChannelGeneric` (for ROOT>=6.23/01 only) 
  1. add method `amplitude` for `Ostap::Math::ChannelBW`
  1. re-remove  generic Breit-wigner channel `Ostap::Math::ChannelGeneric`
  1. change the default `sample` argument for `PDF.generate` from `False` to `True`
  1. tiny fix for the table column alignment 
  1. add `#include <string>` for `NSphere.h` - for  certain configuration it prevents compile error (thanks to Abdul-Kerim Gusseinov) 
  1. Add parameter `accept_fun` for  `toys` - that allows to (re)define the accepance criteria, the default corresponds to `accept_fit` function from `ostap.fitting.toys` module, that checks the    fit status (0) and covariance matrix status ( 3 or -1)   
  1. Add parameter `fit_fun` for  `toys` - that allows to (re)define the default "fit"-policy 
  1. Add parameter `gen_fun` for  `toys` - that allows to (re)define the default "generation"-policy 
  1. `test_plotting_summary_graph.py` : add call for `ROOT.gPad.RedrawAxis` -  thanks to Tom Blake 
  1. add `**kwargs` for all `parallel`-methods, arguments are used for `WorkManager`
  1. remove `evaluateBatch` form all PDFs (folloiimng evolution of ROOT). We need to  gradually introduce `evaluateSpan`
  1. `Ostap::Math::GammaBW3`: use 1/s factor instead of 1/s^3/2. thanks to Misha Mikhasenko! 
  1. Add symmetic   Sinh-Asinh resolution  model `ResoSinhAsinh`
  1. Add symmetic Jonhson's SU resolution  model `ResoJohnsonSU`
  1. fix for the new signature of `TDirectory::CurrentDirectory()` method 
  1. fix/rewrite/improve `ostap.histos.compare` module 
  1. improve reweighting machinery: make it more tunable (and less automatic) and more suitable for multidimensional reweighting. 
  1. Add new context manager `SETPARS` and use it in `PDF.wilks` , `PDF.wilks2` , `PDF.draw_nll`, `PDF.graph_nll`, `PDF.graph_profile` 
  1. Add logistic/"sech-squared" resolution model `ResoLogistic`
  1. Improve `PDG.graph_profile/PDF.graph_nll` : add `draw` argument to draw the graph in progress  

## Backward incompatible changes: 
  1. Reweighting machinery: different signature of `makeWeights` function -  new argument `make_plots`, different meaning of argument `power`,  different return value

## Bug fixes:
  1. fix bugs in `Ostap::Math::BWPS`


# v1.5.0.4

## New features 

  1. add `Ostap::Math::BWPS`, `Ostap::Models::BWPS` and `BWPS_pdf` - function for Breit-Wigner profile, modulated with additional phase-space factors and polynomial degrees of freedom.
 
## Backward incompatible changes: 

  1. Change parameter name `mean` to `m0` for `BreitWigner_pdf`, `BWMC_pdf`, `Voigt_pdf`, `PseudoVoigt_pdf`, ...
 

## Bug fixes:

# v1.5.0.3

## New features 
  1. add argument `callable` for models plotting and their conversion to `TF1`. It allows to draw some derived quantitites 
```
bw = Ostap.Math.BreitWigner( ... )
bw.draw ( xmin = ... , xmax = ... ) ## draw the Breit-Wigner lineshape
bw.draw ( xmin = ... , xmax = ... , callable = lambda x : bw.amp ( x ).real ) ## draw real part of amplitude
bw.draw ( xmin = ... , xmax = ... , callable = lambda x : bw.amp ( x ).imag ) ## draw real part of amplitude
bw.draw ( xmin = ... , xmax = ... , callable = lambda x : cmath.phase ( bw.amp ( x ) ) ## draw the phase 
``` 
  1. add  tiny utilities `lrange` and `log_range` (in addition to existing `vrange`) into  `ostap/utils/utils.py`
```
for x in vrange ( 0.0 , 10.0     , 10 ) : print x  ## "lin-range"
for x in lrange ( 1.0 , 10.0**10 , 10 ) : print x  ## "log-range"
```
  1. add methods `amp_real`, `amp_imag`, `amp_phase` for the Breit-Wigner-like models
  1. add Argand plot for the Breit-Wigner-like models
```
bw = Ostap.Math.BreitWigner(... )
ap = bw.argand ( xmin = ... , xmax = ... , npx = 500 )
ap.draw('alc')
```
  1. add utilities for better visuzalisation of Dalitz densities 
  1. more improvements for Dalitz plot vizualization 
  1. more tweaks for `Ostap::Math::DalitzIntegrator`
  1. more tweaks for `Ostap::Math::BW`
  1. add `bb` ("bounding box") method for `ROOT.TGraph`-like objects.
  1. more  tweaks for `ROOT.RooMinimizer`, in particular better control over printout 
  1. add `PDF.minuit`: add FCN scaling for weighted dataset 
  1. re-add checks for `SumW2/Asymptotic` checks for the `PDF.fitTo` for weighed datasets 
  1. add new test `test_fitting_minuit_weighted`
  1. add two specific  cases for `Ostap::Math::PhaseSpaceNL`

## Backward incompatible changes: 

## Bug fixes:


# v1.5.0.2

## New features 

 1. Simplify interface for `Ostap::Math::Integrator` and `Ostap::Math::DalitzIntegrator` classes: essentially remove large duplication, the tag/label argument for caching is now the last one and always "optional" - no caching is performed if  argument is  zero (default)
 1. add methods to create `Ostap::Math::ChebyshevSum` from `Ostap::Math::ChebyshevApproximation`

## Backward incompatible changes

## Bug fixes:
 1. fix a bit strange "feature" with "derived" variable in `RooDataSet` (thanks  to Alexander Artamonov for reporting).  It happens that `RooDataSet::addColumn`RooDataSet::addColumns` behave a bit differently. the  first one issues the error message and variable behaves weirdly. 
  

# v1.5.0.1

## New features 
 1. `parallel/task.py`  : change master/slave to main/secondary                    (request from Bogdan Popovici)
 1. Modify a bit printout for `Ostap::StatEntity` and `Ostap::WStarEntity` classes (request from Alexey Dzyuba) 
 1. tiny tweak for `ostap.fitting.basic.all_args`
 1. add `all_integers`, `all_numerics` and `all_strings` to `ostap/core/ostap_types.py`

## Backward incompatible changes

## Bug fixes:
 1. fix a bit strange "feature" in the function `make_dataset` from `ostap/fitting/selectors.py` (thanks to Alexander Artamonov for reporting it)

# v1.5.0.0

## New features 
  1. Make`Ostap::Math::Choose` a bit  more  efficient
  1. add `Ostap::Math::choose_array` to get array of binomial coefficients (compile time)
  1. add templated central moments `Ostap::Math::Moment_<N>`
  1. add their python decorators `ostap.stats.moment`
  1. add test for moment-counters `test_stats_moment.py`
  1. add templated weighted moment counters  `Ostap::Math::WMoment_<N>`
  1. large modificatons in `Ostap::Kinematics::Dalitz`
  1. `fitting.basic` : add intermediate mase class `MASSMEAN` that does not hold `sigma`
  1. add `mean_name`, `mena_title` , `sigma_name` and `sigma_title` for `MASSMEAN` and `MASS` base classes : it allows to remove many ugly lines with post-fix for the  variable names
  1. make use of `mean_name`, `meean_title`, `sigma_name`, `sigma_title` for many `PDFs` 
  1. Remove `sigma`(`gamma`) from `Flatte_pdf`
  1. extend interface for `ostap.fitting.simfit.SimFit`, allowinng usage of it for toys
  1. add test `test_fitting_toys_simfit.py`
  1. further extend `Ostap::Math::DaltzIntegrator`
  1. further extend `Ostap::Kinematics::Daltz0` and `Ostap::Kinematics::Daltz` (add more invariants)
  1. extend `Ostap::Math::ChebyshevApproximation` (add scale and bias operators)
  1. add `Ostap::Math::Piecewise` function 
  1. Improve `Ostap::Math::ChebyshevApproximation`
  1. further extend `Ostap::Math::DaltzIntegrator`
  1. extend generic functions, add generic PDFs 
  1. add `binnig` functions to create `RooBinning`
  1. add `ostap.fitting.morphing_pdf` with two PDFs for morphing 
  1. add test for new morphing PDF 
  1. add possiility to use regex for `compressed_shalve.ikeys` method
  1. add integration over s,s2 variables in `Ostap::Math::DalitzIntegrator`
  1. add datetime to the logger format for non-isatty output, e.g. log-files...
  1. add `Ostap::Math::KramersKronig` helper class 
  1. tiny tweaks for `tootshelve`
  1. suppress error prints from  `selectors.valid_formula` 

## Backward incompatible changes
  1. rewrite `Ostap::Math::DalitzIntegrator`
  1. rewrite `Ostap::Math::Integrator`
  1.`Flatte_pdf`: rename arguments and attributes
  1. rename  `ostap.fitting.basic.Resoluton` to `ostap.fitting.basic.CheckMean` and invert its argument 
  1. total re-write of all Breit-Wigner related stuff and in particular temporarily remove all beast like LASS, Bugg, etc...   
  1. fix but with parsing arguments of `PDF.fitTo` : fro certain number of argument the creation of `RooFit::MultiArg` was incorrect (thanks to Pavel Krokovny)

## Bug fixes:
  1. fix missing `hID` in `ostap/fitting/variables.py` (thanks to Alexander Berezhnoy)
  

# v1.4.9.1

## New features 
  1. improve banner
  1. extend `ostap/__init__.py.in`
  1. add new test for splot `test_fitting_splot
  1. extent option `minos`, allow to specify variable name or sequence of names
  ```
  model.fitTo ( ....   , minos = 'S', ...)
  model.fitTo ( ....   , minos = ('S','B') , ... )
  ```
  1. add new test/example `test_fitting_components2.py`
  `  
## Backward incompatible changes
  
## Bug fixes:
  1. fix `truediv` for python3 in several files 


# v1.4.9.0

## New features 
  1. add new cass `P2Quantile` that interfaces P^2 algorithm from GSL for running quatile (approximate)
  1. Add methods `Ostap::Statvar::p2quantile`, `Ostap::StatVar::p2quantiles` and  `Ostap::StatVar::p2interval`. These are much faster (but approximate) versions of `Ostap::StatVar::quantile`, `Ostap::StatVar::quantiles` and `Ostap::StatVar::interval`, using P^2 algorithm. 
  1. update `ostap.stat.statvar` for modified `Ostap::StatVar` methods
  1. allow uisng sqlite3 dbase for compresed shelves
  1. extend compressed shelves to keep some metainformation on database (creating/modification date, versions of ostap, ROOT and python versions
  1. add creation/modification date for the items in compressed shelves   
  
## Backward incompatible changes
  1. change return type from `Ostap::StatVar::quantile` , `Ostap::StatVar::quantiles` and  `Ostap::StatVar::intrval`  method, adding  also number fo events used for quantile/interval estimation.  It allows to judge abotu the precision  
  
## Bug fixes:
  1. fix really stupid bug in `ValueWithError`
  1. Tiny fix in `Ostap.DataFrame.ProgressBar`
  1. fixed in_range option for the case when fit variables are defined as  RooRealVar

# v1.4.8.7

## New features
  1. implement true `read-only mode for `sqldict` (and therefore for `sqlliteshelve`
  1. tune the names of temporary files/directories for `compress_shelve`
  
## Backward incompatible changes

## Bug fixes:
  1. couple of minor fixes in `compressed_shelve` and `dbase`
   
# v1.4.8.6

## New features

  1. improve `compress_shelve` for (much) better treatment of "other' databases, in particular those with several on-disk files 
  1. fix unesessary complains/warnings on redefined varibales 
  1. allow implicit name duplicationn for cloning&copy of `FUNC`/`PDF` objects  (Is it a good idea???)
  1. update `test_fitting_models` 
 
## Backward incompatible changes

## Bug fixes:


# v1.4.8.4

## New features

  1. improvea bot the printout for `compressed_shelve`  
  1. add new module `ostap.math.covtransform` for transformation  of covariance matrices
  1. add tests for `contransfrmm`
  1. extend (part of) linear algebra for `SVectorWithError`
  1. add `io.dbase` module allowing to use `bsddb3` if/when available (for python3)
  `
## Backward incompatible changes

## Bug fixes:

  1. fix bug for `PDF`/`FUNC` evalaution with uncertainty 
  1. fix bug for `Ostap::Math::SMatrixWithError`
  1. fix travis-CI tests   


# v1.4.8.3

## New features

  1. change the inheritance diagrams for PDF/PDF2/PDF FUNC/FUNC2/FUNC4 classes: now PDF3 inherits from PDF2 and FUNC3, FUN2 innherics formm PDF and FUNC2, PDF inherits from FUNC, FUNC3 inherits from FUNC2 and ZVar, FUN2 inherits from FUNC and YVar and FUNC innherits fomr XVar
  1. add new method: derivatives and integration for FUNC/FUNC2/FUNC3 classes 
  1. add the test/example for `graph_summary`  
  1. small optimization for `linalg2'
  1. polishing for `graph_summary` - add `offset` argument
  1. move certainmethdo from PDF to FUNC: `params`, `__contains__` , `parameter`, `parameters` , `load_parameters`
  1. remove usage of `RooAbsReal::getParameter ( None )`
  1. remove usage of `None` as null-pointer
  1. adjust a bit `RooArgList.__contains__`  to use `RooCollection::find` insntead of `RooArgList::index`
  1. add the actual database type to the printut of `compressed_shelve`
  1. add transformation PDF + test/example 
  1. add `M2Q` and `Q2M` transfomoration variables/function
  
## Backward incompatible changes

## Bug fixes:  
  1. fix bugs in `funbasic`, `roofuncs` methods
  1. fix bugs in `PyVAR2`
  1. fix namings in `MakeVar.name`
  1. more fixes in `linalg2/MatrixUtils2/MAtrixUtilsT`
  1. few minor fixes as preparatory for "transform-PDF"
  
# v1.4.8.2

## New features
  1. Modify `Ostap::Models::BWI::evaluate` (temporary action,
     to be properly fixed in the future)
  1. Fix `toys.make_toys` for possible memory leak (thanks to
     Abdul-Kerim Gusseinov for repoting the problem and for solution)
  1. Add `bufstrat` argument for `Convolution` and `Connvolution_pdf`
  1. add more functions `isuint`, `isulong` , `islonglong` , `isulonglong`
  1. more operations with `TMultiGraph`
  1. add `graph_summary`
  1. update `graph_summary` to add colored bands for "averages" `
  1. replace `ROOT.Double` with `ctypes.c_double`
  1. `graph_summary` : add labels and type `Graph`
  1. `graph_summary` : rename classes, remove `TMultiGraph` and add documentation`  
  1. rewrite `ostap.math.linalg` : more functions & mixed operations: S/T-matrices/vectors&numpy
  
## Backward incompatible changes

## Bug fixes:  

  1. fix limits for `right` variable for `PSRight_pdf`(thanks to Tatiana Ovsiannikova for reporting the problem)
  1. fix `pdg_format` for certaint cases
  1. fix missing `name` attribute for `Sum1D/Sum3D` clone  machinery  
  1. fix some bugs in `graphs.py`
  1. more bug fixes in `graphs.py`
  1. fix `Flatte_pdf`
  
  
# v1.4.8.1


## New features:#include "Ostap/MatrixUtils.h"


  1. add new roofit variables  (`RooAbsReal`):  

    useful e.g, for efficiency or phase/amplitude parameterization
     - `Ostap::MoreRooFit::Benrstein`      
     - `Ostap::MoreRooFit::Monotonic` 
     - `Ostap::MoreRooFit::Convex` 
     - `Ostap::MoreRooFit::ConvexOnly` 
  1. More generic `Ostap::MoreRooFit::ShiftAndScale`
  1. make a first try to add `evaluateBatch` for existing PDFs: no large gain is observed :-(
  1. add `roofunc.py` & `funbasic.py`
  1. add `FUNC`, `FUNC2` and `FUNC3` classes. Move some functionality from `PDF`
  1. add `Fun1D`, `Fun2D`  and `Fun3D` wrappers 
  1. update `efficiency.py`
  1. reshuffle a bit the existins strustures
  1. add operations and operators for `FUNC`, `FUNC2` and `FUNC3` objects 
  1. improve baseclasess for functions and PDFs, add more operators/operations 
  1. modify `roofuncs` avoiding the dangling references 
  1. Add `Ostap::MoreRooFit::Id`
  1. make use of `Ostap::MoreRoofit::Id` in `Fun(1,2,3)D`
  1. exclude serialization in `test_fitting_models3_2D` for ROOT 6.18 & python3
  1. remove duplication betwen `ostap/fitting/roofuncs` and `ostap/fitting/variables`
  1. disable `pathos` for python_version<=3.7
  1. change pickling for all `FUNC`-based objects
  
## Backward incompatible changes

## Bug fixes:  

  1. fix `local_roofit.h`


# v1.4.8.0

## New features:

  1. add methods `vars_power`,`vars_exp` and `vars_formula` that allows to create functional variables (and properly store intermediate objects) 
```
pdf = ...
v1  = pdf.vars_power   ( a  ,  2 ) ## a^2
v2  = pdf.vars_power   ( 10 ,  b ) ## 10^b
v3  = pdf.vars_power   ( a  ,  b ) ## a^b 
v4  = pdf.vars.exp     ( a  , -1 ) ## exp(-a)
v5  = pdf.vars.exp     ( a  ,  b ) ## exp(a*b)
v6  = pdf.vars_formula ( '{}*{}/{}'       , vars = ( a , b , c ) ) 
v7  = pdf.vars_formula ( 'x[0]*x[1]/x[2]' , vars = ( a , b , c ) ) 
``` 
  1. add `Asymptotic` (also `AsymptoticErr`, `AsymptoticError`, `AsymptoticErrors`, case-insensitive)
     keyword for fit-related methods, that (for ROOT versions starting from 6.19)
     is decoded to `ROOT.RooFit.AsymptoticError(...)`
  1. add possible file-extensions for   canvas-plots :
```
canvas >> 'plots.tgz' ## make plot in all formats and store them in tar/GZIP  format
canvas >> 'plots.tbz' ## make plot in all formats and store them in tar/BZIP2 format
canvas >> 'plots.txz' ## make plot in all formats and store them in tar/LZMA  format
``` 
  1. Improve a bit the interface for BLUE: Best Linear Unbiased Estimator : 
     combination of correlated measurements 
  1. add test for BLUE  `ostap/stats/tests/test_stats_blue.py`
  1. add SciPy/FFT-based convolution for functions `ostap/math/sp_convolution.py'
  1. add SciPy-based bspline interpolation for functions `ostap/math/sp_interpolation.py'
  1. imporve pseudo-abstract operations `ostap/math/operations.py'
  1. extend interpolation tests
  1. add simple shapes to probe signal/background interference `Ostap::Models::BWI` and `BWI_pdf` 
  1. Fixes for ROOT 6.20/00

## Backward incompatible changes

  1. Remove all `Rotated`-stuff

## Bug fixes:  

  1. fix typo in treatment of `pdf.draw ( ,,, , in_range(1/2/3) = ... , )` 
  1. fix some new typos/errors in `ostap.math`

# v1.4.7.1

## New features:

## Backward incompatible changes

## Bug fixes:  

  1. fix `Ostap::Math::Positive::updateBernstein`. The bug was introduced in 1.4.7.0. Thanks to Tatiana Ovsiannikova for  reporting the problem.

# v1.4.7.0

## New features:

  1. Slight improvements in `ostap.fitting.minuit`
  1. add new test for `ROOT.TMiniut` decorations :  `test_fitting_minuit`
  1. allow `pathos` to be used for paralellization for python version > 3.6  
  1. add test `test_parallel_dill` to check the the problem with combniation 
     of different versions of `ROOT` , `dill` and `python`.    
  1. improve `ostap.fitting.toys`
  1. improve `ostap.parallel.parallel_toys`
  1. add `test_fitting_toys`
  1. add `test_parallel_toys`
  1. add examples/tests for evaluation of significance with toys
  1. improve the output of `ls`-method for for `compressed_shelve` 
  1. few minor tweaks needed for `picalib` 
  1. add `DisplayTree` into `ostap.logger.utils` - useful tool to render the tree-like structures
  1. impove `ls` methods for all shelve-like databases 
  1. add `ls_tree` and `ls_table` method for `ROOT.TDirectory`
```
f = ROOT.TFile ( .. )
f.ls_tree ()
f.ls_table() 
```
  1. speedup construction of Bernstein polynomials from the list of roots 
  1. re-write `PDF.wilks` method to use `ROOT.RooProfileLL`
  1. add methods  `PDF.graph_nll` and `PDF.graph_profile` for 
     a bit more easy and more fast drawing of NLL-scans and profiles 
```
pdf = ...
g1 = pdf.graph_nll     ( 'S' , vrange ( 0 , 20.0 , 100 ) , dataset )
g2 = pdf.graph_profile ( 'S' , vrange ( 0 , 20.0 , 100 ) , dataset , fix = ['gamma','mu'] )
```
  1. allow to suppress certain RooFit message topics from the configuration file, e.g.
```
[RooFit]
RemoveTopics = Plotting            ,
               Caching             ,
               Eval                , 
               Minization          ,
               Integration         ,
               Optimization        ,
               NumericIntegration  , 
               Fitting     
```
  1. Add more flexibility for parallel `Task` : one can specify the additional environment variables, prepend and append path-like environment variables, execution directory, ... variables are expanded.
```
mytask = MyTask  (... ) 
mytask.environment [ 'LD_LIBRARY_PATH' ] = 'some_directory1:some_directory2' 
mytask.prepend_to  [ 'PATH'            ] = 'some_directory1:some_directory2:some_directory3'
mytask.append_to   [ 'PYTHONPATH'      ] = '$HOME/python'  ## will be expanded at remote host
mytask.directory                         = 'some_existing_remote_directory'
mytask.dot_in_path                       = '.' in sys.path 
```
  1. add new  function `random_random` into `ostap/parallel/utils.py`, that sets (hopefully) proper seeds for the `random`, `ROOT.gRandom` and `ROOT.RooRandom` 
  1. add invocation of `random_random` from `initialize_remote` for the parallel toys.
  1. add methods `graph_nll` and `graph_profile` for class `SimFit`
```
pdf = ...
g1 = pdf.graph_nll     ( 'S' , vrange ( 0 , 20.0 , 100 ) , dataset )
g2 = pdf.graph_profile ( 'S' , vrange ( 0 , 20.0 , 100 ) , dataset , fix = ['gamma','mu'] )
```
 
## Backward incompatible changes

## Bug fixes:  

  1. Tiny fix in `ROOT.TMinuit.cor`
  1. Tiny fix in `ROOT.TMinuit.cov`
  1. more fixes in `ostap/fitting/minuit.py`
  1. small fixes in `ostap/fitting/toys.py`
  1. tiny fix in `ostap/math/random_ext.py`
  1. Fix signature of `ds_project` method from `ostap.fitting.dataset.py`
  1. few minor fixes needed for `picalib` 
  1. fix `Ostap::HistoProject::project` for weighted `RooAbsData`, 
     now uncertainties are evaluated correctly, properly accounting the errors in weights  
  1. fix typo in `paralllel_toys2`
  1. fix `random_random` for python3
 

# v1.4.6.0

## New features:

 1. Add new method `dct_params` for `ROOT.RooFitResult`, that gets a dictionary of all parameter values 
 1. `plotting/canvas.py`: Add option to save the canvas congent in tar/tgz/zip-archives:
   ```
   canvas >> 'test.zip'
   canvas >> 'test.tar'
   canvas >> 'test.tgz'
   ```
 1. `ostap/fitting/toys.py`  new module with simple functions to perform fitting toys
```
pdf = ...
results , stats = make_toys ( pdf      ,           ## PDF  to use 
                  1000                 ,           ## number of toys 
                  [ 'mass' ]           ,           ## varibales in dataset 
                  { 'nEvents' : 5000 } ,           ## configuration of pdf.generate
                  { 'ncpus'   : 2    } ,           ## configuration of pdf.fitTo
                  { 'mean' : 0.0 , 'sigma' : 1.0 } ## parameters to use for generation 
                 )
``` 
  1. `ostap/parallel/parallel_toys.py`  a version of the function above for the 
      parallel execution via multiprocessing and/or parallel python
  1. add method `evaluate` for `ROOT.RooFitResult` to evaluate the arbintrary function of 
     fit-parameters with uncertaionties
```
res = ... ## ROOT.RooFitResult object
fun = lambda x,y,z :  x*x+y*y+z*z 
val = res.evaluate ( fun , ('x' , 'y' , 'z') ) 
```
  1. add `EvalNVEcov` and `EvalNVEcor` evaluators into `ostap.math.derivative` module. 
     These objects evaluate the value of N-argument function taking into account 
     the uncertainties and/or correlations 
  1. add method  `max_cor` to `ROOT.RooFitResult` to obtain the maximal correlation coefficient
```
res = ...  ## ROOT.RooFitResult
coefficient , variable = res.max_cor ( 'X' ) ## get the maximal correlation
```
  1. add the column `Global/max correlation` for the table form of `ROOT.RooFitResult`
  1. add `ncpus` argument to `ostap.parallel.parallel_toys.parallel_toys`

## Backward incompatible changes

## Bug fixes:  

  1. Minor fix in `Files.hadd` 
  1. Couple of minor fixes in 
      - `ROOT.RooFitResult.sum`
      - `ROOT.RooFitResult.multiply`
  1. Couple of minor fixes in `ostap.math.derivative`
  
# v1.4.5.0

## New features:

  1. Add new functions:
       - `parameters`
       - `parameter` ( and shortcut as `__getitem__` )
     for the base class `PDF` to allow "easy" access to the values of 
     parameters by name, e.g.  `a = pdf['A']` or `a = pdf.parameter('A')`
  1. New function `mid_point` is added to `PDF`-base class, 
     defined as 0.5*(x_low + x_high), where f(x_low)=f(x_high)=0.5 * f_max.
     x_low and x_high  are the same points used for evaluation of FWHM.
     It characterises the location of the peak, similar to `mode`, `median`,
     `get_mean` and other related functions 
  1. Add `ResoBukin` -  symmetric resolution function based on Bukin-pdf. 
  1. Update `tests_fitting_resolutions.py` 
  1. Add `Losev`-function&pdf - a kind of asymmetric hyperbolic secant function
     - `Ostap::Math::Losev`
     - `Ostap::Models::Losev`
     - `ostap.fitting.signals.Losev_pdf`
  1. tune `RooFitResult.table` printout method
  1. update some fitting examples  
  1. add option/property  `directory` for `ostap.parallel.task.Task` to indicate 
      the directory where the job needs to be executed       
  1. add option/property  `environment` for `ostap.parallel.task.Task` to 
      setup   additional environmental variables (if needed)  
  1. Add  argument `keys`  for the method `clone` for all databases, allowng to 
     copy only certain keys into cloned database. (Default: copy all keys) 
  1. add `json`, `gif` and `jpg` output formats for  default `canvas >> 'aaa'` operator
  1. add `hadd` method for `ostap.trees.data_utils.Files` to merge ROOT file via `hadd` script


## Backward incompatible changes

  1. MASSIVE RENAME/FIX:  Apolonios -> Apollonios
  1. Add `jobid` argument for `Task.process`
      - All existing tasks are updated properly 
      - All new functions and tasks must take this argument into account!


## Bug fixes:  

  1. Fix minor bug in `ResoBukin`-resolution function


# v1.4.4.2

## New  features: 

  1. `table` function for `RooFitResult`: print object in a form of nice table 
  1. `ostap.logger.utils`  set of minor functions for nice printout 
      - `pretty_float` : print float 
      - `pretty_ve`    : print float and error 
      - `pretty_2ve`   : print float with asymmetric errors 
  1. `ostap.logger.table`  : add function `align_column` to aling the column for a given table
 
## Bug fixes: 


# v1.4.4.1

## Bug fixes 
 
  1. Fix for "old" ROOT version in `MoreRooFit.cpp`


# v1.4.4.0 

## New features 

 1. Redesign all *shelve-like data bases:
    - add abstract base class
    - implement concrete databases :
      - zipshelve (ZIP/GZIP compression )
      - bz2shelve (BZIP2 compression) in terms of base class
      - lzshelve  (LZMA-compression) only for python3
 1. impove timing functions
 1. improve `TTree.the_variables`
 1. selector via frames: change order of variables in  snapshot
 1. improve statistics for selectors
 1. improve a bit printout for TTree/TChain (sorted, type,...)
 1. improve a bit printout for DataFrame (sorted)
 1. generate temporary column names for DataFrame using hash instead of random 
      - it allows better debugging (reproducible)
 1. re-enable security key for paralell python servers
 1. Update all tree-collection utilities
      - Files
      - Data
      - Data2
      - DataAndLuimi
 1. Add colorization to progress bars
 1. fix (hope) segfaults for adding new branch to tuple
 1. Improve colorization
 1. make optional use of terminaltables package
 1. add local function to format tables
 1. improve printout of trees and datasets
 1. more improvements in tables&colorization
 1. improve reweighting machinery
 1. add 2D and 3D  moments for the histograms
 1. improve 2D histos comparison functions
 1. improve progress bar
 1. fix addition of new brnaches to `TTree/TChain` and whole `IFuncTree` machinery
 1. change the names for example-tests
 1. add script to check the dependensies
 1. one more attempt to fix the crash for `FuncTH1`
     - `add_branches` - Notifier is not invoked... Invoke it explicitely!
     - `FuncTH1::Notify` : reset/delete formula  instead of Notify...
     - fixes in `Notifier`
     - few more fixes
     - add a test for `add_branches`
 1. make a try to polish a bit the Doxygen/Sphinx machinery
 1. Update `test_tools_reweight2.py`
 1. Add  `Ostap::Math::ChebyshevApproximation`
 1. extend `PDF.nll` to accept all keywords
 1. unify `PDF.nll` keywords with `pdf.fitTo` keywords
 1. tmva & tmva/chopping: fix problem with non-deleted temporary files/directories
 1. `cleanup` :  add concept of "local" trash, to be deleted when `CleanUp` instance is deleted
 1. `mp_pathos` & `parallel/task`: improve the output of the final statistics, add total time and CPU gain due to paralellization
 1. `add_branch`/`add_new_branch`
       - extend current functionlaity  allowing to add several branches at once
 1. Add "clone" methods for all io-databases
 1. Add `dump_root` module/utility/recipe to allow reading of databases 
    created with "old" ROOT versions, with old version of streamers
 1. add/extend math-function for graphs
 1. remove drawing artifacts for `PDF.draw_nll` method
 1. Add `ROOT.TGraph.remove` method
 1. Improve `ROOT.TGraph.filter` method
 1. Add plotting options for `combined_background`, `combined_signal` and `combined_components`
 1. Extend the signature for `Ostap::Math::gauss_pdf`/`gauss_cdf`
 1. Add `PSSmear_pdf` - smeared version for `PhaseSpace`-based PDFs
 1. Add `signals` & `backgrounds` keywords for `Fit1D`-constructor
 1. Add 'args' argument to draw-functions
 1. Add possibility for prefilter of data for `TMVA`/`chopping`.
 1. Allow TMVA/chopping tools to process `RooFit` datasets: (converted internally to `TTree`)
 1. optimize evaluation of polynomials in `Ostap::Math`
 1. Add `Ostap::Math::Clenshaw::term` - evaluate the N-th term of the recursive sequence
 1. Add `Ostap::Math::barrier_factor` -   evaluation of Blatt-Weisskopf angular momentum
       centrifugal form factors for  arbitrary angular  momenta
 1. Fix but in `VE.purity` : Thanks to Alexey Dziuba
 1. Improve configuration of canvas&styles
 1. Add `DalitzIntegrator` for  relatively efficient integration over Dalitz plot
 1. Add true analytical 3-body phase space
 1. Add check for duplicated variable/pdfs names
 1. Add `ProgressBar` action for `DataFrame`, now one can display the progress bar
       during the processing of large frames.
 1. Add context manager to remove/add certain topics to `RooMsgService`

## Bug fixes 

 1. `PDF.draw_nll` : fix bug for the weighted datasets
 1. Few  fixes in `Dalitz`, epsecially for  drawing it
 1. bug fix in `add_branch`
 1. minor fix in `TGraphAsymmErrors.transform`
