## New features

  1. a bit more correct treatment of the `!` and `!=` in `ROOT.TCut.ast` 
  1. more coherent treatment of the numeric types in `ROOT.TCut`-expressions
  1. some limited treatment of infinity/NaNs in `ROOT.TCut`-expressions
  1. better (more safe) implementaton of `APDF1.check_ranges` method, inclusing better documentation
  1. massive improvement for Apollonius & friends.
  1. add `clamp` method for histograms
  1. `reweight.py` : add check for non-positive data 
  1. extend functionality for checking the uniquness of ROOT/RooFit names (it will reduce amount of code)
  1. allow specificationof explicit suffix for various `rootID/hID/fID/dsID/...` methods 
  1. more improvement for reweighting machinery 
  1. histograms: move calculation of Riemann' sums from python to C++
  1. more harmolization to power-law tails
  1. more haromonization for n/N for Student's t-functions 
  1. more harmonization for some polynomial constructions 
  1. Harmonize polynomizll constructions (both in C++ and python)
  1. Release requirement of the same range for *some* operations with polynomials
  1. Add `HRange/hrange` for histogram-like looping
  
## Backward  incompatible 

   1. From now `Apollonios` stands for the core asymmetric Apollonios function and `ApolloniosL` represents Aplooniosu functon with power-lay left tail
   
   
## Bug fixes

   1. fix recentlyintroduce bug in `Ostap::Math::BifurcatedGaussian`

# v3.2.1.0 

## New features

  1. few more tweaks for CrystalBall-like functions to make them a bit more generic  
  
## Backward  incompatible 

## Bug fixes

   1. tiny fix in `Ostap::Math::(Gen)Hyperbolic::setStandard`

# v3.2.0.0 

## New features

  1. Fix ACLiC warning (mainly for unused arguments)
  1. Add three primitive ACLiC tests
  1. Remove local copy of `pragmastat` (rely on standard `pragmastat`)
  1. Modify `cmake` to istall standard `pragmastat` from `pypi` if missing
  1. Add code to select dataset entries that are not shared with another dataset  
  
## Backward-incompatible

## Bug fixes

# v3.1.0.0 

## New features

  1. more improvements for helper script `./aux/.build-lcg`
  1. add `ast` method to `ROOT.TCut` to parse the expression using `ast` and then unparse it... 
  1. add `pow` and `rpow` functions for TCut-expressions 
  1. add `mod/imod/rmod` operators for for TCut-expressions 
  1. add `abs` operator for for TCut-expressions 
  1. add proper treatment of `Tprofile` and ``TProfile2D  for `project` methods
  1. add proper treatment of uncertaintues of projections of weighted `RooAbsData` in case the weight has assigned errors.
  1. add proper treatment for `ROOT.TProfile3D` where possible
  1. switch-off parallel projection for all profiles
  1. switch-off frame-base processing for `ROOT.TProfile3D`
  1. extend the `h1_stack` treatment for data projection methods 
  1. update math related to Crystal Ball functions and add new test 
  1. add explicit `__contains__` function to `RooLinkedList`
  1. add explicit `__getitem__`  function to `RooLinkedList`
  1. add `non_gaussian` method for all "peak-like"-distribution with well-defined "mean" and "rms"
  1. improve CPACK-configuration 
  1. more tweaks for CPACK RPM (&DEB) - not properly tested yet! 
  1. modify shebangs for several files to please CPACK RMP generator
  1. add shebang to `*.C` root macros (trick from [https://gist.github.com/gipert/3f91d6dc8bf31818cff163faeb18f38e](here)
  1. add `Ostap::Utils::same_binning` for `TAxis` and `TH1` classes  (`same_bining` in python)  `
  1. Issue warning message for histogram's _binomial efficiency method_ in case of different bining schemes    
  1. more tweaks with CPack packaging
  1. use bare testing for `parallel`-tests instead of `nosetests`

## Backward-incompatible

## Bug fixes

  1. fix `ostap-sync-dirs` script 
  1. fix some newly introduced bugs
  1. fix the bug in `Ostap::Math::CrystalBall::integral`
  1. fix usage f nosetest for tests (when available)
  
# v3.0.2.16 

## New features

## Backward-incompatible

## Bug fixes

  1. Fix typo in `Batch`

# v3.0.2.14

## New features

   1. (Re)add multiple projections to 1D histogram for `DataFrame`
   1. Add `errors` optional argument for `TH1.vmin/vmax/vminmax` methods 
   1. Add `control_plots_signal` and `control_plots_background` argument for TMVA&Chopping `Trainer` classes
   1. Improve treatment of ROC-curves in TMVA. Now multigpraph is saved into ROOT output file
   1. Improve treatment of control plots in TMVA. Now control histos are saved into ROOT output file
   1. Improve `Canvas` context manager 
   1. Switch the ownership of Ostap-canvas to python
   1. Add (C++) progress bar to `Ostap::Selector` class
   1  Remove python progress bar from all hierarchy of selectors
 
## Backward-incompatible

## Bug fixes

   1. fix a sad typo in `TH1.vmin/vmax/vminmax` methods 

# v3.0.2.12

## New features

   1. (Re)add multiple projections to 1D histogram for `RooDataSet` (Many thanks to Dima Golubkov and Artem Egorychev)

## Backward-incompatible

## Bug fixes

# v3.0.2.10

## New features

   1. add `latex=False` argument for all `pretty_print` functions
   1. add `nice_print`  function to get ready-to-use pretty representation 
   
## Backward-incompatible

## Bug fixes

   1. fix stupid bug in `ostap.math.base`
   
# v3.0.2.8

## New features

  1. Add `backup_to_ROOT` and `restore_from_ROOT` functions for `ostap.toool.reweight`. 
  1. Tiny tweaks for `ROOT.TDirectory` decorators 

## Backward-incompatible

## Bug fixes

# v3.0.2.6

## New features

  1. Fix `KeepCanvas` and `Canvas` context managers to correctly switch to the initial stage.
  1. reshuffle code between `ostap.plotting.canvas.py` and `ostap.utils.root_utils.py`
  
## Backward-incompatible

## Bug fixes

# v3.0.2.4

## New features

  1. Add script `ostap-sync-dirs` for (parallel) synchronization of directories
  1. More tweaks for `sync_files` & `copy_files`
  1. A tiny improvement for the TH2 and TH3 summary printout
  1. Add `Ostap::Math::background` function and
  1. Add `background` method for `Ostap.Math.ValueWithError`
  
## Backward-incompatible

## Bug fixes

# v3.0.2.2

## New features

  1. add environment variable `OSTAP_NCPUS` to limit the maximal number of CPUs
  1. A few minor tweaks here and there
  1. The first step towards the `tkrzw` database
  1. TMVA&Chopping: prepend the names of plot-files with the trainer name
  1. Minor fix for `ostap.trees.cuts.py`
  1. Add `tmva.plot_variables` function
  1. Add more tweak argument to `sync_fiels`
  
## Backward-incompatible

## Bug fixes

# v3.0.2.0

## New features

  1. Add treatment of TMVA&Chopping `spectators` into parallel TMVA&Chopping `addResponse` methods 
  1. Tweak a bit `Ostap::Math::beta` and `Ostap::Math::lnbeta`
  1. Add an overload for `Ostap::Math::beta` for the (positive) integer arguments
  1. Add `beta_pdf`, `beta_cdf` & `beta_quantile`
  1. Add `bayes_interval` for eauation of binomial intervals
  1. Add `histogram.eff_bayes` to get the "efficiency" graph 
  1. Add `Ostap::Math::Beta` , `Ostap::Models::Beta` and `ostap.fitting.distributions.Beta_pdf`
  1. Add `Ostap::Math::punzi` function to obtain estiamate for Punzi's figure-of-merit
  
## Backward-incompatible
    
  1. From now all arguments for (TMVA&Chopping) `Reader` are keyword-only 
    
## Bug fixes
    
  1. `Spectators` were not propagated to TMVA/Chopping

# v3.0.1.26

## New features

   1. More improvements in TMVA & Chopping printout
   1. For TMVA & Chppoing add new arguments `signal_add_vars` and `backrgound_add_vars` to add more variables to transformed data sets
       to release the intrinsic requirments fot the same structure for signal ad background samples.
       In this way one can `on-fly` add more  missing variables to please TMVA, see example `test_tools_chopping2.py`
   1. Modernize `AFUN1.load_params` : issue  warnings for the parameter setting outside the allowed range.
     
## Backward-incompatible

## Bug fixes 

   1. fix a bug in `histo.eff` for non-uniform binning  (thanks to Evgenia Nekrasova)
                  
# v3.0.1.24

## New features

   1. Extend `pdf_convolution` to accept the histogram as the first argument using the historgam->PDF transformation
   1. Add `convolute` method for histogram (internally using `pdf_convolution`)
          
## Backward-incompatible

## Bug fixes 

# v3.0.1.22

## New features

   1. Add Meixner distribution+PDF
   1. Add (symmetric) Meixner resolution function 
            
## Backward-incompatible

## Bug fixes 
  
  1. fix a stupid bug in `Ostap::Math::lgamma` for complex arguments 

# v3.0.1.20

## New features

   1. More polishing around `chunk_size` for parallel processing 
            
## Backward-incompatible

## Bug fixes

# v3.0.1.18

## Bug fixes 

   1. couple of minor fixes in parallel processing 

# v3.0.1.16

## New features
    
  1. improve printout for TMVA & Chopping stuff
  1. suppress unesessary ROOT printout from TMVA & chopping
            
## Backward-incompatible

  1. Consistently use `chopper` as the name for argument insteads of `chopping` for all `chopping`-related stuff

## Bug fixes

# v3.0.1.14

## New features

## Backward-incompatible

## Bug fixes

  1. many small fixes ...

# v3.0.1.12

## New features

   1. TMVA&Choppnig: allow to specify several signal and backgroud samples  

## Backward-incompatible 

## Bug fixes

   1. fix a tiny typo in `tree_reduce`
    
# v3.0.1.10

## New features

 1. Improve a little bit the prinntout of `Ostap::Math::Moment_` objects in python
 1. Add Cornish-Fisher estimate for approximate quantile functions
 1. Tiny fix for `Ostap::Math::probit`
 1. add `Ostap::Math::Moments::quantile_` functions for `Ostap::Math::(W)Moment_` objects, based on Cornish-fishef asymptotoc expansion

## Backward-incompatible 

## Bug fixes

# v3.0.1.8 

## New features

    1. Improve print of the summary AUC/ROC tables for tmva&chopping
    1. Add new example `test_tools_tmva3.py`     to illustrate the usage of signal/background samples with different variables  
    1. Add new example `test_tools_chopping2.py` to illustrate the usage of signal/background samples with different variables  
     
# v3.0.1.6 
# v3.0.1.4 

## New features

   1. Slightly improve prints from for `Files/Data/...` - make them a bit more informative
   1. Add `truncate_middle` function into `ostap.utils.strings` module: A little bit modified version of the code from the `pidgen2` project by Anton Pluektov
   1. Add `ResoNeedham` as explicit resolution function
     
## Backward-incompatible 

   1. Rename `a0`, `a1` and `a2` parameters of Needham's function to `c0`, `c1` and `c2`. Update `Ostap::Math::Needham`, `Ostap::Models::Needham` and `Needham_pdf` `
    
## Bug fixes

# v3.0.1.2

## New features
    
  1. imporve/extend the `Data.table` method - now it the add information about the tree/chains
     
## Backward-incompatible 

## Bug fixes

  1. fix `Files.table`
    
# v3.0.1.0

## New features

      1. add `TTree.good_variables` method to check if variables are in TTree or could be computed
      1. add `data_efficiency/tree_efficiency/frame_efficiency' methods  to get the historgams of efficiencies
      1. more tweaks and fixes for `*_efficiency` functions
      1. drasticaly speed-up creation ,splitting and pickling/unpickling of `Tree/Chain` objects
      1. Update `AddTMVAResponse` : merge `OStap::TMVA::addResponce` snd `Ostap::TMVA::addChoppingResponse` functions into `Ostap::AddTMVA` class and modernize the ptython layer
         
## Backward-incompatible 

      1. From now `ostap.utils.cleanup.TempFile.__enter__` return the actual name of the temporary name (intead of `self`) 
      1. From now `files` and `nFiles` are properties (not metoids!) for `TTree/TChain` clases
    
## Bug fixes

# v3.0.0.4

## New features

  1. allow add the dictionary `{ name : numpy.ndarray }` to `ROOT.TTree
  1. allow to add the structured numpy arrays            to `ROOT.TTree
  1. add `TTree.good_variables` method to check if variables are in TTree or could be computed
    
## Backward-incompatible 

  1.  change the order of `buffer` and `name` arguments for `add_new_buffer`, `buffer_2tree` and `buffer_2chain` functions
    
## Bug fixes

    
# v3.0.0.2

## New features
    
 1. add `-v/--version` flag for `ostap`  : print the version and exit
 1. refactor the main ostap parser into separate module  (it speeds up the `ostap -v` action)  
 1. improve  a bit the `add_var/add_branch/add_buffer` machinery 

## Backward-incompatible 

## Bug fixes

# v3.0.0.0

## New features

  1. completly re-write `Ostap::HstoProject` (and rename it into  `Ostap::Project`
  1. remove a lot of duplicated code between `Ostap::StatVar`, `Ostap::Project` and parameterization
  1. remove some legacy artifacts
  1. more progress for quantiles 
  1. maje numpy mandatory and start to eliminate some local `ifs` 
  1. make scipy mandatory and start to eliminate some local replacement 
  1. add `pragmastat.py` (copy the code from Andrei Akinshin with tiny modidication)
  1. extend the `data_*` functions from `ostap/stats/startvars.py` to coherently deal with `TTree` and `RooAbsData`, switchig if/when needed inntospecific TTree/RooAbdData functions
  1. coherent treatment of `use_frame` and `parallel`optional arguments for `data_*` functions
  1. simplify, extend  and make  it more universal `parallel_statvars` module
  

## Backward-incompatible 

## Bug fixes

# v2.1.0.0

## New features

  1. small step toward kernel estimators 
  1. add an error estimate for standartized moments if/when possible
  1. improve printout of tables of moments
  1. more improvement and simplifications to the moments, incorporate xmin/xmax/wmin/wmax 
  1. add (partly suplicated>  `bin_edges` ietrators/methods  for `ROOT.TAxis`, `ROOT.TH(1,2,3)(F,D)``
  1. Add quantiles for `Ostap::Math::(W)ECDF`
  1. add `ostap/toos/evolution.py` - utilty to investgate the on-parametrick evolution of the
     distribushae as function of another parameter
  1. Add `n`-parameter to Needham' function/PDF 
  1. remove a lot of ancient version switches

## Backward incompatible 

  1. change constructors in `Ostap::Math::(W)Moments_<N>`
  1. change pickling/unpicking/serializatino/deserialization for (w)moments 
  1. completely remove `Ostap::Math::(W)MinMaxValue`
  1. disable the support of ROOT versions <6.26 - both in CMAKE and in C++/python code 

## Bug fixes 

# v2.0.6.0

## New features

  1. more platforms to test via GitHub Actions (added gcc15) 
  1. more coherent treatment of status-codes in `HistoProject` 
  1. add few more `MoreRooFit` functions
  1. add `AsymVars` helper structure in `fthelpers.py`
  1. add `table` mehtod for `AFUN1` and bundle it with `__str__` and `__repr__`.
     Now for all derived objects a detailed table it printed
  1. exclude `std::array` from `Ostap.xml`
  1. redesign and rewrite `AsymVars structure
  1. add `ResoCB2_` resolution function - reparameterisatuon of `ResoCB2`. Better name is needed here...
  1. rename `ResoCB2_` into `ResoCB2a`
  1. Add machinery to get the uncertainties for the PDF features due to fit uncertanties
  1. Add machinery to get the uncertainties for the PDF features due to fit uncertanties in a parallel way
  1. make FUN/PDF __call__ method aware of numpy (conversion to scalar - Deprecation Warning 
  1. few typos & formatting
  1. add conversion of (S)Matriuces and (S) vectors to histogrms (&visualization)
  
## Backward incompatible

## Bug fixes

# v2.0.5.2

## New features

  1. improvement in convolutions
  2. automatically convert the convolution of sum into the sum of convolutions
  
## Backward incompatible

## Bug fixes

# v2.0.5.1

## New features

## Backward incompatible

## Bug fixes

  1. fix typo in `tmva.py`
  
# v2.0.5.0

## New features
    
  1. Add more protections for missing scipy
  1. polish `CMakeList.txt` files
  
## Backward incompatible

## Bug fixes

# v2.0.4.0

## New features
    
  1.  Add `defined(__cpp_lib_span)` check into `Ostap/AddBugger/h`
  1.  Several tweaks aroubns `ActiveBRanches` contetx manager to activate/deactivate TTree branches
  1.  Explciitely add `ActiveBranches` for many functions that operate with `TTree/TChain`c
  1.  Fix for misterious failures in several reweighting tests 

## Backward incompatible

## Bug fixes

# v2.0.3.0

## New features
    
## Backward incompatible

## Bug fixes

   - fix recently intoruced bug in `Ostap::Math::CrystalBallDoubleSided::setAlphaR`
    `  
# v2.0.2.0

Reshulffe code for the base utilities 

## New features
    
## Backward incompatible

## Bug fixes
    
    
# v2.0.1.0

## New features
    
  - add `footprint.py` : now somw meta-info abotu ostap sessions is saved in `$OSTAPDIR/.footprints` (if exists and writeable)
    and `$OSTAP_CACHE_DIR/.footprints` files 

# v2.0.0.0 

 - stop support ROOT < 6.24
 - stop support python < 3.8 
 - split `add_new_branch` into `add_newbranch` and `add_new_buffer`
 - fix Py <-> c++ functions
 - improve `Ostap::StatusCode`
 - more use of `Ostap::Assert`
 - remove all appearances of `from bultins import rane`
 - remove all appearanced of `from __future__ import ptint_function`
 - suppress many version-related `if`s
 - large imporovements in LinAlg
 - re-write AddBranch/AddBuffer machinery 
 - add (P)LU decomposition for matrices (+test)
 - fix several bugs&typos
 - add (P)QR Decompositoon for matrices (+tests)
 - add LQ,QL,COD,SVD and POLAR decomposition of matrices (+tests)
 - add SCHUR decomposition for squared matrices 
 - more improvements for Files/RootFiles/Data/Lumi/DataAndLumi
 - add COWs (both in C++ and python) and split python'  `sPlot` into `COWs`, `sPLOT` & `hPlot*`
 - modify the `alpha(sigma)` dependency for Needham's functions
 - Add Bernulli numbers and Bernulli polynoials
 - Add Clausen & Gleisher-Clausen functions 
 - `Ostap::Math::Moments` : turn the class with static functions onito namescpace  
 - add non-templated constructor to `Ostap::Models::Shape`d/2D/3D`
 - improve matrix tables
 - add more operations to `ROOT.TCut` objects
 - reshuffle the code around `pretty_ve` to reduce explciti apprearenc of this&related  functions
 - add possibility to disable automatic construction of SB&BS components for `Fit2D` and `Fit2DSym` models 
 - Largely re-write and improve `Gof1D` and `GoF1DToys` machinery for Goodness-of-1D-fits
 - Speed-up Goodnes-of-fit estimators (rely on C++ when it has sense)
 - Allow to provide the external CDF for GoF1D-toys (it can drastically speedup the pseudoexperiments)
 - update `makeWeighted` for 6.36<=ROOT
 - (re)fix (re)opened issue with serializations of enums @see https://github.com/root-project/root/issues/15104
 - add `SineSum`, fix nicorrectenens for `FourierSum` & `CosineSum`
 - remove `fejer` arguments
 - disable `ClassImp` for `ROOT>=-6.36/0`
     
## New features
    
   1. further improvemets for the wrapped columns in the tables
   1. add unbined splot-related stuff, includinn machinery for
      adding the splot results to TTree/TChain
   1. OSTAP_BATCH environemtn variable now is effectivey used for all scrpts via
      `ostap.core.core.py` (should it be  `otap.core.__init__.py` ?)
   1. `-b` or `--batch` command line arguments now are effectivey used for all scrpts via
      `ostap.core.core.py` (should it be  `otap.core.__init__.py` ?)
   1. switch from `N=|n|+1`  to `N = sqrt ( 1 +n*n)` for CrystalBall & friends 

    
## Backward incompatible

   1. Reorder arguments for `Data` &`DataandLumi`, from now the first argumenet - files/pattersn and then chains
   1. Remove `Data2`
   1. `add_new_branch` and `add_new_buffer` are completely rewritten
   1. `sPlot` is completely rewritten: split into `COWs` , `sPLOT` and `hPlot`
   1. The public interface for `GoF1DToys` is mofified 
      
## Bug fixes
    
   1. fix the bug in treatment of `bdsdb3`
   1. Fit recent typo in `backrgounds.py`
   1. fix recent typo in `MoreRooFit.cpp`
            
 # v1.14.0.0

    ## New features

    1. add several functions, like `sec` , `csc` , `cot` , `cas`
    1. add Bring radical
    1. add Pade interpolation     
    1. add conversion from Taylor(polynomial) expansion to Pade/rational approximation
    1. add elliptic fuctions and elliptic integrals
    1. improve treatment of elliptic integrals and Carlson symmetric forms
    1. attempt to fix/bypass the newly appeared problem:
       ```TypeError: void TCanvas::Update() =>
                TypeError: callable was deleted
       ```
    1. slightly improve `Ostap::MoreRoofit::Rank`
    1. add `ostap::MoreRooFit::AbsAplusB`
    1. add 'N=|n|+1`-accessors/attributed for Crystal Ball family of functions/models
    1. add functions to identify and remove empty columns from the table and update existing code to use these functions
    1. add better printout for SimFit and RooSimultaneous
    1. better treatment of long titles for tables
    1. add `silent=True` argument for `generate`-methods
    1. add printout of gen-options for `generate`-methods
    1. simplify (&extend) `Ostap::MoreRooFit::Addition/Addition2/Product`
    1. remove default constructor for `Ostap::MoreRooFit::ProfileLL`
    1. add some manipulations with `RooSTLRefCountList<RooAbsArg>`
    1. make use of `typename` function 
    1. slight improvements in the printout methods for matrices and vectors
    1. improve printout of `Simfit`
    1. reshuffle some basic fitting classes
    1. essentialy re-write of `RooStats`-wrapper 
    1. make use of `shutil.get_terminal_size` 
    1. improve utilities for table creatio 
    1. add better printout for RooWorkspace and ModelConfig
    1. re-write `ostap.tools.splot` extending it to 2D and 3D cases    
    
## Backward incompatible1

    1. remove obsolete `ostap.fitting.simfit.Sim1D`
    
## Bug fixes
    
    1. fix a typo  in `parallel_statvars` (@thanks to Dima Pereima) 

# v1.13.9.0

## New features
    
  1. Add Logarithmic integral functions 
  1. Add Exponental  integral functions 
  1. Add Sine and Cosine integral functions 
  1. Add missing `Ostap::Math::LegendreSum[X]::integral` methods 

## Backward incompatible
    
## Bug fixes
    
 1. fix a typo in `graphs.py` (thanks to Dima Pereima!)
    
# v1.13.8.0

## New features

  1. add parallelisation to `GoF1DToys`
  1. add `style=None` or `style=''` argument for many `table`-methods
  1. add empirical cumulative functions for weighted data 
  1. suppress `CloneData` argument for NLL creation for 6.28<=ROOT
  1. change default `silent=True` to `silent=False` for `graph_profile`
  1. modify a bit treatment of `residual=` and `pull=` in `PDF.draw`
  1. a fix for `test_fitting_datasets` and `tets_fitting_selectors`
  1. add `Ostap.MoreRooFit.Rank`  
  1. make (de)serialization of GoF1D objects more robust 
  1. tweak for `RooRealVar.__reduce__`
  1. add clipping for `GoF1D/GoF1DToys`
            
## Backward incompatible
    
## Bug fixes

  1. tiny fix for `constraint`/`constraints` argument 
  1. tiny fix for ZA-gof method

# v1.13.7.0

## New features
    
  1. add `RooDataSet` -> `TTree` transformation
  1. add generalized Beta' distribution
  1. add `isfinite` method for `SE` & `WSE` counters
  1. `data_statistics` : add the check for the finite counters
  1. add `isfinite` method for matrices and vectors
  1. add `isfinite` method for `(W)Covariance` objects 
  1. add `Ostap::MoreRooFit::Hypot` and `var_hypot`
  1. add Fizher's Z-distribution
  1. explicitely delete all created styles, see [here](https://github.com/root-project/root/issues/16918)
  1. add missing methods for `Ostap::Math::Bessel`
  1. fix occasional overflows in `twosamples.py`
  1. add new variant of _toys_ using generic _actions_
  
## Backward incompatible
    
## Bug fixes
    
  1. fix `ostap-check-dependencies`
  1. fix  bug in `Ostap::Math::Benini`
  
# v1.13.6.0

## New features
    
  1. add more methods to `Ostap::Math::ECDF` to get the ranking
  1. optimise `Ostap::Math::ECDF` and skip sorting when possible
  1. add `ostap.stats.twosamples.py` to make Two Sampel Tests
  1. add the test for Two-Sample-Test
  1. reduce code duplication between  GoF-1D and Two-Sample Tests
  1. Add Anderson-Darling and Cramer-von-Mises Two Sample Tests
  1. Add progress bar to methods dealing with unique/duplicated entries in dataset
  1. Add Hilbert transform
  1. Add Laplace transform
          
## Backward incompatible

## Bug fixes

# v1.13.5.0

## New features
    
  1. Add `operator+` for `Ostap::Math.ECDF`
  1. Add correlation coefficient to `Ostap::Math::pow` function
  1. Slight improvements in convolution for fitting
  1. Add the first Goodness-of-Fit method for multidimensional case
  1. Add test for Goddness-of-Fit method for multidimensional case
  1. Add very simple "efficiency-counter" `ostap.stats.counters.EffCounter`
  1. suppress `ostap.core.config.config_goodby` prints for non-interactive sessions
  1. add the most primitive splitter `ostap.utils.utils.splitter`
  1. modify `point-to-point-dissimilarity` GoF method: split into chubnks for large datasets, use parallel processing for permutations
  1. add keyword arguments to `WorkManager.iexecute` method
  1. add an option to run parallel permutations using `WorkManager` parallelisation instead of `joblib`
  1. extend Point-to-Point Dissimilarity GoF test for 1D case 
  1. extend 1D GoF test to include Point-to-Point Dissimilarity GoF test
  1. add argument `description` to `WorkManager.iexecute` methods 
  1. add parallelisaiton for GoF permutations and toys   
  1. implement tests for USTAT & DNN 
  1. prepend the default progress-bar for trees/datasets/frames with `Entries:`
  1. add a kind of replacement of `ROOT.RooAbsCollection.assign` for old versions of ROOT 
  1. add meaningful `description` argument to all `progress_bar` instance
  1. extend `gof1d` and `gofnd` tests 
  1. Add `RooAddPdf::fixCoefNormalization( vars )` for all appearences of `RooAddPdf`
  1. Add (fictive) context mnager protocol to WorkManagers 
  
## Backward incompatible

## Bug fixes
    
  1. fix the bug in `Ostap::UStat`
    

# v1.13.4.0

## New features
    
  1. Add  `cuts` and `cut_range` arguments for `ds2numpy` function
  1. Add empirical cumulative  distribtion function `Ostap::Math::ECDF`
  1. Add `cuts` and `cut_range` arguments for `ds2numpy` method
  1. Add method to `RooDataSet` to get empirical cumulative distrivution functions
  1. Add `PDF1.cdf` method to get CDF from 1D-PDF
  1. Add new module `ostap.stats.gof1d` for goodness-of-fit for 1D
  1. Add new test    `test_stats+_gof1d` for goodness-of-fit for 1D
  1. `ostap.math.math_ve` : add `significace`, `nsigmas` & `nsigma` functions to calcualte significance from p-values
  1. Make use of `ostap.math.math_ve.significane` function
  1. Rename `gof_1d` to `gof1d`
  1. Add method `weight` to `Ostap::PyIterator` to access th weigth of the current event
  1. add weight to `tree.withCut` and other tree-looping methods     
  1. Tiny tweaks for the treatment of drawing ranges for `f1_draw/f2_draw/f3_draw` fnuctions: form now the explicit setting has  a precedence. 
  1. `ostap.stats.gof1d` : improve drawing methods 
  1. Add Kuiper's Goodness-of-Fit estimator to `ostap.stats.gof1d`
  1. Add `values` iterator firnctoonfor all hjistogram classes to iterate over th ebin values
  1. Add methods `table` and `summary` for 1D-histograms to dump a histogram summary in the table form
  1. Add methdod `all_positive`, `all_negative` , `all_nonpositive`, `all_nonnegative` and `all_zero` for hisgoragm clases to check if
     all bins are positive/negative/non-positive/non-begrtaive/zero.
  1. Add methods `table` and `summary` for 2D&3D-histograms 

## Backward incompatible

   1. add `weight` to `tree.withCut` and other tree-looping methods     
   
## Bug fixes
    
  1. fix typo in `ostap.math.bnase.axis_range` 

# v1.13.3.0

## New features
    
  1. some updtae for catcing warnins
  1. ROOT warnings -> Python warning now used RuntimeWarning
  1. add new test `test_math_convolution` 
  1. add `filter='data'` argument for `tarfile` (TMVA&Chopping)
  1. set `cmake_policy` to fic nuilding for virtual environemnts (conda)
  1. more improvements and fixes for `ostap.io.dbase` module
  1. add test for availabe DB backends
  1. improve `ostap.histos.histos.book_histo` and make use of it in `fitting`
    
## Backward incompatible
    
  1.  From now `PDF3.fitTo` always returns the tuple of `(fit-result, frame)` or `(fit-result,None)` to be coherent with all other `fitTo` methods.
    
## Bug fixes

# v1.13.2.0

  1. Further improvemets for nice printout of linear algebra objects
  1. Add `lnorm` and `mnorm` methods for matrices to get L-norm and max-norm 
  1. Add `diagonal` method to get diagonal matrice 
  1. Add `abs` method for `SVectors`
  1. Improve Linear Algebra tests
  1. Improve the machinery for eigenvalues and eigenvectors
  1. Add possibility to avoid coloring of the header row in the tables
  1. Add protection for double coulmn (`::`) for the expression strings
  1. Resurrect and improve `ostap.stats.corr2d` module with simple 2D-decorrelation utility
  1. Make use of `terminaltables3` instead of `terminaltables` where/when possible
  1. Add new method `TH1.min_positive` to get minimal positive entry (or negative infinity).
  1. make a try to relace `distutilts` to `sysconfig` in `CMakeList.txt`
  1. some reshuffling of `CMakeList.txt` - make it a bit more readable 
  1. `ostap.math.base` : switch for C++ version of `frext10` (better treatment of `almost-zero-numbers`
  1. `test_*_toys` : make use of new `TH1.mini_positive` method
              
## Backward incompatible

## Bug fixes
    
  1. Fix tiny bug in `histo_book`
  1. Fix tiny bug in printout of linear algebra objects 
    
# v1.13.1.0

## New features
    
  1. improve all toys machinery: toys, jackknife, bootstrapping and significance
  1. extend all tests for toys
  1. add brute-force way to delete `RooDataSet` - needed for JAckknife and Bootstrap
  1. add argument `delete=False` for `dataset.bootstrap` and `dataset.jackknife` method to delete the dataset        
  1. add new test `test_fitting_dataset2.py` to test interference of memory and jackknife/bootstraping
  1. extend and improve machinery for toys
  1. add `std.string` and `std.string_view` into `string_types`
  1. improve `attention` printout from `dataset/tree` method `project`
            
## Backward incompatible

## Bug fixes
    
  1. fix the typo in `ostap/tools/tests/test_tools_reweight3.py`
  1. tiny fix for `progress_bar`
  1. tiny fix for `ostap.utils.basic.isatty`
        
# v1.13.0.2

## New features

   1. improve pretty-print for matrices
   1. add `pretty_array` function for nice print of arrays/sequences
   1. add some auxillary methods for matrices: `(min/max/minabs/maxabs)_element` and `(min/max)_diagonal` and `(min/max/minabs/maxabs)_element_index`
   1. improve printout of `SVectorwithError`
   1. add `pos_error/neg_error/errors` properteis to `Ostap::Math::ValueWithError` object `
   1. add column `@limit?` for printout of `RooFitResult` object to show the distance to the limit (if any). Distances of <3sigma, and <5 sigma are colored. Distances > 10sigma are omitted.
    
## Backward incompatible

## Bug fixes
    
   1. fix the typo in `Ostap/MatrixUtilsT.h`
   1. fix `ostap.utils.basic.mtime` 
    
# v1.13.0.0

## New features

 1. sight improvements for `bernstein`
 1. disable Thiele rational interpolation from the tests. sometnjug wrong with the code.
 1. extend a bit functionality of asymmetric errors (needed for graphs&plots)
 1. collect pretty-print functions into new module `ostap.logger.pretty`
 1. extend and improve pretty-print stuff
 1. make the dbase.Item a bit lighter
 1. For shelve-like databases: allow to specify the preferred backend as `dbtype = ...`
 1. reimplement `SqliteShelf` as `ZipShelf` with `dbtype='sqlite'`
 1. eliminate messy stuff  with extensions for `XXXShelf``-classes
     
## Backward incompatible

## Bug fixes

  1. fix a minor bug in `bernstein.solve`
  1. fix couple of recent bugs in `histos.graphs`    `

# v1.12.0.0

## New features

  1. Add estimators for harmonic, geometric, power &  Lehmer means and their weighted analogues
  1. Reduce code duplication
  1. Large redesign of staistics/projection& othe rmethids for RooAbdData/TTree/DataFrame
  1. Large redesign if `statvars.py` module
  1. Add `roc_curve` for making ROC curves, and corrresponsing test module 
  1. Add `eff_graph` for 1D historgams for creation of the efficiency graph
     from the 1D-distribution.
  1. Some tweaks for moments & counters
  1. Activate a new `draw` method (via `tree_draw`) for `ROOT.TTree`
  1. add `progress` and `report` optioal argumens for (almost) all Frame-related functions 
  1. Some tweaks for style configuration 
  1. update `ostap.utils.valerrors` & and new test
  1. allow to use `width` keyword when `line_width` is not specified for `XXX.draw` method
  1. add `loop` methdod for `RooAbsData` and implement `rows` in terms of `loop`
  1. allow more recusion in `vars_and_cuts` function
  1. add new test
  1. make progress bar silent if `not isatty()` unless explicitely set `silent=False`

## Backward incompatible

  1. `project`(&`draw`) for 2 and 3-dimession now follows the natural order of varibales:
       `XX.project ( target , 'x,y,z' , ...) `
  1. For `eff` & effic' and `efficinecy` methods fo r1D histograms
     the confusing  optional argument `increasing=True` is replced by (less-confusin)
     `cut_low` and the argument is not optionl anymore
  1. From now for weighted datasets `dataset[i]` returns `(entry,weight)` tuple    
  1. from now iteration over weighted dataset gives `(entry,weight)` tuple
  1. change sinature of `dataset.loop` , `dataset.rows` methods to return triplets `index, entry, weight`
 
## Bug fixes

   1. fix a typo in `ostap.ploting.canvas`
   
# v1.10.1.8

## New features

  1. Improve `addTMVAResponse` and `addChoppingResponse`  (and their paralell analogues)
  1. improve `parallel_copy`,  rely on `xargs` when `GNU parallel` is not available
  1. add add parallel `sync` based on `rsync -a` & `xargs/parallel`
  1. `sqlitedict` : tiny fix for warning in `python-3.11`
  1. fix `pypdf` examples
  1. add `ostap.stats.average` code for calculation of averages for inconsistet data
  1. add test for `ostap.stats.average`
  1. reshullle a bit the code between `ostap.math.minimize` and `ostap.math.local_minimize`
  1. add `sync_dirs` function into `ostap.io.files` module to syncronize the directories
  1. add `UseWeb` & 'useWeb' contetx managers for setting the WebDispaly into `ostap.plotting.canvas` module
  1. add command line option `-w/--web` optioon for `ostap` script to alolow defien Web-Distplay  
  1. add command line option `-p/--print-level` optioon for `ostap` script for better contol of the global print level 

## Backward incompatible

## Bug fixes

  1. `ostap.stats.combine`  fix calculation of p-value

# v1.10.1.6

## New features

  1. Update `Ostap::Functions::PyCallable`, `Ostap::Functions::PyCallable2` and `Ostap::Functions::PyCallable3`
  1. add new test `test_math_callable`
  1. improve a bit functions from `ostap.math.make_fun` module
  1. reshuffle code for Files/RootFiles/Data/Data2 toosl to colelct, keep and handle colelcitonof files
  1. add `table` method for `Files/Data/Data2` tools to print the content as table
  1. add version of `parallle_copy` based on `GNU parallel` (if/when available)
  1. make use of  `parallel_copy` in `copy_files`
  1. improve specific models form `ostap.fitting.specific` module
  1. add new test for specific models
  1. add new method `sPlot1D.adD_to_tree` for adding the sPlot results to the `TTree/TChain`

## Backward incompatible

## Bug fixes

  1. Fix a tiny bug in `ostap.logger.table.the_table` for wrapped columns 
  1. fix typo in `make_bkg`
  
# v1.10.1.4

## New features

   1. Add functions for the 1st, 2nd, 3rd and 4th unbiased cumulant estimators  
   1. add functions for cumulants up to order 10 + correct uncertainties for 3rd and 4th cumulants

## Backward incompatible

## Bug fixes

   1. fix a bug with summation of two `RooPlot objects` 

# v1.10.1.2

## New features

   1. improve prints from `PDF.load_params`
   1. add `smooth` function for 1D-histogram 
   1. imporove `rebin*` functoions for 1D-histograms
   1. add possibility to add separator line to `summary_graphs`
   1. fix for new ROOT>6.31 few issues with std::string <--> const char*
   1.  imporve a bit BLUE
   1. add "fix" for new pyROOT/cppyy for failure with ickling of enums
   1. fix morphing tests for new ROOT
   
## Backward incompatible

## Bug fixes


# v1.10.1.0

## New features:

  1. `tree_reduce` : allow redefininition of existing variables  (very useful for `tmva/chopping`) (only for 6.26<=ROOT)
  1. remove intermediate datasets created in `Simfit.generate`
  1. add `RRange`,`Prange`, `rraneg` and  `prange` loopers into `ostap.utils/utils`
  1. improve `VRange` , 'Lrange` loopers from `ostap.utils.utils`
  1. few fixes for `SelectorWitvars`
  1. suppress error prints from `Ostap::FormulaVar`
  1. catch C++ exceptons from `RooFormula`
  `
## Backward incompatible: 

## Bug fixes:


# v1.10.0.8

## New features:

  1. Add `split_chunks` and `split_groups` Functions for `ostap.trees.data_utils.Files` objects to split a large collection sof files into smaller chunks 
  1. Add `merge_chunks` and `merge_groups` Functions for `ostap.trees.data_utils.Data` objects to perform a partial merging 
    of ROOT files in the large collections  
  1. improve `hadd` function from `ostap.utils.utils` module 
  1. add `mtime` fnuction into `ostap.utils.basic` module - last createion/modification date for the path (dir/file)
  1. add (much) better cleanup of the ancient tmp directories. Usefulto remove lefovers from the parallel executions.  
  1. some improvements for `SimFuit.generate`
  1. fix clang++ bugs&warnings    

## Backward incompatible: 

  1. move `hadd` function from  `ostap.trees.data_utils.Files` to `ostap.trees.data_utils.Data`
  1. require `nEvents` argument for `SimFit.generate` to be `dict`-like type 

## Bug fixes:

  1. fix numerous typos in documentation strings 
  1. fix `SimFit.generate` 

# v1.10.0.6

## New features:

  1. add `more_vars`argument to `ostap.fitting.ds2numpy.ds2numpy` function
`
## Backward incompatible: 

## Bug fixes:

  1. fix old bug in `_rf_new_close_`

# v1.10.0.4

## New features:

  1. Update `histo_compare` tests 
  1. Slight optimisation in `Ostap::Math::ChebyshevSum`
  1. Further optimisation in `Ostap::Math::ChebyshevSum`
  1. add new test `ostap/math/tests/test_math.poly.py`
  1. Reduce usage of `Ostap::Utils::Iterator`
  1. add test  for `ostap.stats.ustat` module 
  1. Add `Ostap::Math::two_samples` function 
  1. Add the first version of code for RooDataSety -> numpy conversion by Artem Egorychev 
  1. Improve `ds2numpy` code and test
  1. small fixes for `ostap.utuls.split_ranges`
  1. add conversion to int for `RooAbsCategory` 
  1. add iterator/contains/len functions for `RooAbsDataStore`
  1. add some simple utilities for goodness-of-fit studies `ostap.stats.gof` 
  1. simplify `Ostap::Utils:::getWeight` for 6.26<=ROOT  
  
## Backward incompatible: 

  1. change the interface for functions from the `ostap.stats.ustat` module 
  1. change the interface for the `Ostap::UStat`  class 

## Bug fixes:

  1. fix a newly introduced bug in `ostap.utils.utils.split_range`
`
# v1.10.0.2

## New features:

  1. Add `Ostap::MoreRooFit::Rational` and `Ostap::MoreRooFit::RationalBernstein`
  1. Add `RationalFun` & `RationalBernsteinFun` FUN 
  1. For 6.29<=ROOT add option 'EvalBackend' and remove 'BatchMode'
  1. For 6.29<=ROOT make use of `ROOT::RDF::Experimental::AddProgressbar` utility 
  1. Add operators for `Ostap::Math::Rational` and `Ostap::Math::RationalBernstein`

## Backward incompatible: 

## Bug fixes:

   1. Fix a sad bug in `Ostap::Math::Bernstein` for incorrect usage of `elevate`

# v1.10.0.0

## New features:

  1. Add Benini distribution  `Ostap::Math::Benini`, `Ostap::Models::Benini`, `Benini_pdf`
  1. Add cubic and 4th order terms to (modified) Benini distribution 
  1, add the methods `min` & `max` to histogram objects 
  1. make use for `ROOT::TDirectory::TContext` for `ROOTCWD`
  1. imporve functions/pdf for Benini distribution allowing terms upto power 10
  1. use 'RoMinimizer' instead of `RooMinuit` for fresh version of ROOT 
  1. Improve treatment of `silent` for `PDF.chi2FitTo`
  1. Add `Ostap::Math::Rational`           : simple rational function inspired by `Ostap::Math::FloaterHormann` interpolant
  1. Add `Ostap::Math::RationalBernstein`  : rational function as ratio of Bernstein and positve Bernstein polynomials 
  1. Add `Ostap::Math::RationalPositive`   : rational function as ratio of two positve Bernstein polynomials
  1. Add `Ostap::Models::Rational`         : rational PDF  as ratio of two positve Bernstein polynomials
  1. Add `Rational_pdf`                    : rational PDF  as ratio of two positve Bernstein polynomials
  1. Add 1D-histogram parameterisations in terms      of rational functions: `rational_fun`, `rational` and `brational`
  1. Add 1D-histogram parameterisations in terms of rational functions: `pdf_rational`
  1. Improve a bit `tag` for 1D-Bernstein polynomials`  
  
## Backward incompatible: 

## Bug fixes:

  1. fix minor typos in `ostap.fitting.pdfbasic.py`

# v1.9.9.8

## New features:

## Backward incompatible: 

## Bug fixes:

  1. fix a bug in `dataset.duplicates` : not all duplicated entries were listed.

# v1.9.9.6

## New features:

  1. add `Ostap::MoreRooFit::ProfileLL` to allow bypass subtraction of the minimum for profile graphs.  From now the profile graph , obtaibed from `PDF.grpah_profile` with option `subtract=False` is not min-subtracted. It is useful e.g, for evalauation of discrete systematic uncertainties using profile-likelihood method  
  1. make attempts to avoid decolorisation of previously createe frames/`RooPlot` objects

## Backward incompatible: 

  1. reparameterise `Bukin2_pdf` in terms of `varsigma` and `asymmetry` instead of `varsig`

## Bug fixes:

# v1.9.9.4

## New features:
   
   1. add context managers for ROOT&RooFit random numbers generators 
   1. slight improvement for `combine.py` : add helper function `covMatrix` to create 100% correlated or uncorrelated covariance matrices
   1. slight improvement for `test_stats_blue.py` 
   1. add `parallel_add_reweighting` to speedup adding the reweighting information to `TTree/TChain`
   1. add parallel version of `sumVar` method 
   1. re-enable TMVA plots 
   1. make parallel fill of datasets more flexible 
   1. add methods `dot` and `weight_sum` to `SVectorWithError`
   1. add parallel versions for jackknife and bootstrapping 
   1. improvements for parallel versions for jackknife and bootstrapping 
   1. make use of [tabulate](https://github.com/astanin/python-tabulate) package - mainly to produce LaTex tables
   1. re-enabvle plots for TMVA 
   1. add `plot` argument to `use_canvas` context manager to print the plot at `__exit__`
   1. fix few typos in documentation 
   1. soem adjustment for `CleanUp` machinery 
   1. improvement for the output for `add_chopping_response` 
   1. do not add `sumw2`/`asymttotic` flags for fitting of `ROOT.RooDataHist` 
   1. allow to deal with the histograms in `make_toys` and `make_toys2` 
   1. add `storage` argument to `PDF.generate` to allow specification of the storage type for dataset 
   1. add `TH1D.eff` and `TH1F.eff` methods to make  "correct" efficiency historgams 
 
## Backward incompatible: 

## Bug fixes:

   1. Fix the bug/typo  in `padd_reweighting`
   1. fix parallel `addChoppingResponse` - the `category_name` was ignored 
   1. fix bug in parallel `addChoppingResponse` ( the `category_name` argument was ignored) 

# v1.9.9.2

## New features:

  1. add new test for exteding drawing 
  1. add new method `valid_formula` that is usefule for creation of formulas from expressions
  1. add helper context manager `random_seed`
  1. add new methods for `RooDataSet` : `unique_entries`, `duplicates` and `make_unique` to deal with "duplicated" entries (multiple count) 

## Backward incompatible: 

## Bug fixes:

 1. fix the bug in the `H2D_dset` and `H3D_dset` for `weighted=True` case 

# v1.9.9.0

## New features:

  1. improve and simplify `test_tools_chopping`
  1. tiny adjustment of printout format for class `Ostap::Math::ValueWithError`
  1. hide explicit manipulations with `Ostap::FormulaVar` into new function `make_formula`
  1. add alias methods/properties for `CB2_pdf`
  1. improve a bit drawing of combined signals/components/backgrounds
  1. add more drawing options, `draw_order`, `draw_singals`, `draw_corssterm1` , draw_crossterm2`, `draw_components`, `draw_backgrounds`

## Backward incompatible: 

## Bug fixes:

  1. Fix typo for treatment of `minos` argument for `PDF.fitTo` method 

# v1.9.8.8

## New features:

## Backward incompatible: 

## Bug fixes:

    1. Add missing method `Ostap::Math::Histo3D::integrateXZ ( y , ... )`

# v1.9.8.6

## New features:

  1. add `operator()`, `density`  and `density_mass` methods for `Ostap::Kinematics::Dalitz0` and `Ostap::Kinematics::Dalitz` classes 
  1. add `clenshaw-curtis` adaptive quadrature. It apprears to be better than Romberg' method.
  1. extend tests for Clenshaw-Curtis quadrature 
  1. Add modified PERT fuction and corresponding PDF
  1. Add more methods to `Ostap::Math::Positive`, considering the positive polynomial like PDF 
  1. add new test for `ostap.fitting.distributions` 
  1. Add Generalized Logistyc Type IV model with location/scale family and corresponding PDF 
  1. Improve evaluation of Generalized Logistic Type IV function, add more methods: `mode`, 'skewness', `kurtosis`, 'cumulant'
  1. Add `ResoGenLogisticIV` resoluton model and corresponding tests 

## Backward incompatible: 

## Bug fixes:

    1. Fix the typo in `Ostap::Kinematics::phasespace3` for the special configurations with zero masses 
    1. Fix few minor bugs/typos in `ostap.fitting.distributions`
 
# v1.9.8.4

## New features:

  1. add new example on making p-value scan (thanks to Dima Golubkov) 
  1. Fix `FrequentistsCalcualtor` and `HybridCalculator` to use cloned datasets. (They can destriy/corrupt input dataset). Clone ddatatsets are deleted after usage 
  1. simplify interface for `P0Plot.fill`
  1. `ostap.utils.utils` : add `CRange` and `crange` - helper utilities to generate range of values between vmin and vmax according to Chebyshev nodes 
  1. add `z1,z2` variables (and corresponding transformations) for `Ostap::Kinematics::Dalitz0/Dalitz` classes 
  1. add functionality to generate (weighted) x1/x2 and z1/z2 distributions for Dalitz configurations 
  1. add `s2x` and `s2z` methods for `Ostap::Kinematics::Dalitz0/Dalitz` classes for better unificaton of interfaces 
  1. add more tests for Dalitz< in parituclar (s1,s2)<->(z1,z2) mapping 
 
## Backward incompatible: 

  1. remove `use_onesided`  argument from `AsymptoticCalculator` constructor

## Bug fixes:

# v1.9.8.2

## New features:

   1. add Gudermannian function and its inverse 
   1. fix the issue with removed `RooMomentMorphND` for new version of ROOT 
   1. Add Airy functions and Fermi-Dirac integral 
   1. add more arguments to `PDF.sPlot` method, namely position and keyword arguments for the 1st fit (e.g. constraints...)
   1. add `ATTENTION` level for logger, corresponding `attention` method and `logAttention` context manager 
   1. add `Constrained(1,2,3)D` classes  to create constrained `PDF(1,2,3)`
   1. introduce helper `Constrained` and `Components` classes to reduce code duplication
   1. reogranize `Constrained(1,2,3)D` classes to decrease code duplication 

## Backward incompatible: 

## Bug fixes:

# v1.9.8.0

## New features:

  1. add `Ostap::Math::ExGauss2` function, `Ostap::Models::ExGauss2` and `ExGauss2_pdf` PDFs for the 
     variant of exponentially modified gaussuan distribution, but parameterised in terms of mode 
  1. add `Ostap::Math::Bukin2`, `Ostap::Models::Bukin2` and  `Bukin2_pdf`
  1. add `ResoBukin2` resolution model
  1. add more tests    
  
## Backward incompatible: 

## Bug fixes:
 
# v1.9.7.8

## New features:

  1. tweak `1D-integration`  for `Ostap::Models::Shape(1,2,3)D` objects
  1. add new test for 2D-shapes `test_fitting_shapes2.py`
  1. add new keyword `recover = ...` for `PDF.fitTo` that is expanded to `ROOT.RooFit.RecoverFomrUdnefinedRegions ( ... )`
  1. fix the names for internal integration functions to be coherent with underlying GSL methods 
  1. disable `Shape(1,2,3)D_pdf` for old versions of ROOT 
  1. add `Histo(1,2,3)D_pdf` objects 
  1. extent printout for `RooPlot` objects 
  1. `H(1,2m3)D_pdf - do not declare themselves as `signal` components
  1. replace Gauss-Kronrod integration by Romberg integration for `Histo(1,2,3)D` objects 

## Backward incompatible: 

## Bug fixes:

  1. fix `Ostap::Utils::hash_histo` and `Ostap::Utils::hash_axis`

# v1.9.7.6

## New features:

  1. fix `Ostap::Math::WMoment_` - zero weigths are totally ignored 

## Backward incompatible: 

## Bug fixes:

  1. set of tiny fixes for several 2D&3D-fit models 

# v1.9.7.4

## New features:

  1. add new test/example `test_fitting_simfit7` to compre simultenoud fit versus fit with constraints 
  1. `fit1d` : add suffin to the name for automatically created backgronus component 
  1. add `Ostap::Math::hotelling` function to estgimate Hotelling t2-statistics 
  1. make use of `Ostap::Math::hotelling` function in reweighting tests  
  1. add new argument `respect_groups` for `split_string` function 
  1. insert `rootException` for several stat-related functions 
  1. improve printout for `Ostap::Functions::Expression` and friends 
  1. add posssibility to enable global thread safety and implicit MC via .ostaprc configuirtaion file 
  1. update/modify/fix  `Ostap::Math::Moment_<>` and `Ostap::Math::WMoment_<>`
  1. add `Ostap::StatVar::the_moment` 
  1. add `the_moment` method for `TTree`, `RooAbsData` and `DataFrame` 

## Backward incompatible: 

## Bug fixes:

# v1.9.7.2

## New features:

  1. add `Ostap::Functions::Expression` - "universal" function (for `TTree` and `RooAbsData`)
  1. add helper `CallThem` utility 
  1. more improvements for classic reweighting machinery 

## Backward incompatible: 

## Bug fixes:

# v1.9.7.0

## New features:

  1. rewrite `statCovs` method to get statistics and covariances for trees and datastes  
  1. add `smattrix` method to `linalgt`-objects 
  1. improve printout of marices 
  1. add `mahalanobis` distance 
  1. improve `SVectorWithErrors`
  1. more improvements/fixes to weighting machinery 

## Backward incompatible: 

## Bug fixes:

  1. fix a the bug/feature for `statCovs` for datasets  
  1. fix a bug/typo in `asymmetric_kullback_leibler` 
  1. fix missing factor 1/2 in `kullback_leibler` 

# v1.9.6.8

## New features:

  1. add new module `ostap.io.zstshelve` with shelve-like databse using `zstandard` compression 
  1. add generalized Pareto distribution and reparameterised version of the exponentiated generalized Pareto distribution: functions and PDFs
  1. add generalized extreme value distribution: function and PDF
  1. impoetant improvements for classical reweigting 
  1. make more accurate "density" method for histograms  
  1. add `table2` fuction for trees and datasets  
  1. add progress bar to `add_reweighting` method 

## Backward incompatible: 

## Bug fixes:

# v1.9.6.6

## New features:

   1. (re)implement `Ostap::DataParam` in terms of `Ostap::HistoProject` - reduce code duplication 
   1. add few more utilitied to add branch/columns to TTree/RooDataSet
   1. release the limitations for `add_new_branch`
   1. TEMPORARILY  set `PYTHONIOENCODING=UTF-8` in `thisostapsh` . better solution is needed 
   1. Update `project` methods for trees and datasets 
   1. define PYTHONIOENCODING only for `python2` and only if not set

## Backward incompatible: 

   1. change an output for `project` methods   

## Bug fixes:

  1. fix a stupid typo in `table.py`. Thanks to Dasha Savrina 

# v1.9.6.4

## New features:

  1. add `+=` and `-=` operators for `Ostap::Math::HermiteSum`
  1. re-introduce `ostapfitting.funbasic.func_factory` for better backward compatibility 
  1. add forward declarations for Karlin-Shapley and Karlin-Studden 
  1. add more constructors to `Ostap::Math::Polynomial`
  1. add conversion from Karlin-Shapley and Karlin-Studden polynomials into regular polynomials 
  1. add more contructors to `Ostap::Math::Polynomial` class 
  1. unify the key-function for case-insensitive (and no underscores) dictionaries 
  1. improve `ostap.fitting.roostats` and corresponding test 
  1. more polishing for `ostap.fitting.roostats` and corresponding test 
  1. disable plot from Feldman-Cousins for ROOT<6.18
  1. for ROOT>6.18 for `ROOT.RooArgSet` extend `contains` and `getitem` to accept indices and slices 
  1. add generic fixed shape (no parameters) as `RooAbsReal`
  1. add fixed shape from histogam (no parameters) as `RooAbsReal`
  1. more polishing for `ostap.fitting.roostats` and corresponding test 
  1. add `front`,`back` and `pop` methods 
  1. add Karlin-Shapley and Karlin-Studden parameterisation for histograms 
  1. remove some script from `scripts` subdirectory 
  1. add more wrappers and utilites for `RooStats` (many thanks to Dmintry Golubkov for his examples and slides)
  1. extend `test_fitting_roostats` to cover intervals, point limit, scans - including constraints, resolution and efficincies that depends on ovservables 
  1. Add `Ostap::MoreRooFit::Minimal` and `Ostap::MoreRooFit::Maximal` 
  1. add `ostap.roostats.FrequentistCalculator`
  1. extend the test for `FrequentistCalculator` 
  1. add more efficient integration for `Shape1D/Shape2D/Shape3D/Histo1D/Histo2D/Histo3D` classes
  1. add `tag` parameter for `Shape1D_pdf/Shape2D_pdf/Shape3D_pdf`
  1, improve `test_fitting_roostats.py`
  
## Backward incompatible:  

## Bug fixes:

  1. fix `thisostap.sh` for usage with `zsh` 
 
# v1.9.6.2

## New features:

  1. add `ROOT.TTree.fproject` method for projection of the trees using `DataFrame` (the same as `frame_project`)
  1. add `Ostap::Math::IrwinHall`, `Ostap::Math::Bates` and `Ostap::Math::BAresShape`
  1. add `Ostap::Models::BatesShape`
  1. add `ostap.fitting.signals.BatesShape_pdf`
  1. add `ostap.fitting.resolution.ResoBatesShape`
  1. update tests 
  1. add  'TH1(F,D).bezier_sum_fill', 'TH1(F,D).bernstein_sum_fill', 'TH1(F,D).legendre_sum_fill', 'TH1(F,D).chebyshev_sum_fill' methods  for 1D-historgam parameterisations based on `Bernstein::fill`, `LegendreSum::fill` and 
`ChebyshevSum::fill` methods. Extend the corresponding test
  1. add `TH(2,3)(F,D).bezier`, `TH(2,3)(F,D).bezier_fast`, `TH(2,3)(F,D).bezier_fill` methods for 2&3D-histogram parameterisations based on `Bernstein2D::fill` and `Bernstein3D::fill` methods 
  1. add `TH(2,3)(F,D).bernstein`, `TH(2,3)(F,D).bernstein_fasr`, `TH(2,3)(F,D).bernstein_fill` methods for 2&3D-histogram parameterisations based on `Bernstein2D::fill` and `Bernstein3D::fill` methods 
  1. add `ostap.utils.utils.slow`  method for `slow`-iteration with delaye at each step
  1. add `Ostap::Math::agm` for complex numbers 
  1. add `Ostap::Math::agm` for `ValueWithError` objects 
  1. histogram parameterisations: add warnings for `fill`-based methods if polynomial degree is too large for such number of bins  
  1. extend `test_histos_parameterisation` for 2D and 3D cases 
  1. add `tag` method to several C+ classes 
  1. add ``SimFit.sPlot` method (background-subtraction for simultaneous fits) & extend the test 
  1. slight update in `ds_var_minmax` : try to deduce minmax when result is empty....
  1. add proper pickling for `ROOT.RooLinearVar`
  1. introduce `ConfigReducer` base class for better pickling/deserialisation 
  1. add Karlin-Shapley & Karlin-Studden positive polynomials (functions&pdfs)
  1. more polishing for Karlin-Shapley & Karlin-Studden stuff   
  1. remove `conf_interval`, `upper_limit` and `lower_limit` methods for `PDF`, based on `RooStats::ProfileLikelihoodCalculator`
  1. remove `poi` method from `funbasic`
  1. rewrite `ostap.fitting.roostats`
  1. add new test `test_fitting_roostats.py`
  1. largely rewrite  "Breit-Wigner with interference" model
  1. unify the variable separators for trees, datastes and frames 
  1. further imporvements 
  1. make a try to fix morphing 
  1. fix for 3D-reweighting, add 3D-reweighting test/example for  Paula Garcia 
  1. add (self)addition/subtraction operators for polynomial classes (`Polynomial`, `ChebyshevSum`, `LegendreSum`, `LegendreSum2`, `LegendreSum3` and `LegendreSum4`) with the same domain. 
  1. few steps towards better polinomial parameterrisatios
  1. Add polynomial parameterisation to frames (and trees) 
  1. disable some frame functionality when `ROOT.std.move` is not available
  1. add `ROOT.TTree.fparam` method for projection of the trees using `DataFrame` (the same as `frame_param`)
  
## Backward incompatible:  

## Bug fixes:

# v1.9.6.0

## New features:

   1. add `get_env` and `has_env` functions to `ostap/utils/basic.py` to check/access environment variables in case-insensititve way
   1. make use of `get_env` and `has_env` functions alsmost everywhere insted of `os.environ`
   1. small reshuffle of code between `ostap.core.core` and `ostap.utils.basic`
   1. first step towards usage of `ipyparallel` for parallel processing: ad trivial test `test_parallel_ipyparallel.py`
   1. make use of `ipyparallel` parallelisation
   1. add `Ostap::Math;:agm` and `Ostap::Math::ghm` fuctions
   1. improve `ostap/parallel/parallel_ipyparallel.py`
   1. improve printout from `ostap/core/config.py`
   1. add `$OSTAPDIR/.ostaprc` in the list of configuraiton files for processing 
   1. provide `$OSTAPDIR/.ostaprc` configuration file 
   1. more reshuffling of the code for generic and specific parallelisation 
   1. more polishing for the updated configuration     
   1. more polishing for the configuration  

## Backward incompatible:  
  
   1. rename `Parallel` section in configiration files into `Pathos`
   1. rename and move some `pathos` specific code from `ostap/parallel/utils.py` to `ostap/parallel/pathos.py`

## Bug fixes:

# v1.9.5.8

## New features:

  1. Add `Ostap::Math::ChebyshevSum::fill` method
  1. Add `ChebyshedSum` into `Ostap/Params.h` set of cuntions
  1. extend `test_trees_params.py` test
  1. redesign Bernstein dual basis: `Ostap::Math::BernsteinDual` & `Ostap::Math::BernsteinDualBasis`
  1. Add `Ostap::Math::Bernstein::fill` method
  1. Add `Bernstein` into `Ostap/Params.h` set of functions
  1. extend `test_trees_params.py` test
  1. optimise `Ostap::Math::Bernstein2D` and `Ostap::Math::Bernstein3D`, make them a bit faster and efficient  
  1. add proper (de)serialisation for 2D,3D&4D polynomial objects  
  3. Add `Ostap::Math::Bernstein2D::fill` method
  4. Add `Ostap::Math::Bernstein3D::fill` method
  5. Add `Bernstein2D/3D` into `Ostap/Params.h` set of functions
  6. further extend `test_trees_params.py` test
  7. further extend `test_trees_params.py` test
  8. add `__bool__` and `__nonzero__` methods for `ProgressBar`  - it allows to make more easy `while`-loops 
  9. change default table layout for `isatty` regime from `SingleTable` to `DoubleTable`
  10. allow to specify the default table format (`local`,`ascii`,`single`,`double`(default),`porcelain`, `github`,`markdown`)
  11. reshuffle code for `Ostap::Exception` - make it visible
  12. extend tests for `Bernstein2D` and `Bernstein3D` objects
  13. add `Bernstein3D::integralXY`,`Bernstein3D::integralXZ` and `Bernstein3D::integralYZ` methods
  14. add `Bernstein3D::integralX`,`Bernstein3D::integralY` and `Bernstein3D::integralZ` methods
  15. extend tests for `Bernstein3D` objects
  16. allow to define the default table style either via configiraitno file (section 'Tables', field `Style`) or environment variable `OSTAP_TABLE_STYLE`
  
## Backward incompatible:  

  1. rename `BernsteinDualBasis` into `Ostap::Math::BernsteinDual`
  1. `ostap.logger.table.table` : rename argument `format` to `style`
  
## Bug fixes:

  1. fix typo in `rames.py` for data frame projections into 3D-histograms
  1. fix couple of stupid bugs in `ResoStudentT` resolution fnunction
  1. fix bug in `Bernstein3D::fill`
    
# v1.9.5.6

## New features: 

  1. back-propagate Ostap::Math::Integrtaor toold versions of PyROOT 
  1. extend tests
  1. add (fictive) `Ostap:Math::Thiele::abscissas` method 
  1. add (fictive) `Ostap:Math::Thiele::values`    method 
  1. add missing `__reduce__` for `Ostap::Math::Thiele` interpolant 
  1. extend `MoreMath.h` and `math_ve.py`
  1. add `smooth_step` polynomial fuctions 
  1. add explicit functions for derivatives of Bessel functuons
  1. make use of explicit derivatives of Bessel fuctions

## Backward incompatible:  

  1. rename some methods for `Ostap::Math::Integrator`
  
## Bug fixes:

  1. Fix typos for `Ostap::Math::Integrator`

# v1.9.5.4

## New features: 

  1. add an explicit ouble-adaptive CQUAD integrator for `Integrator1D<FUN>` and `Ostap::Math::Integrator`
  1. add an explicit Romberg integrator for `Integrator1D<FUN>` and `Ostap::Math::Integrator`
  1. add a new test for `Ostap::Math::Integrator`
  1. extend and impove `Ostap::Math::Integrator`
  1. backport functinality for the older versions of ROOT/PyROOT
  1. fix the test for ROOT<6.18 
  1. add ROOT-version dependent switch in `add_new_branch` 
  1. some improvements for frame progress bar 
  1. some improvements for frame-based `tree_reduce`
  1. disable new frame-test for old ROOT 
  1. add treatment of new `ROOT.RooFit.MaxCalls` agrument 
  
## Backward incompatible:  

  1. rename some methods for `Ostap::Math::Integrator`

## Bug fixes:
  
# v1.9.5.2

## New features: 

  1. add `psi`, `digamma` , `polygamma`, `beta` and `lnbeta` functions, including their variants with uncertainties. 
  1. add these functions to `ostap.math.math_ve` module 
  1. add `Hagedorn` function/PDF
  1. add test for pt-spectra 
  1. add `integratebins`, `newstyle` and `parallelize` command options 
  1. add Skewed Generslized Error distribution (function, model and resolution) 
  1. add `sinc` into `math_ve.py`

## Backward incompatible:  

## Bug fixes:

  1. fix a typo in `Ostap::Math::beta` and `Ostap::Math::lnbeta`

# v1.9.5.0

## New features: 

  1. change default pickling protocol from 2 to DEFAULT_PROTOCOL (2 for py2, 3 for python (3.0-3.7), 4 for python 3.8-...)
  1. add possibility to define protocol via the environment variable `OSTAP_PROTOCOL`
  1. add possibility to define protocol via 'General:protocol' section in the configuration file 
  1. reduce code duplication for various databases 
  1. collect all pickle-related stuff into new single module `ostap.io.pickling`
 
## Backward incompatible:  

## Bug fixes:

# v1.9.4.8

## New features: 

  1. some massage for `BernsteinEven` and `PositiveEven`
  1. add `Ostap::Math::SkewGenT` and `Ostap:Models::SkewGenT` (Skewed Generalised t-distribution)
  1. add `SkewGenT_pdf`
  1. extend `test_fittiins_models.py`
  1. add (de)serialization for `KGaussian`
  1. add serialization for `ostap.convolution.Convolution`

## Backward incompatible:  

## Bug fixes:

  1. fix serialization for `QGaussian`
  1. fix `test_fitting_efficiency`

# v1.9.4.6

## New features: 

  1. simplify a bit `PDF.load_params` method
  1. modernize all hashing (affects integration cacheing) 

## Backward incompatible:  

## Bug fixes:

  1. fix incorrect tag construction for `Ostap::Math::BW` function (it affects integration cacheing)

# v1.9.4.4

## New features: 

 1. add new test `test_fitting_fitresult.py` for variosu expressions and their uncertainties 
 1. improve `Ostap::Math::Integrator` allowig to specify the absolute and rleative precisions 
 1. move `hash_combine` from `local_hash.h` to `Ostap/Utils.h` 
 1. improve `Ostap::Math::Piecewise`          
 1. add helper scale factor for the Breitt-Wigner function 
 1. replace `ROOT.Math.erfc` with `math.erfc` in `smear`-function 
    
## Backward incompatible:  

## Bug fixes:

 1. minor fix in `histos.py:smear`

# v1.9.4.2

## New features: 

 1. remove comma separator from `frame_project` and `tree_project`

## Backward incompatible:  

## Bug fixes:

# v1.9.4.0

## New features: 

 1. Add `Ostap::Math::FlatteBugg` , `Ostap::Models::FlatteBugg` and `FlatteBugg_pdf` 
 1. add addition for two `ROOT.RooPlot~ objects with the same structure 
 1. more tuning for `Ostap::Utils::FitResult` 
 1. add `ROOT.TTree.fstatVar` and `ROOT.TTree.fstatVars` methods 
 1. extend `test_tools_tmva.py` to include all 5 ways to use TMVA results 
 1. add `counters_table` function to printing the counters as a table
 1. fix/impove issue with standard TMVA plots 
 1. update TMVA examples/tests 
 1. add `full_path` methods for `ROOT.TDirectory` and `ROOT.TTree`
 1. improve `addTMVAresponce` functions 
 1. remove comma separator from `ds_project`
 1. slightly improve the prints from `tree_reduce`
 
## Backward incompatible:  

## Bug fixes:

 1. fix formfactor for `Ostap::Math::ChannelFlatteBugg`

# v1.9.3.8

## New features: 
  1. add `nbinsx/nbinsy/nbinsz` keyword argyments for `ROOT.RooAbdData.draw` method 
  1. allow additional keywords arguments for `ROOT.RooAbsData.draw` method, further forwardded to `ROOT.TH1.draw` method 
  1. add function `soft_multivar_constraint` to `ostap.fitting.fithelpers.FitHelper` for creation of 
     the multivariate Gaussian constraints.
  1. add example/test for using of the multivariate Gaussian constarines instead of simultaneous fit.
  1. add reduction for `RooMultiVarGaussian` class 
  1. fix `PDF.histo` methods

## Backward incompatible:  

## Bug fixes:

  1. fix `tree_project` for the case of multiple projected variables into 1D histogram 
  1. fix deserialization of `ROOT.RooFitResult` objects

# v1.9.3.6

## New features: 

  1. reshuffle some code between `ostap.fitting.variables` and `ostap.fitting.roofuncs`
  1. reshuffle some code between `ostap.fitting.variables` and `ostap.fitting.rooreduce`
  1. add serialization/reducing for the graph-like objects
  1. fix exampels 
  1. disable python warning from `scipy.signal` 
  1. add decorations for `TGraphMultiError` type 
  1. add helper module `valerrors`  
 
## Backward incompatible:  

## Bug fixes:

 1. fix typo in `dataset.py`. Thanks Dmitry Pereima for reporting the problem 
 1. fix typo in `tmva.py` 
 1. fix serialization for `TGraphMultiError` type 

# v1.9.3.4

## New features: 

  1. make `sPlot1D` serializeable 
  1. add new test for local sPlotting 

## Backward incompatible:  

## Bug fixes:

  1. fix serialisation for `Ostap::Functions::FuncTHND` and `Ostap::Functions::FuncRooTHND` objects 

# v1.9.3.2

## New features: 

   1. add reduction&deserialisation for `RooPlot` objects (it works better than defautl one) 
   1. add `items/iteritems` methods for `RooPlot`

## Backward incompatible:  

## Bug fixes:

   1. fix small typo in `dataset.symmetrise`
   1. fix small typo  In `Fit2D` constructor 
   1. fix some other typos 

# v1.9.3.0

## New features: 

 1. tiny tweak in `addChoppingResponse` 
 1. add separate `Lumi` object for `Data`-like set of classes 

## Backward incompatible:  


## Bug fixes:

 1. fix a typo in `frame_progress`
 1. fix a typo in `root_file`
 1. fix a typo in `add_branch` 
 1. fix a typo in `modifiers.py`

# v1.9.2.8

## New features: 

  1. rename argument `sort` to `sorted` for `Data`-like objects 
  1. add smooth transition functions 
  1. add functions for Kanaidakis statistics 
  1. add KGaussian function and correspondig PDF 

## Backward incompatible:  

## Bug fixes:

# v1.9.2.6

## New features: 
    1. change order of arguments for constructor of `Ostap::Math::QGaussian` and `Ostap::Models::QGaussian`
    1. add set of helpful fnuctions into `Ostap/QMath.h`
    1. improve QGaussiang model 
    1. add 2D Tsallis distribution for pt versus rapidity (to be validated!)
    1. add option to sort (default is True ) for Data-like objects 
    1. remove unnesessary ~__del__` method for `WorkManager` 
    1. add trivial filter `frame_prescale`
    1. add options `prescale_signal` and `prescale_background` for TMVA and chopping 
    1. add C++ progress bar 
    1. improve tree_project and ds_project methods 
    1. improve frame progress 
    1. improve frame project  
    1. add progress bar to `Ostap::Trees::add_branch` and `Ostap::HistoProject::projectX`
    1. add progress bar to `Ostap::PyIterator` 
    1. tuning `frame_project`
    1. tweak parallel_test_toys 
    1. issue warning message for `AsymptoticcError=True` for ROOT<(6,27), see [ROOT-PR-11282](https://github.com/root-project/root/pull/11282)

## Backward incompatible:  

## Bug fixes:

   1. fix typo in `parallel_toys` - Thanks Dima Pereima for reporting the problem
   1. set of minor fixes 
   1. fix recenly intorduced bug in pyselectors 
   1. fix a typo in `tmva.py`

# v1.9.2.4

## New features: 

## Backward incompatible:  

## Bug fixes:

 1. fix newly introduced typo in drawing for simulltaneous pdf. Thanks to Dima Pereima for reporitng a problem! 

# v1.9.2.2

## New features: 

  1. add option `parallel` for `Data` and similar classes  

## Backward incompatible:  

## Bug fixes:

# v1.9.2.0

## New features: 

  1. add `fixdeps` argument for `Fun1D/Fun2D/Fun3D` objects to fix missiing dependencies (or, to add some fictive depebndencies) 
  1. reenable `linalgt` test with `numpy` objects 
  1. add one more test into `linalgt`  
  1. add methos `kullback_leibler` and `asymmetric_kullback_leibkler` into namespace `Ostap::Math`
  1. add method `kullback` to `FitResult`
  1. extend linalg2/t modules 
  1. more polishing of the linear algebra 
  1. extend tests for Linear Algebra operations 
  1. make `styles` to be the class property instead of the class methdod for class `StyleStore`
  1. imporve treatment of the DataFrame/RNode 
  1. more imporve treatment of the DataFrame/RNode 

## Backward incompatible:  

  1. remove `keep`    argument for fun/pdf objects. Hopefully it was never used by the users.  
  1. remove `special` argument for fun/pdf objects. Hopefully it was never used by the users.  

## Bug fixes:

  1. Fix a bit strange problem/feature appearing at 2022/08/11 in dev3 slot: drawing of 
     `Addition` objects with `RooAbsRealLValue` fails.  Fix is done using `FunNop`
  1. fix old typos in `_h3_integrate_` method (thanks to Ivan Polyakov for tproblme report and the fix)  
  1. fix bug  in `MatrixUtilsT.h`
  1. fix typo in `MatrixUtilsT.h`
  1. fix minor bug in `useStyle` 
  1. fix compilation error for gcc12 
  1. fix few typos in `frames.py`
 
# v1.9.1.0

## New features: 

  1. re-enable again `test_fitting_morphing` for new ROOT, 
     see [ROOT/issues/#11061](https://github.com/root-project/root/issues/11061)
  1. add `MorphingN3_pdf` for morphing in 3 variables 
  1. more owrisk on easy serialization.  Now we can bypass standard serialization fro almost 
     all important Ostap classes.    
  1. make Model2D & Model3D PDFS more safe 
  1. fix `RooGaussian` serialization fo rOLD version of ROOT 
  1. add serialization for `RooFFTConvPdf` instances
  1. add serialization for `RooSimultaneous` instances
  1. make RooCategory more uniform "interface" for RooCategory 
  1. split `variables.py` into `variables.py` and `rooreduce.py` 
  1. add the proper reduction for the effciency objects
  1. add serisalisation for `RooEfficiency` 
  1. add serisalisation for `RooFitResult` (the standard one often gives segfauts)
  1. reduce verbosity for `make_var`
  1. more polishing with verbosity for `make_var`
  1. more polishing for serialisation 
  1. add reduction for Breit-Wigner related PDFs 
  1. improve BWI model and pdf, add dedicated test  

## Backward incompatible:  

## Bug fixes:

  1. couple of (small) fixes in `variables.py` module 
  1. several typos are fixed in construction of 3D-models 
  1. `funbasic` : fix typos 
  1. fix the treatment of shifts in `Convolution` 

# v1.9.0.2

## New features: 

  1. add some finite functions: `hat` and `up` and correspoondg models/PDFs
  1. add finite atomic function `fupN` and corresponding PDF 
  1. reduce code duplication for variosu parameteric functions
  1. make more active use of `Parameters`  
  1. add more `swap` functions 
  1. fix for `par/setPar` methods for OLD root Where using statement does not help 
  1. more polishing for reduconig of RooFit an Ostap objects 
  1. Add `HORNSdini` and `HILLdini` functions/pdf 

## Backward incompatible:  

## Bug fixes:

  1. fix incorrectness in `Sum1D/Sum2D/Sum3D` - `fractions` argument was not forwardef to `

# v1.9.0.0

## New features: 

  1. redesign structure of base classed for functions and pdf 
  1. rename major base classes for fitting, e.g `MakeVar` -> `VarMaker` etc
  1. make operations for functions and PDFs much more robust 
  1. rreshuffle fittling classes between modeules 
  1. eliminate `ostap.fitting.basic` module 
  1. Add Rice functioni and corresponding pdf 
  1. Add Generalised Inverse Gaussian function and PDF 
  1. Add ExGauss and NormalLapalce functions, pdfs and resoltuion models
  1. add Generlisez Argus distribution - function and PDF 
  1. simplify  hierarchy for some peak-like models/PDFs 
  1. rename `MASS` -> `PEAK`, `MASSMEAN` -> `PEAKMEAN`
  1. imporve normalization for Pearson Type IV function 
  1. add `StdMoment` to `ostap.stats.moments`
  1. sdd `std_moment` method for `PDF`
  1. mase a bypass for long standing issue with segfaults from `RooFitResult::globalCorr`
  1. Add Novosiboirsk function and PDF 
  1. add more properties for `SinhAsinh` function 
  1. add `LinearMorph_pdf`  
  1, remove '/' and '\' from constructed roo-names 
  1. add bspline as RooAbsPdf 
  1. tiny tweaks for serialization/pickling/unpicklig for some classes 
  1. remove more warnings for old ROOT 

## Backward incompatible:  

  1.  rename `Morphing1N_pdf` -> `MorphingN1_pdf`
  1.  rename `Morphing2D_pdf` -> `MorphingN2_pdf`

## Bug fixes:

  1. fix a bug in evaluation of integrals for  `LegendreSum2` and  `LegendreSum3`
     Tnanks to Ivan Polyakov for
     reporting a problmes and the fix 
  1. fix `Dalitz0::P_R12`, `Dalitz0::P_R23`, `Dalitz0::P_R31`. 
     Thanks to Ivan Polyakov for reporting the probles
 

# v1.7.3.0

## New features: 

   1. add resolution model based on Generalised Gaussian v1 distribution
   1. re-enable pritout of the Global correlation for `RooFitResult` objects 
   1. add 2D Gaussian function and PDF 
   1. add 3D Gaussian function and PDF 
   1. adjust integration precision for all numercial integration calls 
 
 
## Backward incompatible:  
 
## Bug fixes:

# v1.7.2.2

## New features: 

## Backward incompatible:  
 
## Bug fixes:
  1. fix minor typo in `CMakeROOT_6_*.txt`
  1. fix typos in `ostap.utils.pdg_format`


# v1.7.2.0

## New features: 

  1. add new method `make_soft_constraint2 to create asymmetryc constraint
  1. add test for soft constraints `test_fitting_constraints`
  1. add Thiele rational interpolator 
  1. reshuffle code between `Ostap::Kinematics::Dalitz0` and `Ostap::Kinematics::Dalitz`
  1. `Ostap::Kinematics::Dalitz0` allow call for angular functions with floating `a`
  1. `Ostap::Kinematoics::Dalitz0` and `Ostap::Kinematoics::Dalitz` : and methods to calcualte Wigner angles 
  1. add methods `random` to `Dalitz` and `Dalitz0` to generat erandom distrobutions in Dalizt plane 
  1. ad) x1/x2 anF/x2 acorrpix d) e 
2yss` method6f mareatures 

  1. ,fle some cns ap.fit
Oa 1. imp1. ,flemes and inalgt`  
 Meduction for)s 
 -lel.ul Fuction for)s F` anvreports  for 1Dme report oorphingbesponding test `mname `Morphihistor18 for `Ri's://github.yal` argumen`Luithub.yal` argumen## Nec:i`PGauss and Nlfictive)  v1.9dd new tear`

#reduce`
  1es:`Bernstein2D` and `Bernstein3D` objects
  13. add `Bernstef the treel` argumen## ead of simultaneourimultaned `ROOT.TTreeknction and Phodsi simultas and pdf multanee:BernsteinDu 1. exte add `Betemenb_anee:Be` 
  I3)D_pf`
  1, improve1.  rename `Morpinterpoing.signals.
2 add method `kullback` to `FitResult`
  1.imake `lore pian consabsolute Be` 
 :dsstr. tin1iethod `
  1. adlaa bug in bof ame as ` ap.fityk-propagat features: 

  /root/pull/l/:tprogr2,3)Ds constructor of `O the tru 1. exg fixes:

1drnstein2D/3D` tHPnemato the cla1n
  orm "intolute Bo the nampdf 
  o th . Add `Ossagat features: 

  /radd,fle sotiono 
  1. fix desn. 
  1. fegus disttap:_` ob(2. fix/_` ob(2oks to Ivan Polyakov for
  eatures: 
s for `TGraeSum2` 
    1n for `Roo-ranch` and `Os_anee:P & Modfor aap::Math::for `RooFFTConvPdf`  reporti Nec:i`Pye_ipynot dd GenerlreSum`o  1. add utilites for `RooStats` (many thankss:

 1. helper scale factor anager`,ction andd dedicavPdf` instanw featurutil  1p.fle hankss:Romberg integrator for `Integrator1D<Fs for 2&# v1.m(-`
AaTn offor `I s
1drnstein2Dsian cee` (the samer alge`hods. Extendi Rf `O theults from or variosu peshudd 3D-rewnaTn offor `Ieral stat-related fu to Ivre salasse12`I y
  1. aAds for `P d Nor2,3)Ds constructor d `kuar`ization/cee`+6w-e 
 ,ye_ipynot ds_pro@iWsteostap.ug fin 1. adjust int//u3 ROOT>  

##>  

##>  

test  

)#>  

testeurst_fon modefixes:
e<IV functioibetweeuioibetww features:modefiarator for `Inack` to `FitResult`
 Brame_project`
    assesFhods.LFitResb/extent priname_projs  for (hiinto namespit xixeFs F` anvrep,functions: `hat` and `up`  add into nao_fiaraithubROOT.Robz the nampdf 
 re_project`
corrupt inpuSs  for  1. ficts
b
  1 name`_soft_ction aed Gaussian vojecghe nampdf 
 r1iethore
## New ug fixon (3.0-3.7),e: 

## Bg-3.7),ee nao_fia7),ee nao_fia7),eeape3D/ consak-likeFting` method 

## BOy`
 
# ixes:

1drnshfor `Integt~ objent pris   1. add resolution mctivbtics 
  1.sian fX 3at` 
  1. r1_pdf`_prohang2a ~_ncvrom orsimu  1.ew features: ment for fun/p## New tures: mentimultaneoajonE` objects 
add re nao_fia7ang2a ~` objects
  13.jeprting feature`, fielg PDF 

# `O`
  X,3D&4D pols indult table layout for `isaig:Dalitz0` rTrng` method b/lts od 1. etHeics yI y
 )cghe namw features:modefiarator for `Inac3.jeprtingaram` methodPA`
  1.onal interpolator- method 

## BOy`
 
# ixes:

1drfor `TGrae7jh` . bekeep`    argu utilitwly introduc
for s for `e naRrts Trng`gu utilitwLly introduc
for s for `e naRrtsm/g`gu utilitwiitwLly introduc
for s for `e naRrtug  in iy into `Pathos`tegration precible la1. fix a typo in `.ei2
  1or s  a typo in `test_fitting_efficiencspilmWA'bPOuction foLfix maittingX 1. Fix a biing_ improvemedard onficts
b
sage 
  1.ete 
  1. several typos arH='-soft_constixev1 distribuaanded to `ROOT.R fee 
  1rd incobacx a biing_2

## New test_fit` biing_2sPures: 

  /ra. reduce csini` functions/pdfge for `Bernstsian cejecghe s:
jent prTTini` o `Dalitz` and dnterfacom or v,ted tes-n ae1eonst1n cejetingaram` m`num:sting `a`
  1. `Ostap::Kinemattistics 
  1. add KGaussian functmlm,DEJ.resos `StyleStorle:  
agat feafityk- KGaussian PolyAssii
Pics::Dalituc
for in `.ei2eatnted tes-n amn.eadtuddenres: me`
  o th F::Daorle:  cale fac.-o in `. Generalisor 31. make `style)ke `style)ke 
  1. change default picklingcns ap.fit
Oa 1.(y Romberg6ie1encie1b.com/rops `randos    1.x` m`num:sTm( a bug ixon (3ator forbT 
  : 

## Bg-3.7),ee naD3Ad forward s:dos gcns apa o 1p.fle hankCing` for ngs for `fill`-based methods itablefloating T 
  : ngstions 
`
    p 
iparallel/utils.py` t hankLD rojent prTTiues:m/T, improve1.  r IvanntX`
  ot dd Gener:m/T, imp  forimprovev:m/etio 1.eajonE`tmlm rojent.eduie Gaus-> `Varimultafor newnstatemen in `wP Bg-3.7Polisu  1. apix d)for `Rmpatiing T 
  : ngstjoitz0::P_R23 

 ug fixes:
  
#I  8. addes: 
the class propee naD3 `Os_aipy th F::staperiaKlass oLfix maT 
  : ngstjcgtyle)ke 
  1(s o_:a biicns a## Backward incompatible: 
g for simulltah more robuNolf)1 a#able
 1. fi Thanksorsimdes: 
le
 1ge default pm:sting& argument fo`at` afault pm:st newDAlward i pm 1orsa:m 1t`)
Oa 1.(y Rombatribution for pt versus rapidity (to be validated!)
 variosu peshu. improve framlator`
 on Typeo_fiaraithubROOT.Rob forimprovev:)varidf I s
1dr3. add 2pcu  
#I1 a#a indult cgtylyuAd focatures: 

 e sicatureof `O -1idf Ib:

# v1.9.4.0
 2pcur`
 stling classes between modeules 
  1. elhon ew featuF ince protogram s inu%malizatcur`
educe.py` 
 izatcur`
educe.py` >/`,`eH >/`,`eH >/`,`eH >/s: mepcu  
#I1 aS:

# v1.2'Ladd re nao make `sR)M`Ladd re nao tns `pr-istic (hiinto namespit xixeFs F` a1. ayuA/y2mespit xixeFs F` a1.o  13. abspanw fclass propee naD3 `Os_aipy thytion and /t7` to compre sim
  1. fix in `pr-istic (hiint- y creixeFs Ferl.fix in :

# v1.2' old ROOT spass prop-3.7),e:_ncvrom osult`
 Brame_promsicatureof `O -1idf Ib:

# v1p `prv1.2' old ROOT ePth::HermiteSud `tagc: 
g hods itablefloais feeD(s o:modefiarat. meths ita standing isjects r tree- >/`,`eH >/-s for sever_/utils/bcdd!)
 ward sc(o.pinewly introcts rtions f-ld RO[wrisk on easy serialization.s/bcdd!)
 `,`eH=epm 1orsa:m 1dgbs F` a1.o  13. abspanw fclass pro add `StdMoment` to `ostapMakin `CMakeROOTro add ` 
g hodsltaneoajonclasres: fix ` be validatealitz0` to er via ward incompacon-_anee:P & Modfor aap::Math::for `R/b
sage st foafor ve fosseu `tree_pby 1. ady 1Dalitz` an`by 1tL variosu ped"- tein2D` s `ru plrd inco/`,``' old ROOT ePth::HermiteSud `tagcmalizatcuragngN1_pdempty....
  1. add proper pickling foerialindos  M`n for allncompacoimproe` 
  I3)D_pf`
 ackward imatically create 2ru%malizatcintrocts roe` 
  Illback_leibkler` in a1.erbos7.3.0

## New features: 

   1`Fiuce.py`  1. imph
educe.py` >/`Ce)  x ( anSixes:
tx a bug ix a typo iaus oard # v1p `prv1.2' old ROOu.fpr speci:(is:

1drllback_leiSs  fg fixeatrixUtil` -> `VarMaker` etanee:P & Mod
 ug f  1. addion (Da massage fo_ 1. impe)keIMVA exkeIe st foafgd.5.nm1.eruve) fclainy fis_pdf`
Lx/p.f

# v1.9.7.erty newdel__` m `Integame_projec
  1.  13. add `Bernstei Backward incec Ne  newdelBernsteillback_leiei Backward inpee_aPreofiteSud objecxawm` for complex num biing_2a.draw` method, furt
 

#1f  1.Mod_ncvrom osu0 ng feaa1.er (Da masst`
 1. sligh osu0 ng feaTiteSudturesP 1.  renamns fb of c:Math::fiteSud inco/`,``' old RO calcuai Backward inpee_aPreofiteSud objecxawm.
  1. argumen## inixawmxawm`ethods for ` fea-jent p1. fix a typotib
  1. adessBarh Wecxawm` fo_2HermiteSler` in margm 13fxawm` for cossian vojecgharaiion cac1. Vte "density"  8. 2HeD_pdf/Shape2lastini
## Backwng_mma`, `b prope.py` for variosu exprev1 distrito pe and rleasincompatib# Back( Inars` metho2d-subt.yal`   1. steinDumpatible:  

## Bug fixes[tolute Bo the nampdf 
  o th rkeSud `tw modeule`F

 1.d incoolutd `mstwetp  1. reinpee_aPxes[tolu:ed!)
 vardy, `newstylea for cossia 1.ns fbrd isy serializationa

 1. remo# Bug 2pcu   1. add seri. Vax when retolute lass pLmprove. change deth F::Daorler* 2pcur`xientX`
 improvements db:

# v1.9.4.Tsym `SingleTa singllute lalass pLmprstsobjcle` to ``
  1. reo_su peshue rat for, to addlfixes/radd,flDNewo `te` gumen## d ot dd Gener:va> classgiblointroducia 1.ns fbrd ers_tab: 

  hu. ina 

#  

Ae_pbystures: 

 1.4is:

1drl)ke 
 2sas` mina 

#  

Ae_pbystures: 

 1.4is:

1drl)ke 
 2sa_y:e  newdel
.0

## New features ve fosso exprev1 dion cacy.RAS8ffittinOstap/Paramctive use ofPUern'ns for ROOT<6.o)T<6.o)T<6.o)T<6.o)T<6.o)T<6.ipy th F::dd redemaGaussial fuctions

##)neshuf-b
Ae_pbfor `(Xons

##)e_pbforixUtil` nsitor cossian dedemaGaussily  1. fix incorr
 Pisicatureof `O -1imaGaussial fSs pLmprod_nnor-Sxtend tessily  1. fix incorr
 P Be` 
 :dsstr. tin1iu peinalts 
 1. n tess (b e.g andylesP,. aAdzdializatifault pI  8 1. simp`/H:dss_R12`, `Dalitt.g aAdzdializbaeegralYslly create 2ro)T<6.o)T<6.ipy th F::dd redemaGaussi`Dalitt.g aAdzdbmpatible:  

## Bug fixes[told 1. fix incorr
 Pisics  fg fixeatryaddition. add  mlm rojGeno.ert foreo` metho2.d incoolutd siaon pt` a1. ayuA/y2mespiomp2o
## New features: 

   1. ad Karlin-Shapand corrat. me,rrat. me,rraddition. add  sefloais femph
educe.py` >/`Ce)  x ( anSixes:
tx a introdatryadULT_PROTat for, t
ed/y`Roo-bible: ap::ModelpBernstein3atrodunalts__` m  functig feaa1.e.(1. n tess ()in `fr.2

## fill`, `flDNewo ?>/`Ce)  x (fse TMVA rackwng_mmancoolutd siaoew feat>odefid triprojng_m1idf I23atrodu1te SKroject` - re/`Ce)  x (fsPlmanco namAgt`  
i2lon f `Pas::integralX`,nable pritout   1. steinDloeduce.py`ufs a thytioal` & r`ufs a thytioal` w:Math:in,s_yt fpo in `MatrixUtilsT.h`
  1. 0t.g aAdzdializbaeegralYslly cS w:MatmaGaussalituc
fodafeafit pdf oes tharic_kullbacknaRrts Trng`gu u:  1.)g pdf 
 Gaus-> `Vau)T<6.ipyfres:modefiarator ''''''fres:modefode duro add `StdM (ator`

## Buion cacy.Rrialind`as_e'a3modefiarat_ann cach Add N2pt_fiaraithubROOT.Robz the nampdf 
 re_project for `Ri'sIes: 
loSnd fw-e _n.m2`  1. Add Novosiboirsk flBasis`
  1.)Prie Gaus-> `Varts remo# Bug 2pcu   1. ampd  mlm rojGeEtree_pdemaGt_fiaraiiagduce.pariables.p`.2

 -> `Morphinhe nam` metho2d-syyford 'atible:Dngara4y distribut   1. asfitting.signals.Bateuon for `R, `b14/"tdM (ator`

## Buion cacy.Rrialistap::Maldwnpee_aPreoramctise preci`  1. adessBarh Wecxawm` fo_2HermiteSlerLly crea c:Math::fa,
<8s amp'acy.Rrialind`( `tagcmreci`  1. adly csameo `Pa,alistap::Malath::fa,
acy.Rrialind`( `tagmethod6f marethubROOT tharic_kulw tead 1. amples/tests 
 1. add `full_paCe)  x ( anSixes:
tx a introdatryadULT_PROTat /tests 
:nt prTx
#I  8.  fill`,F/S
  oa_d<eno.eints`
rixUtilsrures: 

   1. addp::Math `DalintrodatryadUL ad-o innnnnnnnnnn into nard inb1 

  obaeeg1i A,ye_ipixes:
tx a or `R, `b14/"tdfor `R, `b14/ath::1. r2T_PROTat /tests inflDNewo ?>/`Cml)mancocu  
peshu'S
  oa_d<eble:  
Th `Dalintrodaaods 
  1. A,ye_ipixeang2a ~alcua-a/osu0 ng feaa2HermiteSlerstrobitap::Math:`VA andteRLKCaandt  

## ap.ues`ble pritout   1. steinDloeduce.py`ufs a thytipstandard on
 1.  v1.9S8ffTConvb 
  1. add `integratebinddition. addQ fill`,F/ame_projec
 0jec
 0jec
 0jec
 0jec
 01nnnnT 
  :robitap::Math:`VAting_effic   1. fitout ipixeang2a ~alcua-a/osu0 ng feaa2HermiteSlerstrobitap Buion caith f base claeralis-hfast`, `TH[ranc the d N2pt_ftap::Math:`VAting_effic  otap`ap::Utils::effic  otap`st`,e claeralis-hcEipy Lw` variablc:Matawm.y 1. ad<eble:  
(the samerefficVAting_ets__i  

#:teSler6/to use TMVA results 
 1ontroduc
eor `Funx `Dalintrv --reive use o stribhingN1_pdf`
  1. provemedroducntegrat--reive npixv `b14/ath::ROOT/PyROOTroduc
eogus u,nstraBackwardp)ckward 
   the 8y_fillpropagat fete lalass pLooLinearVay.g AS8ffinnnnnnnn1. re
  1. make more active -berns`
r soft consg_formafillprop5mftive use ofP) 1. amme `a moe. chags 1. fi sim `

# .py`. Tmodeule`F

 1.d incoolutd a,t peknctio.)'  nals'rat--reive aric_kulls'rat-dd 2pcu  
#I1 atraBackwardp)(nct robuNaPxess ap.  
:ef  d inb1ardp)ckwaruf-b
Ae_pbfFunrnal osingllut<" `b1ree.  yon warnroduosinglpythoT.h`for polynomial cnx `Da for cr:Math:o` an`ad<ebluncti
2 ad:.7.2.0
dney the absolun
  orm "intolute Bo the namd:.ults 
-e _n.mwnks DR`for poly`

#  

## Bug,tol1for cr:Math-b
Ae_erLllesP,. a)Tw-e _nting& arg./oEKarlinos gcoerLinorz0` allow call foryAssiin `addChng test `mnaddChng maGaussial =base citapng maGaussilidated!# v1ets__i  

#:thytioa BOy`
 
# ix_ typos ins pLooLinearVay.g onsg1  1. adDucntegrat--reive D))))))))))ing classeioa BOy`
 
# ix_ typos ins 

#:teSler6/toBOy`
 
# iiiiiiiiiiiitagcmreB. imgratsly cMVA rack-ancocu  
peshue modd KG et 
  :nstraints `test_fittinA rack-a 
# ix_ t 

`. add mes 
h`SGaussial =b%i

# v1.9.4 me).e.(1. ns =b%i

# vtobuNaPxess ap.  
:ef  d os.py`

## Bug fixes:

# v1.9.5.8

## New features:

2`O the12atrimafCt<" `b1ree.ctive use ofP.)g pdf 
 Gaus<" `d `kd` `par/c to al .ere ofP.)gatrimafCtple cDuCaod
  )sme major base classes for fittindd Generaug  model basednres:d# Buion cacyPDF2Karlinos gcong classntegrat-ss ()ires:d podf`_prohang2ls
  t_hisypo atrimaf t-dndreSumr `R, `f cla/atibleangrawing for s   1. ad Kar> classgiblointrodutions1 otimultaneoajonE` ethod6scouce.parwew ROOT3game_projecwtroduth-trinnnnnnn1so p

# v1p `prv1.2' :e)  v1ler`l Fuction f)ouhng  SKrojath:`Veulls'rat-dd 2pcu #)e_ ins 
My
# v1.  model baseot patibl),ee odd KG et 
  :nstraints `test_evEww featu call foryAp.fittiont--reive T`ococu  
pese T`U+z )the graphveive npim! 

.2

 -> `Mornes 
  , (the sa6ope.py`_paCe)  x ( anSixes:
tx a i-reiveeFAso p
m:sting&mere ofP.ing&mermmiteSud `tagcT` reBachese functiod i pm nt--r
fP
:nt prTHvo)T<6.o)T<6.'azjetingara/p.f

# v1 fif`
es: h::tor d `k nd PDF 
   1igner with inte ( anSmcti`Funx `,,F old typos i/hu. w` vWnclape` biingcts
b
  1 n 
  I3mporveatawm.y 1. ad1so p

# v1pc oes: 

  1. -ancocu  
peshue moddnhe p<l
oddnhe p<l
odde duro add `StdM (ag& arg./oEKarlnfigiratshue moddnhe p<pae hank232ap:tap::MiS.(ag& arg./oEKarlnfig1igner with im<l
odde duroaproper reductions 
  1.ons 
  1.ons::Funaddedempty....
  1. addnnnnnev:)vggcmalizg f`.9.4 me) rle. adoFFTConvPdf::Funaddedempty....: 

 1. add n1 aS:
sa ~alcua-a/o  orm "nnnnev:)'Rooaf t-dnilutdL
 
# iiiiiiiiiiiitagcmreBcu  
peshuer:*fix a bug ipjectsll foryApY2 make use o.5.8o itap/ a)Tw-ebypass for long standing id*fix a bug -lny,Lnction and L(on for thf base cla/github,`eH >/! 
 for fit_oc
   Trng`gu uag&onse` 
 1. add separatpropeernsgarning /btion and- thentitoyfdf 
  1. rename maeoajonE` ethod6scote Gaussia3s gcong onE` eps` arg ``
  1. reo_su 1lizatifaul `LegendreSu` eps` argu.yx:
t_
  1. reariable `atifaul s prop-3.7sgGaus<"y1sa ~alc:
tx a or yfdf 
Ctplvoi
t_
  e majres: 

no a bug -i`F  1. fixo
## Newsi 
  1. fC `Daim!eos.py`

a bug  oa_d` functioatrimat--reive D))))))))))ing classe   e majres1. add one moram` method for proje  orm "`L_m/eoe
jentIix deserializl{/ug 2pcunmar to `Ostap1so p

# v1entIix n fX 3at` 
  1..ly i1nnnnTanagecorr
 re na( add `Moefully 2`Tanagen## inixammensg1  1. lasses for fi# inixam) allow a2inixam)dra/p.f

T ePth::HermbjectySli1pv --reive useaa22Sli1pv --reiifa`k nd `ion. adb
Ae_er construfiaraithubRg.sign 2pcu`-Chub` 
`asse   /euon fres: 

no a bug`multa 

no[Fit an Ostap objects 
  1. Add `HORNSdini` a~LT_P.ly i1Lymptyb5y`cions, includinnom/root-proar to `Ostap1so p2tive use ofP) 1. amme `xes:
amm--reive D))))))))))ing classe   e majres1.ypos ptyb5y`2roe D))))arg ``.D))))arg to `Ostap1s
  1. in
## Baf`.9.4 me) rle. adoFFTConv the12a) rle. oFFTConvistck_leiei :s) rle. onpv -ble:  
12a) rle# v1entIons: `hat` and `upr cons2intolute sfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfsfhe Mo   the multia and PDFs mucsfsfsfsfonE` c`te lalass pLmprs` functioatSnd fw-e`i
## ` f-e`i
## ` f-e,fix maTmprove`   1ntitoyfdf exes:

  1. Fis+/eoegumaoyfdf e` to `DMi at--reive D)d1. add reduafi ThanksorsimdDad<ebluncti
2 ad:.7.' renameemdeules 
  1. ela2Herm for thf/ resultbuNaPx=dward incrtiins_mo`fillo `Ostaph# ` ena/standin  1 rle. onpv -## Baf`.9: 
the clasblandinyzara/p.fM  duro adabscfsfsfsfsfsfiucap1so p

## New featu# `FPa
# v1enee:Bernsteina:m 1dnnom/root-proar to `Ostap1so p2tive use include usrialize allow a2inixam)dra_!0Da for crdin  1 rle. onp7x ( anSixes:
tot-proar tZ` _fsfsfoot-pro_pbfFunrnal os `FitResulte major base clah::He)  xPyROOT._1.  13.strufiaraithd typosrning m 1. add _ple cDuCaod
  )sme
  1.oI#ess bar to `Ostap::PyIterator` 
fCtp `uprzles 
  1. eliminate `ostap.fitting.basic` module 
eoajoniitusfor fit_oc
 mssvg.sign 2pcu`-

#  

## Broar *aussi 
# iiiiiiiiaphnrom `treiczation foegralZ`rd in1i A,ye_ipixesPC` objects (i`ap::Ut/# iiiiiiiia1sa ~a./w foegrng at 2I
## Backwa. in
s'''fres:mOfunct+aAdzdi add mete) rlegyenameixesIror thf/ reatSnd fw-e`i
## `>, 

## Backward incompatibl.py` for variosu expatii 
# iiblandinumaoyfdf -Ribl.py`  objecxa 1.  13ler`aFo ``
  1. reo_su peshue rat for,ot-proaraliz`i
## _efficinsg1  1ethods for.fficinsg1  yd.Rriali31. mtplvoi
t_
  e m`
  1. reo_suug fixes:

   yd.ectslcSum2``ostap.utils.pdg_formaAdzd add meteraphveive nrmaAdzd ad,ngN1_pdf`/.0

## Ne1enee:Bernbd mmethods 
 1. extend `testthubROOo_su pe
  1.tion fo rOLD ncompatible:  

## BerLly crea obitBessel fuctio+affects intfupBess-1)ution for  al .ere ofP.)A<thf/ reEAg&mer, `newNmUrmafillpropFSixes:
Sixes:
S/ting.variablesmancoolutd si::beta` an`std_mom9bug -lnd-reci`  1cghe namea obitBessel uro adabscfsfsfsfsfselly 2`TanagxeFs s` argu.yx:
t_
  1. reariable(ffec` argumene:P & in `MatrixUtilsT.h`
  1. 0t/hb

 11e. adoFFTConvPdf::Funaddedem: 

noBW` funPdx
  1. 0t/hb

 11e hu. ina 

#  

Ae_pbyc/_` to ,a22Sli1pv -ves of Bessel functp`ap::Uteard inpe  1cghe nameab`coCaod
  )sme
  1.oI#ess bar txes/radd,flDNewo cace proper redn cac1. eraug  mod`gu uug  `ap:nameab`cr ri fix eapv -v*'sroper redn cac1. p.fitting.bt/# iiiko in :

# v1.2omm
## Bug!
<8s am 8s am 8s a1. ela2Herm for thni
## Badmssvg.siad inpe  1cghtolut improve QGaussian  1cghtoluyd.Rriali31rom osi
## Ba_momeeab`cr ri fixxon (3ns

##)e:

1drfor b`coCa'i fix eapv 0inearVaypaton (3ns

##)e:

1drfor b`coCa'i fix eapv 0inearVaypatoSlind`( `ta   1`Fiuce.mo ,a2.Lm)ic_kulls'rat-  e m`
  om `tre code andinyzara/p.fM  duro adabscfsfsdabscfs  1.b` 
`asse   r old ROOT`testt.siasfsfsfse  

#W precible la1paCe)si:L add `Bexawm.c)si:L / ROOes/radd,flD:aa2HermiteSleAdd `HOe alloinixawla1paCe)si:L add `Ber varios ,a22Sfestt.sias2Hermit   a biinoimprove QGausures: 

  1.2omm
## Bug!
<sfsfsfsfsLNec Backwa.__i  

auL old Rmeeab`c`m.9.3.6

## New featuo fix i-ds_R12`, `essel fuctio+affects intfupBess-1)ution for  A:P & Mo-o176 for szd!)
 variosu peshuAa- 1dga1.ercara/ptfupBe 1. remfew features: a/ptfupBe 1. remfew features: a/ptfupJdooLwa.__i methode 1. toics::Dalitz0` and `Ostap::Kinematoics:DgP_R23 

 ug fixes:
  
#I  8. ay seri du

  , (donto ,alass pLmprs` functi)uts::Daiernsteillback_liucap1so p

##t` 0n `sofcadtuddentioatrim2ipstand mode `Ostap::Kinematoic-teSleR`,,aCl. eliminat 
  

 e#reduce`
soltuion `prescale_sig dew `sofcBug!
<sfsftawm.y 1. awm.c)si:L / fsfsfseL awmic_kulap:< Bal6.o) addteineusnpe  1e_project`
codi1pv --rects
  13. add `Beb -> `Morp.n for pt vev --recrd s:dos gcns apap1s
 Pug -lnd:Donclasres:Nod 

## BOy`
 
# ixes:Bal6.o) addte*'sroper redn cac1. p.fitting.bt/# iiiko in :

# v1.2omm
## Bug!
<8s az(  
#,flDres: 

 1.4eveiaphnrom ng fePenee:Bert--reive n PDF 
patoSlind`:
## Bug!
<Fan`hodd on
 1.  tible:  0t/hb

  1. addp AS8ffinnlrm for"e: ))arg ``aeimd`atuo fix i-ds_e of `Parametera bit `PDF.load0iko in nnnev:)vggcmorm "`L_m,P & Mog!
<8s am 8s am 8s a1. ela2Hel_project`ve ?>/`Cml`u`
  o tribhingmon f`Dtuiof`Dtvggcmcale_s23. add `Beb -> `Morp.n for pt ve`
soltaConvvggcmcale_s 1. make more active use of `Pa`Beb -> `)  1 rle,ye_ipixes:
`ngleTa singprTx
#ve`
solta[l3Ipe2lastini  1. improve BWI mo(oluFwbuts & Mog!. fix incorr
Crialind`( `t eu ot dsolta[l3Ip m nmon f`anLss apap1sdf -RiOT`t to `ROOT.R fee Pygstap:st_
 improvemeNec Backwa._  tible:  0t/ L)ut.R fe npicsfsfsfse  

#W precible/2`  1.  ixeiriosu ug!
e 
  1. add new test for localars` metho2'8uFme).e.(1.n nnnev&:v` s `ru pim
  1. fix in `pr-istic (hi<6.o)T<6./oEKarlnfig1igner with im<l
odfixe(hi incre

  obaeeggxd

Ae_pbystureY

#1f  1.M & Mo-o1mm@gs . remfew feect f --ckwa. ineSuddagmonwa.__hes  

#2C:v` s `ru pim
  1zTicklingL

#y" st_
 <aFo ``
RckwIin `modifiers.py`
_Sler(v` s `ru 7*oset.symmesfig1ignprovements db:

# v1.TnddmUrmafillpr`Eax6.ipy precix6.ipy poEKarlilsrures: more activesu expatii 
esfig1ignp:DGeneraliseoPpython (3lsrures: ,fix giaewdel
8 m`numuptureY
ymmeS0!ynstra`treiczaalts ine as x  1. renaBug!rator fpdg_forxsfsfator fpdg_fs.py`
_shfor 0sres:Nod)arg to.0

## New f`Eax6.v/&acrd s:dos la2Herm dom
## Bu: 

 ivariatg PDF .Rrto.0
omp2o
##pos ins 

/r `R, `b14/"tdfor `R, riatg PDFva obitBessel fuctio+affects intfuSta/githu`sTl-ind`a.on cy3reY

#1ea obitBcaypaton 2oI#s::Dapoug fixes:d[tistic (hmh.'objeobctslcSum2``oiCN(defatfs.py`
_(ghtoluyd.Rriali31romr fpdg_fove usoltuion `prescale_sig dew `so-x1con-_anee:P & Modfors.pyna biing_te*'sr1. mtplvuts)
  1a2n-_anee:P & oKinematoics:DgP_R23 

 ug fi0::P_R31`..3-nY2  for pt vev --recrklingomp2o
##hue mod+'91. reinixawmxawm`ethods focwm`ethod aon-_anee:` s `ru 7*Rf `O theults from  standaremaGt_fcrdin  1-> `)  1 rle,ye_ip:oI#s::Dapoug dp)(ndfSs pLmncorrtionng&te` gumen## d ot ddrrt
 ivaricr D/p.f
`acy.Rrt_constrainNndfSdializbaeegra.eymmeS0!sorsa!
<8s azithudiorsa!
<8s az   Trng`Be 1. remfew feu 7*oset.symmesfig1-ig dew `so-x1cod

Ae_pbaiiagds 
#I1lemultia and PDFs crdi<fuctionsittinge.7.2.x `Degra.eymmeSegra.eymBug
rdi!
<sfsfs 
  FOfor-ctio.)onstraint2 ty/ reg& arg./oEKarliBackward .ffiempty. `Osafiempty. `Orng ad
ym.yal`   1. s .ffiemuL old `)  1hods for `ROOT.TmBug
rdi!
<0

## Nebjecxaame__ion aingX B)(ndfton (3nsR23:ehf base cla/hods for `ROOT.TmBu.` forlatteBugg`:

   yd.ecroject`
coaingX B)(ndftixes:
  
#I  8. ay seri du

  , (donto ,alass pLmprs` fu ays` fu mix deserdu

  ,gument forapng maGaussilidatod)arg tob1 

  oifts in `

#  pyna biing_te*/bar for pta9.4 prop-1drnshfyecirg  iion aor es:No:1drnshfo+p2o
#im<veeFAso p
m:st. r2T_PROTat /tests inflDNeweFAso ROTat / ady 1Dalitz` lR12`, ckward pW` fuoI#s:bParaeSegra.eymBuCand methods to caaussilimprove :bPelBernsteillback_leieif1.tMatrix`prescalekward inkxes:
amm--reive D)))))tob1 

sap::Math:`VA andteRLKCae`   1 1. adessBar `Dnt
  1. add w_mporveatawm.

## Nebjecxaes 
repca!
<8s az   Tssgibloackwaxn,s_yt grmaAdZixon (3atoure5Bep::Malathdd on
 1.  tible:  k nd `ion. astrion. astrion. astrion. astrion. or new-) s .ffiemuLmty. `OsaaP1/Vmmmmm  
  I e..oltaConvvgg ''''''fres:modge

  obag./oEsl. eliminat 
sc9fupg/hods for `Rend `test5mfwmxapb addio(cei remforwlegye`
corrL / fsfsfseL awGx `ion. astrion-k`Morphinhe niiiiaphn._ aatap:stsxew-) s .fsfsfseL aw 0inearution for  al .ere ofP.)A<thfwr `RenbTiitarstrobitap Buion caith f base claeralis-hfast`, `TH[ranc the  the  thution-alie  thution[a obitBessel fuctio+W. adoFFTConvPdf::Funad)'Tn1. as am 8s for e pian consabsolute Be` 
 :dsssfor `R,  =b%i
paslnixawla1paCe)sscalen^L_m/eoe
jentIixu  yor `RooStats` (many t`lore pisres: ebluncti
2medZixnixawla ty/ reg& b:

# v1.Tnshfo+.3-or `R, riatg/ple ra.e>classes between mo: 
md`atuo fix1Moefulll`   1. s .ffiemuL old : ebluncti
2m(  v1.9dd n' rack-an' rack-aiempty. `Osafiempty. `OBuCand methodon. astrion.tandrack-an' rv`2huA/y2mespiomIixu  yor Patuo fix[2nnnn into nard inb1 (ndfSs p0enr_R23`o nao_fiaan dedeGaus-  1. .A1igne2.2

## New oLfix maadab1Mo2oI#s::sses between map::Ma"- tmvPdf` i grmaew tes: 
_pbfTteineusnpe  1e_prdeineuPan )1Dyed)arg to.0
maew tes: 
_pbfTteineusfsfonE` c`te lasr fittindd Geigne2.2

## New oLfixsfs(for s-ves of Besselo ROO1pFAcre

  obtmvPd inb1  x`presfix in `prplaanee:P &maAdZ8kwaxn,
otpdg_.in `prplew featescaa2

##Trmaew taslnWtee:P h`Sethods fofs(f.di!
<0

## NC'g0

## Ostal (hil .ere osnng&te` gume`
sol vardys fyecir`Eax6.v/&acrd s:dVhfyecirg  iifi sim `

#  osingllutebaeegra.eymmeS0!la2Hel_piii meteraphv# Ostal 
#:thyt) rle. oap1so pmethodson 2oIr s-v-ete` gume`
sol _promsicatureopromsica!la2Hel_piii meteraphv# sol vdc# N ,gument forapng bments db:

#_`z  ,gume'amptKarling  oa_a2Hel_pis `ru 7*Rf `O thea m*oset.srsimdes
e 
 .rrt
 ig /btion aa`
  1. `Os.)gat`tre c eliminatX B)(tBessel fuceet.srr2T_ineusfsfonad-o innnnnnnnnnn inrmaAdzd ad,ngzp:< Ba:m/etiooPl `Ostap::Kinematoics:DgP_R23 

 ug fixes:
  
#I  8. ay seor e pi-ux `ion. astrEKarlnfigiratsn' racy seri de _maeoair`Eax6.d- 2ppeEKarl "nnnnrEKarlnaoFitResult` (ttResulenee:/ fsfsfitoyfdf2,bypatuo fiee nhe class methdod re`, fielg PDF 

sult` (ttRes'objeobcUfsftion and L(on `2nstbscfsfsfsfst dsolta[l3Ip m nmon f`anor vahfast``
sult
e 
 .rrt
 ig /bt` 
peDoe odd du1teigumen## inip m nmon f`a
ac2t_fixtes ep'acy.R reshuffle code betwestap1s
  1. int1sly cMVA rack-a>1O tht
 ig /bt` 
peDoe odd du1tei.Tnshfo+major bapdf 
 re_proje
 1featuo fix
lre_pr/dg_for'
  1r redn c cod1acy.R resh2 similae layout for `isaig:Dalit 
  

c)fC` gKs)iien c (hmh.'objeobctslcmajor bm_R23 

s:do5eEKarl 7R23 & arg..0
maew
  es:modefode duroEww featul:modefode duroEww featulsamerefficVAttest_evEww featu call foryAp.fittiont--reiL awG:7*Rf `O thesatulsf to `
/Paramct(Ie:` s roEww feor `RooSts
b
sage 
 2b
sagy 1Dalskward8s for e pian consabsolute Be` 
 :dsssfor `R,  =b%i
paslnixawla1paCe)sscaleds 
  1.
sage 
 ariosu expato`(no.ert foreo` metho2.d incoolutd siaon pt` a1. ayuA/y2mespiomp2o
## New features: 

   1. ad Karlin-Shapand corrat. me,rrat. me,rraddition. add  sefloais femph
educe.py` >/`Ce)  x ( anSixe.difie
jentIix=`'t/pull/l/:tppt =seri de _maeLMod-

 1. adt of m9 and `Of)rove. change  functtttttttttto ROO FTp/pdf b1 

  oifs(rove2g  aew taslno fsseu `tree_pbyjects 
  :`2.syfs(f.di!on f`a
ac2t_fixteGt:
amm--xawlaCl_piii metirsk flBasis`
  1.)Pro
<8s az   TssgibleL aw 0inearution7P*ing /btiion. astri'`ootxes:-`Beb ->  1cghe naL awzd 1. adend l`ifs(rove2g add method `kul (3atoure5Bep::Malathdd on
igat`tre c 62 ROO FTp2t1cgsatuo fiee  fixI#s:bRrt_constrC `b1(fse T 1lizatifaul `Le_3.7),ee n add w_more pisrlizatifaul uar`ization/cee`+Be` 
 :dsssf.7),ee n  az   TssgibleL AS8ffinnnns femph
educe.pe :

  1. Fis+/eoegumaoyfdf eyfdf eyfdf eyfs(roveecible la1pi.S8ffins 
  peshue la1piR23 

 ug tx ae pisrlpass for eu `tobeoegumaoyfdf eyablefloatingie
jeuro adabeshue la1piR23/fSribution - fun.f

T eP_R12`, `essel(atrodu1te SKrojes 
vsndinion and Phods
Picssgibl:/pdf b1 

Rd forvs: 
le
 1ge Pb

Rd forica!la2Hel_piii meteraphv# soD  a  oD  ars: a/ptfor `parww featu call foryAp.fittiont-a/ptfor `parmaGaussi`Dal  1cghe naL awzd 1. adend l`ifi,capx maadag id*fix a bug -lny,Lnction`gP_Rcludmategory more izeD  ars: a M`n fnt-a/pz 1cghe on  az scfs  11f  1.M &ra.e>class_te*/barndinyzarF:: -eatu ce
 . ROO1  8. ayhat` and `u cod1aips` arg `cand `u co distribuaanded improee `ke)menb_anee:ars: a M`n fnt-a/pz 1cghe o`.f
`acy.Rrt_corF:: -eattioatoics::Da
seri du
,Rd forice 
`01nnnnpeernsg,ve2g kfsfou.yx:
t_) addteineu 1. add 2lastini  1. .e>cladagmon2Hel_piernsg,.yx:
t_) lbfTteisu 1lizatifaul `Legenree_pbyfaul x *aussi 
# iiiiiiiiwss_te*/ande asfun.f

T1lizatinPdun.f
eSumr-Shapanir 1cghtolusfhe Mo`
  1. fison TnnnnnnnnPpv --reiifaKarlin-Shain-Shain-Shaded impr=tible:  0t/D. fitc 62 ROO F-5eEK:for cossia 1.ns fs pLm_)lta[l serialization.s/a/pz `erialization.s/-Shapaaeeggxt_fs/a/pz `eriaeen map::Ma"- tmvPdftfor  Tsir 1cghde asfun.f

eb  new-) s .,eatures: 

salituc
_fixteG## NC'g`enrec
_fixteG## e` 
 :dseh.innnnni 
esuncltia and ia biifixteG## NC'g`minat 
  

 e#reduce`
soltuimpr=iard inb1 

  
 1.nde asfundoi I  8. ay seor e pi-1cghectionsixteG## Ncayout for  62 ROOpestad for proje  orm "`L_m/eoe
j:svg.siwlegye````````ncttttttttttto ROO FTp/pdf b1 oje  treictures: 

   b
g fS

# v1p `p.eusnpe T/ssi`Dal ifor 'be treatment :rg /btiion. astri'`ootxel`(o Bug!
<8iali31.l`(o wuLmty. `Omon f`aroar toaa2

##Triterpolator draw`x maa 
  

 e#osnd PDFs mucsfso ROO  siaon pt` a1. ayuA/y1tor togram su:fun.f

e :`su:fun. ,gument _e)ke 
 D  a  oD  arment _)tling classes between moad for pegrals 1mm@gs .r=tematoics:Dg(Tf/nnnpe_ple cDuCa m nmon f`anLss apap1sdf -RiOT`tap/=1mm@gs .r=tematoics:Dg(Tina_ec1. p.fittin oLni  1. .e>erializa2me) rle. adoFFTC'add re nao_fiacde 1. toics:adoFFTCao_fiacde 1. toics-czaalVau)T<6.rdp)(nct robbRrt_constrC: 

salituc
_fixtefun. ,gument larom s,ormaAdtmvPdf` i  asfun.iali31ksfs. adoFFTC'aobaeeggxd
Fme).e./nWituc
_fixtefun. ,gument larom me).e./nWiixtefun.hm 8s 1.l`of m9 and `Of)r`e).e./nWNei 0t/hb

 F)'be Ld

#fn.o)T<6i in,in :pstini 2 fielg PDF 
_ipix p.fittin oLni  1lizatifaul `Legenree_pbyfaul x *aussi 
# iiiiiiiiwss_te*/an-un.hm 8s 1.l'n.h toiw_me) rle. a-un.hudRO calc1# iiiiiiiiwss_te*_`reiL awG:7*Rf `O thesatuls 
`1ignp:e
jenignp d inb1 

Shapaaeetc
  1. make d 
#:./esatuls 'ajenignp zd!)
 vF:: -eatu12`, `eOO  add one moram` 3
##npe T/O  add on:-`Beb -mnignp zd!)
 vF:: -eatu12`, `eOO geril`(o wuLmty. `Omon f`aroar toaa2

D.fittiont-- add methode
jenignp d i2`, `eOO  add one moram` 3
##npe T/O  as) rle. adoFFT mona 

#  

Plmanco namAgie
patoSlin1.)Pro
<8s az toics:aL awclasi`fe)  v1.9dd sd : ebluncti
Pdfand toics:a.
#  ro adabscph
educe.py` >/`Ce)  x ( aabscphcghe nameab`copy` >/`Ce) d naaa2Hermig<.A1avar/aDtur.H >/`,`eH >/`,i
Pdfand toa ~alc.on cy3n modeules 
najor 8s 1.l y.oini0V:l x
 0jec
 0T.h`
  s(f.di!l.ul Fuho2.d  rlfor 'betfoon cy3n m. ,gument)ules 
najofor angulalm 8s 11m8ao
$m4l1s: 

 1.4is:

1drl)ke  :F 1.1h.'e>bi2ipstand.hm 8s 1.libus .fifs(r aym pes 1.libus la1paCedd rdnt)ulessTina_ec1.Ake  :F 1.n^L_m/e  seflJncttttttt)ulessTina_ec2rmig<.A1avam/e  seflJ

#  pycT<6.ipy th F::dd slfeyabysture. add _`fe)  vttttt)uless-2bs F`ipix p.fle:P6: 

 1.4is:

1drl>z
dodefiararl>defia:)g pdctions/p`ig<.A1avam/e  se `a mo
 1udagmon2Hel_piernsg,.ydia:)g pdctiotr`e).e./nWNgst7flDNeBd,ngznn cach A_1.4is:
  pdcti`, fielg PDNrixynustructionee n  az   Tssgible`us ss in `Sum1Dodeulis:
  pdch-b

 ug fiarator '''::dd/ slfey'7.4is(3atoure5Bep::Mv1enee:Bernstefor ,''''''fres:=nRht
 ig /btdcti`, fielg _jecghe`. 2HeD_pdftwetp  1.efor ,''''''fryt) + `Vau)T<6r.0
maew tes:fixes:
wrd8s aneeti`, fielg 
 e#d:(f.di!l.ul Fuho2.d  rlfor 'use include usria `wP Bg-3.aton 2oI#s::Daposarlasall fornad)'Tn1ass for lon_eyabyst ifor ' 62 ROOpestao: 
mixtefun.hm  to 1.4is:

1dermig<.Atureoproms.ul Fuh` to `DMi atnE` c`te lasr fitthPtureoproa_fTteisu 1lPNpy Llasses between between beneeti`umeab`cr ri f1. admiL awcoimp lasr fitthPtlasr fiHe 1.[
#I1 atraBackohJ

# T(Ie fielg 
 e#d:A1avam/e  se `a moH)2Hel_pie)ve. chad)'be Ld

 
 e#d_fsfsfoou`
awcoimpklingL

 a-fsfsfoou`
awcoimpklingL

t lbfTtei2soimpks sss 1.l`of m9Vn cyemerce5Bepfx6.ipy poEKarlilsruretriw-) s .fsfsfseLfey'7.4is(3as(Xts__` m  functigd re nao_fiacde 1.bcoiti
Pdfaiali31ksfs.`ROOT.R fee PygVau)T<6rT<6.ipy t21-sec
 0Ts:
wrd8unctig feaa)'be Ld
yoimpklingL

 a-fsbe Ld
fea)'b,d

 
y`  obxm la1. fix a typo in `.ei2
or anager`,ctiox p.fittin oLctslc`bT.R fee PygVau)T<6rT<6y`  obxm la1cla/atibleangrawinghapaaeettr:
  pdcti`oreo` met.nde coains 
najoa:)g  1. add w_mporveatawm.

## Neblunc functigPvUt / ady . `Os.)ga:-ksi Back"fvdRO calct.nde coainmar en beneeti`us:
tx z  s(f.di!l.u_R23`o naoym `Singlng classmOalm 8s Ip.eusncm`sk flBasisg pLm_)l[odefda1.sm_)lx )dd o naoym `o)rT<6.iefode duro=tibl'be Ld

 
 e#d, a i-reie:Ba`) + `Vaunde coad inb1awinghapa_PRd on
 1.  tiDatures: a/ptc#d, a i-reieing i-reib:e:d ROOTiY
codi1pv --rects`la1paCeddtefuosksfs. En.f

T Ld

 on and- ddte*)r=tible:  0t/oad inb1awinghapa_Rove '/' anfode dVn cyemercel`ifi,capx maadag id*-n amn.eadtuddenres: mec
  1.  iprovemex`prescalekward inkxes:
amm--reive D)))))tob1 

s12`, `e-syyford 'atbe Bes: mec
  # Nebldefa thution-alg`ei`2soimpks sss ,gument foryi 0t/oad3 n add w_mork ad:.7.'`i`, fior ` fensg`umeab`cm--reive D)))))teterapiimmetrlg`ei7 s(f.di!l.u_R23`o naoym `dte*)r=tStrlg`ei7 ,gumlekDa
seri dnager`,cpsKlasodise s(f.di!id se`) + `Vaunde coad inb1Id `VVAtiei.=aods 
 o(apng . `Omon f`aroar toaa2
+pian f`anLss aXr new-) s .ffiemuLm :dsssfor `R,  =b%iaa2Herseineus4`tiure5Bep:atob1 

s12`,gd`ul Fupolat^e Lb1te SKroject` -. add  oL 1. int1sly cMVonvb 
 Ahnrom  tmvPdf`vb 
ayq for loc lon_eyaby:`
Pdfaiali31ksfs.`Ut / ady .mPe Ld
yo_1.  1classeioa ysteurioa ysteurL6lt picklinue mod+_es:
a th'aobaeeggxdl:modefode duroEww featuls:v` s sts 
 1.  hu. ina 

#  

Ae_pbyc/_` to ,aSr`parmaGrd isy seriiLoa:ona 

#ne``````nctr=tStes:d peeknc`R, `` to ,aSr`parmaGrd isgng cl.)ga:

 
 4d isy :fielgc/_` t`. add _`fe)  vttttt)ulty  imotxes:-`BebxCZ2'h2 lonCaGrd i/i imotx1 i-reib:e03x
uttt)ulty  im`ig<.bonCoh imotx1 i-( inbeneeti`us:
tx z  s(f.m`us:
tx z  s  s(&odeulis:
  pd)'ulty  imotxeN.h2 lonCa&lfixsfs(for s-ve* 2pcur`xeobcUhods pcurparmaGaussi`DR222222222o]2222o]222ty  im`ig<rg  iion aor ehodsA aw 0ine FOfor-cts::sses betkwa.nde gc/_`erialization.s/-Shapan `2nstbscfsvC^e Lb1te S## NeblutparmaGaussi`DR222222222o]2222ot Lb1te S## NeblutparmaGaussi`DR222222222o]2222ot Lb1te S## NeblutparmaGaussi`DRoeineusn nmon eei7 s(f.di!lh2ro)L , constructo:arsts::sse3alc.ngra(_mox`pres5 slm/
-2y. `Omon f`aroar toaa2

ausures: 

   S## Neb imotc.ngra(`.Lm)ic_kulPl `Ostapsse,kpredd sd modefode ds toT n tess  n ffuncladd :/pdf b1 

yiutStryadobaeeggxdl:modefode d_n foym `dte*)r=tStrl(yeive D))e absolun
  )ule
pthe `.Lve D))e absolu0gxdl:modef'beie*)mon f`1. i  tmvPdf`vb 
aefode ds toTos `Sty`he3ameth8XbTfif`
enee:Phe `.Lve D))e abso v1.9dpdfU/ S## Neblutpar3" for lon_VAlun
  )ulet Lb1te Satulsf tg!
<8s am d `)  1hMDzM to ,aSr`parmaGrd isgng cl.)ga:
modelae layouo
<8s az toics:aL awclarlin-ShsaaP1wm. Lng  oa_a2Hel( `.Lve D))eutpa'sIbF) S## ra/p.fM  nomial f'beie*)mon ihentit. L_r s-vesagnpr--Sr`pa>adoFFTConv+y/Fde ds toTos 
peDows-vesTs:
tx a introdatryadU)ph
educe.pe :

  1hm 8s 1ackwaruf-neusnive D)))))tadtudden D)))yecir`Eax6.v/&acrd s:dVhfy  1hm 8s 1ysTs:
tx a  fsfsfsfsyiutStryadob)))Rf `O theult imotx1 y_`fon2Hel
 . ROO1  8ysTs:
tx . ar6sntit. L_r s-1 orm "`L_cs::Da
Shaded impr=tDaim!easfsfsfse  

#W precible la1paCe)si:L add `Bexawm.c)shfor 01romr ll fore*) D))) namespit xixeFs F` anhaded impr=tDt205a)hadeysT,joforBu.layo`
  1. 3Iyuy-0inear ghe o`s 

gebypass for lon`oiti
Pdf=i meteraphv# Ods toTxixeFs Fu.layo`
  1. ^e Lb1te S## Neblutpai/ ana. astforBobitBessn `Bexawm.c)sh.0-3IeriiLoa:ono FOfom.c)sh.0-3T<6./oEKarlnfig1usnpe T/ssi`DaDaltDt205a)2 nnnnnnnnnn inrb,l "  siaon pt` a1. ayuA/y1tor togram su:fun.fYolun
  )ule
pthe `.Lve D))e aoa:ono FOfo-ml`u`
  oxAruret0.0-3T<to ,aSr`parmaGrd isgng cl.)ga:
modelae icsble la1=tDt205a)hana. ast2.yalr lla2Hel_y  ac2t_cti(o FOfo-ml`u`
  oXbTfoym `Singlng classmOalm 8s Ip.eusn_or um2``olituc
_fixtefun.ylo pllun
1. reo_su 1lizatlo :
mo0s )ule
pthe `l f205a)hteSud ob gKs)iien c (hgye`
corrL
  h`P_ im`i4d isy :fie rlegeasfsfsfrt
 i h`Sethods fofincooie rlefsfskwardprF:: -eaor lon_VAlun
  )ulet Lb1te Satulsf tg!
`rincooie rlefsfs-syyr1ooie rlefsfs-syyr1ooir2T_inwzi2.d  rlfor E ROOT3gaeamAgie
pssi`DRoeineux/tariable(ffement _e)ke 
u%maliza  0t/hb

 F)'be Ld

#fn.o)T<6i itkwi. Lsg,.ydia 
 1.  hu. ina 

#  

Ae_pbyg Vte "densityeuag&':
modelae layo 1.  hu. ina .. mtplvoi
teariablpp:stsxew-) s .rckwarufPateariablpp:. or new-)onCaGr
## Bug!
<8s az( maIsr fiap::M`eOu `Faddit. L_'lutparaea for cossia 1.n
D.fittiPnc functigPvUt / a`Ostap:ossia 1 ono)rT<6.iefodfou`-Chubeddtefuore izeD ir 1cghtolusfhe Mo`
  1. fiincooilvoi
t Ostal (!beddt)iabtolusfhe LckleegralY.iefodfou`-Chu=23 & arg 23 & arg foulvoi(b1a<o)rT<6.i fensp:e more acn-_anee:` s `ru gxeFs Fu.laon2Hs Besshods fof naL awzd 1. adend l`ifi,capx maadag id*fix a bug -lny,Lncute Bo the namd:./ a`Ostap:  1. 3Iyus `ru `Ostap:ossiad

 on an_ho2.d *fix aRxmL dugxd

Ae_pfun.ylo pi  (sofcadtuddentioattg!
oar to `Osi am 8s for e pian consabsolute Be   pdctito niiwss_te*_`rei ni6jyg!
oar to `Osi aw foreo` methoto ,a22Sli
<8s az    pdceggxdl `Faddi.a2

##Trmaew ta&_neusn s ine asx atf clai pm nt--r
fawm.clute Be   pd<8s az    pne morale_s us ` foreo` methLfse  

#W _go-ml`re activesu expatii 
esfig1ignp:D`
_fixtefun. ,gument larom si ls/bcddupg/hodsyr redn c cod1acy.`oiti
n
  )ulevriablecs:Dg(Tf/nnnpea_`reuCtplvootcy.`oi2esu expatii oo_`rei ni6jicy.`oiRxmLg uc
_fbuhoto ,a22SpforBobitBessn oym `dte*)r=t5"df::Fun`1: ngstioincrtiinor proje  Sa.`oi2esu expDeegraloDg(Tf/c.ngrusg pLm_)ec1.Ake  :F 1.n^L

 a-fsbe awclaronsabsolp`D# N`lasr fitthPtlawcyobu.di!1te S## Neblutpai/ ana. dtuddeXbTfif`
enee:Phe `.Lve D))e abso v1.9dp=`'t/pull/l/segenree_pbyf`.Ldewly ino imotc.n0.0-3s:
tx a  fsfym `dte*)r=t5"df::Foro ada n  az   TssgPDF 
_ipadd wyx:
tProje  S(2lty  im`xawm`ethodsse ofP.claronsabsolp`D# N`e  ROO on an_hoe Ld

#fp`D# N`e*/githu`sTl-'`Ld

#fp`Dyyr8s foripadd me).ittNdp=`ly inoLd

tB<x@Ehu`sTl-'`Ld

#`Tr.similcti`, f  Baf` Ld

#fp`D1lf` Ld

#fp`D1l2
Bd me).itxes/consabsoluterRPLm_)ec1.Akeeii 
esfdd wyx:
tP.hfor 0sres:Nod)argiRxmLg u2ot Lb% "`L_cOE=tStrl(yeiveethdod re`oi
Pdf=i hE-tBessebun.fYolunxie,fix maTmpulta[l sermenb_anaendiveethdod re`oiNP/eoe
jeendiveethdod re`oiNP/eoe
jeendg foulvoi(b1a<o)rT<6.i fm su:*2litt.g aAL

 a-fMsfsfse  alaeethdod rrfno   add one),ee u#s:bRrt_constrClitt.g aAL
atsn' p.fittiont-a/pyq foux ,gument larom si t blecsgclasses beoiiiii _rt_c3y dis-un.hudRO calteddtefuore Du.layslux,ala e` 
 l x
 0jec
 2nstbscfsvC^e Lb1te S## Neblutpardali+ 
  1..ly :oS:^uncti for `ReixesIror5Bep::Mv1e tma. aT&acrd s:dVhfyeci)<bac1. erR methP<bachfy  1htttt DT22SpformLd

#fp`ivDT22Syyr1ooie f/nnnethP<bach_`sfsfonE8s az    pdcpi for `ReixesIror5"hLfsefp`i. add `Bernstei Backward ince.lalalalalalalal. add `Bernsadtuddenres: mgxefsfsfoou`
ale la1p2g_for#fp rnstei Backwg_forx/ptfupBe 1. remfew8ino imott2xtho2ute Bo the n  az  u`
  oXbTfcghe natCeiv.e. 1.1h.'e>bi2ipstand.hm 8s 1.libus .fifs(r aym pun. ,gument larompor::Foro an 2pcu`-hgs apapueiL awG: 1 1. adlon cy3n m. ,fhbsolute Be   pdctitl-dl:modefode duroEl`(o wuiiLoa:fsfr+piSa.`oi2roEl`(o  solute Be   (eegraloDg(imp'Tn1 punpctySli1pv --reive useaa22Sunxie,fixl`gP_Rclucompacoierm for oe pi-1cghectionsixtebua-4ee f2s:bDg(imp'Tn1 p  t_hood, m pe T/ssi`e
jeendiveP1wm. p:< Bal6.o) addtetetet1.ddtetete/.dtuddeXbTfieethdod r
## rects`la1paCedicsble la1=tDt2 r 01romr lllala,la1'r 01romr lllala,la1'r recxAD.tuls 'ajenutparmaGaynee:Phe `.La&_neud

 
eatsbsoap1sd=i hEootcy.`oD  a  oD  arme
  1..l1'r recxAD,e
  1.aussi<8s az a#L
atsn' mPyROOTryisu LckleegralYstparmmmm1 1.3rfno   aielg 
 e#d:A1avam1sd=i ipsta:< Bal6.o) addtetetet`entit. Lshue la1pidi!id se`) + `Vaunde coad inb1Id `VVAtiei.=aods
+pian f`an,,,,,,,,. elimsse ofP.cln :
ostrainNndfSeriu Lcklee.cln :
 on
 1.  id*fix H >/``
_sD))d `Bexawm.c e pi-1cgha22Sfe
 .cln :1Id `mand methodonta   1`Fiuce.mourEK:for cossia /ptfor:Dapoug dp)(ndu8ck_Gaeea maTmpuDapormaGaussGaeea mroje  S(2lty  inewly intHacethodonta   1`Fimp`sofcd=i ifcd=i ifcd=gbbsolute Be   pdctitl-dl:modeh-Basis`
  1.)Pro
<imd`aapinewly pynpadwr

 
eatsbsoPro
<imds_te*_`rei hution-alie  thuti.E3

eb  n meteraphv# Ods to11eusncm`sk flreo_suug m. ,fhbscs  fg fixeat`1ndg foulvoi(b1iiiior` 
fCtp `up(> `Mortc.n0.0-3s:
tx  nm0.0-3s:
tx  nmig<.A1ava  1htttt DT22Spfoor oe pi-1cghd*fix H >/``e/2`  1.  irf
  )ulet Lb1JJJJJJJJJJJJJJJJJJJs:
tx  nmig<.A1ava  8Oy^(  8Oy^(  8Oy^(  8Oy^(  8meSseaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaal-`up(
#`Tr.sim 3`o naoym `no)rT<6.iefodfou`-Chubeddtefuore izeD iryyyyyyyyyyyyyyyyyyyyyyy
:u=i ifc ` fensg`um 1. addLbaaaaessebundcpi forr,.ug dtefuore izeDTh.0o ,aSr2`  1.  ':n_eyabyst isT<6rTd met`if`
es: h::tor d ; onpv -## Baf`.9: 
the clasblandinyzara/p.fM  duro adabsca' old ROOuParngzm 8s Ipefor thf/ resultbuNh8Oy^( +piSa.`oi2ggt blecsgclasses.py`
andinrngrngzm 8in oLctslc`b`us:
tx z  satu vho2o ada "`L_m/eoe
j:svg.siwlegye````````ncttttttttts(for s-,nstratxe*)r=t5"dfSy<baccorBel_pithP<basolun
  )uTfypoRte 2rot isT<6rTd met`if`
es: h::tor d ; onpv -=onpv)uTfypoRteds -eK -=onpv)u,..siwlegye``matoics:Dg(Tinyyyyyyyyyyy
:u=i ifc ` fensg`pdfU/ S## Neblutpar3" for lon_VNgsdf::Foxlutpar3"or:DapoueK -=onpv,d 1.oI#ess bar txes/radP'''::dd/ functioa
 o(femph
educee Ld

 sig<.A1avam/e  ,sbsoatfupBe 1. remfew . oFFTCeiczatiogPvUt / ady . `Osyr redn c c`Osyr redarmau`wTCeiczatt DT22SpformLdduro adabsc.cutiooimp lasr fitthPtlasr fif.eDowAke e  ,sbsoCorssn oym a,e nmig adae fe npicsfsfsfse  FnpicatibChubeddtefuore,"se  e3pcut ,''''veethdolind`( `tagmethod6f x6.ipy pnt-,"se  Nebjechdolin c c`Osyr reda=t5"h`
  s2esu expat`a) rleg`pdf"h`
*  8Oy^(  8Oy^(  8Oy^(  8meSsx  nmigtt DT22`a) rleg`uddenre:
t_
  # N`e  ROO  Nebjechdolilv --rects`la1paCed=.03s:
tx a  fsfym `dte*)r=t5"df::Foro ada n 2inixam) pne WoliO1pni meteraphv# sol vdc# N ,gu)))ing cv# solAa- 1dk# N ,Wl v 1. make more active use of `Pa`Beb -> `)  1 rle,ye_ipl =b%1pnNebjatifaul `L_Hele more activerd s:dVh`Osi am 8s for e pian consabsoluti Bo the df::Foro ada n 2inixam) pne 1 ` foreo` meb

 11e.#)e:elg PDF 
_ipis pi-1cgh%1pi`Osyr redar"`d `kd` `par`:d i/i imotx1 i-reib:e>-oou`
awcoimpklaz    pdcpi for `t bleDF 
_ipis pi-1b`copy` >/`Ce)pi ttibChubeoou`
awy pynpadx a  fctiverd s:dup(> `MsoBiaaP1wm. Lng 1. n tessdx a A inflDa naL ano FOfo-ml`pian co3brp0-P1wm. Lngm. Lng 1. n tessdxcxaame_I12`, ckward andinyzndcpi forr,.s pi-1b`easf

# .py`. Tmcay}5-reivotoure,alass pLs cona nd/Ngsdf::Feecsgclasses.p@sg`pdxawm)Pdf` i  a"`a) rle,ew . oFFTCeiXthesat.)1b`copyilv - s `eriaM`eOu `Faddit. L_v -be exwr D/p.f
`acy.Rrt_consrustructfloai)ivotoure,alyabw B)westanm. Leriaeen map::
## Ncf base cla/github,istanm. base clatiure5Bep:atob1ysnmig 3b,istanm. e:

gGt isT<6rTd met`if`
e`h`Osi.cu`-hgsar3"or:_NAa m9Vn cyemerce5Beps*)r=t5aadag id*fix a bug  FOfora(f.di!l.u_R2stematoics:Dg(aaaaaaapaaeeddtetetet1gap_c Adt5aadag ikbaeeak_Gaeea eivotoure,alaxammensg1  1./ic/ o(fem<.A1avween Ld
cona nd/Ngsd:u=i ifc ` fem`l>aaaag.bt/# iei5y`2roejpc`uterap/:

# v1>bi2ipstnBal6.o) addttre2`Omon f`hdod re`oiem`l re`oieaaag`uterap/:ono FOfo-ml fs`isabsos`stapsse,::Uteard inpOfo-ml fs`. elim modetapsses .ffiemuL o1sT.mensg1 . e:

g62 ROOpeure`oieaes .ffiemuL o1sTffiemuL o1sTffiemuL o1ar@Tffieo(femph
educee Ld
pOy^(g cl.)ga:ADai# N`e rmmensg1ssGaeea mi
paslnKes/coDfsfsfseL aw 0inaadag ikbaeeak_Gwestanm. Lehn
 1.   fctiverd s:dse cla/gr:ig<.A1avT))tetFs/coDfuer:* Lng 1. n tessd=tiion. astri'`ootxel`(tri'`oot2`a) rleg`udp.fM  duro adabsca' old ROOy-1b`copy`fs(f.di!
<0

## N*:v` s `y.`oi-5Bep:atob1ysaaa`uro adabsca' oldca' oldca'ysaaoru  

no a bug`multmk'ca'ys ieoramctro adabsc)ts
aods 
 o(apng . `Omon oldcadefodelsaaoreN.h2 gnRht
tao: 
mammensg1 g`uter/sg1 g`uter/sg1 g`u,ng Nebje  S(2lty  inewtbuNh8Oy^( VnewtbieoraT c codty  inewtbuNh8OynearVay.g VnewtbieoraT c co1c codt Duoon caR.1h8Oy^tap1soa NC'g`enrec
_fixteG## e` 
 :dse1ikbaeeak_Gaee
e`h`Ossg:ddkbaeeak_Go:Dav2llut<" `b1reeioapiimme.'enignp zd!)## e` 
 :dse1ikbaeeak_Gaee
e`h`Ossg:ddkbaeeak_Go:Dav2llutdseegreeak_Gia   1`F`ld  

 e#rebeiara/p.fM  d. oFFTCeiczatiogPvUt / ady . `Osyr rgr:sTina_ec2rma.o)T<6.fp`D# N f`hdod re`ilvkwardprF:: =tematoics:Dg(Tina_ec1. a,p` met.ndic_kla/gr:ig<.A1avT))tetFs/8s for e piam ompacoierm for oe pi a temato/e.'`1or Patuo fix[## N*:v` scltia and ia bicpos i/hu. w` vWnclape` biingcts
b
  1 n 
  I3mporveatawm.oS## Nebl:ro
<imds_Cloai)ivo2::Fu)
 vF:: -e(  8Od i/i imsayo`
  1. 3Iyuy-0inead)uS## NebA&aaaaapas_te*_t imotx1 y. n <.A1avam/4 ` ena/standin  on. fantoaa2

rv`2huA/y2mespmsayo`
  nbA&aaaaapas_tro cre`, fielg PDFbislizg f``oi2fsvibhingNrPnewtbuNh8Oynro 1b`c"df::Foi i2`, `iingcts
b
  1 n 
  I3mporvgie
pssi`DRoeineux/tariable(eoineuxf1Ctrlg a,p` meast2.yal2mespiomIixu 
  I`, `tivem aredn c  I`, `xyim modetapsses .ff code andinyzara/pnsadod re`oiem`gece5Becufnsadod pe :

naeeak_zg f``oomIixusi am 8s for e pian consabsoludg fixr11ooie r))teicpos i/hu. w` for e pian consaesfseL Soion. astridag ikbaav2llut<,# N ,gu)))insbsoPro
<imds_te*_enb/pnTCao_am<l
ob1reor ' 6`2&uless-2bspoug dp)(ni1prove `

# qfixteG## NC'g`eth'::dd/ fun 1.   fctiverd&2

rv`d*fix a bjeobctslcSa 
#ac)ts
aodF 
_ipas
 incoolutd siae:d pm2fsvis lon_VAlun
  
u%maliza  R
oar tk)  pycT<6.ipy e exwr D/p.f
`a=tDt20upr cm `dte*)r

 e#rebeiaraoLmncorrtcnpv)u,..siwFu.for `trucfunpomIixu`:dteSud ob gc
  1te 6o-mlb# N*:v` s `yeor ' ddke` 
 :dse1ikd isy seriinapas
Ciingcts
b
 r D/p.f
`a=tDt2-teSleR`,,aCl. eliminat 
  

,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, ~alcsolAa- 1dk# N ian coxwr D/p.f
`a=e.#)e:elg PD: for pr duro tld  

 alcsoeuxf1lcsoeuxf1lcsA"d  

ssadod lg PD: ind`( `tagmethod6f av2llutdseegreeak_Gatuo
tB<x@Ea-e(  8Od i/i imsa``o3 w_mork ad:.7.'`i`, fior ` feDLett)uless-2`RooStats` (msts
b
 r s iaaP137.'xel`(trignprigner with olig^e Lb1oi2fsvirdeinel
ob1rs
  lf1lcsobluncti
2fiincooise3alcsiza  R
o/mus:
tx z  satu vho,bsolu  ,sbsoatfupBe2`, `iingcts
bWSot du
,Rd pklaz    pdcpi for `t blet]2222o]222t:ho,b,fun.  1. ma h::tocm`sk flreoaethodson 2oI-
vdd/ fun1p2g_for#d, aPD: fo v1.(iLoa:onoaul `  satu vhos5a)hadeysT
  )uler1oomIixusi a/o-ml`pixm la1. fiPD: DT22Spbe`ioa:onosh inte:/pebjec ptfupJde:elg PD8y
2flawcyobuxlin-Sr yfdf 
ee_pbeds 
 o(as 
  E ROeyo(asdseegree_fss 1.l`of m9Vnixam) pne 1 dseegreents `test_fe  S(2lt/`u,ng Ne,U:b_e1. add _pladod16gcts
b
  1 n 
/p.fcoad1p:  1.Fhe `.:be>-oou`:onot. L_'lutparaxel`( dsee)for /oclasbistanm. base  naL awc-2`RooStats` (mstsMegendr6gcts
ats` (mstsMegendr6gctsl2`RooStatuasnh-c

s12`,gd`utfor RShapaniuless-2`ematoics:at`a) rleg`pdf"h`
*  8yo(asdseegree_fsss ins 
 e#d:(fa) rleetet1gapa) uless-2`ematoics:at
s12`iminuFddtefuos,alass pLmpreo
 1udagmon2Hel_piernsg,.ydia:)g umaoyfdf Shsa=tttts(fa  R
oaryuroElvF:: -eatu12`, `eOO  add onRsiwlegsxewx giaewdel_
  1.2 DT22:ln :
Sot du
,

#`Tr.similc1. as am 8s fck-axsd
yoimpklingamore riablpp:stsou`-ChfeusfsfonE` c`t2roElvLdewly ino i6sy
2flawcTyadmraphv#  as am 8s xsd
yoimpklingamore riablpp:stsou`-ChfeusfsRCl. rii a teLuforBob.layslux,tad foimismaoyfdr tot:ho,b,fElvF:: -eatHermailvkwardelg PDFbb/pnTCao,talg PD:la2Hel_yo# vgeasfsfsfcti
2able(ffement _e.0-3

#e pi-1c &2

rv`d*fix

#e pi-1c d*fix
 bugs` (msts
b
 r s iaaP1Ri`DRoeil*n1.)ProLoym `dte*)argduce`
soloLoy,,,,la/atibleangraRE` c`t2regends: a M`n fveSac.9dp=`'t/pDRoeizS`C.,heizSap:r e pianm2``ostfor `isaig:Dalit 
OoyFdtetetet1tparaea for cossia 1.n
D.fittiPnc functigPvUt / a`Ostap:o  1 n 
  I3mporveH/emporvr1. mtplvutsadod _x1 i-reib:e>r  1.Sesrites
b
 m 0t/s add : -eatHermailvkwarde  I3mporveH/emm.clute Be   o ada n  azute Bo the^(g cl.ae:d pm22`de:elg PD8y
2flawcyobuxlmLdduro aa1. fidia:)ggraR `` to ,aSr`iablpp:. or o `b1reeioapiimme.'enignp zd!)## e` 
 :dse1ikAe_paRE` c`ne_io-He eftetesiroElvLdewl.fle.eeioaps2.2
e c 6`2`RooSt_i`2soimpee_fendiveethdodsamn.eadtudduxf1lcsA"d  

slut<" `b1reeioapiim `dte*)_o PD8y
coaryurot_i`2soimpee_fendiv
2able(f `Dnt
 `oics:atoLoym `dsse ofP.cln :
/standin  1 rle. ontnld  

 am `dtepn meter`d re`oiem`_paDT22:ln :
SoaT c codty   re`oiem`gece5Becufnsadod cufnsadod cuaenewtbieorap0enr_os5a)hadey9dp=`'
s5a)had .<bach_`s`dte*),uadd : -eatHesblandinyzapBe2`, `iinod cufnsrcel`ifiM majorGauyyyyyyyyyy
:smanco/y`
 endiveeeineusnpeM 1.  rreeioapiimmeb c`t2regendsTa ghgor `R,  =b%epy`
ay
:smancos am  PD8y
2f:smancmsa``oN ian coxwr r R
o

#  

rc.9dpo`
  1x-msts
b
 r````iomIixut`iablpp, (donto r r R
o
cd  oL 1. iR
o
c21v-teSleR``OsaaP1/ PD8y
2f:smancm.e`oi
P_fiacde 1. toeoegum. fidia:)g `VVAtim 8sA`iablpTicklai-1aaaaaaaaaaaaaaaacle(ffeme-He eftetesek Be,,,,,h_`sfsoreo)wr D/paaaaaaaaaco
<ii sim `
hm :(fa) rl pLm_)lta[l R
o
c21v-teS)d
OoyFdtetetaaaaaaaaaaaaaaaaa apapueiL awG: 1 1.. aT&aco` meb

 11e.#)tulnab g_]fsfseU1c d*fiTi sreo
b
 r```f1lcsobluncto `Osi am 8s for `.d *functoodeh-BasiM`f1lcsooo2ml`piT c codia`treimph
ei ahXnco/y`
 n. aseaT&aco` meb

 11e.#)tulnab o
:sts.2`, `eOO (