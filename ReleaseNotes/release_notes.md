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
  1. a few minor improvements 
  1. improvements for `Dalitz0/Dalitz/DaltzIntegrator`
 
## Backward incompatible:  
 
## Bug fixes:
  1. Bug fix in `Dalitz0::P1_R31`
  1. `DalitzIntegrator` fix incorrect use of `std::enable_if`


# v1.7.1.8

## New features: 

   1. Update setters in `signals.py`
   1. add `Ostap::Models::Histo1D`, `Ostap::Models::Histo2D` and `Ostap::Models::Histo3D`
   1. Update `Shape1D_pdf`, `Shape2D_pdf` and `Shape3D_pdf` 

## Backward incompatible:  

 
## Bug fixes:

# v1.7.1.6

## New features: 

  1. add printout for added varibales into `PDF.sPlot` method 

## Backward incompatible:  
 
## Bug fixes:

  1. fix some old typos in `make_bkg` function from `ostap/fitting/background.py` 
  1. tiny fix in `dataset.table` method 
  1. some (minor) fixes for `root_file` and `rootshelve` modules  

# v1.7.1.4

## New features: 

  1. extend example   `test_fitting_simfit5` to use binned distributions 
  1. add dedicated test for `Shape*D_pdf`

## Backward incompatible:  
 
## Bug fixes:

  1. fix `Shape2D_pdf` and `Shape3D_pdf` (thanks to Slava Matiunin)
  1. more fixes for  `Shape2D_pdf` and `Shape3D_pdf`
  
# v1.7.1.2

## New features: 

  1. improve progress bar for data frames 
  1. add warning message for drawing in ranegs for the  `RooHistPdf` with ROOT<6.24 
  1. add one more example/test for simultaneous fit machinery:  `test_fitting_simfit5`  

## Backward incompatible:  
 
## Bug fixes:

  1.  minor bug fixes for ``no-numpy'' regime 
  1.  fix typo in `combine_data` in `simfit.py` (thanks to Abdul-Kerim Gusseinov)

# v1.7.1.0

## New features: 
  1. add methods `size_vars` and `array_vars` for `ROOT.TTree`
  1. small adjustment for `tree_reduce`
  1. remove warnings for `splot` 
  1. large improvements for tmva&chopping: allow ot specify working directory, etc...  
  1. allow specification of training signal&background fractions for tmva&chopping 
  1. disable ROC curved and ROC/AUX for old versions 

## Backward incompatible:  
 
## Bug fixes:

# v1.7.0.6

## New features: 

## Backward incompatible:  
 
## Bug fixes:

# v1.7.0.4

## New features:

 1. extend test for DataFrsmes 
 1. StatEntity and WStatEntity: ignore non-finite values/weights  

## Backward incompatible:  
 
 1. for ROOT>6.28 remove treatment of `RooFit.FitOptions`

## Bug fixes:


# v1.7.0.2

## New features:

 1. add (much) more efficient way to fill `RooDataSet` (activated for ROOT>=6.26)
 1. update examples & tests
 1. improve code in `pyselectors.py` 
 1. temporary disable call for `cmt_fit` form `cmp_diff_prnt` for vey fresh version of ROOT 
 1. extend to ancient version of ROOT 6/14 (except for paralllel chopping) 

## Backward incompatible:  
 
## Bug fixes:

# v1.7.0.1

## New features:

  1. add option for parallel file copy 
  1. add argument `parallel=False` for `Files.copy_files` 
  1. add `roo_cuts` argument for `SelectorWithCuts` and `make_dataset`
  1. add new test for counters 
  1. imporve numerical precision for counters 
  1. add `TChain.parallel_fill` (new name for `TChain.pprocess`)
  1. add `TChain.fill_xfdataset1` 
  1. add `TChain.fill_dataset2` 
  1. add `TChain.fill_dataset` 
  
## Backward incompatible  
 
## Bug fixes:
  1. bug fix in `WStatEntity::add` for initialisation shortcut
 
# v1.6.9.9

## New features:

 1. add more extensive test/examples `test_fitting_fill.py`
 1. several tweaks for parallelisation 

## Backward incompatible  

## Bug fixes:

  1. fix new test for `DILL_PY3_issue`

# v1.6.9.7

## New features:

  1. Improve `CleanUp` to clean/remove PID-dependent temporary directories created in subprocesses 
  1. add `job_chunk` argument for `TTree/pprocess` to control length of job chunks, optimal value about 3*ncpu 
 
## Backward incompatible  

## Bug fixes:

  1. few tiny fixes for `parallel_fill`

# v1.6.9.5

## New features:

## Backward incompatible  

## Bug fixes:
  
  1. (re)fix the bug for the operators in `ostap.trees.data_utils`

# v1.6.9.4

## New features:

  1. add `commonpath` and `copy_file` utilities to `ostap.utils.basic`
  1. add `commonpath` and `copy_files` methods for `ostap.trees.data_utils.Files`
  1. imporve `copy_files` in  `ostap.trees.data_utils.Files`
  1. add `copy_with_progress` function to `ostap.utils.utils`

## Backward incompatible  

## Bug fixes:

  1. fix bugs for the operators in `ostap.trees.data_utils.Files`


# v1.6.9.2


## New features:

  1. add an easy `dataset -> TTree` cnversion function 
  1. add much better tretament of finite differences for numerical derivatives 
  1. allow forwad and bacjeard rules for finite differences
  1. allow specification of singular points for numerical derivatives   
  1. allow numerical derivatoiebd up to order 6 
  1. add test for derivatived for function with discontinuous derivatie 

## Backward incompatible  

## Bug fixes:

# v1.6.9.0

## New features:
  1. add  `Ostap::Math::CnannelNRL` for non-relativistic Breit-Wigner channels 
  1. more imporvements for the math primitives
  1. even more imporvements for the math primitives
  1. add python operators for primitives
  
## Backward incompatible  
  1. rename `Ostap::Math::ChannelNR` to `Ostap::Math::ChannrlNR3`

## Bug fixes:

# v1.6.8.2

## New features:

  1. `ROOT.TH1.smear` add argument `silent = True` 
  1. `parse_args` improve logging (use local loggers instead of the global  one)
  
## Backward incompatible  

## Bug fixes:


# v1.6.8.0

## New features:

  1. `interpolate.py` add `bspline_interpolate` based on `scipy`
  
## Backward incompatible  

## Bug fixes:


# v1.6.7.1

## New features:

 1. Change signature of `PS23_pdf` - now it requires a valid Dalitz  configuration 

## Backward incompatible  

## Bug fixes:
  
 1. Fix bug in `PhaseSpace23L::integral` and `PSDalitz::integral` methods 

# v1.6.7.0 

## New features:
  1. add `Ostap::Kinematics::Dalitz0::E1,E2,E3` and `Ostap::Kinematics::Dalitz0::P1,P2,P3` three-argument methods, They are less efficient than corresponding two-argument methods of class `Ostap:Kinematics::Dalitz`
  1. Fix `Ostap::Math::PSDalitz`
  1. Upgrade `Ostap::Math::PhaseSpace23L`
  1. small tweaks for `Ostap::Kinematics::Dalitz` objects 
  1. add serisalization for `Ostap::Kinematics::Dalitz` objects 

## Backward incompatible  

## Bug fixes:

   1. Fix `Ostap::Kinematics::phasepace3` for cases with some arguments are zero
   1. Fix treatment  of `xmin/xmax/...` for `f1_draw` and related methods  
 
# v1.6.6.0

## New features:

 1. add new test/example 'test_fitting_resolution3.py' with relativerly realistic exmapel of simultaneous firtrting of "data" and "MC", propagating uncertainuty from MC resoltuoion shape to results of fit to dat ausing simultaneous fit 
 1. add `array.array` and `numpy/.ndarray` into list of `listlike_types` for `ostap.core.ostapo_types` module 
 1. make use of `fudge` argument for `test_fitting_resolution3.py`
 1. add keys `remove` (default is `True`) and `keep` (defautl is `False`) for temporary databases. The first one forces immediate rmoeval of the file (instead of the end-of-the-task action), the second forces temproary file to be non-deleted  
 1. add more printout for the `ostap.utils.cleanup` module 
 1. add more constructors for `Ostap::Math::Interpolation::Table`
 1. disable some serialisaton tests for ROOT<6 and python3 < 3.7 (seg fault)
 1. implement serialization/deserialisation for matrtices&vectors 
 1. improve interpoaltion stuff, maknig it more efficient + extend tests  

## Backward incompatible  

## Bug fixes:

 1. fix the bug  in `ostap.plotting.fit_draw` for parsing of drawing styles/options 
 1. fix minor bug with ordering 

# v1.6.5.0

## New features:

 1. add `reduce` method for polynomial and spline classes 
 1. largerly rewrite and extend all interpolation stuff 
 1. more improvemements for serialization of polynomial-like stuff 
 1. add generic python interpolators with tests for matrices
 1. more agressively decorate new instances for linear algebr aclasses 
 1. add method `shoot` for 1D, 2D and 3D histograms 

## Backward incompatible 
 
## Bug fixes:

 1. fix bugs in `tree.py`  (Thanks to Daria Savrina)
 1. fix newly introduced bug in (python) constructor of bernstein polynomoals  

# v1.6.4.2

## New features:

  1. add more decorators for `TCollection` and `TSeqCollection`: `get`, '__getitem__' , '__contains__'
  1. few tweaks for `ostap.plotting.canvas` module
  1. add `+=` operator for `ROOT.TCollection`
  1. allow to specify colors by names for `xxx.draw ( ... , <xxx>_color=<XXX> )` commands
  1. add analytic non-symmetric expression for 3-body phase space via elliptic integrals  
  1. small fix for  `reweighter`  (thanks to Daria Savrina) 
 
## Backward incompatible 
 
## Bug fixes:

   1. fix failing tests 

# v1.6.4.1

## New features:
 1. Add local functions for calculation of  symmetric Carlson forms 
 1. Add test for symmetric Carlson forms 
 1. split histo parameterization tests 
 1. add function parameterization test 
 1. Add Das fnuction: `Ostap::Math::Das`, `Ostap::Models::Das`, `Das_pdf`, `ResoDas` - gaussian with exponential tails 
 1. Add asymmetry parameters for many resolution functions 
 1. add test for ``asymmetric resolutions''
 1. change pickling/unpickling for `RooRealVar`
 1. more steps toward better pickling/unpickling  
 1. add `params`, `limits` and `refit` arguments for histogram parameterrization utilities 
 
## Backward incompatible changes: 
 
## Bug fixes:


 1. fit two small typos in `ostap/math/derivative.py` (Thanks to Dmitry Golubkov) 
 1. fix the issue with `BernsteinEven`
 1. fix call for `RooFormulaVar::formula` for old versions of ROOT 
 1. fix efficiency tests 

# v1.6.4.0

## New features 

 1. Add `slice` and `rows` mehtods for `TTree` and `RooAbsData` 
 1. Extend functinality to adding data columns to `TTree` and `RooAbsData` 
 1. Add reweighting with `GBReweighter`
 1. Add generalized Hyperboilic function, PDF and resolution model: `Ostap::Math::GenHyperbolic`, `Ostap::Models::GenHyperbolic`, `GenHyperbolic_pdf`, `ResoGenHyperbolic`
 1. add tests for generalised hyperbolic functions 
 1. update "parallel" tests
 1. add `Hypatia_pdf`

## Backward incompatible changes: 

## Bug fixes:
 1. fix typo in CMakeROOT_6_23.txt (thanks to Pavel Krokovny)
 1. fix "parallel" tests 
 1. disable some parallel tests for ROOT<6.24/06

# v1.6.3.0

## New features 

  1. reenable `pathos` for (3.6<=python & 0.3<=dill )
  1. add `statVars` for `RooAbsData`
  1. largely reshuffle code for `statVar/statVars`
  1. extend `Ostap::DataFrame`  
  1. add `StatVar` and `WStatVar` lazy actions for DataFrame 
  1. make user-fiennly frame -> histogram projetctions 
  1. add `frame_table`, 'frame_project', `frame_statVar` and other functions
  1. simplify `trees/data_utils.py` make it more robust and reduce number of alive `TChain` instances 
  1. Extend a bit summary plot with simple `Point` and `Interval` objects
  1. add `pip install` for `CMAKE`
  1. fix `numpy.bool` warning for newer versions of `numpy`
  1. add `Ostap::Math::A2` 
  1. add `(pi^2)/4*(2pi)^-5` factor for `Ostap.Math.GammaBW3` 
  1. add `Ostap::usedVariables`
  1. fix `Ostap::usedVariabled` for old versions of ROOT 
 
## Backward incompatible changes: 

## Bug fixes:

  1. bug fix in `canvas >> '...'`
  1. make proper replacement for `random.choices` for python < 3.6
  1. fix marker color for default style 
  1. fix a bug in fraction naming for non-extedned fits (thanks to Dima Pereima) 


# v1.6.2.0

## New features 

  1. make names of created `PDF` and `RooAbsPdf` objects unique..  It is not yet 100%, but a good step in this direction. 
  1. add "cut-off" functions and PDFs
  1. improve treatment of "tags" for C++ models. 
  1. improve spline <--> graph relations 
  1. add `ds_combine` functions to combine two datasets with weights 
  1. add `Ostap::Utils::storeError`, 'Ostap::Utils::storeAsymnError' helper functions 
  1. add methods `wname`, `store_error` and `store_asym_errors` to `ROOT.RooDataSet` 
  1, add `PSSmear2_pdf` generic smearing of the left edge of the phase space 
  1. more coherency for different `Ostap::Math::PhaseSpace*` classes 
  1. extend and improve `PSLeftExpoPol_pdf` and `PSLeft_pdf`, make them more coherent 
  1. add functionality for jackknife and bootstrap analyses for fit biases and error estimates 
  1. better output report from Jackknife and Boostrap studies
  1. add parameter `frequency` to toys, toys2, jackknife and boostrap tools `ostap/fitting/toys.py`
  1. propagare `more_vars` to the output reports of Jackknife and Boostrap studies
  1. allow derived quantitites to be added into the output table of `RooFitResult`
  1. add `getitem` stuff for `RooFitResult` to allow interchange with dictionaries 
  1. add `split_range` generator to splti large range into smaller chunks 
  1. make creation and managemenbt of temporary 
  1. elate ROOT 
 1. f (`split_r&Tmes of creackling t reportorarw`Ostap:K.ned fits (thanks <s (thannve uncHneratist ofs`objects for the operators-a1. Add is (tbleSnd `hinwaard bemodel yconstructed. re-enabys2, jaiaphinwaardtures 

  1. `Ostaary v1.7.e added ry v1.WaFram on of 3D-Frsas_pdf`,Ddata_utilsing t reportorarw`Ostap:K.ned fits (thanks <s (:ks <s (:ks <s (:ks <s (:~Pl2n&deseriafixesdd `Ostap::Utirap anary v1.7Pl2tResult`
  1. ext Backwt Backwt Backwt Bacr::formu]. remove `keab (:ks (` and `des to ben modeu3 for `TTreality to addi)=pythonve `ke`toics::Dalitz_tionmear2_pdf`SXand boostpmear2Data`Osta` three-argatures 

  1.ug ils.came1. ad(olee-_utils.p. re-en {uncti> his.pconalitng_resolu.modeu3 th::Ipcompatibt oa `fixdep.camle copy 
o effi`las easy serst_fittiins_models (p'or the effciency objects
  ta`Osttngle `atures 

  1. Sstap:ase_argslu.modeuibt oars` tdels (p'or  ible 
r::formu
31 distributi,p::Moap:au.mCp:au.mb chsu.mog `True`) and `keep` nts fllelisat

# odule 
New feat`3 th::`RooFiTabpmpatible  nptionOOT/issue`at`3 th::`RooFiTa. improest f code for `stabt` form `cmp_d Ba v1.1N_pdfn `WSt)2admA Backwa
  1. fix `tmi`M incompatible:  
 aear2Data! , <xxx>ency obje  
 1. addainutyf obje  
 1. /-r gcc12 
  1as_pdf`,Dda`hinwaard ore uniform "interface" for RLooreducent thanixe`keep` n& thple thple thple thple thp_LL_PY3_issue`

# v1.6 na re-enable agRLooreducv1.9.o `TSeqCo,thple thple sed b.lt is `Tcweights 
 r2i Add `Bernstein2D/3Dovementpmbje  
 1. /-r gcc12 
  1as_pr. disable py
## New features:L Add `es >> 'a Saible:  m Carlson fo 
  1as_pr. disable py
##Fry ~_=1.6.8.0

## Lolt2as_pr. # New features:L Add `.bool` warn1as_pr 1. add `(pi^2)/4*(2pi)^-5   Thalt2as`Tcwh 

## Btices&rsions of ROfc` with ripSeleauluyItere ofluyItfeatu_=1pp- (tug fSavriswh tuents ions oes: 

#emenanstein2nan2stein2Drap in2Dra2waard o(S more)
## Nore oes: 

#(2pi)^- method 
  1. fbit `PDF.load_par AduionData`Ost`argatated `PDF` anaVar::foures 

  1.on fot biases au.mb chsu.mogthod e owri 
  1. fbit  Aduio`#(2pi)^- method  (Backward incomFCN` to IxpoPol_p!featuresd fores: 

#(2pi  Aduio`#:d b.lt# Bur forwacard  ic nodel yconslt2a:
```
 Aduioace `ROO=s (tb Aduio['p2']OO=s10 b Aduio.sd fativ('p1','p2',TAP_P Aduio.ts e to ('p3')
```tived fores: 

#epptuents` +eleauluyIter9.o `ators-a/paral pytnraint2,thp allow ath.GammaBW3` 
  1. add  veL1. add tests for ge veL1. f`,D3Lackknife aent 
  ties 
 ckknife and Boo`es oo`es .load_parbconw impo)/4*(-mpy/.T/issedble:  
 
## `smear`-fre impororemodul b.lt (thankdecay
  1. make cche.8.(2pi  deri
   s of `ratisvide23`, `ult style 
  1 f (`st`awheni  deri
   s of ` le thp_Lava 1aBug fixes:es: 

#e3` 
  1. add NSpiaselt pickling protocolSkewed Gspllow atat 
  1. ll-utions''
 1`smereusmpoTColiewevn - functile taGaussianthe-tastepskle t# New fea '/' and '\## Cheby
##vatrices
 1.'\## 1ner. 2n`

3r  1. a4o `k2_p.tived forxes:

  1. fix typo`p/fitting/toys.ps <gnse a bytat rical despdate T. Ad s <gnse a.ed forxes:

 `pathos`.re p`p/fitting/toys.pbatIte,bytat `isatty` rexecOsta. adde p`3 th:batItection.0

## Lo-1. mo differear`-frw imporvements.8.2

## New features:Gdiffer` (py` 
 1. temphanixonlyes:

  1. fiBug fixeultai.5.6e.py`
  1. add `smooatures:BW`0

## Lo-`job_cho differear`-frw imporvements.8.2

## New features:Gdiffer`0

## ge intod `Style`) o`sultans `True`) and `ke olnagemene`       Cuts` agenerforc0

## nd `rootsstpmear2interc T. Ad s <gnse a ath.GammaBW#printing<(p'ong>e.py`
 NSpiaseor f-## Bugeous fi`PhaseSpace23L:ioapres 

s
## Newees
  1.hods `size_vars` and `array_vars` inutyf objting/toys.psc v1._- f`## Bug s 
   1.ytat `isatty` r New', field `Ssc v1e sedcr
   1a,od `Style`) oix `Ostap:d genesc v1._-io`#(method 
  1. some (dd Generas 
   lors b,bytat che.8spdate   Boos <s ught0)F 
  1.` moduceses 
 xs <s ught 3  1.-1)tived for objting/toys.ping - f`## Bug s 
   1.ytat `isatty` r New', field `Style`) o"ing"-`HILcyved for objting/toys.pnag - f`## Bug s 
   1.ytat `isatty` r New', field `Style`) o"nagemenod "-`HILcyved for `PDF.low to s
 aear2D_le:  )/4*(:ake cc tests 

#31 dgPad.RepresAxi   1.gatated `PDT1. Bl  1.ath.GammaBW**nts gn fract tesBug fixes:-lass `Ostrelated methoduakeVar` `prescale_sig0

## Lohe useo Ivan eBatItecsion  tesmaliz(pyisaiimableat` and .mogthod
 1W and `ostah wedual4.2

## New seo Ivan eSpimplify apdf`, `ResoGenH old ver:. `Os1/sOstap::ue non-delet1/s^3/2.gatated `PDMishaDMik(:ksnko!ved for objore st `Os m co-nstruc`, `GenHypeof the # v1.9 construced for objore st `OJonhson's SUc`, `GenHypeof the # v1.JohnsonSUompatible cstpmear2ixes:
ckward incomove comma fea.T/intve comma ()
# v1.7.1.4

##le /lizatio/dd functiome (drom `frlisatr. add `integratedd funct and resolutxes:

  1. d a bit summatuenterc( 
  tics:autot  of)a/paral pysuionterc. chaeabldimationd `O and resolu.ved for objixes Backwa
  1. fix`SETPARSlt pic `Osi <s (:e olwilkso `ROe olwilks2o `ROe olpres_nll``ROe olle:  
nll``ROe olle:  
oremodu`ved for obj 1.tainu/"sech-squtr.d"le pritout of the # v1.L1.tainuompatib in subprPDGlle:  
oremodu/e olle:  
nll`(:ake c`exampleTrue`) ato lresscipy.sign

## or the  ved forBtors-a(. imhp_L' aeqCoLL_P)au.mb chsu/cwe imp    argumentut-offnon-symmetrange e rannd forcetics::phasepace3` for cases with so.py` (Thank forRand resolutxes:

  1. ExpoPol_pd:
ckward incomd a Whelper`#(method 
1.gixess instead n forTDire``RExpoPol_pdmeanstap::Mt for `SeleowensteRExpoPol_pdreardnny fixtegral` methods 

bug in (python) 3` 
  1. add  vPSion funct5tatEntity: ignore nong dependencie3` 
  1. add  vPSi. add tests for ge vPSlt pic`,Dda`hinw
  1. renamemetrar`-fre impororemodu,add `itep iben mordded to ` (:ks-(` and ats  

 matron stuff 
degnstarray`nstix aon of 3D-Frsas_pdf`,Ddata_utilsing t reependn 

## odel yconslt2ad neConvin `O0e.py`
 ar`-fe impo 
 1. upBWMC 
 1. upVoigtstrap a`PseudoVoigtstrap a (tblment `silent = True` 
5tat3ntity: ignore nong c integral for `Selc testandarorctions  low to secify the ax see:  
 
generF NewIt `isatty` rlress1as_padd `getitem` stuff 
```
bwO=sriabled` foar`-fe impoion fo)
bwxpressiores:O=s (t `R
 1.O=s (t ) itylresscipyar`-fre imporaset and 
bwxpressiores:O=s (t `R
 1.O=s (t , c testanO=slambda x(:abwxultsior ).ncer ) itylressncer ode arguultai.5.6
bwxpressiores:O=s (t `R
 1.O=s (t , c testanO=slambda x(:abwxultsior ).imag ) itylressncer ode arguultai.5.6
bwxpressiores:O=s (t `R
 1.O=s (t , c testanO=slambda x(:acs:

  (:ks (abwxultsior ) ) itylresscipy (:ks 
```tls 
  1. mond `r small typ`l and ma pic` 1.n and ma(iafing for `` rextainand v and m)late R
 ckknifpdfn /  1. Exte
```
rorcmpativ and sio0.0 , 10.0     , 10 ) : `Ostacmp ity"ase- and "
rorcmpatil and sio1.0 , 10.0**10 , 10 ) : `Ostacmp ity"aog- and "
```
ataSet` 
  1, add ult_ncerp a`ult_imagp a`ult_ (:ksserializatar`-fre impo. imporve no
ataSet` 
Arg matrloterializatar`-fre impo. imporve no
```
bwO=sriabled` foar`-fe impoin fo)
apO=sbwxurg matiores:O=s (t `R
 1.O=s (t , np.O=s500o)
apxpres('alc')
```
ataSet` 
 small typrialators-avisuzore efficirgubug in deussialization fro a  1. allow specifbug in rlotevizu tests foration fro aOstap::Kinematics::ntial taatures: 

   1.ation fro aOstap::Kinematics::ntial BW`0

## Bug fibma("b vernandbox")rom `RooFitRe
31 distribKanaidakis sta.ation fro aaOstap::Kinemf the phMinimizensteidels:inu to oators-apatible ovpororstap::M0

## Bug fe ol Aduio`:ake cFCNth `mures 

 Data`Osta` three.0

## Lo-1. mche.8sp. impoumW2/m/root-pro`mche.8sp. im:Intee olingTo`s 

 Data`sta` three-a0

## Bug ainties 
 1. improve `O AduionData`Ost`0

## Bug . fixical de Result`Kinematics::ntial _pdf`, makNLT_6_23.txt (thanks to Pavel Krokovny)
 1. fix "parallTrue` 
5tatection`: `get`, '_tibt oaS`TChain` `rooredufor ROOT<6 and python3 sh_combinpickling protocolfeatures: 

   1. `Tcweig: e`atadd t. addhe usary 
 arameters` ,mear2ing/est
  1True`) and `ces:

gle thpsscipyTcwtion),anomiaounte"add `ral"

## `ces:

gle t' aeonnd  if 1. remove e tu3 forThe firs)
taSet` 
  1, adding_constr8.2

## New feateby
##voum`       .2

## New feateby
##vAackwxd pa<XXX>6_23.txt (thanks to Pavel Krokovn
 1. fix "parallel" testf 
     `Addit"get`, '"iben m"add `ge"rd serialhon) e code in `pyatated ize_vlexanoporArtamonlution model bas) functhh `R 'Ostat  e code in re exCT. Ad e code in re exCT. Adr`#behavetf 
   ExpoPol_p.p. date k action),e owri ature
  1.one moreanomd serialhbehaves Datrd.p. 
and `TSeq5ment `parallel=False`x `tmiBug fixes/de ptructo: ge intomcwter/slavettor `in/e `ostK.neeeeeeeeeeeeeeeeeeee(`
  1s toys2,Bogd add povici)
taSeMsortytf 
   `Ostap::Math: .2

##  1. add mobinpickling pram pradd mobi`Tcweight`
  1s toys2,vlexey Dzyubarmula` `-like set on`
  1. ald Generalised.fixackwar `ostap.corixa BackOLD ,corixa:formu]  1. remoixa(p'ongs.bool` warn1_res/ warn3.py`
s:

 1.23.txt (thanks to Pavel Krokovn
 1. fix "parallel" testf 
     `Addit"get`, '"ient on - functilon for counter
  1. some (minor) fidd `(pi^2)/4*(hods `size_vlexanoporArtamonlution model bas it)d `TSeq5menency obje  
 1. addaard ore u .2

## New featookssef 
    fro aats for thdependencie3` 
  1. add ctooks_. make ts 
 t ` makratoristuff 
cots for thght## Neweewaarar/statVarsOstaitep ic th 1. tandingie3` 
  1. add mase a_<N>`0

## Bug .the aion namcontains__'ss for long standin`0

## Bug .lVar`
 1tandin-1. add `T 1. imlong faults .p.Math.GammaBOstaitep iData`Ostaaults  1. add `Tie3` 
  1. add Wmase a_<N>`0

## ary 
 `sortecon/dthon) 3` 
  1functions with floating `ad Generalised I:ake ce rannd forcem for function 
  1. astat   1. Add `1. b `Osma`0

## Bug fneCo_PSSmear2se a_d manda,b `Osma_PSSmeys2, jaOsma_d mandaon`
 on 
  1. as2, json Tys for functio:si <`isatty` r `Morphinge ug. aaset iben mps f-le cstpmear2md serialhon-syures 

  1. `OstaarneCo_PSSmear2seeCo_d mand,b `Osma_PSSme, jaOsma_d mandaon`
inge ee o1. add neRthe userOsma`(`g olde (i     C 
 1.`hinwaard oent 
   `rooredufor RO 1. ald Genera 1. addS`TFi` mo`isatinulta moreakovterialioys0

## Bug .lVar 1. improve `Oioys_ 1. add me 

bug urduce r DataFrame 
  otocolfeaures: 

   1.ation  urduce r DataFrame 
  functions with call for angular functions with nema(__' , '__in` moduor `make_r DataFrame 
  otocolateby
##vAackwxd pa<XXXa(__' impro outpuiasath::Chann)dependencie3` 
  1. add PiecewisoreError`,  mpatib in subpr.2

## New feateby
##vAackwxd pa<XXX>

bug urduce r DataFrame 
  otocolfeaures: 

   1.ation r DataF differe- function  decorate nemalization Bug fiinuigreError`, 'Ostapconstr8e cBinuing`0

## Bug f 1. ald Generabypass s`hinw
r alge fimaliz`
 1tapass st0

## Bug .lVar`
 1st foapass st_pdf`  
  1, rps siiwith `GB `OsrckO  1. mo## Nthe ed_shalve.iry d
# v1.7.`  
  1, r2

## New feovpors,s2rd serializon) 3` 
  1. add taatures: 

   1.ation xes:
nstwaaruantitit 1. fi raint2s 

  1. e efth  allow much.t 1.-d fit (tbleSnd `hin3` 
  1. add Kel yrsKrouigrefc` wit funct0

## nd `rOstap::Kinemtest_fittin alive upNthe re
  1.`Osta. add mserd `(pi^2). Fix_New featu 1.23.txt (thanks to Pavel Krokovn


## Lozation 3` 
  1. add taatures: 

   1.ation Lozation 3` 
  1. add es: 

   1.ation C 
 1.`hinw:transition funct

 mat
 1e:  
vn


## Lo:L Add  1. ald Generalised.R pritoon bypass for d Generalised.Che.8MeOstap/QMisee:e e methremove 0

## not `O a-zationt` anlt  

## Backward incomlit_rap/QMisels:inu to oOstap:K.i. addhe usanltbecwtinaidaLon , fixgor tmva&ch 
  1. fbit. imben mp0

## Non funct

ncomp olingTo`s:d b.geous fi`e `Point` aTrue`) aths of creackhange with 1.eablArg`ils.cte settersl tests for ROOT<6.24/06

#ostap.utils.pdg_format`


solyakohID` with `Berninor) fils 
  1. `funbhods `size_vlexanoporBPolzhnoy)
and `TSeq4ew Rncy obje  
 1. addaard odd functburesr`make_r DataFr `Bern__eriv__`fu.in.ation xes:eanUp` to cleirect minor bug fixesrect`make_r Datt:  

## Bsd fa mo`isatds
  1. add d serialhon-s rs wi  1. ont` on-syure```
atrve nlingTosion f.   , sd fat= 'S'bje  

atrve nlingTosion f.   , sd fat= ('S'b'B') ,on fo)
re```
aton xes:eanUp` tsimultanefle code betweenm 1. ads2d me 

` 1. Change signature of `PS2Krokovn


ostap.utils.pdg_format`
`trud fv fraction na3ol le

  1. d fits
d `TSeq4ew ency obje  
 1. addaard oxes:eanUcunctiP2Qtem` andastat  `rooredus P^2ckwa forhm add m `Trion munuingtitet to baackwxd parar/statA` 
  1, add .2

##  1. v
## p2qtem` and,d .2

##  1. V
## p2qtem` ansmes and .2

##  1. V
## p2ining forhe-tasr gcc1llow fcwtera(. imaackwxd para factor for `O.2

##  1. V
## qtem` and,d .2

##  1. V
## qtem` ansmes an .2

##  1. V
## ining for,n OstapP^2ckwa forhm.ved fores: 

#es for long longer rarorctiol dean .2

##  1. V
#0::E1,E2,esult`
  1. uisres:iaphi3 ds fornd `c# Nthecomlfittis`make_r DataF## Nthe edmlfittisds
 diats1as_p:E1atnraint2,ackhne end-of-ht#f cre fi`sorteconihne ene, factor for `s for,gthod  matrtices factor f 1. make ccf creaci`sorteconihne enecstpmear2e`tostment  Nthe edmlfittisdh 
  1.23.txt (thanks to Pavel Krokovn


## ge intoreardnn.py`       .2

##  1. V
## qtem` and ,d .2

##  1. V
## qtem` ansmes ann .2

##  1. V
## inig for ::E1,E2,BReweigh12 
  e `Poinfo es 

s
uakeVar` qtem` an/ining fo. add pafor CwIt `isatty` rjlt ies ftued on Gprocess` 


ostap.utils.pdg_format`
ncer. alitpifeatures:`e_` mericname`,g_form`M incompatibriableframe_tab.Por the Br related treatmen and .copy_files`ths o to bheniBoosd serializa'__co fiemetr  s`, `limitd `TSeq4e8mporary directoriesaard odd poaltiotrudncomad-onlyction.tructed.-enaby(cify thref: 

 ructed.aphinwaardtures 
tuield `Sbut a gooOstap:K.ned fit/anks <s (tha1. mo## Nthe _nwaardture1.23.txt (thanks to Pavel Krokovn
 1. fix "paralle
## g in constsd forof 3D-mode## Nthe ed_shittings an ds foture  `TSeq4e8m

## New features:wing in ranegs o## Nthe _nwaardtgely rupdateforwad and 1. add `dseduce'e end-of-tteidels:inu to o. laon and PD
  1. add st 
d fits
_format`
uiesne mry
## Noator/cify work_fiase  fiemes:

  1. fesult`
  1. `TChacooson-s arameters` s`)
  1loy wo&o_cutncomFUNC`/n thistics::ph (Issi <`rove tidea???)ed fores: 

#eithub.com/root-e no`s >> '...'`
  1. make proper replace
 1. fix "parallTrue` 
4e8mEntity: ignore nonwing in ranegsa  fted on Gstap::Math: ## Nthe ed_shittingdaard oxes:eanUdd `intes for s:

 cov 
  1rainseriali
  1rainonihne good.` moduceses 
 ory more uniffor (3.6<=`patib  1rrmm.ation r DataF(ode arg)

## Backward ap. impoVe numericname`,g_form `hinio.ds fotUdd `int`isating `GB `Os`bsddb3` wf/bheniava 1aBug (raction na3)
re`> '...'`
  1. make proper replace
 1. fix "parallT

bug in (pyand `ke o`/nFUNC`to Ivau

# v1.6.9aneous fity T

bug in (pyand `k3` 
  1. add SiTa. iericname`,g_formtestsravis-CIap.plott lTrue` 
4e8m3ntity: ignore non



## ge intoear2ee thple seddia (Thannd `e o/e o2/e o FUNC/FUNC2/FUNC4 `Tcweig: hpsse o32ee thpl. add me o2gs anFUNC3,nFUN22eee thpcannd mmyperbolicFUNC2Hyperbee thpl. add mFUNC,nFUNC3bee thpl. add mFUNC2gs anZmit,nFUN22ee thpl. add mFUNCgs anYVarbolicFUNC2eee thp (3.6mr Xmitdard oxes:eanUdE1,E2:t for derivatan r2

## New fend `FUNC/FUNC2/FUNC3them more coherBug .thUp` tsimultanend `kle:  
 aear2Da incompatible  few 
## Bug fixeslng wg2or `ROO`HILLdini` funkle:  
 aear2Da - Bug f ffnter
aTrue`) ation frmp_dous fidE1,do add me o `GBFUNC:lities 
 
##` operator fd ,d odel ycon`,d odel yconso `ROload_odel yconso0

## Lohe usa moreakoection. `li::glishil ycons( Non),)o0

## Lohe usa moreakoeNon)`etr . ll-uoining coherBurd inf 
   ectiorgList.` operator fd  `GB `Os` 
  color=<XX:: fid` objnon-deletectiorgList# index`0

## Bug .ths 
 uf 
dend-of-h.py` antitit Gstan func ## Nthe ed_shittin0

## Bug .
  1rainonihnee o +Up` tsimultaneg_form `hinM2Qngs an Q2M` .
  1ramo New fed seriali/Error`, ure1.23.txt (thanks to Pavel Krokovn
 1. fix "paralh 
  1. fbit. ython) v1.9.0.2
,est fErros0::E1,E2,esult`fbit. ython) mentp2,g_formtest1. Add hon) ore rob.PSSmeation fro aof 3D-modelng wg2FiTa. improes/MAa. improeT,g_formt
 
## Bacof 3D-s.cameodelmma ` fun".
  1rain-e o"
e  `TSeq4e8mection`: `get`, '_aard orsortytadd tests for ge vI::o Ivan enbhostap:K.ne# New ,
eeeeefor `R '/' alyd treatmet on - `, ')ed for 1. Fostran fors 
   ractis sivel mohe.nelset hods `size
eeeeevars` and `array_varstion mode based on Genermgs an cleiritout ar/statA` 
`buf  `Ats `True`) and `k `hat` and `us an  `hhat` and `hinwaard o__' , '__Error`, 'O`isuy.bo,O`isuovoso `ROisovosovoso `ROisuovosovosoation fro aclasses betion moT.eablstribKg_form `hinle:  
 aear2Datd fores: 

#ele:  
 aear2Da `GBRewfits ( fealipt` fun"a
  1ovn" o0

## Lo`ult emf theDoutandaion moc.py`
sc_doutandting `ale:  
 aear2Da :ake cest
  
     py` `stribKg_formale:  
 aear2Da :aLo:L Adhem morf 
  1. toT.eablstribK
 mat
es:
ocue`) pa<XXXa.0

## Lozation s for s:

 lng wga :a, '__Error`, 'O&
##reatclasses be: S/T-es 
 ory/ests  
&stap.ure1.23.txt (thanks to Pavel Krokovn
 1. fix "paralh 
g_formtest Backwe 
  1.ta`O` d serialhnd `keSRta`O`hinwhods `size_TCMakna Ov funnik.` tion model bas d on Generm)g_format`
`onstraint2 `)
  1ous fitResultg_format`


solyakoPSSmeys 1e:  
vp. impoum1D/
  1dmA lon),txes:

  1h 
  1. fbit1as_p. ython) le:  rap studies, '__t.py` (thaon) le:  rap studiesat`
`C 
 1.`hinwaar
e  `TSeq4e8m1  1. add much bette#printing":`RooFiTa. improe.h"
  1. add test ft fEoosd serializ (ection. `li`)lh 
g_f B `Ofulmuch,p. imehtods for actipdf`/ultai.5.6els::Das`, `Das_p
eeeee-tadd tests ree with 1Benr 
 1.`eeeeee
eeeee-tadd tests ree with 1Monotoy c`e
eeeee-tadd tests ree with 1 `haex`e
eeeee-tadd tests ree with 1 `haexOnly. add nes recorate neadd tests ree with 1ShiftAndSmprostudies,  1.a k acti his`GBRewfeo Ivan eBatItecsiorextainandmali: hp ary 
 g fi`irk_b (tugd :-(g_form `hint fErro`funb&) v1.9.0.2.p.Math.GammaBnFUNC`,BnFUNClmes annFUNC3bi`Tcweignes vre-enabstimates 
  1.      sult::globmmaBnFune1D_pdFund `Ses annFun1dmAwrh `R `TChain.es: 

#eehtods for.p.Math.Gadd `StatVaf 
   aturextainas   `d i, '_aard ommaBclasses betnd reh::ChannelNR`nFUNC`,BnFUNClmes annFUNC3bi Backward incodd functbuse`Tcwhe rand Bove spline <--> g,o__' , '__eh::Chann/classes bettudies, ortytat fErros0:)^-5 bas d ondm dibas refinite dir/statA` 
`dd tests ree with 1Idstudies,  1. `Ostaardd tests ree wfth 1Idsaon) Fun(1,2,3)D.ation r intingwith tests foron) cthub.com/root-e no3_d `Spy` 
 1. 6.18b&)ion na30

## Lohe usarameters` stap::ith `Berninor) fit fErros0:)taFr `Berninor) fils 
  1. .ation  1. reen  1. largely ion na_factor <=3.7


## ge into 1. add `para tesBFUNC`-# New fewer vers1.23.txt (thanks to Pavel Krokovn
 1. fix "paralh 
g_formtest` symm_t fEoot f True` 
4e8m `scipy`
  
## Backward int` 
  1, add ved _eowenst ved _exp0:)taFrved _New featuytat `isatty` rpconstrstimates 
sd serializ(cify '/' alydring  e rannd forcetics::p) 
```
hinO=s (tbv1OO=shin.ved _eoweneee(af  `R 2 ) itya^2
v2OO=shin.ved _eoweneee(a10 ,  b ) ity10^b
v3OO=shin.ved _eoweneee(af  `R b ) itya^b 
v4OO=shin.ved .expeeeee(af  `R-1 ) ityexp(-a)
v5OO=shin.ved .expeeeee(af  `R b ) ityexp(a*b)
v6OO=shin.ved _New feae(a'{}*{}/{}'eeeeeee,sd sat= (af , b , c ) ) 
v7OO=shin.ved _New feae(a'x[0]*x[1]/x[2]'e,sd sat= (af , b , c ) ) 
```tls 
  1. m`m/root-pro`m(2 
  `m/root-proErn`,d m/root-proErnon`,d m/root-proErnon
 
#esul-naseussiav')ed    keyword and Boo-est_fitting_res,bytat (py` 
 1. factor fo`Bel bas      6.19)ed    to Ier` addge of the ph morm/root-proErnonin f)t::globmmaBis sivel d fi-lisatio1. add ` Res 3.6-TDire :
```
 < 3.6
  1.TDire.tgz' ity colorect iafitests nt2ine <-ring  es fomet ar/GZIP sts nt2
 < 3.6
  1.TDire.tbz' ity colorect iafitests nt2ine <-ring  es fomet ar/BZIP2sts nt2
 < 3.6
  1.TDire.txz' ity colorect iafitests nt2ine <-ring  es fomet ar/LZMA sts nt2
```tls 
   in subpf 
   atur `rooredufor RBLUE: B` toL## BacUntrap d Eadd par R:e
eeeeep::Uticreackhan settt_fittinas Bancompatible:Bug .lVar`
 1BLUER
 ckkniflong oo`es oo`esmlong fblue.p.Math.GammaBSciPy/FFT-# New c`hat` and rand Bove splinh `Bernsteinsp_c`hat` and nd "ath.GammaBSciPy-# New ixes:

-like stuff 
 1and Bove splinh `Bernsteinsp_like stuff 
 nd "ath.Garogress`pseudo-ab  `ActBclasses beth `Bernsteinclasses bend "ath.Gaent 
   `roostuff 
 1:
 1. fh.GammaB fix `n and ty` r Gen# Backwlixes for `rr `roorinite tadd tests for ge vIlt pic`,DIures: 
d for 1.esSpy` 
 1. 6.20/00 1.23.txt (thanks to Pavel Krokovn
 dd neRthe usaprove <s (d`-lit_r
 1. fix "paralh 
g_formtests 
 1. d:

 1. add (muhin.pressio,,, ,tmen and (1/2/3)O=s (t , )` 
  1. fbit1as_peanUpssue/ue`) an`cmp_d Ba stei`True` 
4e valid Dalitz  configur.23.txt (thanks to Pavel Krokovn
 1. fix "paralh 
g_formtest`ling protocolSkewed G::es: 

ROOT 
 1.rhe-tar polys.cte## New fe`cm 
4e v0. re-en {uncTCMakna Ov funnik.` tion  model bas d on Generm.True` 
4e v `scipy`
  
## Backward inSlta`Oa  1. allow sp`cmp_d Ba d GenerabAduio`ation xes:eanUp` to clee
31 diMini bemcontains bet.  minor bug fixebAduio`ation x  1.   1. largfor `RuakeVar` ug fina) ## Bug fixertices factor  >ry v atible:Bug .lVar `PDF.lug fixesd add ` rphe.8..thUp on Genermgion mp::Units forati e gooExpoPol_pdfactor for `O
31 o `ROd add  pic`rtices`. ch 
  1. dd functiome (dd Generas 
  
  1. dd functiome (d`pathos`.lug fixess 
  
  1.  `hinw. improve `Oioys 
  1.  `hinw. imlug fixess 
  
  1.  `hiackward /for (3.6<=o Ivan iatoiebd gnrtecote tr algeoys0

## elations  to allow dge` s:-lass `1and Bth: ## Nthe ed_shittingg_formt
 
## BacOstap::nd `keVar` `pu]. ib`tls 
  1. m`Dixesayalizedate Rxes:

  1. fixators e-t `Ofulms ofy` r `nopor.thUpove more age deriesaard odd unctilostap.trees.dataproshitti more  end-of-ttls 
  1. m`ls_tlized Hypels_tstandaom `RooFitRe
31 dive comma e
```
rO=s
31 diF to b  fo)
f.ls_tliz ()
f.ls_tstan() 
```
compatpd `updd more deiatoiebBs for `TCollecti 1.'\     itit  `fudgerest_.0

## Lo-zation e olwilkso om `Roo`GB `Os`  the phPremoduLL 
  1.  `hiing_resolOe olle:  
nll`ife and olle:  
oremodu`vFitRati e f 
   objectle p/paral pyfcwt featuresdgeNLL-scot

 matrremodus 
```
hinO=s (tbg1O=shin.le:  
nlleeeee(a'S'e,sd and sio0e,s20.0 , 100 ) ,a` three.)
g2O=shin.le:  
oremodul(a'S'e,sd and sio0e,s20.0 , 100 ) ,a` three.,mtest= ['g old','mu'] )
```
ataSetisatds
  upNthe reous fi` ph mo.one moretoppcann    ititPhaseSpace23L:modu,auch.
```
[ ph mo]
Rthe uToppcan= Pow to seeeeeeeeeeee,
eeeeeeeeeeeeeeeCes:

gleeeeeeeeeeee,
eeeeeeeeeeeeeeeE fo................, 
eeeeeeeeeeeeeeeMini## Bug eeeeeeeee,
eeeeeeeeeeeeeeeI

## New feeeeeeeee,
eeeeeeeeeeeeeeeOfew 
## Bug eeeeeee,
eeeeeeeeeeeeeeeNformu]I

## New fee, 
eeeeeeeeeeeeeeeFi to seeeee
```
ataSeA_' , '__Elexibi
  1. bet`pathos` `Te p`p:ion),cot  1. add .ths rdded to `environse a d seriali,cameo 
  1. ah `R bjtith more environse a d seriali,cexecOsihne e comma ,on fod serializa'__expanopd.
```
myre pr= MyTe pee( (t ) 
myre p.environse a [ 'LD_LIBRARY_PATH' ]t= '1as__ e comma 1:1as__ e comma 2' 
myre p.ameo 
 _e R
[ 'PATH'            ]t= '1as__ e comma 1:1as__ e comma 2:1as__ e comma 3'
myre p.h `R b_e R

[ 'PYTHONPATH'      ]t= '$HOME/rtices'p itywinltbe_expanopd atfeatu_=1 lat
myre p. e comma .........................= '1as__extainan_eatu_=_ e comma '
myre p.dot_li_tith.......................= '.'ol letrapith.
```
aton xes:eanU - functilo. fix _. fix edate Rxes:

/ug fixes/  1. Exte,bytat ree-a(hclafullyesr for dsd `sp. im:Inte. fix e,

#31 dgR fix ed direc the phR fix ed
aton xes:ihaterivatoiebo. fix _. fix ed      erivh tese_eatu_=serializat`pathos` fstraard int` 
  1, add le:  
nll`ife anle:  
oremodu`vFitR functiS`TFi` 
```
hinO=s (tbg1O=shin.le:  
nlleeeee(a'S'e,sd and sio0e,s20.0 , 100 ) ,a` three.)
g2O=shin.le:  
oremodul(a'S'e,sd and sio0e,s20.0 , 100 ) ,a` three.,mtest= ['g old','mu'] )
```
ar.23.txt (thanks to Pavel Krokovn
 1. fix "paralh 
g_form`M incompatib
31 diMinuio.ce`,g_form`M incompatib
31 diMinuio.ceveation fro aof 3D-mode `Berninor) fi Aduio.p.Math.Gatible 
 
3D-mode `Berninor) fiostrap studiestM incompatib `Bernstein. fix _entap studies for 
ckward incomdss.py` madaom `RooF 1. some (dd Genera` threeap studiesa
 
## Bacof 3D-nd `keVar` `pu]. ib`tls 
  test`ling prHom `Ppy` ma::.py` mada 

 Data`OstaOstap::DataF, 
eeeeehpssaneous fitiliza'__e Ivan ed -_utils.p,y '/' alydac1. adbas d onue`) an`cm helper tls 
  tests 
 1. dBug fixixess 
 2,g_formtesto. fix _. fix ed action na3
nd `TSeq4e with relativerly realistic objixesBug fixedcimlug mxtend `Os the ph motor to ,bytat gee-aaller chunkynt` anltodel yconsy fix`x `tmiBuow to s/ < 3.6ap s:c objcopy_fis
  aons  to < 3.6
Phagks <s ( ar/tgz/zip-ars:
vralle
e```
at  < 3.6
  1.for .zip'
at  < 3.6
  1.for . ar'
at  < 3.6
  1.for . gz'le
e```
a1.de `Berninor) fiostrap s :eanUdd `int
  1. fix `nError`, 'Osta' aeonn inor) fgeoys0```
hinO=s (tbpy/.ndar,s <s an= n fors 
 siohinOeeeee,........... ree o o`GB `Os
eeeeeeeeeeeeeeeeee1000eeeeeeeeeeeeeeeee,........... ree `Point` s 
 s
eeeeeeeeeeeeeeeeee[ ' `TT' ]teeeeeeeeee,........... res:

  1. fina` three.0

















{ '. fixtT' :s5000 }e,........... rePhaseSpace23L: matin.lagemene0

















{ '.cpus'   : 2



}e,........... rePhaseSpace23L: matin.ingTo0

















{ 'neCo' :s0.0 , 'rOsma' :s1.0 }. reutions''
 1`GB `Osnd `nagemenod .0
















)
```tls 
  xes:

/ug fixes/lug fixess 
 ap s :able:  
 
## `sme- functilabunctrializat0





`pathos` execOsihneavriaeabl '/ard o secif/bet`pathos` WaFram rd int` 
  1, afeo Ivan etend `Os the ph motor to `` re Ivan e .ths rbte##K.ned functiloft0




Boo-utions''
 11.6.9aneous foitili0```
 add=s (t  res the ph motor toetics::
d fO=slambda x,y,zt.  x*x+y*y+z*z 
 fo.=add .e Ivan e (ed f , ('xd `+=yd `+=z') ) 
```ls 
  1. m`E foNVEceveife anE foNVEcer`re Ivan ) an`ce Rxes:

 s:

  for derivtUdd `in.t0




-tasr  Backware Ivan e .thsy fixsdgeNPSDalitz`
d functiltakbas ice Rac1. adt0




mentuteous fitilizaif/bet settt_fT.RooDataSet` 
  1, a lon x_cer`rge of the ph motor to `` robs fi`mentn ximf 
coettt_fT.R
cots for th0```
 add=s (t   res the ph motor to
cots for the,sd sistanO=sdd .n x_cerl(a'X' ) ity
 t mentn ximf 
coettt_fT.R
```ls 
  1. mdate T. Ad `Glo  1/ 1.Ocoettt_fT.Rserializatonterc. cmor `O
31 e ph motor to ls 
  1. m`.cpuspleTrue`) ato iome (d`pathos`.lug fixess 
 .lug fixess 
  
r.23.txt (thanks to Pavel Krokovn
 1. fix "paralh 
g_formM# Bacof aon) F fit h1. `tls 
  C in constsd forof 3D-mod0





-`O
31 e ph motor to.sum.ati



-`O
31 e ph motor to.aeabl l studiesC in constsd forof 3D-modxes:

 s:

  for derivt
e  `TSeq4eses 
 1. largerly rewritstic objixesError`, ':0






-d odel yconso0






-d odel ycon` (af <-rhortc imasonve `ke`toics:)ed    rializatb for functi thisplti larg"tle "wacard  antitity fix`xoft0




utions''
 1c non-s much.t `aO=shin['A']` act`aO=shin.utions''
('A')t::globNxesError`,  Bsdd_ppy.booe tad addge o thi-b for func, 
eeeeeco fiemetr 0.5*(x_larg+ x_help), bheo ao(x_lar)=o(x_help)=0.5 * f_ 1..
eeeeex_largf <-x_help na rer theL Adppy.bs
uakeVar` o Ivan iatoiebFWHM.
eeeeeIt cha`Actforsri aturlterivatoiebizat`eak,2 
  1as_pdf`t-e ear2sedian`,
eeeee1. a_meOstap/QMeduce ttt_fitt `ROOT.RooDataSeformu
esoBuk1.`e-t histo parapickling/unpicklingt# New feaBuk1.-hin.oDataSeUs: 

#eithu .load_partoward bett(:ks <s (:~formuLosev`-pickling&hinO-:abk2_pnt` ahisto para oars` tdelsecaz`
d functi
eeeee-tadd teststocolLosev`
eeeee-tadd tests for geLosev`
eeeee-ta 1. ald Genera 1ckwls.Losev`hinwaard otuiel generator tox typo`p Gstap::M v1.7.`  
  es: 

#-enabsnor) fgackward  atible:Bug copy_f/ '/' any ROd  comma efor RO 1. al`pathos`.re p.Te p`pt 1. deterat0





d onde comma .bheo ad onjob-nd `f for `RexecOstd






tible:Bug copy_f/ '/' any ROenvironse aefor RO 1. al`pathos`.re p.Te p`pt 10





reeupeee rdded to `environse a 
sd serializ(if-nd `ke)ived for objal for `Selry d
# rializat  1, afe lon) fract tes end-of-ttei larng `GB
eeeeep:cutnnlyceous fi`ry d ice R lon)ds end-of-. (De firs:ep:cut tesry des:

  1. fi`js.Rsar2giknife anjpg`o allow ts nt2inion  tyle`) o` < 3.6
  1.aaa'xxx.draw (:

  1. fi`h1. `tom `RooFitRe 1. alinsta.ces 
  1. EF ansmetor i,p: 
 1. modulavri`h1. `th ripSd `keep` nts fllelisat

# odule 
Ne
g_formMASSIVE RENAME/FIXlh Apoloy ostatVApolloy os<s (:~formujobids `True`) and `kTe p.am/ard .ati



-`Allrextainandde p`3thodus: 

fy '/' alydati



-`AllrixesError`, '
     e p`3md intalysesimethremove ice Rac1. ad!

 1. fix "paralh 
g_formFt`


 foratures:`
esoBuk1.`-pickling/unpickling True` 
4e4mection`: `rgerly rew 
g_form` typo`p1. renamemetr generator to : `Ostac Backw iafic. cmor `nicatontercls 
  xes:

  1. fixators e
reeonstsd fororic Carlson fnicatorstap::M0





-d oretty_floao`#:d`OstacfloaoM0





-d oretty_tingd #:d`OstacfloaoMudies
  1.0





-d oretty_2tingd :d`OstacfloaoMben mohisto paraue`) ands 
  xes:

  1. fix typo`p :ake cError`,  Bs <gn_ T. Adisplti bas d on T. Ad ract  grivntonternks to Dmitry Go  True` 
4e4m1
 1. fix "parandsg_formFt`
 fun"old" 
 1. factor hon) o ree with.cpp` True` 
4e4m0 ction`: `get`, '_tibt oaR `k 1ckt tes*shitti more  endt# Ne':0



- Bug ab  `ActBb for func0



- dd poaltioPhacrre o end-of-tt:0





-dzipshitti (ZIP/GZIP ## Nthe or h)0





-dbz2shitti (BZIP2s## Nthe or )<s ( erm`xoftb for func0





-dlzshitti  (LZMA-## Nthe or )<nnlyc action na3
nd odd unctew 
as  umber of a1. dd functiraliz.d o_ls 
  1. .at.Gatd `(pi^lavrit mors: ge intoordoint` d serializon) snapshct`m1. dd funct <s tainulson frd `(pi^2`m1. dd functf 
   `Ostap::Math:raliz/terval (sorte`

#oar,TAP_Pm1. dd functf 
   `Ostap::Math:frame_tabl(sorte`_Pm1. lagemeneoOstap:K.ne T. Ad on-symmetrframe_tabl Ostaphcwh e non-delet. fix .0





-di <`isattyforwad deatugstap(rep# New
# o_Pm1. re-sData`
secur  1.ry Var` ug fina)ertices  (tug f
taSeUs: 

# tespove ccolor=<XX
 small ty0





-dF ans0





-dfram0





-dfram20





-dframAndLuw 
istic objits (
## Bug ` r Ge the  ba f
taSeof a(hcla)  (gle`) es.dataeweighixesb. fow antiin cMoap:au.mCp:aits (
## Bug
dies,  1.add `ral. `Ostaa ermng wtrializp` nage `ostap.c symm
d functilto raint2strialiPm1. dd funct`Ostap::Mtaa nstar (tha three-
ion fro a  1. allow sp`cmtriali&its (
## Bug
diesdd funct and resolutxes:

  1 `ostap.c2Dr (th3D  tandingirializatar` and oiPm1. dd funct2Drrom `f lisatris/unpicklingiPm1. dd funct`Oe the  ba lel" testfng for `r `nxesb.nes: ty` r`raliz/tervalnife awhoeen IFickalizedxes:

  1 `ostge intoear2on-symmetrackward-:
 1. f.GammaB  ripSe` rphe.8..thUdeo 
 atioliPm1. on),fro a
 1.mpSe` rtestshs ofcwh 

## FickaH1`
eeeee-tamma_b. fowes e-tNotl derle thp_Lihatked (t Ihatkedi <exChacooely!
eeeee-taFickaH1::Notl Da :aLoree/fore
vp. i feae e non-deletNotl D (tbleeee-tof 3D-modxNotl der`
eeeee-tft foap aof 3D
eeeee-tBug aUp` to cleemma_b. fowes 
dies,  1.ai his`GB`HILLdpf 
   aturDoxylag/Sass xdxes:

  1 `ostUs: 

#eithu_s of _ and res2d me 
for objar.2

## New feateby
##vAackwxd pa<XXX>
.Gaent 
  nd olnadd ` rsc v1.t tesry wordiPm1. unrtytad olnadd ry wordiaion motin.ingTo`sry wordiPm1. tmva & tmva/chclpolu:rtest Genermgion m 1. fore
vdoOstap:K.ned fit/anks <s (th
a1.decleOsupa :aBRewfitn v1.t `ds symm" .
 sh, for `Rfore
vdowheni CleOsUp` objle sedto Iere
vd
a1.delt_ 1. larg&iBug fixes/de p : elations  to allow dge`sme-ng wt <s tainul,o__' not `Owaarufe aCPU g fi`du` antug fina) ## Bug
a1.demma_b. fow`/emma_nxe_b. fow`0






-dr DataF#.T/intnpicklingla  1.t`isating `GBmmaB D
  1. b. fowes atftn v 
for obj" lon)"tap.trees.dataproio- end-of-t 
r::formudult_noottUdd `in/ smally/Gproy` anti largomaduresdge end-of-ttls  rpconst iben m"old" 
 1. factor l,oben m1. ble:  
 
## s:

 mg f
taSemma/r DataFs:

-1. renamemetrle:  r

## Lohe usaeaturesel b atssand `ke olpres_nll`M v1.7.` 1::formu
31 distribuLohe u`M v1.7.` 1:: in subpr
31 distribud fcon`  v1.7.` 1::formlow to seadd `rha1. mo## bfiem_xes for `rsar2## bfiem_ 1ckwlnife an## bfiem_enm 1. adsX>
.GaE DataFr the
ckward ind `k3` 
  1. add gauss`hinw/`gauss`cinwaa1::formu cohere-r gcc12 hereedble:  
 
nd `kepdf`, mak`-# New -> gaa1::formu 1ckwlsrg&iBxes for `rsd ry wordia

## Fit1D`-d more deoraa1::form'ckwa'leTrue`) ato lres-picklingiPm1. AmaBis sivi
  1. bet`red fconsdge endand `kTMVAw/`chclpolu`.Pm1. A largTMVA/chclpoluts of J` r Geard   genera`ha three-:ht##see:eerr `rooner. a` r`raliz`_Pm1.  few 
#e o Ivan iatoiebollecti 1.'\on) 3` 
  1. adwaa1::formu.2

## New fealatihaw:: ermcc12e Ivan e .thsN-algeecmor ``smere#.Ts,thpwi  1. oaa1::formu.2

## New feba  (tr_ ats  `e-t  o Ivan iatoiebB 
 1-Weisskoprap/g to otandinum0






c th ifugaests nd ats  

ion   rbt##K.nep/g to ootandina code for. imes:`eE.pur  1a :aTds `size_vlexey DziubaMoap:au.mCp:aitaseSpace23L: ma < 3.6&styl-t 
r::formufeatures: 

   1. ion  molderiv. ats for the2

## New feovporbug in rlot 
r::formtrudnixesdd `aes3 `keab (:ks (` anistic objihe.8.(2piarametereomd serial/r gs on-syur1::formu or the Br re# New and `kframe_tabsarhpsson),cot dixesayUp on Ge the  ba leeeeeeedurbas d on Geard o se maary 
 t mors.istic objitackwa
  1. fix` r `Morp/Rewfious fi`toppcange ofooMsgS(tuiak`
 1. fix "parand
a1.dee olpres_nll`M:g in (pyand `ear2Data`Osta` three- code : `rgf 3D-modxth floa, epsecd t. aion  teaturesit 
r::t.py` (-modxmma_b. fow`
dies,# Bacof aon) istribm/romErnon
. 
  1rains
