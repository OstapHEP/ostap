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
    
## Backward incompatible1

    1. remove obsolete `ostap.fitting.simfit.Sim1D`
    
## Bug fixes

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
  1. make creation and managemenbt of temporary files and directories more robust, probably more efficient, use better namings, ...
  1. add `timeout` parameter for `sqlitedict` and `sqliteshelve`
  1. make use of `berkeleydb` for 3.6<=python

 
## Backward incompatible changes: 

## Bug fixes:

  1. fix assertion statement in `dalitz.py` 
  1. fix a bug in analytic  three-body phase space for cases with zero masses 
  1. fix a bug in numerical three-body phase space for case with all zero masses 
  1. fix a bug in `__getitem__` for range/slice/index sequecne for the weighted datasets - the event weigth was propagates incorrectly. Thanks to Dmitry Pereima.

# v1.6.1.0

## New features 

 1. Add `Ostap::Math::Hyperbolic` hyperbolic distribition 
 1. Add `Ostap::Models::Hyperbolic` hyperbolic distribition 
 1. Add `Hyperbilic_pdf` hyperbolic distribition 
 1. Add `ROOT.TGraph.merge`
 1. Improve treatment of GSL errors 


## Backward incompatible changes: 

## Bug fixes:

 1. Fix some tiny incorrectnesses in `Ostap/MatrixUtils2.h`
 1. Fix small problem in `ostap.utils.utils.KeepCWD` context manager 
 1. Tiny fix in `graph_summary`

# v1.6.0.0

## New features 

  1. Add `FlattePS_pdf` - similar to `BWPS_pdf`
  1. Make few important steps towards ROOT 6.23/01 ("New PYROOT"). Full  adaptation is not yet achieved, there   are some pending problems with effective inheritance from C++ classes (namely `TSelector`, `PyPDF`, `PyVAR`, ...). There are also some puzzling  crashes... 
  1. rename  tests, make test selection more transparent and easy to navigate back 
  1. rename `ostap/fitting/selectors.py` to `ostap/fitting/pyselectors.py` to avoid the name clash for `python3` 
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
  1. Better (but not perfect yet) treatment/assignement of the unique names for many intermediate objects 

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
