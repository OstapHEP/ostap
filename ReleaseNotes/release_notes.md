
## New  features: 

  1. `table` function for `RooFitResult`: print object in a form of nice table 
  1. `ostap.logger.utils`  set of minor functions for nice printout 
      - `pretty_float` : print float 
      - `pretty_ve`    : print float and error 
      - `pretty_2ve`   : print float with asymmetric errors 
  1. `ostap.logger.table`  : add function `align_column` to aling the column for a given table
 
## Bug fixes: 


======================================================
# v1.4.4.1
======================================================


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
 1. add 2D and 3D  moments for the historgams
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
