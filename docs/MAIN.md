Ostap Project                    {#mainpage}
=============

Nowadays [ROOT] and [PyROOT] are de-facto standard tools for performing physics analysis. The Ostap project is a community-driven initiative aiming to provide more user friendly and more intuitive interface to [ROOT], [PyROOT] and [RooFit] and extending the existing functionality.

Project started in 2009 from the private collections of python functions used in [Kali] framework for calibration of LHCb electromagnetic calorimeter. A lot of functionality is picked from [Bender] project, the python based physics analysis environment used in LHCb experiment. Till Autumn 2016 the project was a part of LHCb software suit and with great success has been used for preparation of approximately 30 physics papers. A standalone, LHCb independent version, has appeared at start of 2017

Key features include
--------------------

-   Very easy manipulations with [ROOT] and [RooFit] objects: histograms, trees, datasets, etc
-   Very easy interface to [RooFit] machinery
-   Extended set of `models.PDFs` for [RooFit]
-   Powerful, pickle-based persistency for object
-   Interactive `ostap` analysis environment

Dependencies
------------

- _mandatory_: [ROOT], [RooFit]
- _highly desirable_: [numpy]
   - mandatory for Fast Fourier Transform, used in histogram/function parameterization;
   - optional for some other issues, in particular for the prime number treatment;  
- _optional_: [scipy]
   - numerical integration (quadratures, cubatures), root finding, minimization; 
   - `ostap` offers home-made replacements, but the native methods from `scipy` are more efficient;
- _optional_: [pathos], [dill], [multiprocess] and [ppt]
   - needed for parallel processing; 
   - `ostap` offers a [multiprocessing]-based replacement with reduced functionality; 
- _optional_: [terminaltables]
   - nice format of tables 
         - in particular for nice printout for `TTree`, `TChain`, `RooDataSet`, ... ;
   - `ostap` offers a home-made replacement with a bit reduced functionalty.  


[ROOT]: http://root.cern.ch
[PyROOT]:http://root.cern.ch/drupal/content/pyroot
[RooFit]: https://root.cern.ch/roofit 
[numpy]: https://numpy.org 
[scipy]: https://www.scipy.org 
[pathos]: https://github.com/uqfoundation/pathos 
[dill]: https://github.com/uqfoundation/dill
[multiprocess]: https://github.com/uqfoundation/multiprocess
[ppt]: https://github.com/uqfoundation/ppft
[multiprocessing]:https://docs.python.org/2/library/multiprocessing.html
[terminaltables]: https://pypi.org/project/terminaltables
[Kali]:http://inspirehep.net/record/1111459
[Bender]:http://lhcb-release-area.web.cern.ch/LHCb-release-area/DOC/bender


