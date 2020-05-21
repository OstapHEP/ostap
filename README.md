Ostap Project                   
=============
[![Build Status](https://travis-ci.org/OstapHEP/ostap.svg?branch=master)](https://travis-ci.org/OstapHEP/ostap)
[![Coverage Status](https://coveralls.io/repos/github/OstapHEP/ostap/badge.svg?branch=master)](https://coveralls.io/github/OstapHEP/ostap?branch=master)
[![Build Status](https://dev.azure.com/OstapHep/OstapHep/_apis/build/status/OstapHEP.ostap?branchName=master)](https://dev.azure.com/OstapHep/OstapHep/_build/latest?definitionId=5&branchName=master)
[![pipeline status](https://gitlab.cern.ch/ostapHep/ostaphep/badges/master/pipeline.svg)](https://gitlab.cern.ch/ostapHep/ostaphep/commits/master)
[![Join the chat at https://gitter.im/OstapHEP/ostap](https://badges.gitter.im/OstapHEP/ostap.svg)](https://gitter.im/OstapHEP/ostap?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/81464356.svg)](https://zenodo.org/badge/latestdoi/81464356)
<!--[![build status](https://gitlab.cern.ch/amazurov/ostap/badges/master/build.svg)](https://gitlab.cern.ch/amazurov/ostap/commits/master)-->

Nowadays [ROOT](http://root.cern.ch/) and [PyROOT](http://root.cern.ch/drupal/content/pyroot) are de-facto standard tools for performing physics analysis. The Ostap project is a community-driven initiative aiming to provide more user friendly and more intuitive interface to [ROOT](http://root.cern.ch/) and [PyROOT](http://root.cern.ch/drupal/content/pyroot) and extending the existing functionality.

Project started in 2009 from the private collections of python functions used in [Kali](http://inspirehep.net/record/1111459) - framework for calibration of LHCb electromagnetic calorimeter. A lot of functionality is picked from [Bender](http://lhcb-release-area.web.cern.ch/LHCb-release-area/DOC/bender/) project - python based physics analysis environment used in LHCb experiemnt. Till Autumn 2016 the project was a part of LHCb software suit and with great success has been used for preparation of approximately 30 physics papers. A standalone, LHCb independent version, has appeared at start of 2017

Key features include:

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
- _optional_ (only for python3) : [bsddb3]
   - python interface to Berkeley DB (`libdb` needs to be installed!)
   
Setup
-----

There are several possibilities to start working with Ostap, you can build Ostap on the Linux or lxplus/7, run the docker container and now the Ostap is available on conda and SWAN.   

The possible  setup options are described  [here](INSTALL.md)

[ROOT]: http://root.cern.ch
[RooFit]: https://root.cern.ch/roofit 
[numpy]: https://numpy.org 
[scipy]: https://www.scipy.org 
[pathos]: https://github.com/uqfoundation/pathos 
[dill]: https://github.com/uqfoundation/dill
[multiprocess]: https://github.com/uqfoundation/multiprocess
[ppt]: https://github.com/uqfoundation/ppft
[multiprocessing]:https://docs.python.org/2/library/multiprocessing.html
[terminaltables]: https://pypi.org/project/terminaltables
[bsddb3]: https://pypi.org/project/bsddb3/
