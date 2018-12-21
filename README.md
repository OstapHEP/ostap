Ostap Project
=============
<!--[![Build Status](https://travis-ci.org/OstapHEP/ostap.svg?branch=master)](https://travis-ci.org/OstapHEP/ostap)-->
[![Build Status](https://dev.azure.com/OstapHep/OstapHep/_apis/build/status/OstapHEP.ostap)](https://dev.azure.com/OstapHep/OstapHep/_build/latest?definitionId=5)
[![Join the chat at https://gitter.im/OstapHEP/ostap](https://badges.gitter.im/OstapHEP/ostap.svg)](https://gitter.im/OstapHEP/ostap?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![DOI](https://zenodo.org/badge/81464356.svg)](https://zenodo.org/badge/latestdoi/81464356)

<!--[![build status](https://gitlab.cern.ch/amazurov/ostap/badges/master/build.svg)](https://gitlab.cern.ch/amazurov/ostap/commits/master)-->

Nowadays [ROOT](http://root.cern.ch/) and [PyROOT](http://root.cern.ch/drupal/content/pyroot) are de-facto standard tools for performing physics analysis. The Ostap project is a community-driven initiative aiming to provide more user friendly and more intuitive interface to [ROOT](http://root.cern.ch/) and [PyROOT](http://root.cern.ch/drupal/content/pyroot) and extending the existing functionality.

Project started in 2009 from the private collections of python functions used in [Kali](http://inspirehep.net/record/1111459) - framework for calibration of LHCb electromagnetic calorimeter. A lot of functionality is picked from [Bender](http://lhcb-release-area.web.cern.ch/LHCb-release-area/DOC/bender/) project - python based physics analysis environment used in LHCb experiemnt. Till Autumn 2016 the project was a part of LHCb software suit and with great success has been used for preparation of approximately 30 physics papers. A standalone, LHCb independent version, has appeared at start of 2017

Key features include:

-   Very easy manipulations with [ROOT](http://root.cern.ch/)  and [RooFit](https://root.cern.ch/roofit) objects: histograms, trees, datasets, etc
-   Very easy interface to [RooFit](https://root.cern.ch/roofit) machinery
-   Extended set of models.PDFs for [RooFit](https://root.cern.ch/roofit)
-   Powerful, pickle-based persistency for object
-   Interactive `ostap` analysis environment

Setup
-----

You need to provide ROOTSYS environment for building the library. Library should be built with the same compiler version as was used to build ROOT.

    source /path/to/root/bin/thisroot.sh

or, if you are working at the LHCb environment

    lb-run ROOT bash

e.g. at lxplus/7 one can do 

    lb-run --ext Python --ext pytools --ext pyanalysis --ext ROOT LCG/93 bash --norc
    git clone git://github.com/OstapHEP/ostap.git
    cd ostap
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<INSTALL_DIRECTORY>
    make -j8
    make install
    source <INSTALL_DIRECTORY>/thisostap.sh 

