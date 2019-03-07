Ostap Project
=============
<!--[![Build Status](https://travis-ci.org/OstapHEP/ostap.svg?branch=master)](https://travis-ci.org/OstapHEP/ostap)-->
[![Build Status](https://dev.azure.com/OstapHep/OstapHep/_apis/build/status/OstapHEP.ostap%20(1)?branchName=master)](https://dev.azure.com/OstapHep/OstapHep/_build/latest?definitionId=7&branchName=master)
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

There are several possibilities to start working with Ostap, you can build Ostap on the Linux or lxplus/7 and also run the docker container.  

Linux
-----
Ostap requires the python2 version >2.7. 
You need to provide ROOTSYS environment for building the library. The library should be built with the same compiler version and C++ standard as was used to build ROOT. Set ROOT environment with following commands

    source /path/to/root/bin/thisroot.sh

or, if you are working at the LHCb environment

    lb-run ROOT bash

e.g. at lxplus/7 one can do 

    lb-run --ext Python --ext pytools --ext pyanalysis --ext ROOT LCG/95 bash --norc
then clone the repo and build Ostap package 

    git clone git://github.com/OstapHEP/ostap.git
    cd ostap
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<INSTALL_DIRECTORY>
    make -j8
    make install
    source <INSTALL_DIRECTORY>/thisostap.sh 
Docker
-----
We also provided Dockerfile to build the OstapHep image.  You can run Ostap interactively using the command line or via Docker Desktop which is available for MacOS and Windows. To create the docker image from the Ostap directory run:

    sudo docker build --network host -t <dockerID>/ostaphep:latest .
Run the image iteractively:

    sudo docker  run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -v ${WORKDIR_PATH}/work_dir:/work_dir  -it <dockerID>/ostaphep:latest
To know more about docker, please check the documentation: https://docs.docker.com/.
