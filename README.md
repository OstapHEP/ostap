Ostap Project
=============
[![Build Status](https://travis-ci.org/OstapHEP/ostap.svg?branch=master)](https://travis-ci.org/OstapHEP/ostap)
[![Build Status](https://dev.azure.com/OstapHep/OstapHep/_apis/build/status/OstapHEP.ostap?branchName=master)](https://dev.azure.com/OstapHep/OstapHep/_build/latest?definitionId=5&branchName=master)
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


Dependencies
------------
- _mandatory_: [`ROOT/PyROOT`](https://root.cern.ch)
- _highly desirable_: [`numpy`](https://numpy.org)
   - mandatory for Fast Fourier Transform, used in histogram/function parameterization;
   - optional for some other issues, in particulat for prime number treatment;  
- _optional_: [`scipy`](https://www.scipy.org)
   - numerical integration (quadratures, cubatures), root finding, minimization; 
   - `ostap` offers home-made replacements, but the native methods from `scipy` are more efficient;
- _optional_: [`pathos`](https://github.com/uqfoundation/pathos), [`dill`](https://github.com/uqfoundation/dill), [`multiprocess`](https://github.com/uqfoundation/multiprocess) and [`ppt`](https://github.com/uqfoundation/ppft)
   - needed for parallel processing; 
   - `ostap` offers [`multiprocessing`](https://docs.python.org/2/library/multiprocessing.html)-based replacement with reduced functionality; 
- _optional_: [`terminaltables`](https://pypi.org/project/terminaltables) 
   - nice format of tables (in particular for nice printout for `ROOT.TTree`, `ROOT.TChain`, `ROOT.RooDataSet`, ...);
   - `ostap` offers home-made replacement with a bit reduced functionalty.  

Setup
-----

There are several possibilities to start working with Ostap, you can build Ostap on the Linux or lxplus/7 and also run the docker container.  

Linux
-----
Ostap requires the python2 version >2.7 or python3. 
You need to provide ROOTSYS environment for building the library. The library should be built with the same compiler version and C++ standard as was used to build ROOT. 

At lxplus/7 you can do it with several LCG  versions  (94,95). Check on which platform is the preferred version of LCG located and run LbLogin, for instance for LCG 95 and  x86_64-centos7-gcc8-opt:

    LbLogin -c x86_64-centos7-gcc8-opt
    source /cvmfs/sft.cern.ch/lcg/views/LCG_95/${CMTCONFIG}/setup.sh

After setting the enviroments clone the latest released version and build Ostap package 

    git clone —-branch <latest tag> git://github.com/OstapHEP/ostap.git
    cd ostap
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<INSTALL_DIRECTORY>
    make -j8
    make install
    source <INSTALL_DIRECTORY>/thisostap.sh 
For the latest tag check the page https://github.com/OstapHEP/ostap/releases
To update the package to latest version use following command:

    git pull origin <latest tag>
or to get the head version use:

    git pull origin master

Docker
-----
We also provided Dockerfile to build the OstapHep image.  You can run Ostap interactively using the command line or via Docker Desktop which is available for MacOS and Windows. To create the docker image from the Ostap directory run:

    sudo docker build --network host -t <dockerID>/ostaphep:latest .
Run the image iteractively

    sudo docker  run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -v ${WORKDIR_PATH}/work_dir:/work_dir  -it <dockerID>/ostaphep:latest
To know more about docker, please check the documentation: https://docs.docker.com/.

Сonda
-----
OstapHep is now available on the conda-forge channel https://github.com/conda-forge/ostaphep-feedstock. You can get ostap with conda using the following steps:

Install mininconda

    wget http://repo.continuum.io.miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh
    
Set the conda environments

    source path-to-miniconda/etc/profile.d/conda.sh
    
Add the conda-forge to your channels

    conda config --add channels conda-forge
    
Check which versions are available

    conda search ostaphep --channel conda-forge

Create an environment with a specific version of python

    conda create --name ostap-env ostaphep python=3.7

To activate or deactivate the ostap environment use the following command

    conda activate ostap-env 
    conda deactivate  
The list of available version you can find here: https://anaconda.org/conda-forge/ostaphep/files.
To know more about conda-forge, please visit conda page: https://conda-forge.org.

