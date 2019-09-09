.. include:: ./references.rst

.. _installation:

Installation
============


Setup
-----

There are several possibilities to start working with `ostap`, you can build `ostap` 
on the Linux or lxplus/7 and also run the docker container.  

Linux
-----
`ostap` requires the `python2` version >2.7 or `python3`. 
You need to provide `$ROOTSYS` environment for building the library. 
The library should be built with the same compiler version 
and C++ standard as was used to build ROOT_. 

At CERN, at `lxplus/7` you can do it with several LCG versions  (94,95,97,...). 
Check on which platform is the preferred version of LCG located 
and run `LbLogin`, for instance for LCG 95 and  `x86_64-centos7-gcc8-opt`:

    LbLogin -c x86_64-centos7-gcc8-opt
    source /cvmfs/sft.cern.ch/lcg/views/LCG_95/${CMTCONFIG}/setup.sh

After setting the enviroments clone the latest released version and build Ostap package 

    git clone  git://github.com/OstapHEP/ostap.git
    cd ostap
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<INSTALL_DIRECTORY>
    make -j4
    make install
    source <INSTALL_DIRECTORY>/thisostap.sh 

For the latest tag check `this page <https://github.com/OstapHEP/ostap/releases>`_
To update the package to latest version use following command:

    git pull origin <latest tag>

or to get the head version use:

    git pull origin master

Docker
------
We also provided `Dockerfile` to build `ostap` image.  
You can run `ostap` interactively using the command line or 
via `Docker Desktop` which is also available for `MacOS` and `Windows`. 
To create the docker image from the Ostap directory run:

    sudo docker build --network host -t <dockerID>/ostaphep:latest .

Run the image iteractively

    sudo docker  run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -v ${WORKDIR_PATH}/work_dir:/work_dir  -it <dockerID>/ostaphep:latest

To know more about docker, please check the `documentation. <https://docs.docker.com>`_

Conda
-----
`ostap` is now available on the `conda-forge channel. <https://github.com/conda-forge/ostaphep-feedstock>`_ You can get `ostap` with conda using the following steps:

Install mininconda

    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
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

The list of available version you can find `here. <https://anaconda.org/conda-forge/ostaphep/files>`_
To know more about conda-forge, please visit `conda page. <https://conda-forge.org>`
