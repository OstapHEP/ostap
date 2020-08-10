Setup
-----

There are several options to start working with ostap, you can build ostap on the Linux or lxplus/7, run the docker container and futhermore now the ostap  is available on conda and SWAN services.  

Linux
-----
ostap requires the python2 version >2.7 or python3. 
You need to provide ROOTSYS environment for building the package library. The library should be build with the same compiler version and C++ standard as was the one used to build ROOT. 

After setting the enviroments clone and build ostap package: 

    git clone  git://github.com/OstapHEP/ostap.git
    cd ostap
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<INSTALL_DIRECTORY>
    make -j8
    make install
    source <INSTALL_DIRECTORY>/thisostap.sh 
    
On lxplus/7 you can do it with several LCG versions (95,96,97), using the scripts/setup.sh. Check the location of the preffered LCG version. For instance for LCG 97 and  x86_64-centos7-gcc8-opt:

    lb-set-platform x86_64-centos7-gcc8-opt
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97/${CMTCONFIG}/setup.sh
    ./scripts/setup.sh
    source LCG_$LCG_VERSION/INSTALL/thisostap.sh


To update the package to latest version:

    git pull origin <latest tag>
or to get the head version:

    git pull origin master
For the latest tag check the page https://github.com/OstapHEP/ostap/releases

Docker
-----
We also provide Dockerfile to build the ostap image. You can run Ostap interactively using the command line or via Docker Desktop which is available on MacOS and Windows. To create the docker image from the Ostap directory run:

    sudo docker build --network host -t <dockerID>/ostaphep:latest .
Run the image iteractively:

    sudo docker  run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix  -v ${WORKDIR_PATH}/work_dir:/work_dir  -it <dockerID>/ostaphep:latest
For the further info on docker, please check the documentation: https://docs.docker.com/.

Сonda
-----
ostap is now available on the conda-forge channel https://github.com/conda-forge/ostaphep-feedstock. To get ostap with conda:

Install mininconda:

    wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh
    
Set conda environments:

    source path-to-miniconda/etc/profile.d/conda.sh
    
Add conda-forge to your channels:

    conda config --add channels conda-forge
    
Check available versions: 

    conda search ostaphep --channel conda-forge

Create an environment with a specific version of python:

    conda create --name ostap-env ostaphep python=3.7

To activate or deactivate ostap's environment:

    conda activate ostap-env 
    conda deactivate  
To get the latest version of the ostap package:

    conda install -n ostap-env  ostaphep=version
The list of available version you can find [here](https://anaconda.org/conda-forge/ostaphep/files).
To know more about conda-forge, please visit [conda page](https://conda-forge.org).

Additionally you can install the master version of the ostap package using dependancy packages which are also availible on  conda:

Set conda environments:

    source path-to-miniconda/etc/profile.d/conda.sh
    
Add conda-forge to your channels:

    conda config --add channels conda-forge
Create an environment with required dependancies:

    conda create -n ostap-req-env root_base root-binaries root-dependencies gsl  future configparser  numpy scipy pathos dill multiprocess ppft terminaltables binutils-meta c-compiler compilers cxx-compiler fortran-compiler python ipython cmake

If you computer has installed Berkeley-DB (`libdb`) and you are using `python3`, it is desirable also  to add `bsddb3` in the list

Activate the  environment  with requirement packages:

    conda activate ostap-req-env 

Clone and build ostap package:

    git clone  git://github.com/OstapHEP/ostap.git
    cd ostap
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<INSTALL_DIRECTORY>
    make -j8
    make install
    source <INSTALL_DIRECTORY>/thisostap.sh 

Note that the latest command should be used every time when the ostap-req-env is activated again.

SWAN
-----
The ostap installation uses the SWAN CERN service (described  [here](SWAN.md).)
