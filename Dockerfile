FROM gitlab-registry.cern.ch/lhcb-docker/os-base/centos7-devel:latest
MAINTAINER tatiana.ovsiannikova <tatiana.ovsiannikova@cern.ch>
LABEL description="ostap HEP framework"

RUN #!/bin/bash
RUN  yum  install -y git wget
RUN wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /root/miniconda
ENV PATH="/root/miniconda/bin:${PATH}"
RUN echo $PATH
RUN conda config --set always_yes yes --set changeps1 no
RUN conda config --add channels conda-forge
RUN conda create -q -n ostapenv python=3.7 coverage coveralls cmake nose ninja ipython root_base=6.20 root-binaries root-dependencies gsl future configparser numpy scipy pathos dill multiprocess ppft terminaltables gdbm libdb bsddb3 psutil more-itertools 
ADD . /ostap
WORKDIR /ostap

ENV PATH="/root/miniconda/envs/ostapenv/bin:${PATH}"

RUN  mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=./INSTALL/ && make -j12 && make install && echo "source build/INSTALL/thisostap.sh" >> ~/.bashrc

CMD /bin/bash
