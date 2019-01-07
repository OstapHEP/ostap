FROM tovsiann/ostap-base:latest 
MAINTAINER tatiana.ovsiannikova <tatiana.ovsiannikova@cern.ch>
LABEL description="ostap HEP framework"

RUN yum -y install doxygen && yum -y clean all

ADD . /ostap
WORKDIR /ostap

RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py" && alias python="/usr/local/bin/python2.7" && python get-pip.py && pip install ipython 

RUN alias python="/usr/local/bin/python2.7" && export LD_LIBRARY_PATH=/usr/local/lib64  && source /usr/src/root-6.14.06/builddir/bin/thisroot.sh && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=./INSTALL/ && make -j12 && make install 

RUN echo "alias python="/usr/local/bin/python2.7"" >> ~/.bashrc && echo "export LD_LIBRARY_PATH=/usr/local/lib64" >> ~/.bashrc && echo "source /usr/src/root-6.14.06/builddir/bin/thisroot.sh" >> ~/.bashrc && echo "source build/INSTALL/thisostap.sh" >> ~/.bashrc

CMD /bin/bash
