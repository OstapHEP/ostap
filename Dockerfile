FROM mazurov/cern-root:latest
MAINTAINER alexander.mazurov <alexander.mazurov@gmail.com>
LABEL description="ostap HEP framework"

RUN yum -y install doxygen && yum -y clean all

ADD . /ostap
WORKDIR /ostap

RUN ./scripts/bootstrap

ENV LD_LIBRARY_PATH "/opt/ostap/build/lib:$LD_LIBRARY_PATH"
ENV PYTHONPATH      "/opt/ostap:$PYTHONPATH"

CMD /bin/bash