FROM centos:7

RUN yum groupinstall -y "Development Tools"
RUN curl https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.bz2 > openmpi-2.1.0.tar.bz2
RUN tar -xvjf openmpi-2.1.0.tar.bz2
RUN cd openmpi-2.1.0
RUN ./configure --prefix=/usr/local
RUN make
RUN make install