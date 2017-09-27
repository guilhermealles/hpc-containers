# This is the BASE Image to be used when building container clusters
FROM centos:7

# Install development tools
RUN yum groupinstall -y "Development Tools"
RUN yum install -y compat-gcc-44-gfortran

# Install OpenMPI 2.1.0
RUN curl https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.0.tar.bz2 > openmpi-2.1.0.tar.bz2
RUN tar -xvjf openmpi-2.1.0.tar.bz2
RUN cd openmpi-2.1.0 && ./configure --prefix=/usr/local && make && make install