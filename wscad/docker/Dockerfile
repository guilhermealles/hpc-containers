FROM debian:9.4

RUN apt-get update && apt-get upgrade -y
RUN apt-get install build-essential gfortran curl wget sudo -y

RUN adduser --disabled-password --gecos "" alles \
  && echo "alles   ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

WORKDIR /home/alles

# Install MPICH
ARG MPICH_VERSION="3.2"
RUN mkdir mpich
WORKDIR /home/alles/mpich
RUN wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
  && tar xvzf mpich-${MPICH_VERSION}.tar.gz  \
  && cd mpich-${MPICH_VERSION}  \
  && ./configure \
  && make && sudo make install \
  && rm -rf /home/alles/mpich

USER alles
RUN mkdir /home/alles/nas
WORKDIR /home/alles/nas
RUN curl https://www.nas.nasa.gov/assets/npb/NPB3.3.1.tar.gz > NPB-3.3.1.tar.gz
RUN tar -xvzf NPB-3.3.1.tar.gz
RUN cd NPB3.3.1/NPB3.3-MPI/ && cp config/make.def.template config/make.def 
RUN cd NPB3.3.1/NPB3.3-MPI/ && sed -i 's/f77/mpif77/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=1
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=4
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=16

RUN cd NPB3.3.1/NPB3.3-OMP/ && cp config/make.def.template config/make.def
RUN cd NPB3.3.1/NPB3.3-OMP/ && sed -i 's/f77/gfortran/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-OMP/ && sed -i 's/^FFLAGS\s=\s-O$/FFLAGS  = -O -fopenmp/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-OMP/ && sed -i 's/^FLINKFLAGS\s=\s-O$/FLINKFLAGS  = -O -fopenmp/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-OMP/ && sed -i 's/^CFLAGS\s=\s-O$/CFLAGS  = -O -fopenmp/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-OMP/ && sed -i 's/^CLINKFLAGS\s=\s-O$/CLINKFLAGS  = -O -fopenmp/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-OMP/ && make ep CLASS=B