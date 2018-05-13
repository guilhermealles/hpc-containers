FROM debian:9.4

RUN apt-get install gfortran curl wget

# Install MPICH
ARG MPICH_VERSION="3.2"
RUN mkdir /mpich
WORKDIR /mpich
RUN wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
  && tar xvzf mpich-${MPICH_VERSION}.tar.gz  \
  && cd mpich-${MPICH_VERSION}  \
  && ./configure \
  && make && make install \
  && rm -rf /mpich