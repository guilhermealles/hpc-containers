FROM debian:9.4

RUN apt-get update && apt-get upgrade -y
RUN apt-get install build-essential gfortran curl wget sudo -y

RUN adduser -D alles \
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