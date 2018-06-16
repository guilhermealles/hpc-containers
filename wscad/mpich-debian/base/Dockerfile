FROM debian:9.4

RUN apt-get update && apt-get upgrade -y

# Install basic development tools
RUN apt-get install build-essential gfortran curl wget openssh-server sudo net-tools -y

# Install MPICH
ARG MPICH_VERSION="3.2"
RUN mkdir /mpich
WORKDIR /mpich
RUN wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz \
    && tar xvzf mpich-${MPICH_VERSION}.tar.gz  \
    && cd mpich-${MPICH_VERSION}  \
    && ./configure ${MPICH_CONFIGURE_OPTIONS}  \
    && make ${MPICH_MAKE_OPTIONS} && make install \
    && rm -rf /mpich

# Clean up 
WORKDIR /
RUN rm -rf /tmp/*

# Add default user
ARG USER=mpi
ENV USER=${USER}
RUN adduser --disabled-password --gecos "" ${USER} \
    && echo "${USER}   ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
ENV USER_HOME /home/${USER}
RUN chown -R ${USER}:${USER} ${USER_HOME}

# Create working directory for user
ARG WORKDIR=/project
ENV WORKDIR ${WORKDIR}
RUN mkdir ${WORKDIR}
RUN chown -R ${USER}:${USER} ${WORKDIR}

WORKDIR ${WORKDIR}
USER ${USER}
