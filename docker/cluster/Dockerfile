FROM guilhermealles/centos-openmpi:base

RUN yum install -y net-tools
RUN yum -y install openssh-server openssh-clients

#### ADD DEFAULT USER ####
ARG USER=mpi
ENV USER ${USER}
RUN adduser ${USER} \
    && echo "${USER}   ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

ENV USER_HOME /home/${USER}
RUN chown -R ${USER}:${USER} ${USER_HOME}

#### CREATE WORKING DIRECTORY FOR USER ####
ARG WORKDIR=/project
ENV WORKDIR ${WORKDIR}
RUN mkdir ${WORKDIR}
RUN chown -R ${USER}:${USER} ${WORKDIR}

WORKDIR ${WORKDIR}

# # ------------------------------------------------------------
# # Utility shell scripts
# # ------------------------------------------------------------

COPY mpi_bootstrap /usr/local/bin/mpi_bootstrap
RUN chmod +x /usr/local/bin/mpi_bootstrap

COPY get_hosts /usr/local/bin/get_hosts
RUN chmod +x /usr/local/bin/get_hosts

COPY auto_update_hosts /usr/local/bin/auto_update_hosts
RUN chmod +x /usr/local/bin/auto_update_hosts

# # ------------------------------------------------------------
# # Miscellaneous setup for better user experience
# # ------------------------------------------------------------

# Set welcome message to display when user ssh login
COPY welcome.txt /etc/motd

# Default hostfile location for mpirun. This file will be updated automatically.
ENV HYDRA_HOST_FILE /etc/opt/hosts
RUN echo "export HYDRA_HOST_FILE=${HYDRA_HOST_FILE}" >> /etc/profile

RUN touch ${HYDRA_HOST_FILE}
RUN chown ${USER}:${USER} ${HYDRA_HOST_FILE}

# Auto go to default working directory when user ssh login
RUN echo "cd $WORKDIR" >> ${USER_HOME}/.profile

# # ------------------------------------------------------------
# # Set up SSH Server
# # ------------------------------------------------------------
# Add host keys
RUN cd /etc/ssh/ && ssh-keygen -A -N ''

# Config SSH Daemon
RUN  sed -i "s/#PasswordAuthentication.*/PasswordAuthentication no/g" /etc/ssh/sshd_config \
    && sed -i "s/#PermitRootLogin.*/PermitRootLogin no/g" /etc/ssh/sshd_config \
    && sed -i "s/#AuthorizedKeysFile/AuthorizedKeysFile/g" /etc/ssh/sshd_config

# Unlock non-password USER to enable SSH login
RUN passwd -f -u ${USER}

# Set up user's public and private keys
ENV SSHDIR ${USER_HOME}/.ssh
RUN mkdir -p ${SSHDIR}

# Default ssh config file that skips (yes/no) question when first login to the host
RUN echo "StrictHostKeyChecking no" > ${SSHDIR}/config
# This file can be overwritten by the following onbuild step if ssh/ directory has config file

USER root

COPY ssh/ ${SSHDIR}/

RUN cat ${SSHDIR}/*.pub >> ${SSHDIR}/authorized_keys
RUN chmod -R 600 ${SSHDIR}/* \
    && chown -R ${USER}:${USER} ${SSHDIR}

# Switch back to default user when continue the build process
USER ${USER}

# # ------------------------------------------------------------
# # Build MPI project
# # ------------------------------------------------------------

# Put all build steps and additional package installation here

# Note: the current directory is ${WORKDIR:=/project}, which is
# also the default directory where ${USER:=mpi} will SSH login to

# Install NAS Benchmarks
RUN curl https://www.nas.nasa.gov/assets/npb/NPB3.3.1.tar.gz > NPB-3.3.1.tar.gz
RUN tar -xvzf NPB-3.3.1.tar.gz
RUN cp NPB3.3.1/NPB3.3-MPI/config/make.def.template NPB3.3.1/NPB3.3-MPI/config/make.def
RUN sed -i 's/f77/mpif77/g' NPB3.3.1/NPB3.3-MPI/config/make.def
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=C NPROCS=16 && cp bin/ep.C.16 /project/ep.C.16