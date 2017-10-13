FROM alpine:3.4

# Install basic development tools
RUN apk add --no-cache sudo build-base gfortran curl

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