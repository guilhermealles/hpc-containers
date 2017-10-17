FROM guilhermealles/alpine-mpi:onbuild

# # ------------------------------------------------------------
# # Build MPI project
# # ------------------------------------------------------------

# Put all build steps and additional package installation here

# Note: the current directory is ${WORKDIR:=/project}, which is
# also the default directory where ${USER:=mpi} will SSH login to

# Copy the content of `project` directory in the host machine to 
# the current working directory in this Docker image
COPY project/ .

# Normal build command
USER root
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
RUN chown -R ${USER}:${USER} *
USER mpi

# ####################
# For Docker beginner:

# After Docker syntax `RUN`, you can execute any command available in 
# the current shell of the image.

# To switch to root:    USER root
# To switch back to default user: USER ${USER}
