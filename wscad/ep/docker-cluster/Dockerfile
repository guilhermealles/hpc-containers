FROM guilhermealles/mpich-debian:onbuild

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
USER mpi
RUN curl https://www.nas.nasa.gov/assets/npb/NPB3.3.1.tar.gz > NPB-3.3.1.tar.gz
RUN tar -xvzf NPB-3.3.1.tar.gz
RUN cd NPB3.3.1/NPB3.3-MPI/ && cp config/make.def.template config/make.def 
RUN cd NPB3.3.1/NPB3.3-MPI/ && sed -i 's/f77/mpif77/g' config/make.def
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=1
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=4
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=8
RUN cd NPB3.3.1/NPB3.3-MPI/ && make ep CLASS=B NPROCS=16

# ####################
# For Docker beginner:

# After Docker syntax `RUN`, you can execute any command available in 
# the current shell of the image.

# To switch to root:    USER root
# To switch back to default user: USER ${USER}
