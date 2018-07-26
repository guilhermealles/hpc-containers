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
RUN sudo apt-get update -y && sudo apt-get install -y libxml2-dev libxslt1-dev git
RUN git clone https://gitlab.inria.fr/simgrid/platform-calibration.git
RUN cd platform-calibration/src/calibration && make
RUN cp ./alles.xml platform-calibration/src/calibration/examples 
USER mpi

# ####################
# For Docker beginner:

# After Docker syntax `RUN`, you can execute any command available in 
# the current shell of the image.

# To switch to root:    USER root
# To switch back to default user: USER ${USER}
