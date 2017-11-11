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
RUN rm ./SRC/options.h && cp ./ESSAI-XML/options.h ./SRC
RUN sed -i 's/^\(TESTFLAGS.*+=.*-DMPI\)$/#\1/' ./SRC/Makefile
RUN cd SRC/ && make clean && make
RUN mv ondes3d ondes3d-omp
RUN sed -i 's/^#\(TESTFLAGS.*+=.*-DMPI\)$/\1/' ./SRC/Makefile
RUN sed -i 's/^\(TESTFLAGS.*+=.*-DOMP\)$/#\1/' ./SRC/Makefile
RUN sed -i 's/^\(OMP_FLAGS.*+=.*-fopenmp\)$/#\1/' ./SRC/Makefile
RUN cd SRC/ && make clean && make
RUN mkdir ESSAI-OUTPUT
RUN mv ondes3d ondes3d-mpi
RUN chown -R mpi:mpi *
USER mpi

# ####################
# For Docker beginner:

# After Docker syntax `RUN`, you can execute any command available in 
# the current shell of the image.

# To switch to root:    USER root
# To switch back to default user: USER ${USER}
