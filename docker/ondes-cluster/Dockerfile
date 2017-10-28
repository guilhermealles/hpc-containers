FROM guilhermealles/alpine-mpi:onbuild

# # ------------------------------------------------------------
# # Build MPI project
# # ------------------------------------------------------------

# Put all build steps and additional package installation here
COPY project/ .
# Note: the current directory is ${WORKDIR:=/project}, which is
# also the default directory where ${USER:=mpi} will SSH login to

# Copy the content of `project` directory in the host machine to 
# the current working directory in this Docker image

# Normal build command
USER root
RUN rm ./SRC/options.h && cp ./ESSAI-XML/options.h ./SRC
RUN sed -i 's/^\(TESTFLAGS.*+=.*-DMPI\)$/#\1/' ./SRC/Makefile
RUN cd SRC/ && make clean && make
RUN cp ondes3d ondes3d-omp
RUN sed -i 's/^#\(TESTFLAGS.*+=.*-DMPI\)$/\1/' ./SRC/Makefile
RUN sed -i 's/^\(TESTFLAGS.*+=.*-DOMP\)$/#\1/' ./SRC/Makefile
RUN sed -i 's/^\(OMP_FLAGS.*+=.*-fopenmp\)$/#\1/' ./SRC/Makefile
RUN cd SRC/ && make clean && make
RUN mkdir ESSAI-OUTPUT
RUN cp ondes3d ondes3d-mpi
USER mpi