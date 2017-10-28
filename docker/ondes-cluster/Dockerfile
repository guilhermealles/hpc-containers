FROM guilhermealles/alpine-mpi:onbuild

# # ------------------------------------------------------------
# # Build MPI project
# # ------------------------------------------------------------

# Put all build steps and additional package installation here

# Note: the current directory is ${WORKDIR:=/project}, which is
# also the default directory where ${USER:=mpi} will SSH login to

# Copy the content of `project` directory in the host machine to 
# the current working directory in this Docker image
RUN mkdir /project
WORKDIR /project

# Normal build command
USER mpi
RUN git clone https://bitbucket.org/fdupros/ondes3d.git
RUN cd ondes3d && cp ./ESSAI-XML/options.h ./SRC
RUN sed -i 's/^\(TESTFLAGS.*+=.*-DMPI\)$/#\1/' ./ondes3d/SRC/Makefile
RUN cd ondes3d/SRC/ && make clean && make
RUN cp ondes3d/ondes3d ondes3d/ondes3d-omp
RUN sed -i 's/^#\(TESTFLAGS.*+=.*-DMPI\)$/\1/' ./ondes3d/SRC/Makefile
RUN sed -i 's/^\(TESTFLAGS.*+=.*-DOMP\)$/#\1/'
RUN sed -i 's/^\(OMP_FLAGS.*+=.*-fopenmp\)$/#\1/'
RUN cd ondes3d/SRC/ && make clean && make
RUN cp ondes3d/ondes3d ondes3d/ondes3d-mpi