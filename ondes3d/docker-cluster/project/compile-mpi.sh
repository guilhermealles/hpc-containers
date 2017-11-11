#!/bin/sh

rm ./SRC/options.h && cp ./ESSAI-XML/options.h ./SRC
sed -i 's/^#\(TESTFLAGS.*+=.*-DMPI\)$/\1/' ./SRC/Makefile
sed -i 's/^\(TESTFLAGS.*+=.*-DOMP\)$/#\1/' ./SRC/Makefile
sed -i 's/^\(.*-fopenmp\)$/#\1/' ./SRC/Makefile
cd SRC/ && make clean && make
cd .. && mv ondes3d ondes3d-mpi