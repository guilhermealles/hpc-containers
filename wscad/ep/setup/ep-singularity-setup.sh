#!/bin/bash

if [ -z $1 ]; then
    echo "Please provide a directory to create the images in"
    exit
fi

SINGULARITY_IMAGES_DIRECTORY=$1

sudo singularity build --writable "$SINGULARITY_IMAGES_DIRECTORY/ep-16.img" ../singularity/mpi-ep-16.bootstrap
sudo singularity build --writable "$SINGULARITY_IMAGES_DIRECTORY/ep-8.img" ../singularity/mpi-ep-8.bootstrap
sudo singularity build --writable "$SINGULARITY_IMAGES_DIRECTORY/ep-4.img" ../singularity/mpi-ep-4.bootstrap
sudo singularity build --writable "$SINGULARITY_IMAGES_DIRECTORY/ep-1.img" ../singularity/mpi-ep-1.bootstrap