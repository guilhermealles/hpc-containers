#!/bin/bash

if [ -z $1 ]; then
    echo "Please provide a directory to create the images in"
    exit
fi

SINGULARITY_IMAGES_DIRECTORY=$1

singularity create -s 250 "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ep-16.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ep-16.img" ../singularity/alpine-mpi-ep-16.bootstrap

singularity create -s 250 "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ep-4.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ep-4.img" ../singularity/alpine-mpi-ep-4.bootstrap

singularity create -s 250 "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ep-1.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ep-1.img" ../singularity/alpine-mpi-ep-1.bootstrap

singularity create -s 250 "$SINGULARITY_IMAGES_DIRECTORY/alpine-omp-ep.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-omp-ep.img" ../singularity/alpine-omp-ep.bootstrap
