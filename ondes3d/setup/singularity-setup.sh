#!/bin/bash

if [ -z $1 ]; then
    echo "Please provide a directory to create the images in"
    exit
fi

SINGULARITY_IMAGES_DIRECTORY=$1

[ ! -d $SINGULARITY_IMAGES_DIRECTORY/ondes3d-mpi ] && mkdir "$SINGULARITY_IMAGES_DIRECTORY/ondes3d_mpi"
cp -R ../docker-cluster/project/* "$SINGULARITY_IMAGES_DIRECTORY/ondes3d_mpi"
singularity create -s 300 "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ondes3d.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ondes3d.img" ../singularity/alpine-mpi-ondes3d.bootstrap

[ ! -d $SINGULARITY_IMAGES_DIRECTORY/ondes3d-omp ] && mkdir "$SINGULARITY_IMAGES_DIRECTORY/ondes3d_omp"
cp -R ../docker-cluster/project/* "$SINGULARITY_IMAGES_DIRECTORY/ondes3d_omp"
singularity create -s 300 "$SINGULARITY_IMAGES_DIRECTORY/alpine-omp-ondes3d.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-omp-ondes3d.img" ../singularity/alpine-omp-ondes3d.bootstrap