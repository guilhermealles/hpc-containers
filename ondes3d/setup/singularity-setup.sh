#!/bin/bash

if [ -z $1 ]; then
    echo "Please provide a directory to create the images in"
    exit
fi
SINGULARITY_IMAGES_DIRECTORY=$1
SCRIPT_DIRECTORY="$(dirname "$0")"

cd $SCRIPT_DIRECTORY
[ -d $SINGULARITY_IMAGES_DIRECTORY/ondes3d-mpi ] && rm -rf $SINGULARITY_IMAGES_DIRECTORY/ondes3d-mpi
cp -r ../docker-cluster/project "$SINGULARITY_IMAGES_DIRECTORY/ondes3d-mpi"
singularity create -s 300 "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ondes3d.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-mpi-ondes3d.img" ../singularity/alpine-mpi-ondes3d.bootstrap
cd $SINGULARITY_IMAGES_DIRECTORY
chmod +x ./ondes3d-mpi/compile-mpi.sh
singularity exec alpine-mpi-ondes3d.img ./ondes3d-mpi/compile-mpi.sh

cd $SCRIPT_DIRECTORY
[ -d $SINGULARITY_IMAGES_DIRECTORY/ondes3d-omp ] && rm -rf $SINGULARITY_IMAGES_DIRECTORY/ondes3d-omp
cp -r ../docker-cluster/project "$SINGULARITY_IMAGES_DIRECTORY/ondes3d-omp"
singularity create -s 300 "$SINGULARITY_IMAGES_DIRECTORY/alpine-omp-ondes3d.img"
sudo singularity bootstrap "$SINGULARITY_IMAGES_DIRECTORY/alpine-omp-ondes3d.img" ../singularity/alpine-omp-ondes3d.bootstrap
cd $SINGULARITY_IMAGES_DIRECTORY
chmod +x ./ondes3d-omp/compile-omp.sh
singularity exec alpine-omp-ondes3d.img ./ondes3d-omp/compile-omp.sh