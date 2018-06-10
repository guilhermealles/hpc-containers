#!/bin/bash

if [ -z $1 ]; then
    echo "Please provide a directory to create the images in"
    exit
fi
SINGULARITY_IMAGES_DIRECTORY=$1

[ -d $SINGULARITY_IMAGES_DIRECTORY/ondes3d-mpi ] && rm -rf $SINGULARITY_IMAGES_DIRECTORY/ondes3d
cp -r ../../../software/ondes3 "$SINGULARITY_IMAGES_DIRECTORY/ondes3d"
sudo singularity build --writable "$SINGULARITY_IMAGES_DIRECTORY/ondes3d.img" ../singularity/ondes3d.bootstrap
echo "$(cd $SINGULARITY_IMAGES_DIRECTORY && singularity exec ./ondes3d.img ./ondes3d/compile-ondes3d.sh)"