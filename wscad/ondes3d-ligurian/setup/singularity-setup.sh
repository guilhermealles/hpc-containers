#!/bin/bash

if [ -z $1 ]; then
    echo "Please provide a directory to create the images in"
    exit
fi
SINGULARITY_IMAGES_DIRECTORY=$1
CURRENT_DIRECTORY=$PWD

[ -d $SINGULARITY_IMAGES_DIRECTORY/ondes3d ] && rm -rf $SINGULARITY_IMAGES_DIRECTORY/ondes3d
cd "$SINGULARITY_IMAGES_DIRECTORY"
git clone https://guilhermealles@bitbucket.org/schnorr/ondes3d-cases.git
mv ondes3d-cases ondes3d
cp "$CURRENT_DIRECTORY/compile-ondes3d.sh" ondes3d/
sudo singularity build --writable "$SINGULARITY_IMAGES_DIRECTORY/ondes3d.img" ../singularity/ondes3d.bootstrap
echo "$(cd $SINGULARITY_IMAGES_DIRECTORY && singularity exec ./ondes3d.img ./ondes3d/compile-ondes3d.sh)"