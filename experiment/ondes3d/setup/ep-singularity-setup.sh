#!/bin/bash
mkdir ./images

singularity create -s 500 ./images/ondes-mpi-16.img
sudo singularity bootstrap ./images/ondes-mpi-16.img ../../singularity/bootstrap/alpine-mpi-ondes3d-16.bootstrap

singularity create -s 500 ./images/ondes-mpi-4.img
sudo singularity bootstrap ./images/ondes-mpi-4.img ../../singularity/bootstrap/alpine-mpi-ondes3d-4.bootstrap

singularity create -s 500 ./images/ondes-mpi-1.img
sudo singularity bootstrap ./images/ondes-mpi-1.img ../../singularity/bootstrap/alpine-mpi-ondes3d-1.bootstrap

singularity create -s 500 ./images/ondes-omp.img
sudo singularity bootstrap ./images/ondes-omp.img ../../singularity/bootstrap/alpine-omp-ondes3d.bootstrap
