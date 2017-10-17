#!/bin/bash
mkdir ./images

singularity create -s 500 ./images/mpi-c-16.img
sudo singularity bootstrap ./images/mpi-c-16.img ../../singularity/bootstrap/alpine-mpi-ep-16.bootstrap

singularity create -s 500 ./images/mpi-c-4.img
sudo singularity bootstrap ./images/mpi-c-4.img ../../singularity/bootstrap/alpine-mpi-ep-4.bootstrap

singularity create -s 500 ./images/mpi-c-1.img
sudo singularity bootstrap ./images/mpi-c-1.img ../../singularity/bootstrap/alpine-mpi-ep-1.bootstrap

singularity create -s 500 ./images/omp-c.img
sudo singularity bootstrap ./images/omp-c.img ../../singularity/bootstrap/alpine-omp-ep.bootstrap
