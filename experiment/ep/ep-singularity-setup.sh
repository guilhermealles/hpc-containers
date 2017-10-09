#!/bin/bash
mkdir ./images

singularity create -s 1000 ./images/mpi-c-16.img
sudo singularity bootstrap ./images/mpi-c-16.img ../../singularity/bootstrap/mpi-ep-c-16.bootstrap

singularity create -s 1000 ./images/mpi-c-4.img
sudo singularity bootstrap ./images/mpi-c-4.img ../../singularity/bootstrap/mpi-ep-c-4.bootstrap

singularity create -s 1000 ./images/mpi-c-1.img
sudo singularity bootstrap ./images/mpi-c-1.img ../../singularity/bootstrap/mpi-ep-c-1.bootstrap

singularity create -s 1000 ./images/omp-c.img
sudo singularity bootstrap ./images/omp-c.img ../../singularity/bootstrap/omp-ep-c.bootstrap