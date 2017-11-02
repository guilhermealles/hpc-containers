# hpc-containers
Experiments with HPC applications running inside Docker and Singularity containers.

This repository contains code to build a Docker cluster infrastructure. Highly influenced by [NLNKguyen's work](https://github.com/NLKNguyen/alpine-mpich): Docker images and the multi-host orchestration are kept mostly intact, with some minor tweaks.

Also, based on the same Docker images, there are some Singularity definition files to create singularity containers for performance analysis.

An automated script to run experiments (configured in a .csv file) is also included. Results will be included in the future.
