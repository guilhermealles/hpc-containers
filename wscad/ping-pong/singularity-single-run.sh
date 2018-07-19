#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

SINGULARITY_IMAGES_DIR="$HOME/singularity_images"
EXPERIMENT_HOME_DIR="$PWD"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/pingpong-results.csv"
NAME=$1
ENVIRONMENT=$2
SIZEEXP=$3

EXEC_TIME=$(mpirun --hostfile $HOSTFILE -np 2 $SINGULARITY_IMAGES_DIR/ping-pong-sing.img $SIZEEXP)

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME" >> $RESULTS_FILE