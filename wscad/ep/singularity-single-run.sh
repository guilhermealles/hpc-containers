#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

SINGULARITY_IMAGES_DIR="$HOME/singularity_images"
EXPERIMENT_HOME_DIR="$PWD"
SOFTWARE_UTILS_DIR="$EXPERIMENT_HOME_DIR/../software/utils"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/ep-results.csv"
NAME=$1
ENVIRONMENT=$2
CONTEXT=$3
PARALLELISM=$4

if [ $CONTEXT = "mpi" ]; then
    EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE -np $PARALLELISM $SINGULARITY_IMAGES_DIR/mpi-ep-$PARALLELISM.img)
fi

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$CONTEXT,$PARALLELISM,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$CONTEXT,$PARALLELISM,$EXEC_TIME" >> $RESULTS_FILE