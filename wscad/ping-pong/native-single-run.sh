#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

LOCAL_PINGPONG_BUILD="$HOME"
EXPERIMENT_HOME_DIR="$PWD"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/pingpong-results.csv"
NAME=$1
ENVIRONMENT=$2
SIZEEXP=$3

EXEC_TIME=$(mpirun --hostfile $HOSTFILE -np 2 $LOCAL_PINGPONG_BUILD/ping-pong-mpi $SIZEEXP)

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME" >> $RESULTS_FILE