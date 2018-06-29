#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

LOCAL_ONDES3D_BUILD="$HOME/ondes3d"
EXPERIMENT_HOME_DIR="$PWD"
SOFTWARE_UTILS_DIR="$EXPERIMENT_HOME_DIR/../../software/utils"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/doe-results.csv"
NAME=$1
ENVIRONMENT=$2
PARALLELISM=$3

cd "$LOCAL_ONDES3D_BUILD"
EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE -np $PARALLELISM ./ondes3d)

cd "$EXPERIMENT_HOME_DIR"

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$PARALLELISM,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$PARALLELISM,$EXEC_TIME" >> $RESULTS_FILE
