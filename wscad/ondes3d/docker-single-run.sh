#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

DOCKER_CLUSTER_DIR="$PWD/docker-cluster/"
EXPERIMENT_HOME_DIR="$PWD"
SOFTWARE_UTILS_DIR="$EXPERIMENT_HOME_DIR/../../software/utils"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
SWARM_HOSTFILE="$EXPERIMENT_HOME_DIR/setup/swarm_hosts.txt"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/doe-results.csv"
NAME=$1
ENVIRONMENT=$2
PARALLELISM=$3

if [ $PARALLELISM = "4" ]; then
    MPI_NODES="1"
    THREADS_COUNT="4"
elif [ $PARALLELISM = "8" ]; then 
    MPI_NODES="2"
    THREADS_COUNT="4"
elif [ $PARALLELISM = "16" ]; then
    MPI_NODES="4"
    THREADS_COUNT="4"
else # parallelism = 1
    MPI_NODES="1"
    THREADS_COUNT="1"
fi



cd $DOCKER_CLUSTER_DIR
EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec mpirun -np $MPI_NODES ./ondes3d $THREADS_COUNT)

cd "$EXPERIMENT_HOME_DIR"

sleep 10

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$PARALLELISM,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$PARALLELISM,$EXEC_TIME" >> $RESULTS_FILE
