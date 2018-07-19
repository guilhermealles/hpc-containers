#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

DOCKER_CLUSTER_DIR="$PWD/docker-cluster/"
EXPERIMENT_HOME_DIR="$PWD"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/pingpong-results.csv"
NAME=$1
ENVIRONMENT=$2
SIZEEXP=$3

cd $DOCKER_CLUSTER_DIR
echo $(./swarm.sh up size=3)
sleep 5
EXEC_TIME=$(./swarm.sh exec mpirun -np 2 /project/ping-pong-mpi $SIZEEXP)
echo $(./swarm.sh down)
sleep 10

cd "$EXPERIMENT_HOME_DIR"
echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME" >> $RESULTS_FILE