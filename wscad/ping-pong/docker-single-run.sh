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

echo $(./setup/ep-start-docker-cluster.sh up 3)
cd $DOCKER_CLUSTER_DIR
EXEC_TIME=$(./swarm.sh exec mpirun -np 2 /project/ping-pong-mpi $SIZEEXP)

cd "$EXPERIMENT_HOME_DIR"
echo $(./setup/ep-start-docker-cluster.sh down 16)

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$SIZEEXP,$EXEC_TIME" >> $RESULTS_FILE