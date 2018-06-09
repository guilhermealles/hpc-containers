#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

DOCKER_CLUSTER_DIR="$PWD/docker-cluster/"
EXPERIMENT_HOME_DIR="$PWD"
SOFTWARE_UTILS_DIR="$EXPERIMENT_HOME_DIR/../software/utils"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
SWARM_HOSTFILE="$EXPERIMENT_HOME_DIR/setup/swarm_hosts.txt"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/doe-results.csv"
NAME=$1
ENVIRONMENT=$2
PARALLELISM=$3

echo $(./setup/assemble-swarm.sh create $SWARM_HOSTFILE)
echo $(./setup/start-docker-cluster.sh up $PARALLELISM)
cd $DOCKER_CLUSTER_DIR
EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec mpirun -np $PARALLELISM ./ondes3d-mpi)

cd "$EXPERIMENT_HOME_DIR"
echo $(./setup/start-docker-cluster.sh down 16)
echo $(./setup/assemble-swarm.sh destroy $SWARM_HOSTFILE)

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$PARALLELISM,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$PARALLELISM,$EXEC_TIME" >> $RESULTS_FILE