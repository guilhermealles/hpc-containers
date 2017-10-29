#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

DOCKER_CLUSTER_DIR="$PWD/docker-cluster/"
EXPERIMENT_HOME_DIR="$PWD"
SOFTWARE_UTILS_DIR="$EXPERIMENT_HOME_DIR/../software/utils"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
HOSTFILE_FORCE_COMM="$EXPERIMENT_HOME_DIR/config/hosts-force-comm.txt";
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/ep-results.csv"
NAME=$1
ENVIRONMENT=$2
CONTEXT=$3
PARALLELISM=$4

if [ $CONTEXT = "mpi" ]; then
    ./setup/assemble-swarm.sh create $HOSTFILE
    ./setup/ep-start-docker-cluster.sh up $PARALLELISM
    cd $DOCKER_CLUSTER_DIR
    EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec mpirun -np $PARALLELISM NPB3.3.1/NPB3.3-MPI/bin/ep.B.$PARALLELISM)
elif [ $CONTEXT = "mpi-high-comm" ]; then
    ./setup/assemble-swarm.sh create $HOSTFILE_FORCE_COMM
    ./setup/ep-start-docker-cluster.sh up $PARALLELISM
    cd $DOCKER_CLUSTER_DIR
	EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec mpirun -np $PARALLELISM NPB3.3.1/NPB3.3-MPI/bin/ep.B.$PARALLELISM)
else # OpenMP
    ./setup/assemble-swarm.sh create $HOSTFILE
    ./setup/ep-start-docker-cluster.sh up 1
    cd $DOCKER_CLUSTER_DIR
    EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec OMP_NUM_THREADS=$PARALLELISM NPB3.3.1/NPB3.3-OMP/bin/ep.B.x)
fi

cd "$EXPERIMENT_HOME_DIR"
./setup/ep-start-docker-cluster.sh down 16
./setup/assemble-swarm.sh destroy $HOSTFILE

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$CONTEXT,$PARALLELISM,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$CONTEXT,$PARALLELISM,$EXEC_TIME" >> $RESULTS_FILE