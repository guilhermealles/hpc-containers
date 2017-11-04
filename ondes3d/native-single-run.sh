#!/bin/bash

if [ -z $1 ] || [ -z $2 ] || [ -z $3 ] || [ -z $4 ]; then
    echo "Invalid parameters for $0. Please double check."
    exit
fi

LOCAL_NAS_BUILD="$HOME/NAS"
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
    EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE -np $PARALLELISM $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$PARALLELISM)
elif [ $CONTEXT = "mpi-high-comm" ]; then
    EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE_FORCE_COMM -np $PARALLELISM $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$PARALLELISM)
else # OpenMP
    export OMP_NUM_THREADS=$PARALLELISM
    EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-OMP/bin/ep.B.x)
    export OMP_NUM_THREADS=''
fi

echo ">>>>>>>>>> $NAME,$ENVIRONMENT,$CONTEXT,$PARALLELISM,$EXEC_TIME <<<<<<<<<<"
echo "$NAME,$ENVIRONMENT,$CONTEXT,$PARALLELISM,$EXEC_TIME" >> $RESULTS_FILE