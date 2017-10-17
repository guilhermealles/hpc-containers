#!/bin/bash
INPUT_FILE='./ep-configuration.csv'
RESULTS_FILE='./ep-results.csv'
LOCAL_NAS_BUILD="$HOME/NAS"

echo "Running experiments from input CSV file..."
OLDIFS=$IFS
IFS=","
while read -u 11 name order number rp environment context parallelism block exectime
do
    HOSTFILE="./hosts.txt"
    HOSTFILE_FORCE_COMM="./hosts-force-comm.txt";
    if [ $environment = "singularity" ]; then
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE -np $parallelism ./images/mpi-c-$parallelism.img'"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism ./images/mpi-c-$parallelism.img'"
        else #OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh 'OMP_NUM_THREADS=$parallelism ./images/omp-c.img'"
        fi
    elif [ $environment = "docker" ]; then
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh '../../docker/cluster/swarm.sh exec mpirun -np $parallelism NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        elif [ $context = "mpi-high-comm" ]; then
	    EXEC_COMMAND="../../software/utils/ms-time.sh '../../docker/cluster/swarm.sh exec mpirun -np $parallelism NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        else # OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh '../../docker/cluster/swarm.sh exec OMP_NUM_THREADS=$parallelism NPB3.3.1/NPB3.3-OMP/bin/ep.B.x"
        fi
    else # native
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE -np $parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        else # OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh 'OMP_NUM_THREADS=$parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-OMP/bin/ep.B.x'"
        fi
    fi

    # Exclude first line
    if [ $name != "name" ]; then
        echo "$EXEC_COMMAND";
    fi
    sleep 5
done 11< $INPUT_FILE
IFS=$OLDIFS

echo "Done."
