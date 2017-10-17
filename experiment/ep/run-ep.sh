#!/bin/bash
INPUT_FILE='./ep-configuration.csv'

echo "Creating singularity containers..."
echo "$(./ep-singularity-setup)"
echo "Singularity containers created successfully."


echo "Running experiments from input CSV file..."
OLDIFS=$IFS
IFS=","
while read name order number rp environment context parallelism block exectime
do
    HOSTFILE="./hosts.txt"
    HOSTFILE_FORCE_COMM="./hosts-force-comm.txt";
    if [ $environment = "singularity" ]; then
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE -np $parallelism ./images/mpi-c-$parallelism.img > /dev/null 2&>1'"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism ./images/mpi-c-$parallelism.img > /dev/null 2&>1'"
        else #OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh 'OMP_NUM_THREADS=$parallelism ./images/omp-c.img > /dev/null 2>&1'"
        fi
    elif [ $environment = "docker" ]; then
        if [ $context = "mpi" ]; then
            ./setup/assemble-swarm.sh create $HOSTFILE
            ./setup/ep-start-docker-cluster.sh up size=$parallelism
            EXEC_COMMAND="../../software/utils/ms-time.sh '../../docker/cluster/swarm.sh exec run_ep_$parallelism.sh'"
        elif [ $context = "mpi-high-comm" ]; then
            ./setup/assemble-swarm.sh create $HOSTFILE_FORCE_COMM
            ./setup/ep-start-docker-cluster.sh up size=$parallelism
            EXEC_COMMAND="../../software/utils/ms-time.sh docker mpi-high-comm ep c 16"
        else # OpenMP
            ./setup/assemble-swarm.sh create $HOSTFILE
            ./setup/ep-start-docker-cluster.sh up size=1
            EXEC_COMMAND="../../software/utils/ms-time.sh docker openmp ep c 16"
        fi
    else # native
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE -np $parallelism ep.C.$parallelism > /dev/null 2&>1'"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism ep.C.$parallelism > /dev/null 2&>1'"
        else # OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh 'OMP_NUM_THREADS=$parallelism ./ep.C.x > /dev/null 2&>1'"
        fi
    fi

    # Exclude first line
    if [ $name != "name" ]; then
        echo "$EXEC_COMMAND";

        # --- EXECUTION
        EXEC_TIME=$($EXEC_COMMAND)
        if [ $environment = "docker" ]; then
            ./setup/ep-start-docker-cluster.sh down size=16
            ./setup/assemble-swarm destroy $HOSTFILE
        fi
    fi
done < $INPUT_FILE
IFS=$OLDIFS

echo "Done."
