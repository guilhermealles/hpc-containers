#!/bin/bash
INPUT_FILE='./ep-configuration.csv'
RESULTS_FILE='./ep-results.csv'
LOCAL_NAS_BUILD="$HOME/NAS"
RUN_EP_PWD=$(pwd)

echo "Creating singularity containers..."
echo "$(./setup/ep-singularity-setup.sh)"

echo "Setting up docker keys..."
chmod 600 ../../docker/cluster/ssh/*

echo "Running experiments from input CSV file..."
sleep 5
OLDIFS=$IFS
IFS=","
while read -u 11 name order number rp environment context parallelism block exectime
do
    cd "$RUN_EP_PWD"
    HOSTFILE="./hosts.txt"
    HOSTFILE_FORCE_COMM="./hosts-force-comm.txt";

    if [ $context = "openmp" ]; then
        export OMP_NUM_THREADS=$parallelism
    fi

    if [ $environment = "singularity" ]; then
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE -np $parallelism ./images/mpi-c-$parallelism.img'"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism ./images/mpi-c-$parallelism.img'"
        else #OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh './images/omp-c.img'"
        fi
    elif [ $environment = "docker" ]; then
        if [ $context = "mpi" ]; then
            ./setup/assemble-swarm.sh create $HOSTFILE
            ./setup/ep-start-docker-cluster.sh up $parallelism
            EXEC_COMMAND="../../software/utils/ms-time.sh './swarm.sh exec mpirun -np $parallelism NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        elif [ $context = "mpi-high-comm" ]; then
            ./setup/assemble-swarm.sh create $HOSTFILE_FORCE_COMM
            ./setup/ep-start-docker-cluster.sh up $parallelism
	        EXEC_COMMAND="../../software/utils/ms-time.sh './swarm.sh exec mpirun -np $parallelism NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        else # OpenMP
            ./setup/assemble-swarm.sh create $HOSTFILE
            ./setup/ep-start-docker-cluster.sh up 1
            EXEC_COMMAND="../../software/utils/ms-time.sh './swarm.sh exec NPB3.3.1/NPB3.3-OMP/bin/ep.B.x"
        fi
        cd ../../docker/cluster/
    else # native
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE -np $parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh 'mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism'"
        else # OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh '$LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-OMP/bin/ep.B.x'"
        fi
    fi

    # Exclude first line
    if [ $name != "name" ]; then
        echo "$EXEC_COMMAND";

        # --- EXECUTION
        EXEC_TIME=$($EXEC_COMMAND)
        echo $EXEC_TIME
        echo "$name,$environment,$context,$parallelism,$EXEC_TIME" > $RESULTS_FILE
        if [ $environment = "docker" ]; then
            ./setup/ep-start-docker-cluster.sh down 16
            ./setup/assemble-swarm destroy $HOSTFILE
        fi
    fi
    sleep 5
done 11< $INPUT_FILE
IFS=$OLDIFS

echo "Done."
