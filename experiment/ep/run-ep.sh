#!/bin/bash
LOCAL_NAS_BUILD="$HOME/NAS"
SINGULARITY_IMAGES_DIR="$HOME/singularity_images"
DOCKER_CLUSTER_DIR="$HOME/singularity-mpi/docker/cluster"
EXPERIMENT_HOME_DIR="$HOME/singularity-mpi/experiment/ep"
SOFTWARE_UTILS_DIR="$HOME/singularity-mpi/software/utils"
INPUT_FILE="$EXPERIMENT_HOME_DIR/ep-configuration.csv"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/ep-results.csv"

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
    if [ $name = "name" ]; then
        continue
    fi
    cd "$EXPERIMENT_HOME_DIR"
    HOSTFILE="./hosts.txt"
    HOSTFILE_FORCE_COMM="./hosts-force-comm.txt";

    if [ $context = "openmp" ]; then
        export OMP_NUM_THREADS=$parallelism
    fi

    if [ $environment = "singularity" ]; then
        if [ $context = "mpi" ]; then
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE -np $parallelism $SINGULARITY_IMAGES_DIR/mpi-c-$parallelism.img)
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism $SINGULARITY_IMAGES_DIR/mpi-c-$parallelism.img)
        else #OpenMP
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh $SINGULARITY_IMAGES_DIR/omp-c.img)
        fi
    elif [ $environment = "docker" ]; then
        if [ $context = "mpi" ]; then
            ./setup/assemble-swarm.sh create $HOSTFILE
            ./setup/ep-start-docker-cluster.sh up $parallelism
            cd $DOCKER_CLUSTER_DIR
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec mpirun -np $parallelism NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism)
        elif [ $context = "mpi-high-comm" ]; then
            ./setup/assemble-swarm.sh create $HOSTFILE_FORCE_COMM
            ./setup/ep-start-docker-cluster.sh up $parallelism
            cd $DOCKER_CLUSTER_DIR
	        EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec mpirun -np $parallelism NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism)
        else # OpenMP
            ./setup/assemble-swarm.sh create $HOSTFILE
            ./setup/ep-start-docker-cluster.sh up 1
            cd $DOCKER_CLUSTER_DIR
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh ./swarm.sh exec NPB3.3.1/NPB3.3-OMP/bin/ep.B.x)
        fi
        cd "$EXPERIMENT_HOME_DIR"
        ./setup/ep-start-docker-cluster.sh down 16
        ./setup/assemble-swarm.sh destroy $HOSTFILE
    else # native
        if [ $context = "mpi" ]; then
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE -np $parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism)
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh mpirun --hostfile $HOSTFILE_FORCE_COMM -np $parallelism $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-MPI/bin/ep.B.$parallelism)
        else # OpenMP
            EXEC_TIME=$($SOFTWARE_UTILS_DIR/ms-time.sh $LOCAL_NAS_BUILD/NPB3.3.1/NPB3.3-OMP/bin/ep.B.x)
        fi
    fi

    echo ">>>>>>>>>> $name,$environment,$context,$parallelism,$EXEC_TIME <<<<<<<<<<"
    echo "$name,$environment,$context,$parallelism,$EXEC_TIME" >> $RESULTS_FILE
    sleep 10
done 11< $INPUT_FILE
IFS=$OLDIFS

echo "Done."
