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
    if [ $environment = "singularity" ]; then
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh mpirun ./images/mpi-c-$parallelism.img"
        elif [ $context = "mpi-high-comm" ]; then
            # TODO
            EXEC_COMMAND="../../software/utils/ms-time.sh high_comm mpirun ./images/mpi-c-$parallelism.img"
        else #OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh OMP_NUM_THREADS=$parallelism ./images/omp-c.img"
        fi
    elif [ $environment = "docker" ]; then
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh docker mpi ep c 16"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh docker mpi-high-comm ep c 16"
        else # OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh docker openmp ep c 16"
        fi
    else # native
        if [ $context = "mpi" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh mpirun --hostfile file ep.C.$parallelism"
        elif [ $context = "mpi-high-comm" ]; then
            EXEC_COMMAND="../../software/utils/ms-time.sh mpirun --hostfile file2 ep.C.$parallelism"
        else # OpenMP
            EXEC_COMMAND="../../software/utils/ms-time.sh OMP_NUM_THREADS=$parallelism ./ep.C.x"
        fi
    fi

    # Exclude first line
    if [ $name != "name" ]; then
        echo "$EXEC_COMMAND";
    fi
done < $INPUT_FILE
IFS=$OLDIFS

echo "Done."