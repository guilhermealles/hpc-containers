#!/bin/bash
LOCAL_NAS_BUILD="$HOME/NAS"
DOCKER_CLUSTER_DIR="$PWD/docker-cluster/"
EXPERIMENT_HOME_DIR="$PWD"
SOFTWARE_UTILS_DIR="$EXPERIMENT_HOME_DIR/../software/utils"
INPUT_FILE="$EXPERIMENT_HOME_DIR/config/ep-configuration.csv"
RESULTS_FILE="$EXPERIMENT_HOME_DIR/results/ep-results.csv"
HOSTFILE="$EXPERIMENT_HOME_DIR/config/hosts.txt"
HOSTFILE_FORCE_COMM="$EXPERIMENT_HOME_DIR/config/hosts-force-comm.txt";

echo "Setting up docker keys..."
echo $(./setup/fix-docker-ssh-permissions.sh $DOCKER_CLUSTER_DIR/ssh/id_rsa)

echo "Running experiments from input CSV file..."
sleep 2
OLDIFS=$IFS
IFS=","
while read -u 11 name order number rp environment context parallelism block exectime
do
    if [ $name = "name" ]; then
        continue
    fi
    cd "$EXPERIMENT_HOME_DIR"

    if [ $environment = "singularity" ]; then
        echo $(./singularity-single-run.sh $name $environment $context $parallelism)
    elif [ $environment = "docker" ]; then
        echo $(./docker-single-run.sh $name $environment $context $parallelism)
    else # native
        echo $(./native-single-run.sh $name $environment $context $parallelism)
    fi
    sleep 2
done 11< $INPUT_FILE
IFS=$OLDIFS

echo "Done."
