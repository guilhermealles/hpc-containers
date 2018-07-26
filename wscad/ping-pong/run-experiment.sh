#!/bin/bash
EXPERIMENT_HOME_DIR="$PWD"
INPUT_FILE="$EXPERIMENT_HOME_DIR/config/ping_pong_doe.csv"

echo "Running experiments from input CSV file..."
sleep 2
OLDIFS=$IFS
IFS=","
while read -u 11 name order number rp environment sizeexp block exectime
do
    if [ $name = "name" ]; then
        continue
    fi
    cd "$EXPERIMENT_HOME_DIR"

    if [ $environment = "singularity" ]; then
        echo $(./singularity-single-run.sh $name $environment $sizeexp)
    elif [ $environment = "docker" ]; then
        echo $(./docker-single-run.sh $name $environment $sizeexp)
    else # native
        echo $(./native-single-run.sh $name $environment $sizeexp)
    fi
    sleep 1
done 11< $INPUT_FILE
IFS=$OLDIFS

echo "Done."
