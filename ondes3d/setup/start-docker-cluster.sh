#!/bin/bash

if [ -z $1 ] || [ -z $2 ] ; then
    echo "USAGE: $0 <up/down> <cluster-size>"
    exit
fi

START_COMMAND=$1
PARALLELISM=$2
CLUSTER_SIZE=$(($PARALLELISM+1))

cd docker-cluster
./swarm.sh "$START_COMMAND" size="$CLUSTER_SIZE"

# Sleep for 20 seconds in order to allow the containers to spin up
sleep 20