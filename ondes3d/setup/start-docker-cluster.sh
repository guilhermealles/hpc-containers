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

if [ "$START_COMMAND" = 'up' ]; then
    ./swarm.sh exec wait_worker_connections "$PARALLELISM"
elif [ "$START_COMMAND" = 'down' ]; then
    sleep 10
fi