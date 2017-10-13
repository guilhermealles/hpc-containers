#!/bin/bash

if [ -z $1 ] ; then
    echo "USAGE: $0 <cluster-size>"
    return;
fi

CLUSTER_SIZE=$1

cd ../../docker/cluster
./swarm.sh size="$CLUSTER_SIZE"

# Sleep for 15 seconds in order to allow the containers to spin up
sleep 15