#!/bin/bash

# TODO arrumar as combinações de parâmetros (down com size não faz sentido)
return
if [ -z $1 ] || [ -z $2 ] ; then
    echo "USAGE: $0 <up/down> <cluster-size>"
    return;
fi

START_COMMAND=$1
CLUSTER_SIZE=$2

cd ../../docker/cluster
./swarm.sh "$START_COMMAND" size="$CLUSTER_SIZE"

# Sleep for 20 seconds in order to allow the containers to spin up
sleep 20
