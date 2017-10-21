#!/bin/bash

if [ -z $1 ] || [ -z $2 ]; then
    echo "Usage: $0 <create/destroy> <hostfile>"
    exit
fi

COMMAND=$1
HOSTFILE=$2

if [ $COMMAND = "create" ]; then
    SWARM_COMMAND=$(docker swarm init | tr '\\\n' ' ' | sed 's/Swarm initialized.*command:\(.*\)To add a manager.*/\1/')
    echo $SWARM_COMMAND
elif [ $COMMAND = "destroy" ]; then
    SWARM_COMMAND="docker swarm leave --force"
fi

while read -u 10 LINE
do
    HOST=$(echo $LINE | sed 's/:*$//')
    echo "Processing $COMMAND on host $HOST"
    echo $(ssh -o StrictHostKeyChecking=no alles@$HOST "$SWARM_COMMAND")
    echo "Host $HOST processed."
done 10< $HOSTFILE
