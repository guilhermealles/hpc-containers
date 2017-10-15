#!/bin/bash

if [ -z $1 ] || [ -z $2 ]; then
    echo "Usage: $0 <create/destroy> <hostfile>"
    return
fi

COMMAND=$1
HOSTFILE=$2

if [ $COMMAND = "create" ]; then
    SWARM_COMMAND=$(docker swarm init | grep 'docker swarm join' | sed 's/To\sadd\sa\smanager.*$//')
elif [ $COMMAND = "destroy" ]; then
    SWARM_COMMAND="docker swarm leave"
fi

OLDIFS=$IFS
IFS='\n'
while read LINE
do
    HOST=$(echo $LINE | sed 's/ .*$//')
    ssh alles@$HOST $SWARM_COMMAND
done < $HOSTFILE
IFS=$OLDIFS

if [ $COMMAND = "destroy" ]; then
    docker swarm leave --force
fi