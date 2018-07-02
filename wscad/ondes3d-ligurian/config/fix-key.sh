#!/bin/bash

hosts=$(cat /home/alles/hpc-containers/wscad/ondes3d-ligurian/config/hosts.txt)
for host in $hosts; do 
	echo $host
	ssh "alles@$host" ls
done < './hosts.txt'
