#!/bin/bash

if [ -z $1 ]; then
    echo "You must provide a ssh directory."
    exit
fi

echo "Fixing permissions on directory $1..."
chmod 600 "$1"