#!/bin/bash
PWD="$(pwd)";
docker build -t doe-base .
docker run -v "$PWD/files:/files" -it doe-base /bin/bash
