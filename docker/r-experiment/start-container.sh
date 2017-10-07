#!/bin/bash
PWD="$(pwd)";
sudo docker build -t doe-base .
sudo docker run -v "$PWD/files:/files" -it doe-base /bin/bash
