#!/bin/bash
ts=$(date +%s%N) ; $@ > /dev/null 2>&1 ; tt=$((($(date +%s%N) - $ts)/1000000)) ; echo "$tt"
