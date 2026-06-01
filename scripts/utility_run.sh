#!/bin/bash

#For situations where running manual scripts is needed (reruns or manual edis to jsons)
DATE=$1
#Appropriate modes: (BCLConvert, Cellbender, Cumulus) 
MODE=$2

for f in scripts/$DATE/alto_terra_$MODE*.sh ; do bash "$f"; done
