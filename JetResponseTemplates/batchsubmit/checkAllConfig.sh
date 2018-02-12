#/bin/bash

condor_dir=$1

for condorfile in $condor_dir/*; do
    python checkConfig.py $condorfile
done
