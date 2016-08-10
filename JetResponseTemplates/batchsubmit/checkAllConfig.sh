#/bin/bash

condor_dir=$1

for condorfile in $condor_dir/*; do
    # echo "RUNNING " $condorfile
    if [[ $condorfile != *"resubmit.cmd"* ]]; then
        ./makeResubmitConfig.sh $condorfile
    fi
done
