DIR=$1
for FILE in `ls ${DIR}/condor*resubmit.cmd`; do
    condor_submit $FILE
done