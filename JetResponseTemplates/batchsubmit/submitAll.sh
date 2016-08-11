DIR=$1

for FILE in `ls ${DIR}/condor*.cmd`; do
    condor_submit $FILE
done