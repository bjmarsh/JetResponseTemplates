#! /bin/bash

INDIR=/nfs-7/userdata/bemarsh/JRTbabies/94x_v4/

for FILE in `ls ${INDIR}/*.root`; do
    bn=`basename ${FILE}`
    fn=${bn%.*}
    echo JRTlooper $fn $FILE
    nohup nice -n 10 JRTlooper $fn $FILE &> logs/log_${fn}.txt &
done

