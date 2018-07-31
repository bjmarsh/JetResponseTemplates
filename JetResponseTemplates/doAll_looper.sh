#! /bin/bash

INDIR=/nfs-7/userdata/bemarsh/JRTbabies/94x_v5/
HADOOP=/hadoop/cms/store/user/bemarsh/JRTbabies/94x_v5

mkdir -p logs

# for FILE in `ls ${INDIR}/*.root`; do
#     bn=`basename ${FILE}`
#     fn=${bn%.*}
#     echo JRTlooper $fn $FILE
#     nohup nice -n 10 JRTlooper $fn $FILE &> logs/log_${fn}.txt &
# done

for SAMP in `ls ${HADOOP}`; do
    bn=`basename ${SAMP}`
    echo JRTlooper $bn ${INDIR}/${bn}*
    nohup nice -n 10 JRTlooper $bn ${INDIR}/${bn}* &> logs/log_${bn}.txt &
done

