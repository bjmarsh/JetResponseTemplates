#! /bin/bash

INDIR=/nfs-7/userdata/bemarsh/JRTbabies/94x_v7_newSamp_newJEC
HADOOP=/hadoop/cms/store/user/bemarsh/JRTbabies/94x_v7_newSamp_newJEC

mkdir -p logs

# for FILE in `ls ${INDIR}/*.root`; do
#     bn=`basename ${FILE}`
#     fn=${bn%.*}
#     echo JRTlooper $fn $FILE
#     nohup nice -n 10 JRTlooper $fn $FILE &> logs/log_${fn}.txt &
# done

for SAMP in `ls ${HADOOP}`; do
    bn=`basename ${SAMP}`
    if [[ $bn == *"ext"* ]]; then
        continue
    fi
    echo JRTlooper $bn ${INDIR}/${bn}*
    nohup nice -n 10 JRTlooper $bn ${INDIR}/${bn}* &> logs/log_${bn}.txt &
done

