#!/bin/bash

function run () {
    echo root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/$1/\",\"${OUTPUTDIR}/$1.root\"\)
    nohup nice -n 19 root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/$1/\",\"${OUTPUTDIR}/$1.root\"\) >& ${LOGDIR}/log_merge_$1.txt &
}

TAG=v3

HADOOPDIR=/hadoop/cms/store/user/${USER}/JRTbabies/${TAG}
OUTPUTDIR=/nfs-7/userdata/bemarsh/JRTbabies/${TAG}
LOGDIR=mergeLogs

mkdir -p $OUTPUTDIR
mkdir -p $LOGDIR
chmod -R a+wrx $OUTPUTDIR

# run qcd_pt15to30
# run qcd_pt30to50
# run qcd_pt50to80
# run qcd_pt80to120
# run qcd_pt120to170
run qcd_pt170to300
# run qcd_pt300to470
# run qcd_pt470to600
# run qcd_pt600to800
# run qcd_pt800to1000
# run qcd_pt1000to1400
# run qcd_pt1400to1800
# run qcd_pt1800to2400
# run qcd_pt2400to3200
# run qcd_pt3200toInf
