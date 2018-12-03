#!/bin/bash

function run () {
    echo root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/$1/\",\"${OUTPUTDIR}/$1.root\"\)
    nohup nice -n 10 root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/$1/\",\"${OUTPUTDIR}/$1.root\"\) >& ${LOGDIR}/log_merge_$1.txt &
}

TAG=94x_v7_newSamp_newJEC

HADOOPDIR=/hadoop/cms/store/user/${USER}/JRTbabies/${TAG}
OUTPUTDIR=/nfs-7/userdata/bemarsh/JRTbabies/${TAG}
LOGDIR=mergeLogs

mkdir -p $OUTPUTDIR
mkdir -p $LOGDIR
chmod -R a+wrx $OUTPUTDIR

run qcd_pt15to30
run qcd_pt15to30
run qcd_pt30to50
run qcd_pt30to50_ext1
run qcd_pt50to80
run qcd_pt50to80_ext1
run qcd_pt80to120
run qcd_pt80to120_ext1
run qcd_pt120to170
run qcd_pt170to300
run qcd_pt170to300_ext1
run qcd_pt300to470
run qcd_pt300to470_ext1
run qcd_pt470to600
run qcd_pt600to800
run qcd_pt600to800_ext1
run qcd_pt800to1000
run qcd_pt800to1000_ext1
run qcd_pt1000to1400
run qcd_pt1000to1400_ext1
run qcd_pt1400to1800
run qcd_pt1400to1800_ext1
run qcd_pt1800to2400
run qcd_pt2400to3200
run qcd_pt3200toInf

# run qcd_ht200to300
# run qcd_ht300to500
# run qcd_ht500to700
# run qcd_ht700to1000
# run qcd_ht1000to1500
# run qcd_ht1500to2000
# run qcd_ht2000toInf
