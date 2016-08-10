#!/bin/bash

#
# args
#

FILEID=$1
FILE=$2
COPYDIR=$3

echo "[wrapper] FILEID    = " ${FILEID}
echo "[wrapper] FILE      = " ${FILE}
echo "[wrapper] COPYDIR   = " ${COPYDIR}

echo "[wrapper] printing env"
printenv
echo 

echo "[wrapper] hostname  = " `hostname`
echo "[wrapper] date      = " `date`
echo "[wrapper] linux timestamp = " `date +%s`

#
# untar input sandbox
#

echo "[wrapper] extracting input sandbox"
tar -zxf input.tar.gz

cd CMSSW_8_0_11/src/JetResponseTemplates/JetResponseTemplates
echo "[wrapper] in directory: " ${PWD}
echo "[wrapper] attempting to build"
scramv1 runtime -sh
scramv1 b ProjectRename
scram b

cd test
sed -i -- "s#PUTFILENAMEHERE#${FILE}#g" JRT_cfg.py

echo "[wrapper] running: JRTbabymaker JRT_cfg.py"

JRTbabymaker JRT_cfg.py


#
# do something with output
#

echo "[wrapper] output is"
ls

#
# clean up
#

echo "[wrapper] copying file"
OUTPUT=out.root
echo "[wrapper] OUTPUT = " ${OUTPUT}

if [ ! -d "${COPYDIR}" ]; then
    echo "creating output directory " ${COPYDIR}
    mkdir ${COPYDIR}
fi

gfal-copy -p -f -t 4200 --verbose file:`pwd`/${OUTPUT} srm://bsrm-3.t2.ucsd.edu:8443/srm/v2/server?SFN=${COPYDIR}/${FILEID}.root

echo "[wrapper] cleaning up"
for FILE in `find . -not -name "*stderr" -not -name "*stdout"`; do rm -rf $FILE; done
echo "[wrapper] cleaned up"
ls
