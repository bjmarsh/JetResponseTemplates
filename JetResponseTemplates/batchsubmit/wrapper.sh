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

export SCRAM_ARCH=slc6_amd64_gcc630
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd CMSSW_9_4_1/src/JetResponseTemplates/JetResponseTemplates
echo "[wrapper] in directory: " ${PWD}
echo "[wrapper] attempting to build"
eval `scramv1 runtime -sh`
scramv1 b ProjectRename
scram b
eval `scramv1 runtime -sh`

echo "PATH: " $PATH
echo "LD_LIBRARY_PATH: " $LD_LIBRARY_PATH

ls ../../../bin/slc6_amd64_gcc530

cd test
sed -i -- "s#PUTFILENAMEHERE#${FILE}#g" condor_template_cfg.py

echo "[wrapper] running: JRTbabymaker condor_template_cfg.py"

cmsRun -n4 condor_template_cfg.py

if [ ! -f out.root ]; then
    # if it doesn't exist, try again
    echo "[wrapper] Output file not produced, trying once more..."
    cmsRun -n4 condor_template_cfg.py
else
    SIZE=`stat --printf="%s\n" out.root`
    if [ "$SIZE" -lt 10000 ]; then
        #file produced, but it is empty (usually xrootd error)
        echo "[wrapper] Output does not seem to be valid, trying once more..."
        rm out.root
        cmsRun -n4 condor_template_cfg.py
    fi
fi

if [ -f out.root ]; then
    SIZE=`stat --printf="%s\n" out.root`
    if [ "$SIZE" -lt 10000 ]; then
        echo "[wrapper] still invalid output. quitting."
        rm out.root
    fi
    RES=`python -c "import ROOT as r;f=r.TFile(\"out.root\");print 1 if f.IsZombie() else 0"`
    if [ "$RES" -eq 1 ]; then
        echo "[wrapper] file is a zombie. quitting."
        rm out.root
    fi
else
    echo "[wrapper] still no output produced. quitting."
fi


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

export LD_PRELOAD=/usr/lib64/gfal2-plugins/libgfal_plugin_xrootd.so
gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUT} gsiftp://gftp.t2.ucsd.edu${COPYDIR}/${FILEID}.root

echo "[wrapper] cleaning up"
for FILE in `find . -not -name "*stderr" -not -name "*stdout"`; do rm -rf $FILE; done
echo "[wrapper] cleaned up"
ls
