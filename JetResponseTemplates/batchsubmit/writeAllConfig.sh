rm ./input.tar.gz
DIR=$PWD
tar -hcf input.tar --exclude='.git' --exclude='*.root' --exclude='PhysicsTools' --exclude='batchsubmit' ../../../../../CMSSW_8_0_11
gzip ${DIR}/input.tar

while read LINE; do
    echo $LINE
    ./writeConfig.sh $LINE
done <samples_aod.txt
