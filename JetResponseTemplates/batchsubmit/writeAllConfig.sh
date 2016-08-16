rm ./input.tar.gz
tar -hcf input.tar --exclude='.git' --exclude='*.root' --exclude='PhysicsTools' --exclude='batchsubmit' ../../../../../CMSSW_8_0_11
gzip ./input.tar

TAG=v3

while read LINE; do
    echo $LINE $TAG
    ./writeConfig.sh $LINE $TAG
done <samples_aod.txt
