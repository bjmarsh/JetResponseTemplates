# rm ./input.tar.gz
# tar -hcf input.tar --exclude='.git' --exclude='*.root' --exclude='PhysicsTools' --exclude='batchsubmit' ../../../../../CMSSW_9_4_1
# gzip ./input.tar

TAG=94x_v4

while read LINE; do
    echo $LINE $TAG
    ./writeConfig.sh $LINE $TAG
done < samples_ptbin_aod.txt
