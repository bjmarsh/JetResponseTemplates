# rm ./input.tar.gz
# tar -hcf input.tar --exclude='.git' --exclude='*.root' --exclude='PhysicsTools' --exclude='batchsubmit' ../../../../../CMSSW_9_4_7
# gzip ./input.tar

TAG=102x_v1

while read LINE; do
    echo $LINE $TAG
    ./writeConfig.sh $LINE $TAG
# done < samples_ptbin_aod.txt
done < samples_ptbin_aod_102X.txt
