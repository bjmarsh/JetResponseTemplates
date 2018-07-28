if [ -f input.tar.gz ]; then
    rm ./input.tar.gz
fi
tar -hcf input.tar --exclude='.git' --exclude='*.root' --exclude='PhysicsTools' --exclude='batchsubmit' ../../../../../CMSSW_9_4_7
gzip ./input.tar
