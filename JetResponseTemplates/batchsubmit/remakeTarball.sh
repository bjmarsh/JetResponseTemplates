if [ -f input.tar.gz ]; then
    rm ./input.tar.gz
fi
tar -hcf input.tar --exclude='.git' --exclude='*.root' --exclude='PhysicsTools' --exclude='batchsubmit' ../../../../../$CMSSW_VERSION
gzip ./input.tar
