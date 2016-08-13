import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        'root://cmsxrootd.fnal.gov//store/mc/RunIISpring16DR80/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/60000/20FD8A42-1D01-E611-85CD-5065F381A2F1.root',
        # 'root://cmsxrootd.fnal.gov//store/mc/RunIISpring16DR80/TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/50000/BCD12251-8702-E611-A467-001E67D8A423.root',
        # 'file:../sample_files/ttbar_MINIAOD.root',
        # 'file:../sample_files/ttbar_AOD.root',
        ),
    maxEvents   = cms.int32(500),                             ## optional
    outputEvery = cms.uint32(100),                            ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('test.root')
)
