import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        # 'root://cmsxrootd.fnal.gov//store/mc/RunIISpring16DR80/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/60000/20FD8A42-1D01-E611-85CD-5065F381A2F1.root',
        'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17DRPremix/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/AODSIM/94X_mc2017_realistic_v10-v1/40000/2E4A341E-35D5-E711-97BC-5065F37D4131.root',
        # 'file:../sample_files/ttbar_MINIAOD.root',
        # 'file:../sample_files/ttbar_AOD.root',
        ),
    maxEvents   = cms.int32(10),                             ## optional
    outputEvery = cms.uint32(100),                            ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('test.root')
)
