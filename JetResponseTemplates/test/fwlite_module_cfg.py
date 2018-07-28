import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        # 'root://cmsxrootd.fnal.gov//store/mc/RunIISpring16DR80/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/60000/20FD8A42-1D01-E611-85CD-5065F381A2F1.root',
        # 'file:/nfs-7/userdata/bemarsh/JRTbabies/145DDAC0-ACE4-E711-8E2E-24BE05CECBD1.root',
        # 'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17DRPremix/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/AODSIM/94X_mc2017_realistic_v10-v2/60000/124C02A8-67E2-E711-9381-0CC47A4C8E2A.root',
        'file:../scratch/pickevents.root',
        # 'file:../sample_files/ttbar_AOD.root',
        ),
    maxEvents   = cms.int32(1),                             ## optional
    outputEvery = cms.uint32(100),                            ## optional
)
    
process.startEvent = cms.int32(13548)
process.recoJetIdx = cms.int32(1)
process.genJetIdx = cms.int32(1)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('test.root')
)
