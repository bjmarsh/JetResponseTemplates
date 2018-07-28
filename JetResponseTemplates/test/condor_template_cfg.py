import FWCore.ParameterSet.Config as cms

process = cms.Process("JRT")

process.load("FWCore.MessageService.MessageLogger_cfi")
# load event level configurations
process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
# process.load('Configuration.Geometry.GeometryRecoDB_cff')
# process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v10', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:../scratch/pickevents_split.root'
        'root://cmsxrootd.fnal.gov/PUTFILENAMEHERE',
    )
)

process.JRT = cms.EDAnalyzer('JRTbabymaker',
    genjets      = cms.InputTag("ak4GenJets"),
    pfjets       = cms.InputTag("ak4PFJetsCHS"),
    genparticles = cms.InputTag("genParticles"),
    pfcands      = cms.InputTag("particleFlow"),
    vertices     = cms.InputTag("offlinePrimaryVertices"),
    pfmet        = cms.InputTag("pfMet"),
    genmet       = cms.InputTag("genMetTrue"),
    muons        = cms.InputTag("muons"),
    fixedGridRho = cms.InputTag("fixedGridRhoFastjetAll"),
    ecalTP       = cms.InputTag("EcalDeadCellTriggerPrimitiveFilter"),
    hbheNoise    = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"),
    hbheNoiseIso = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),
    eeBadSC      = cms.InputTag("eeBadScFilter"),
    badPFMuon    = cms.InputTag("BadPFMuonFilter"),
    ecalBadCalib = cms.InputTag("ecalBadCalibFilter"),
    badChargedCandidate = cms.InputTag("BadChargedCandidateFilter"),
    globalTightHalo2016 = cms.InputTag("globalTightHalo2016Filter"),

    outFile      = cms.string("out.root")
)

from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter, primaryVertexFilter, eeBadScFilter, HBHENoiseFilterResultProducer, BadPFMuonFilter, BadChargedCandidateFilter, globalTightHalo2016Filter, ecalBadCalibFilter
process.EcalDeadCellTriggerPrimitiveFilter = EcalDeadCellTriggerPrimitiveFilter.clone()
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
process.eeBadScFilter = eeBadScFilter.clone()
process.eeBadScFilter.taggingMode = cms.bool(True)
process.primaryVertexFilter = primaryVertexFilter.clone()
process.HBHENoiseFilterResultProducer = HBHENoiseFilterResultProducer.clone()
process.BadPFMuonFilter = BadPFMuonFilter.clone()
process.BadPFMuonFilter.taggingMode = cms.bool(True)
process.BadChargedCandidateFilter = BadChargedCandidateFilter.clone()
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)
process.globalTightHalo2016Filter = globalTightHalo2016Filter.clone()
process.globalTightHalo2016Filter.taggingMode = cms.bool(True)
process.ecalBadCalibFilter = ecalBadCalibFilter.clone()
process.ecalBadCalibFilter.taggingMode = cms.bool(True)

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("out.root"),
                               outputCommands = cms.untracked.vstring("keep *"),
)

process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter * 
                                                           process.eeBadScFilter * 
                                                           # process.primaryVertexFilter * 
                                                           process.HBHENoiseFilterResultProducer * 
                                                           process.BadPFMuonFilter * 
                                                           process.BadChargedCandidateFilter * 
                                                           process.globalTightHalo2016Filter * 
                                                           process.ecalBadCalibFilter * 
                                                           process.JRT)
# process.p = cms.EndPath(process.out)

# from FWCore.ParameterSet.Utilities import convertToUnscheduled
# process=convertToUnscheduled(process)
