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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        # 'file:../scratch/pickevents.root'
        'root://cmsxrootd.fnal.gov/PUTFILENAMEHERE'
    )
)

process.JRT = cms.EDAnalyzer('JRTbabymaker2',
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

    outFile      = cms.string("out.root")
)

from RecoMET.METFilters.metFilters_cff import EcalDeadCellTriggerPrimitiveFilter, primaryVertexFilter, eeBadScFilter, HBHENoiseFilterResultProducer
process.EcalDeadCellTriggerPrimitiveFilter = EcalDeadCellTriggerPrimitiveFilter.clone()
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
process.eeBadScFilter = eeBadScFilter.clone()
process.eeBadScFilter.taggingMode = cms.bool(True)
process.primaryVertexFilter = primaryVertexFilter.clone()
process.HBHENoiseFilterResultProducer = HBHENoiseFilterResultProducer.clone()

# process.out = cms.OutputModule("PoolOutputModule",
#                                fileName = cms.untracked.string("out.root"),
#                                outputCommands = cms.untracked.vstring("keep *"),
# )

process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter * 
                                                           process.eeBadScFilter * 
                                                           # process.primaryVertexFilter * 
                                                           process.HBHENoiseFilterResultProducer * 
                                                           process.JRT)
#process.p = cms.EndPath(process.out)

# from FWCore.ParameterSet.Utilities import convertToUnscheduled
# process=convertToUnscheduled(process)
