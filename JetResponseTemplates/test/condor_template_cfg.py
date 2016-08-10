import FWCore.ParameterSet.Config as cms

process = cms.PSet()

process.fwliteInput = cms.PSet(
    fileNames   = cms.vstring(
        'root://cmsxrootd.fnal.gov/PUTFILENAMEHERE',
        ),
    maxEvents   = cms.int32(-1),                             ## optional
    outputEvery = cms.uint32(100),                            ## optional
)

process.fwliteOutput = cms.PSet(
    fileName  = cms.string('out.root')
)
