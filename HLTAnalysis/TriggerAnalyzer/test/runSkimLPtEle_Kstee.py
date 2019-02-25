import FWCore.ParameterSet.Config as cms


IsData = False
globaltag = '102X_upgrade2018_realistic_v15'
HLTsave = False ; 
if IsData:
   print "We have established we Run on data"
   globaltag = '102X_dataRun2_Prompt_v11'
   HLTsave = True ; 
else:
   print "We have established we Run on MC"
print "Run parameters ",globaltag," HLT ",HLTsave

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag,'')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#"file:../testFiles/step3_inMINIAODSIM.root",
#"file:../testFiles/step3.root",
#"file:/vols/cms/amartell/BParking/testFiles/test2_RAW2DIGI_L1Reco_RECO.root",
"file:/vols/cms/vc1116/low_pt_electrons/step3_mc_Kstee.root"
                            ),
                            secondaryFileNames=cms.untracked.vstring(
                            )
                         )


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.AOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = [
   'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff', 
   'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff', 
   #   'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff', 
]
for idmod in my_id_modules:
   setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
   

#For Run A
#print "Obj ","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q"

#print "path ","HLT_Mu9_IP6_part","HLT_Mu8p5_IP3p5","HLT_Mu10p5_IP3p5","HLT_Mu8_IP3","empty","empty"

#foir Run B
#print "Obj ","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q","empty"

#print "path ","HLT_Mu9_IP6_part","HLT_Mu9_IP5","HLT_Mu7_IP4","HLT_Mu8_IP3","HLT_Mu12_IP6","empty"

process.demo = cms.EDAnalyzer('SkimAnalyzer',
                              electrons    = cms.InputTag("gedGsfElectrons"),
                              muons = cms.InputTag("muons::RECO"),
                              tracks = cms.InputTag("generalTracks::RECO"),
                              vertices     = cms.InputTag("offlinePrimaryVertices"),
                              beamSpot = cms.InputTag('offlineBeamSpot'),
                              conversions  = cms.InputTag('allConversions'),
                              gen = cms.InputTag("genParticles::HLT"),
                              clusters = cms.InputTag("particleFlowClusterECAL"),
                              lowptgsfTracks = cms.InputTag("lowPtGsfEleGsfTracks::RECO"),
                              mvaSeeds = cms.VInputTag( cms.InputTag("lowPtGsfElectronSeedValueMaps:unbiased"),
                                                        cms.InputTag("lowPtGsfElectronSeedValueMaps:ptbiased") ),

                              eleIdMapVeto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                              eleIdMapSoft = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wpLoose"),
                              eleIdMapMedium = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
                              eleIdMapTight = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
                              eleIdMapValue = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),

                              triggerresults = cms.InputTag("TriggerResults::HLT"),
                              triggerobjects = cms.InputTag('hltTriggerSummaryAOD','','HLT'),                   
                              
                              #muonIdMap = cms.InputTag("muons","muidTM2DCompatibilityLoose","RECO"),

                              HLTFilter=cms.vstring("hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q","hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q","empty"),                          
                              HLTPath=cms.vstring("HLT_Mu9_IP6_part","HLT_Mu9_IP5","HLT_Mu7_IP4","HLT_Mu8_IP3","HLT_Mu12_IP6","empty"),
                              #HLTPath=cms.vstring("HLT_Mu9_IP6_part","HLT_Mu9_IP5","HLT_Mu7_IP4","HLT_Mu8_IP3","HLT_Mu8p5_IP3p5","_HLT_Mu10p5_IP3p5"),
               
                              RunParameters = cms.PSet(                                 
                                 Data = cms.bool(IsData), SaveHLT = cms.bool(HLTsave),
                                 LeptonFinalStateID = cms.int32(11), # or 13
                                 IsKll = cms.bool(False),
                                 runOnPfEle = cms.bool(False),

                                 ##not really used
                                 #weight = cms.string("HLTAnalysis/TriggerAnalyzer/data/trk_only.xml"),
                                 #MinMVA_Cut = cms.double(-10.4), MaxMVA_Cut = cms.double(-10.9),
                                 
                                 MuTrgMatchCone = cms.double(0.02),
                                 
                                 #track selection
                                 PtTrack_Cut = cms.double(0), 
                                 EtaTrack_Cut = cms.double(2.5), 
                                 MinChi2Track_Cut = cms.double(0), 
                                 MaxChi2Track_Cut = cms.double(1000),
                                 TrackSdxy_Cut = cms.double(0),

                                 #lep lep pair 
                                 ObjPtLargerThanTrack = cms.bool(True), 
                                 MaxMee_Cut = cms.double(100), 
                                 MinMee_Cut = cms.double(0), 
                                 Probee_Cut = cms.double(1.e-34),
                                 EpairZvtx_Cut = cms.double(1000), #llvtx_z wrt trigger vtx_z
                                 Cosee_Cut = cms.double(-10.99),

                                 ##triplet
                                 MaxMB_Cut = cms.double(100), 
                                 MinMB_Cut = cms.double(0),
                                 PtB_Cut = cms.double(-5.0),
                                 SLxy_Cut = cms.double(-4.0), 
                                 ProbeeK_Cut = cms.double(-0.001), 
                                 CoseeK_Cut = cms.double(-90.9),
                                 Ksdxy_Cut = cms.double(0),   #3rd track wrt PV

                                 MinKst_Cut  = cms.double(0),
                                 MaxKst_Cut  = cms.double(100.),
                                 PtKst_Cut = cms.double(0.),
                                 ProbKst_Cut = cms.double(1.e-34),
                                 CosKst_Cut = cms.double(-90.9),
                                 SLxyKst_Cut = cms.double(-4.0),

                                 MaxBeeKst_Cut  = cms.double(100),
                                 MinBeeKst_Cut  = cms.double(0),
                                 PtBeeKst_Cut = cms.double(-5.0),
                                 ProbBeeKst_Cut = cms.double(1.e-34),
                                 CosBeeKst_Cut = cms.double(-90.9),
                                 SLxyBeeKst_Cut = cms.double(-4.0),

                                 PtMu_Cut = cms.double(0.0), 
                                 QualMu_Cut = cms.double(0), 
                                 PtEl_Cut = cms.double(0), 
                                 PtKTrack_Cut = cms.double(0.0),
                                 
                                 MuTrgMuDz_Cut = cms.double(0.3), 
                                 ElTrgMuDz_Cut = cms.double(0.3),
                                 TrackMuDz_Cut = cms.double(0.7),    #wrt reco muon matched to trigger muon

                                 MuTrgExclusionCone = cms.double(0.4), 
                                 ElTrgExclusionCone = cms.double(0.4),
                                 TrkTrgExclusionCone = cms.double(0.4), #wrt trigger muon

                                 #for triplet objects and tracks
                                 TrkObjExclusionCone = cms.double(0.02),
                                 #MuTrkMinDR_Cut = cms.double(0),                                  
                                 #TrkTrkMinDR_Cut = cms.double(0.02), 

                                 #configuration
                                 SaveOnlyTracks = cms.bool(False),
                                 SaveOnlyEPairTracks = cms.bool(False),                                 
                                 UseOnlyBKeeMCForTriplets = cms.bool(False), 
                                 EarlyStop = cms.bool(False),
                                 SkipIfNoMuMatch = cms.bool(True)
                              )
                           )

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   #SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( "outputBSkim_KsteeVince_test_lowPtele.root" )
)
process.fevt = cms.OutputModule("PoolOutputModule",
                                # SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
                                outputCommands = cms.untracked.vstring(#"drop *",
                                   # "keep *_*_*_Demo",
                                   #  "keep *_offlineBeamSpot_*_*",
                                   #  "keep *_offlinePrimaryVertices_*_*",
                                   #  "keep *_offlinePrimaryVerticesWithBS_*_*",
                                ),
                                fileName = cms.untracked.string("edm_output.root")
                             )

#process.p = cms.Path(process.egmGsfElectronIDSequence)#* process.demo)
#process.endjob=cms.EndPath(process.fevt)
process.p = cms.Path(
   process.egmGsfElectronIDSequence   
   +process.demo
)   
   
#process.p = cms.Path(process.demo)


