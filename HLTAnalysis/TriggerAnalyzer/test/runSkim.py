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

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,globaltag,'')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(

#data
#'/store/data/Run2018B/ParkingBPH5/AOD/PromptReco-v1/000/317/650/00000/E216E20E-E06E-E811-8F13-FA163E15F06E.root'
#mc
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FF7DE4ED-3ED2-1040-BD7F-87729F17AE92.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/F9FCE450-3579-FB4A-B5B6-3B792929EA6C.root', 
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FC68B2C9-79D1-7348-859C-3F0E4C362A4C.root', 
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/FF7DE4ED-3ED2-1040-BD7F-87729F17AE92.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/E8E7B542-E3FC-2740-81A3-E4272BBBB006.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/E8EE9916-18D4-6547-883C-055A8345C330.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/EAF93F92-4D7A-FD4C-B726-4EDA797DD991.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/ECD4A3A2-57FC-BE4C-A354-2F83B59210E6.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/CBF34D3D-411B-7041-809C-A074911D8977.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/CD12E151-0572-9F47-8A1D-A66CF5A1DA97.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/CF2FA866-2269-B343-9D1B-6D83EAFF4906.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/E10A7C57-D2F7-3249-A263-BF8B080AD3B3.root',
'/store/mc/RunIIAutumn18DR/BuToK_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/90000/C484C9AE-FCDA-1345-BA8A-0CE3144316E4.root'
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
                                 IsKll = cms.bool(True),

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
                                 SkipIfNoMuMatch = cms.bool(False)
                              )
                           )

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
   #SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( "outputResHPt_MC10k.root" )
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




'''
/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_1.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_2.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_3.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_4.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_5.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_6.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_7.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_8.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_9.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_10.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_11.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_12.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_13.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_14.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_15.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_16.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_17.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_18.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_19.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_20.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_21.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_22.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_23.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_24.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_25.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_26.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_27.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_28.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_29.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_30.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_31.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_32.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_33.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_34.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_35.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_36.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_37.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_38.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_39.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_40.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_41.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_42.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_43.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_44.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_45.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_46.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_47.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_48.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_49.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_50.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_51.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_52.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_53.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_54.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_55.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_56.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_57.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_58.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_59.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_60.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_61.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_62.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_63.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_64.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_65.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_66.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_67.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_68.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_69.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_70.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_71.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_72.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_73.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_74.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_75.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_76.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_77.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_78.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_79.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_80.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_81.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_82.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_83.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_84.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_85.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_86.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_87.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_88.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_89.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_90.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_91.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_92.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_93.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_94.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_95.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_96.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_97.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_98.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_99.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_100.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_101.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_102.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_103.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_104.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_105.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_106.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_107.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_108.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_109.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_110.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_111.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_112.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_113.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_114.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_115.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_116.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_117.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_118.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_119.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_120.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_121.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_122.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_123.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_124.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_125.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_126.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_127.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_128.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_129.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_130.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_131.root',

'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_133.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_134.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_135.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_136.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_137.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_138.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_139.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_140.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_141.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_142.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_143.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_144.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_145.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_146.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_147.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_148.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_149.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_150.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_151.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_152.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_153.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_154.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_155.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_156.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_157.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_158.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_159.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_160.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_161.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_162.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_163.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_164.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_165.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_166.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_167.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_168.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_169.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_170.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_171.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_172.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_173.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_174.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_175.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_176.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_177.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_178.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_179.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_180.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_181.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_182.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_183.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_184.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_185.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_186.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_187.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_188.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_189.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_190.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_191.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_192.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_193.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_194.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_195.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_196.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_197.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_198.root',
'/store/user/tstreble/BToKee_Pythia/BToKee_Pythia_AODSIM_18_03_22/180321_162718/0000/BToKee_AODSIM_199.root'
'''
