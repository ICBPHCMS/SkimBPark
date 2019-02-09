#ifndef SkimAnalyzer_h
#define SkimAnalyzer_h


// Package:    HLTAnalysis/SkimAnalyzer
// Class:      TSkimanalyzer
// 
/*class SkimAnalyzer SkimAnalyzer.cc
 Description: [one line class summary]
*/
//
// Original Author:
//         george karathanasis, georgios.karathanasis@cern.ch
//         Created:  Thu, 5 Nov 2018 17:40:23 GMT
//
//


// system include files
#include <memory>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"


#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "HLTrigger/Egamma/plugins/HLTGenericFilter.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "L1Trigger/L1TNtuples/interface/MuonID.h"
#include <vector>
#include "TTree.h"
#include <string>
#include <iostream>
#include "TMath.h"
#include "DataFormats/Common/interface/Ref.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "TLorentzVector.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/CombinedKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackPointingKinematicConstraint.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "Epropagation.h"
#include "MVAReader.h"




class SkimAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

public:
  explicit SkimAnalyzer(const edm::ParameterSet&);
  ~SkimAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:  
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  float Dphi(float phi1,float phi2);
  float DR(float eta1,float phi1,float eta2, float phi2);
  void Init();

  std::pair<std::vector<float>, std::vector<std::vector<std::vector<float> > > > HLTAnalyze(const edm::Event& iEvent, 
											    const edm::EventSetup& iSetup,
											    std::vector<std::string> HLTPath,
											    std::vector<std::string> Seed );

  void genAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  
  std::vector<float> SelectTrg_Object(std::vector<std::vector<float> > &tr1, std::vector<std::vector<float> > &tr2, 
				      std::vector<std::vector<float> > &tr3, std::vector<std::vector<float> > &tr4,
				      std::vector<std::vector<float> > &tr5, std::vector<std::vector<float> > &tr6);

  std::vector<float> SimulateTrigger(std::vector<float> & genMu_pt, std::vector<float> & genMu_eta, 
   				     std::vector<float> & genMu_phi, std::vector<float> & genMu_ch);
  
  std::vector<std::vector<float> > genMuAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  

  bool debugCOUT;

  edm::EDGetToken electronsToken_;
  edm::EDGetToken muonsToken_;
  edm::EDGetToken Tracks_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapVetoToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapSoftToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapMediumToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapTightToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > eleIdMapValueToken_;
  edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_; 
  edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsToken_;
  vector<string> HLTFilter_;
  vector<string> HLTPath_;
  edm::EDGetToken GenToken_;
  edm::EDGetTokenT<reco::PFClusterCollection> clusters_;

  TTree * t1;
  edm::Service<TFileService> fs;

  int nmuons;
  std::vector<float> muon_pt;
  std::vector<float> muon_eta;
  std::vector<float> muon_phi;
  std::vector<float> muon_qual;
  std::vector<float> muon_charge;
  std::vector<float> muon_dxy;
  std::vector<float> muon_dz;
  std::vector<float> muon_edxy;
  std::vector<float> muon_edz;
  std::vector<float> muon_d0;
  std::vector<float> muon_ed0;
  std::vector<float> muon_vx;
  std::vector<float> muon_vz;
  std::vector<float> muon_vy;
  std::vector<float> muon_trkpt;
  std::vector<float> muon_trketa;
  std::vector<float> muon_trkphi;
  std::vector<float> muon_iso;
  std::vector<bool> muon_medium;
  std::vector<bool> muon_loose;
  std::vector<bool> muon_tight;
  std::vector<bool> muon_soft;


  int nelectron;
  std::vector<float> el_pt;
  std::vector<float> el_eta;
  std::vector<float> el_phi;
  std::vector<float> el_charge;
  std::vector<float> el_vx;
  std::vector<float> el_vy;
  std::vector<float> el_vz;
  std::vector<float> el_dxy;
  std::vector<float> el_dz;
  std::vector<float> el_edxy;
  std::vector<float> el_edz;
  std::vector<float> el_trkpt;
  std::vector<float> el_trketa;
  std::vector<float> el_trkphi;
  std::vector<float> el_mva_out;
  std::vector<float> el_mva_iso;
  std::vector<float> el_iso;
  std::vector<float> el_mva_map_value;
  std::vector<bool> el_veto;
  std::vector<bool> el_soft;
  std::vector<bool> el_medium;
  std::vector<bool> el_tight; 
  
  std::vector<float> vertex_x;
  std::vector<float> vertex_y;
  std::vector<float> vertex_z;
  float beam_x;
  float beam_y;
  float beam_z;


  int ntracks;
  std::vector<float> track_pt;
  std::vector<float> track_eta;
  std::vector<float> track_phi;
  std::vector<float> track_norm_chi2;
  std::vector<float> track_charge;
  std::vector<float> track_dxy;
  std::vector<float> track_dz;
  std::vector<float> track_edxy;
  std::vector<float> track_edz;
  std::vector<float> track_mva;
  std::vector<bool> track_MuCleaned;
  std::vector<float> track_vx;
  std::vector<float> track_vy;
  std::vector<float> track_vz;

  //lepton pair 
  std::vector<std::vector<float> > Epair_PtEtaPhiM;
  std::vector<std::vector<float> > Epair_XYZ; 
  std::vector<float> Epair_cos;
  std::vector<float> Epair_chi_prob;
  std::vector<float> Epair_ObjIndex;
  std::vector<float> Epair_TrkIndex;
  std::vector<float> Epair_Lxy;
  std::vector<float> Epair_eLxy;
  std::vector<unsigned int> Epair_ObjId;
  //track triplet
  std::vector<std::vector<float> > TTrack_PtEtaPhiM;
  std::vector<std::vector<float> > TTrack_XYZ;
  std::vector<float> TTrack_chi_prob;
  std::vector<float> TTrack_ObjIndex;
  std::vector<float> TTrack_TrkIndex;
  std::vector<float> TTrack_kid;
  std::vector<float> TTrack_piid;
  std::vector<float> TTrack_mll;
  std::vector<float> TTrack_mKst;
  std::vector<float> TTrack_cos;
  std::vector<float> TTrack_Lxy;
  std::vector<float> TTrack_eLxy;
  std::vector<unsigned int> TTrack_ObjId;

  //Kst
  std::vector<std::vector<float> > Kstpair_PtEtaPhiM;
  std::vector<float> Kstpair_cos;
  std::vector<float> Kstpair_chi_prob;
  std::vector<float> Kstpair_Lxy;
  std::vector<float> Kstpair_eLxy;


  int genpart_B_index;
  int genpart_lep1FromB_index;
  int genpart_lep2FromB_index;
  int genpart_KFromB_index;
  int genpart_KstFromB_index;
  int genpart_KFromKst_index;
  int genpart_PiFromKst_index;

  int genpart_B_pdg;
  int genpart_lep1FromB_pdg;
  int genpart_lep2FromB_pdg;
  int genpart_KFromB_pdg;
  int genpart_KstFromB_pdg;
  int genpart_KFromKst_pdg;
  int genpart_PiFromKst_pdg;

  std::vector<float> genpart_B_PtEtaPhiM;
  std::vector<float> genpart_lep1_PtEtaPhiM;
  std::vector<float> genpart_lep2_PtEtaPhiM;
  std::vector<float> genpart_K_PtEtaPhiM;
  std::vector<float> genpart_Kst_PtEtaPhiM;
  std::vector<float> genpart_Pi_PtEtaPhiM;

  std::vector<float> genMu_pt;
  std::vector<float> genMu_eta;
  std::vector<float> genMu_phi;
  std::vector<float> genMu_ch;
  std::vector<float> genMu_motherId;
  std::vector<float> genMu_gmotherId;

  std::vector<std::vector<float> > TrgObj1_PtEtaPhiCharge;
  std::vector<std::vector<float> > TrgObj2_PtEtaPhiCharge;
  std::vector<std::vector<float> > TrgObj3_PtEtaPhiCharge;
  std::vector<std::vector<float> > TrgObj4_PtEtaPhiCharge;
  std::vector<std::vector<float> > TrgObj5_PtEtaPhiCharge;
  std::vector<std::vector<float> > TrgObj6_PtEtaPhiCharge;
  std::vector<float> SelectedTrgObj_PtEtaPhiCharge; 

  int SelectedMu_index;  
  float SelectedMu_DR, MuTrgMatchCone;

  int trigger1, trigger2, trigger3, trigger4, trigger5, trigger6;
  unsigned int event, run_number, ls;

  //options
  bool IsData, SaveHLT;
  int LeptonFinalStateID;
  bool isKll;

  float PtTrack_Cut, EtaTrack_Cut;
  float MinChi2Track_Cut, MaxChi2Track_Cut;
  float TrackSdxy_Cut;
 
  //leplep
  float MaxMee_Cut, MinMee_Cut, Probee_Cut, Cosee_Cut;
  float EpairZvtx_Cut;
  //triplet
  float MaxMB_Cut, MinMB_Cut, PtB_Cut, SLxy_Cut, ProbeeK_Cut, CoseeK_Cut;
  float Ksdxy_Cut;

  float MinKst_Cut, MaxKst_Cut, PtKst_Cut, ProbKst_Cut, CosKst_Cut, SLxyKst_Cut;
  float MaxBeeKst_Cut, MinBeeKst_Cut, PtBeeKst_Cut, ProbBeeKst_Cut, CosBeeKst_Cut, SLxyBeeKst_Cut;

  float PtMu_Cut, QualMu_Cut;
  float PtEl_Cut;
  float PtKTrack_Cut;

  float MuTrgMuDz_Cut, ElTrgMuDz_Cut, TrkTrgMuDz_Cut; // was TrackMuDz_Cut
  float MuTrgExclusionCone, ElTrgExclusionCone, TrkTrgExclusionCone;  //was TrgExclusionCone
  float TrkObjExclusionCone; //objects non overlapping
  float TrkTrkMinDR_Cut;  //same s above
  //  float MuTrkMinDR_Cut; // same as above

  bool ObjPtLargerThanTrack;
  bool SaveOnlyTracks, SaveOnlyEPairTracks, UseOnlyBKeeMCForTriplets, EarlyStop, SkipIfNoMuMatch;

  /*
  MVAReader bdt;
  std::string weights;
  double MaxMVA_Cut, MinMVA_Cut;
  */

  //bdt
  /*
  float trk_pt; float trk_eta; float trk_phi; float trk_p; float trk_charge;
  float trk_nhits; float trk_high_purity; float trk_inp; float trk_outp;
  float trk_chi2red; float preid_trk_ecal_Deta; float preid_trk_ecal_Dphi;
  float preid_e_over_p; 
  */
 };


#endif
