#include <memory>
#include <iostream>
#include "HLTAnalysis/TriggerAnalyzer/plugins/SkimAnalyzer.h"


float BChargedMass_ = 5.279;
float B0Mass_ = 5.279;
float KaonMass_ = 0.493677;
float Kaon0StarMass_ = 0.89581; //nominal K*(892) mass => need to change
//float Kaon0StarMass_ = 0.89176; //nominal K*(892) mass => need to change
float PionMass_ = 0.139570;
float MuonMass_ = 0.10565837;
float ElectronMass_ = 0.5109989e-3;


SkimAnalyzer::SkimAnalyzer(const edm::ParameterSet& iConfig):  
  electronsToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>  ("electrons"))),
  muonsToken_(consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
  Tracks_(consumes<std::vector<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  conversionsToken_(consumes< reco::ConversionCollection > (iConfig.getParameter<edm::InputTag> ("conversions"))),
  eleIdMapVetoToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapVeto"))),
  eleIdMapSoftToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapSoft"))),
  eleIdMapMediumToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapMedium"))),
  eleIdMapTightToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleIdMapTight"))),
  eleIdMapValueToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("eleIdMapValue"))),
  trgresultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag> ("triggerresults"))),
  trigobjectsToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag> ("triggerobjects"))),
  HLTFilter_(iConfig.getParameter<vector<string> >("HLTFilter")),
  HLTPath_(iConfig.getParameter<vector<string> >("HLTPath")),
  GenToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("gen"))),
  clusters_ (consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("clusters") )),
  lowPtGsfTracks_(consumes<std::vector<reco::GsfTrack> >(iConfig.getParameter<edm::InputTag>("lowptgsfTracks"))),
  mvaSeedTags_(iConfig.getParameter< std::vector<edm::InputTag> >("mvaSeeds"))
  //runParameters(iConfig.getParameter<edm::ParameterSet>("RunParameters"))
{
  debugCOUT = false;

  for ( const auto& tag : mvaSeedTags_ ){ 
    mvaSeeds_.push_back( consumes<edm::ValueMap<float> >(tag) ); 
  }


  edm::ParameterSet runParameters = iConfig.getParameter<edm::ParameterSet>("RunParameters");     
  IsData = runParameters.getParameter<bool>("Data");
  SaveHLT = runParameters.getParameter<bool>("SaveHLT");
  LeptonFinalStateID = runParameters.getParameter<int>("LeptonFinalStateID");
  isKll = runParameters.getParameter<bool>("IsKll");
  lookAtpfEle = runParameters.getParameter<bool>("runOnPfEle");

  MuTrgMatchCone = runParameters.getParameter<double>("MuTrgMatchCone");

  PtTrack_Cut = runParameters.getParameter<double>("PtTrack_Cut");
  EtaTrack_Cut = runParameters.getParameter<double>("EtaTrack_Cut");
  MaxChi2Track_Cut = runParameters.getParameter<double>("MaxChi2Track_Cut");
  MinChi2Track_Cut = runParameters.getParameter<double>("MinChi2Track_Cut");
  TrackSdxy_Cut = runParameters.getParameter<double>("TrackSdxy_Cut");

  ObjPtLargerThanTrack = runParameters.getParameter<bool>("ObjPtLargerThanTrack"); 
  MaxMee_Cut = runParameters.getParameter<double>("MaxMee_Cut");
  MinMee_Cut = runParameters.getParameter<double>("MinMee_Cut");
  Probee_Cut = runParameters.getParameter<double>("Probee_Cut");
  Cosee_Cut = runParameters.getParameter<double>("Cosee_Cut");
  EpairZvtx_Cut = runParameters.getParameter<double>("EpairZvtx_Cut");
  
  MaxMB_Cut = runParameters.getParameter<double>("MaxMB_Cut");
  MinMB_Cut = runParameters.getParameter<double>("MinMB_Cut");
  PtB_Cut = runParameters.getParameter<double>("PtB_Cut");
  SLxy_Cut = runParameters.getParameter<double>("SLxy_Cut");
  ProbeeK_Cut = runParameters.getParameter<double>("ProbeeK_Cut");
  CoseeK_Cut = runParameters.getParameter<double>("CoseeK_Cut");
  Ksdxy_Cut = runParameters.getParameter<double>("Ksdxy_Cut");

  MinKst_Cut = runParameters.getParameter<double>("MinKst_Cut");
  MaxKst_Cut = runParameters.getParameter<double>("MaxKst_Cut");
  PtKst_Cut = runParameters.getParameter<double>("PtKst_Cut");
  ProbKst_Cut = runParameters.getParameter<double>("ProbKst_Cut");
  CosKst_Cut = runParameters.getParameter<double>("CosKst_Cut");
  SLxyKst_Cut = runParameters.getParameter<double>("SLxyKst_Cut");

  MinBeeKst_Cut = runParameters.getParameter<double>("MinBeeKst_Cut");
  MaxBeeKst_Cut = runParameters.getParameter<double>("MaxBeeKst_Cut");
  PtBeeKst_Cut = runParameters.getParameter<double>("PtBeeKst_Cut");
  ProbBeeKst_Cut = runParameters.getParameter<double>("ProbBeeKst_Cut");
  CosBeeKst_Cut = runParameters.getParameter<double>("CosBeeKst_Cut");
  SLxyBeeKst_Cut = runParameters.getParameter<double>("SLxyBeeKst_Cut");


  PtMu_Cut = runParameters.getParameter<double>("PtMu_Cut");
  QualMu_Cut = runParameters.getParameter<double>("QualMu_Cut");
  PtEl_Cut = runParameters.getParameter<double>("PtEl_Cut");
  PtKTrack_Cut = runParameters.getParameter<double>("PtKTrack_Cut");


  MuTrgMuDz_Cut = runParameters.getParameter<double>("MuTrgMuDz_Cut");
  ElTrgMuDz_Cut = runParameters.getParameter<double>("ElTrgMuDz_Cut"); 
  TrkTrgMuDz_Cut = runParameters.getParameter<double>("TrackMuDz_Cut");

  MuTrgExclusionCone = runParameters.getParameter<double>("MuTrgExclusionCone");
  ElTrgExclusionCone = runParameters.getParameter<double>("ElTrgExclusionCone");
  TrkTrgExclusionCone = runParameters.getParameter<double>("TrkTrgExclusionCone");

  TrkObjExclusionCone = runParameters.getParameter<double>("TrkObjExclusionCone"); 
  //   MuTrkMinDR_Cut = runParameters.getParameter<double>("MuTrkMinDR_Cut");
  //   TrkTrkMinDR_Cut = runParameters.getParameter<double>("TrkTrkMinDR_Cut");



  SaveOnlyTracks = runParameters.getParameter<bool>("SaveOnlyTracks");
  SaveOnlyEPairTracks = runParameters.getParameter<bool>("SaveOnlyEPairTracks");
  UseOnlyBKeeMCForTriplets = runParameters.getParameter<bool>("UseOnlyBKeeMCForTriplets");
  EarlyStop = runParameters.getParameter<bool>("EarlyStop");
  SkipIfNoMuMatch = runParameters.getParameter<bool>("SkipIfNoMuMatch");

  /*
  weights=runParameters.getParameter<string>("weight");
  MaxMVA_Cut=runParameters.getParameter<double>("MaxMVA_Cut");
  MinMVA_Cut=runParameters.getParameter<double>("MinMVA_Cut");

  bdt.SetWeight(weights);  
  */
  
  t1 = fs->make<TTree>("mytree","mytree"); 
  t1->Branch("event",&event);
  t1->Branch("run_number",&run_number);
  t1->Branch("ls",&ls);
  
  t1->Branch("vertex_x",&vertex_x); 
  t1->Branch("vertex_y",&vertex_y);
  t1->Branch("vertex_z",&vertex_z); 
  t1->Branch("beam_x",&beam_x);
  t1->Branch("beam_y",&beam_y); 
  t1->Branch("beam_z",&beam_z);

  t1->Branch("HLT_path1",&trigger1); t1->Branch("HLT_path2",&trigger2);
  t1->Branch("HLT_path3",&trigger3); t1->Branch("HLT_path4",&trigger4);
  t1->Branch("HLT_path5",&trigger5); t1->Branch("HLT_path6",&trigger6);
  t1->Branch("TrgObj1_PtEtaPhiCharge",&TrgObj1_PtEtaPhiCharge);
  t1->Branch("TrgObj2_PtEtaPhiCharge",&TrgObj2_PtEtaPhiCharge);
  t1->Branch("TrgObj3_PtEtaPhiCharge",&TrgObj3_PtEtaPhiCharge);
  t1->Branch("TrgObj4_PtEtaPhiCharge",&TrgObj4_PtEtaPhiCharge);
  t1->Branch("TrgObj5_PtEtaPhiCharge",&TrgObj5_PtEtaPhiCharge);
  t1->Branch("TrgObj6_PtEtaPhiCharge",&TrgObj6_PtEtaPhiCharge);
  t1->Branch("SelectedTrgObj_PtEtaPhiCharge",&SelectedTrgObj_PtEtaPhiCharge);

  t1->Branch("SelectedMu_index",&SelectedMu_index);
  t1->Branch("SelectedMu_DR",&SelectedMu_DR);
   
  t1->Branch("nmuon",&nmuons); t1->Branch("muon_pt",&muon_pt);
  t1->Branch("muon_eta",&muon_eta); t1->Branch("muon_phi",&muon_phi);
  t1->Branch("muon_charge",&muon_charge); t1->Branch("muon_dxy",&muon_dxy);
  t1->Branch("muon_dz",&muon_dz); t1->Branch("muon_edxy",&muon_edxy);
  t1->Branch("muon_edz",&muon_edz); t1->Branch("muon_d0",&muon_d0);
  t1->Branch("muon_ed0",&muon_ed0); t1->Branch("muon_vx",&muon_vx);
  t1->Branch("muon_vy",&muon_vy); t1->Branch("muon_vz",&muon_vz);
  t1->Branch("muon_iso",&muon_iso); t1->Branch("muon_soft",&muon_soft);
  t1->Branch("muon_loose",&muon_loose); t1->Branch("muon_medium",&muon_medium);
  t1->Branch("muon_tight",&muon_tight);t1->Branch("muon_trkpt",&muon_trkpt); 
  t1->Branch("muon_trketa",&muon_trketa);t1->Branch("muon_trkphi",&muon_trkphi);
  
  t1->Branch("nelectron",&nelectron); t1->Branch("el_pt",&el_pt);
  t1->Branch("el_eta",&el_eta); t1->Branch("el_phi",&el_phi);
  t1->Branch("el_charge",&el_charge); t1->Branch("el_dxy",&el_dxy);
  t1->Branch("el_dz",&el_dz); t1->Branch("el_edxy",&el_edxy);
  t1->Branch("el_edz",&el_edz); t1->Branch("el_vx",&el_vx);
  t1->Branch("el_vy",&el_vy); t1->Branch("el_vz",&el_vz);
  t1->Branch("el_mva_out",&el_mva_out); t1->Branch("el_mva_iso",&el_mva_iso);
  t1->Branch("el_iso",&el_iso); t1->Branch("el_veto",&el_veto);
  t1->Branch("el_soft",&el_soft); t1->Branch("el_medium",&el_medium);
  t1->Branch("el_tight",&el_tight); t1->Branch("el_mva_map_value",&el_mva_map_value);
  t1->Branch("el_trkpt",&el_trkpt); t1->Branch("el_trketa",&el_trketa); 
  t1->Branch("el_trkphi",&el_trkphi);


  t1->Branch("ngsfTracks", &ngsfTracks);
  t1->Branch("gsfTrk_pt", &gsfTrk_pt);
  t1->Branch("gsfTrk_eta", &gsfTrk_eta);
  t1->Branch("gsfTrk_phi", &gsfTrk_phi);
  t1->Branch("gsfTrk_charge", &gsfTrk_charge);
  t1->Branch("gsfTrk_seedBDTunb", &gsfTrk_seedBDTunb);
  t1->Branch("gsfTrk_seedBDTbiased", &gsfTrk_seedBDTbiased);

  t1->Branch("genpart_B_index", &genpart_B_index);
  t1->Branch("genpart_lep1FromB_index", &genpart_lep1FromB_index);
  t1->Branch("genpart_lep2FromB_index", &genpart_lep2FromB_index);
  t1->Branch("genpart_KFromB_index", &genpart_KFromB_index);
  t1->Branch("genpart_KstFromB_index", &genpart_KstFromB_index);
  t1->Branch("genpart_KFromKst_index", &genpart_KFromKst_index);
  t1->Branch("genpart_PiFromKst_index", &genpart_PiFromKst_index);

  t1->Branch("genpart_B_pdg", &genpart_B_pdg);
  t1->Branch("genpart_lep1FromB_pdg", &genpart_lep1FromB_pdg);
  t1->Branch("genpart_lep2FromB_pdg", &genpart_lep2FromB_pdg);
  t1->Branch("genpart_KFromB_pdg", &genpart_KFromB_pdg);
  t1->Branch("genpart_KstFromB_pdg", &genpart_KstFromB_pdg);
  t1->Branch("genpart_KFromKst_pdg", &genpart_KFromKst_pdg);
  t1->Branch("genpart_PiFromKst_pdg", &genpart_PiFromKst_pdg);

  t1->Branch("genpart_B_PtEtaPhiM", &genpart_B_PtEtaPhiM);
  t1->Branch("genpart_lep1_PtEtaPhiM", &genpart_lep1_PtEtaPhiM);
  t1->Branch("genpart_lep2_PtEtaPhiM", &genpart_lep2_PtEtaPhiM);
  t1->Branch("genpart_K_PtEtaPhiM", &genpart_K_PtEtaPhiM);
  t1->Branch("genpart_Kst_PtEtaPhiM", &genpart_Kst_PtEtaPhiM);
  t1->Branch("genpart_Pi_PtEtaPhiM", &genpart_Pi_PtEtaPhiM);
  
  t1->Branch("genMu_pt",&genMu_pt);  t1->Branch("genMu_eta",&genMu_eta);
  t1->Branch("genMu_phi",&genMu_phi);  t1->Branch("genMu_ch",&genMu_ch);
  t1->Branch("genMu_motherId",&genMu_motherId);  t1->Branch("genMu_gmotherId",&genMu_gmotherId);
   
  t1->Branch("ntracks",&ntracks); t1->Branch("track_pt",&track_pt);
  t1->Branch("track_eta",&track_eta); t1->Branch("track_phi",&track_phi);
  t1->Branch("track_norm_chi2",&track_norm_chi2);
  t1->Branch("track_charge",&track_charge);
  t1->Branch("track_dxy",&track_dxy); t1->Branch("track_dz",&track_dz);
  t1->Branch("track_edxy",&track_edxy); t1->Branch("track_edz",&track_edz);
  t1->Branch("track_vx",&track_vx); t1->Branch("track_vy",&track_vy);
  t1->Branch("track_vz",&track_vz); t1->Branch("track_mva",&track_mva);

  t1->Branch("Epair_PtEtaPhiM",&Epair_PtEtaPhiM);
  t1->Branch("Epair_XYZ",&Epair_XYZ); t1->Branch("Epair_cos",&Epair_cos);
  t1->Branch("Epair_chi_prob",&Epair_chi_prob);
  t1->Branch("Epair_ObjIndex",&Epair_ObjIndex); t1->Branch("Epair_TrkIndex",&Epair_TrkIndex);
  t1->Branch("Epair_Lxy",&Epair_Lxy); t1->Branch("Epair_eLxy",&Epair_eLxy);
  t1->Branch("Epair_ObjId",&Epair_ObjId);

  t1->Branch("TTrack_PtEtaPhiM",&TTrack_PtEtaPhiM);
  t1->Branch("TTrack_chi_prob",&TTrack_chi_prob);
  // t1->Branch("TTrack_min2trk_prob",&TTrack_min2trk_prob);
  t1->Branch("TTrack_XYZ",&TTrack_XYZ); t1->Branch("TTrack_ObjIndex",&TTrack_ObjIndex);
  t1->Branch("TTrack_TrkIndex",&TTrack_TrkIndex); 
  t1->Branch("TTrack_KstarIndex",&TTrack_KstarIndex); 
  t1->Branch("TTrack_kid",&TTrack_kid);
  t1->Branch("TTrack_piid",&TTrack_piid);
  t1->Branch("TTrack_mll",&TTrack_mll); 
  t1->Branch("TTrack_mKst",&TTrack_mKst); 
  t1->Branch("TTrack_cos",&TTrack_cos);
  t1->Branch("TTrack_Lxy",&TTrack_Lxy); t1->Branch("TTrack_eLxy",&TTrack_eLxy);
  t1->Branch("TTrack_ObjId",&TTrack_ObjId);

  t1->Branch("Kstpair_PtEtaPhiM", &Kstpair_PtEtaPhiM);
  t1->Branch("Kstpair_cos", &Kstpair_cos);
  t1->Branch("Kstpair_chi_prob", &Kstpair_chi_prob);
  t1->Branch("Kstpair_Lxy", &Kstpair_Lxy);
  t1->Branch("Kstpair_eLxy", &Kstpair_eLxy);

}


SkimAnalyzer::~SkimAnalyzer()
{
  // cout<<"total "<<trg_counter<<" fires "<<fire_counter<<" l3 "<<l3_counter<<endl;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  //  cout<<"total trg2(Mu17) "<<totb1<<" trg3(Mu20) "<<totb2<<" trg5(Mu27) "<<totb3<<endl;
  
  
}


void SkimAnalyzer::genAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  edm::Handle<edm::View<reco::GenParticle> > genPart;
  iEvent.getByToken(GenToken_,genPart);

  //LeptonFinalStateID
  int nGenPart = genPart->size();
  
  if(debugCOUT) std::cout << " \n new event " << std::endl;


  for(int i_Bu=0; i_Bu<nGenPart; ++i_Bu){
    genpart_B_index = -1;
    genpart_lep1FromB_index = -1;
    genpart_lep2FromB_index = -1;
    genpart_KFromB_index = -1;
    genpart_KstFromB_index = -1;
    genpart_KFromKst_index = -1;
    genpart_PiFromKst_index = -1;
    

    if(isKll){
      if(abs((*genPart)[i_Bu].pdgId()) == 521){
	int  nD = (*genPart)[i_Bu].numberOfDaughters();
	
	if(debugCOUT) std::cout << " found a B " << " isKll = " << isKll << " nDaug = " << nD << std::endl;
	
	if(nD < 3) continue;
	genpart_B_index = i_Bu;

	for(auto iD=0; iD < nD; ++iD){
	  const Candidate* daug = (*genPart)[i_Bu].daughter(iD);
	  
	  int pdgId = daug->pdgId();
	  float partPt = daug->pt();
	  float partEta = daug->eta();
      
	  if(debugCOUT) std::cout << " daug " << iD << " pdgID = " << pdgId << " pt = " << partPt << " eta = " << partEta << std::endl;

	  if(partPt < 0.2) continue;
	  if(std::abs(partEta) > 2.6) continue;
	  
	  if(abs(pdgId) == LeptonFinalStateID && genpart_lep1FromB_index < 0)
	    genpart_lep1FromB_index = iD;
	  else if(abs(pdgId) == LeptonFinalStateID)
	    genpart_lep2FromB_index = iD;
	  else if(abs(pdgId) == 321)
	    genpart_KFromB_index = iD;
	  else if (abs(pdgId) == 22) continue;
	  else {
	    genpart_lep1FromB_index = -1;
	    genpart_lep2FromB_index = -1;
	    genpart_KFromB_index = -1;
	    break;
	  }
	}//daug
	if(genpart_lep1FromB_index >= 0 && genpart_lep2FromB_index >= 0 && genpart_KFromB_index >= 0) {
	  if(debugCOUT) std::cout << " found B -> K decay " << std::endl;
	  break;
	}
	else{
	  genpart_B_index = -1;
	  genpart_lep1FromB_index = -1;
	  genpart_lep2FromB_index = -1;
	  genpart_KFromB_index = -1;
	}
      }//B+
    }//isKll
    if(!isKll){
      if(abs((*genPart)[i_Bu].pdgId()) == 511){
	int  nD = (*genPart)[i_Bu].numberOfDaughters();
	
	if(debugCOUT) std::cout << " found a B " << " isK*ll = " << isKll << " nDaug = " << nD << std::endl;

	if(nD < 3) continue;
	genpart_B_index = i_Bu;
	
	for(auto iD=0; iD < nD; ++iD){
	  const Candidate* daug = (*genPart)[i_Bu].daughter(iD);
	  
	  int pdgId = daug->pdgId();
	  float partPt = daug->pt();
	  float partEta = daug->eta();

	  if(debugCOUT) std::cout << " daug " << iD << " pdgID = " << pdgId << " pt = " << partPt << " eta = " << partEta << std::endl;	  

	  if(partPt < 0.2) continue;
	  if(std::abs(partEta) > 2.6) continue;
	  
	  if(abs(pdgId) == LeptonFinalStateID && genpart_lep1FromB_index < 0)
	    genpart_lep1FromB_index = iD;
	  else if(abs(pdgId) == LeptonFinalStateID)
	    genpart_lep2FromB_index = iD;
	  else if(abs(pdgId) == 313){
	    int  ngD = daug->numberOfDaughters();
	    
	    if(ngD < 2) continue;
	    genpart_KstFromB_index = iD;
	    
	    for(auto igD=0; igD < ngD; ++igD){
	      const Candidate* gDaug = daug->daughter(igD);
	      
	      int pdgId_gd = gDaug->pdgId();
	      float partPt_gd = gDaug->pt();
	      float partEta_gd = gDaug->eta();
	      
	      if(debugCOUT) std::cout << " gdaug " << iD << " pdgID = " << pdgId_gd << " pt = " << partPt_gd << " eta = " << partEta_gd << std::endl;

	      if(partPt_gd < 0.2) continue;
	      if(std::abs(partEta_gd) > 2.6) continue;
	      
	      if(abs(pdgId_gd) == 321) 
	      genpart_KFromKst_index = igD;
	      else if(abs(pdgId_gd) == 211)
		genpart_PiFromKst_index = igD;
	      else if (abs(pdgId_gd) == 22) continue;
	      else{
		genpart_KstFromB_index = -1;
		genpart_KFromKst_index = -1;
		genpart_PiFromKst_index = -1;
		break;
	      }
	    } //loop over  K*
	  }
	  else if (abs(pdgId) == 22) continue;
	  else {
	    genpart_lep1FromB_index = -1;
	    genpart_lep2FromB_index = -1;
	    genpart_KstFromB_index = -1;
	    genpart_KFromKst_index = -1;
	    genpart_PiFromKst_index = -1;
	    break;
	  }
	}//daug
	if(genpart_KFromKst_index != -1 && genpart_PiFromKst_index != -1 && genpart_lep1FromB_index != -1 && genpart_lep2FromB_index != -1) {
	  if(debugCOUT) std::cout << " found B -> K* decay " << std::endl;
	  break;
	}
	else{
	  genpart_B_index = -1;
	  genpart_lep1FromB_index = -1;
	  genpart_lep2FromB_index = -1;
	  genpart_KstFromB_index = -1;
	  genpart_KFromKst_index = -1;
	  genpart_PiFromKst_index = -1;
	}
      }// found B*
    }// K*ll
  }// loop over gen
  
  
  if(isKll && genpart_KFromB_index < 0) return;
  if(!isKll && (genpart_KstFromB_index < 0 || genpart_PiFromKst_index < 0 || genpart_KFromKst_index < 0)) return;
  if(genpart_lep1FromB_index < 0 || genpart_lep2FromB_index < 0) return;
  
  if(debugCOUT) std::cout << " now filling vectors  " << std::endl;
  
  if((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->pt() > (*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->pt()){
    genpart_lep1_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->pt()));
    genpart_lep1_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->eta()));
    genpart_lep1_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->phi()));
    genpart_lep1_PtEtaPhiM.push_back(float((LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_));
    
    genpart_lep2_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->pt()));
    genpart_lep2_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->eta()));
    genpart_lep2_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->phi()));
    genpart_lep2_PtEtaPhiM.push_back(float((LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_));

    genpart_B_pdg = (*genPart)[genpart_B_index].pdgId();
    genpart_lep1FromB_pdg = (*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->pdgId();
    genpart_lep2FromB_pdg = (*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->pdgId();
  }
  else{
    genpart_lep2_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->pt()));
    genpart_lep2_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->eta()));
    genpart_lep2_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->phi()));
    genpart_lep2_PtEtaPhiM.push_back(float((LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_));
    
    genpart_lep1_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->pt()));
    genpart_lep1_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->eta()));
    genpart_lep1_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->phi()));
    genpart_lep1_PtEtaPhiM.push_back(float((LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_));

    genpart_B_pdg = (*genPart)[genpart_B_index].pdgId();
    genpart_lep1FromB_pdg = (*genPart)[genpart_B_index].daughter(genpart_lep2FromB_index)->pdgId();
    genpart_lep2FromB_pdg = (*genPart)[genpart_B_index].daughter(genpart_lep1FromB_index)->pdgId();    
  }
  
  if(isKll){
    genpart_K_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KFromB_index)->pt()));
    genpart_K_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KFromB_index)->eta()));
    genpart_K_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KFromB_index)->phi()));
    genpart_K_PtEtaPhiM.push_back(KaonMass_);

    genpart_KFromB_pdg = (*genPart)[genpart_B_index].daughter(genpart_KFromB_index)->pdgId();
  }
  else{
    genpart_Kst_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->pt()));
    genpart_Kst_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->eta()));
    genpart_Kst_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->phi()));
    genpart_Kst_PtEtaPhiM.push_back(Kaon0StarMass_);

    genpart_K_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_KFromKst_index)->pt()));
    genpart_K_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_KFromKst_index)->eta()));
    genpart_K_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_KFromKst_index)->phi()));
    genpart_K_PtEtaPhiM.push_back(KaonMass_);

    genpart_Pi_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_PiFromKst_index)->pt()));
    genpart_Pi_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_PiFromKst_index)->eta()));
    genpart_Pi_PtEtaPhiM.push_back(float((*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_PiFromKst_index)->phi()));
    genpart_Pi_PtEtaPhiM.push_back(PionMass_);

    genpart_KstFromB_pdg = (*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->pdgId();
    genpart_KFromKst_pdg = (*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_KFromKst_index)->pdgId();
    genpart_PiFromKst_pdg = (*genPart)[genpart_B_index].daughter(genpart_KstFromB_index)->daughter(genpart_PiFromKst_index)->pdgId();
  }
  
 
  return;
}



std::vector<std::vector<float> > SkimAnalyzer::genMuAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ 
  using namespace std;
  using namespace edm;
  using namespace reco;

  edm::Handle<edm::View<reco::GenParticle>> genPart;
  iEvent.getByToken(GenToken_,genPart);

  std::vector<float> genMu_pt,genMu_eta,genMu_phi,genMu_ch,genMu_motherId,genMu_gmotherId;
  std::vector<std::vector<float>> genMu;

  for(typename edm::View<reco::GenParticle>::const_iterator gen=genPart->begin(); gen !=genPart->end(); gen++){ 
    if (fabs(gen->pdgId())!=13) continue;
    const Candidate * mo= gen->mother();
    if (fabs(mo->pdgId()) == 13) continue; 
    const Candidate * mo2= mo->mother(); //const Candidate * mo3= mo2->mother();
    genMu_pt.push_back(gen->pt()); genMu_eta.push_back(gen->eta()); 
    genMu_phi.push_back(gen->phi()); genMu_ch.push_back(gen->charge()); 
    genMu_motherId.push_back(mo->pdgId()); genMu_gmotherId.push_back(mo2->pdgId());

    // std::cout << " in genMuAnalyze found 1 muon pt = " << gen->pt() << std::endl;
  }
  genMu.push_back(genMu_pt); genMu.push_back(genMu_eta); genMu.push_back(genMu_phi);
  genMu.push_back(genMu_ch); genMu.push_back(genMu_motherId); genMu.push_back(genMu_gmotherId);
 

  return genMu;
}



std::pair<std::vector<float>,std::vector<std::vector<std::vector<float> > > >  SkimAnalyzer::HLTAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup,
													std::vector<std::string> HLTPath, std::vector<std::string> Seed ){

  using namespace std;  using namespace edm;  using namespace reco;
  using namespace trigger;
 
  edm::Handle<trigger::TriggerEvent> triggerObjectsSummary;
  iEvent.getByToken(trigobjectsToken_ ,triggerObjectsSummary);
  
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);
  
  trigger::TriggerObjectCollection selectedObjects;
  std::vector<float> fires;
  std::vector<std::vector<std::vector<float>>> trg_event; 

  for (unsigned int ipath=0; ipath<Seed.size(); ipath++){ 
    std::vector<std::vector<float>> tot_tr_obj_pt_eta_phi;
    if (triggerObjectsSummary.isValid()) {  
      size_t filterIndex = (*triggerObjectsSummary).filterIndex(InputTag(Seed[ipath],"","HLT"));
      trigger::TriggerObjectCollection allTriggerObjects = triggerObjectsSummary->getObjects();     

      if (filterIndex < (*triggerObjectsSummary).sizeFilters()) { 
	const trigger::Keys &keys = (*triggerObjectsSummary).filterKeys(filterIndex);
	for (size_t j = 0; j < keys.size(); j++) {
	  trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
	  std::vector<float> tr_obj_pt_eta_phi;
	  if (fabs(foundObject.id()) != 13) continue;
	  tr_obj_pt_eta_phi.push_back(foundObject.pt());
	  tr_obj_pt_eta_phi.push_back(foundObject.eta());
	  tr_obj_pt_eta_phi.push_back(foundObject.phi());
	  tr_obj_pt_eta_phi.push_back(foundObject.id()/fabs(foundObject.id()));
	  tot_tr_obj_pt_eta_phi.push_back( tr_obj_pt_eta_phi);
	}
      }
 
    }    
    trg_event.push_back(tot_tr_obj_pt_eta_phi);
  }

  //paths
  float fire0 = 0, fire1 = 0, fire2 = 0, fire3 = 0, fire4 = 0, fire5 = 0;
    if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);    
    //cout << "new" << endl;
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      TString TrigPath = trigName.triggerName(i_Trig);    
      if (TrigPath.Contains(HLTPath[0]) && trigResults->accept(i_Trig)){
	fire0 = 1;}
      if(TrigPath.Contains(HLTPath[1])  && trigResults->accept(i_Trig) ){
          fire1 = 1;}
      if (TrigPath.Contains(HLTPath[2]) && trigResults->accept(i_Trig)){ 
          fire2 = 1;}
      if (TrigPath.Contains(HLTPath[3]) && trigResults->accept(i_Trig)){
	  fire3 = 1;}
      if (TrigPath.Contains(HLTPath[4]) && trigResults->accept(i_Trig)){
          fire4 = 1;}
      if (TrigPath.Contains(HLTPath[5]) && trigResults->accept(i_Trig)){ 
	  fire5 = 1;}
    } 
    }
    
    fires.push_back(fire0);  fires.push_back(fire1);   fires.push_back(fire2);
    fires.push_back(fire3); fires.push_back(fire4); fires.push_back(fire5);
    return std::make_pair(fires,trg_event);
}



std::vector<float> SkimAnalyzer::SelectTrg_Object(std::vector<std::vector<float> > &tr1, std::vector<std::vector<float> > &tr2, 
						  std::vector<std::vector<float> > &tr3, std::vector<std::vector<float> > &tr4,
						  std::vector<std::vector<float> > &tr5, std::vector<std::vector<float> > &tr6){
  
  //each vector contains pt, eta, phi, charge
  std::vector<std::vector<float> > max1;
  for (auto & vec: tr1) max1.push_back(vec);
  for (auto & vec: tr2) max1.push_back(vec);
  for (auto & vec: tr3) max1.push_back(vec);
  for (auto & vec: tr4) max1.push_back(vec);
  for (auto & vec: tr5) max1.push_back(vec);
  for (auto & vec: tr6) max1.push_back(vec);
 
  //vector of pt sorted trigger objects
  std::sort(max1.begin(), max1.end(),
	    [](const std::vector<float>& a, const std::vector<float>& b) {
	      return a[0] > b[0];
	    });
  return max1[0];
}


std::vector<float> SkimAnalyzer::SimulateTrigger(std::vector<float> & genMu_pt, std::vector<float> & genMu_eta, 
						     std::vector<float> & genMu_phi,std::vector<float> & genMu_ch){
  std::vector<std::vector<float>> max1;
  for(unsigned int igen=0; igen<genMu_pt.size(); igen++){

    if (genMu_pt.at(igen) < 9 || fabs(genMu_eta.at(igen)) > 1.5) continue;
    std::vector<float> temp; 
    temp.push_back(genMu_pt.at(igen));
    temp.push_back(genMu_eta.at(igen));  temp.push_back(genMu_phi.at(igen));
    temp.push_back(genMu_ch.at(igen));
    max1.push_back(temp);
  }
  if (max1.size() > 0)
    //sort by decreasing pt   
    std::sort(max1.begin(), max1.end(),
	      [](const std::vector<float>& a, const std::vector<float>& b) {
		return a[0] > b[0];
	      });
  else {
    std::vector<float> temp; temp.push_back(-99); temp.push_back(-99); temp.push_back(-99); temp.push_back(-99);
    max1.push_back(temp);
  }
  return max1[0];
}



float SkimAnalyzer::Dphi(float phi1,float phi2){
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}


float SkimAnalyzer::DR(float eta1,float phi1,float eta2, float phi2){
  return TMath::Sqrt((eta1-eta2)*(eta1-eta2)+Dphi(phi1,phi2)*Dphi(phi1,phi2));
}


int SkimAnalyzer::sharedHits(reco::Track& track1, reco::Track& track2) const{

  int match = 0;

  for (trackingRecHit_iterator hit1 = track1.recHitsBegin(); hit1 != track1.recHitsEnd(); ++hit1){
    if ( !(*hit1)->isValid() ) continue;
    DetId id1 = (*hit1)->geographicalId();
    //if ( id1.det() != DetId::Muon ) continue; //ONLY MUON

    //LogTrace(category_)<<"first ID "<<id1.rawId()<<" "<<(*hit1)->localPosition()<<endl;
    GlobalPoint pos1 = (*hit1)->globalPosition();

    for (trackingRecHit_iterator hit2 = track2.recHitsBegin(); hit2 != track2.recHitsEnd(); ++hit2) {
      if ( !(*hit2)->isValid() ) continue;

      DetId id2 = (*hit2)->geographicalId();
      //if ( id2.det() != DetId::Muon ) continue; //ONLY MUON
      //          LogTrace(category_)<<"second ID "<<id2.rawId()<< (*hit2)->localPosition()<<endl;

      if (id2.rawId() != id1.rawId() ) continue;

      GlobalPoint pos2 = (*hit2)->globalPosition();
      if ( ( pos1 - pos2 ).mag()< 10e-5 ) match++;
    }
  }
  return match;
}



void SkimAnalyzer::Init(){

  vertex_x.clear();  vertex_y.clear();  vertex_z.clear();  
  beam_x = -9999, beam_y = -9999, beam_z = -99999;

  nmuons=0; nelectron=0; ntracks=0;
  muon_pt.clear(); muon_eta.clear(); muon_phi.clear(); muon_charge.clear();
  muon_dz.clear(); muon_dxy.clear(); muon_edz.clear(); muon_edxy.clear();
  muon_d0.clear(); muon_ed0.clear(); muon_vx.clear(); muon_vy.clear();
  muon_vz.clear(); muon_iso.clear(); muon_soft.clear(); muon_loose.clear();
  muon_medium.clear(); muon_tight.clear(); muon_trkpt.clear();
  muon_trketa.clear(); muon_trkphi.clear();

  el_pt.clear(); el_eta.clear(); el_phi.clear(); el_charge.clear();
  el_vx.clear(); el_vy.clear(); el_vz.clear(); el_dxy.clear(); el_mva_out.clear();
  el_dz.clear();el_edxy.clear(); el_edz.clear(); el_mva_out.clear(); 
  el_mva_iso.clear(); el_iso.clear();   el_mva_map_value.clear();
  el_veto.clear(); el_soft.clear(); el_medium.clear(); el_tight.clear();
  el_trkpt.clear(); el_trketa.clear(); el_trkphi.clear();
     
  ngsfTracks = 0;
  gsfTrk_pt.clear();
  gsfTrk_eta.clear();
  gsfTrk_phi.clear();
  gsfTrk_charge.clear();
  gsfTrk_seedBDTunb.clear();
  gsfTrk_seedBDTbiased.clear();

  TrgObj1_PtEtaPhiCharge.clear(); TrgObj2_PtEtaPhiCharge.clear();
  TrgObj3_PtEtaPhiCharge.clear(); TrgObj4_PtEtaPhiCharge.clear();
  TrgObj5_PtEtaPhiCharge.clear(); TrgObj6_PtEtaPhiCharge.clear();
  trigger1=0,trigger2=0,trigger3=0,trigger4=0,trigger5=0,trigger6=0;
  SelectedTrgObj_PtEtaPhiCharge.clear(); 

  SelectedMu_index = -1;
  SelectedMu_DR = 1000;

  genpart_B_index = -1;
  genpart_lep1FromB_index = -1;
  genpart_lep2FromB_index = -1;
  genpart_KFromB_index = -1;
  genpart_KstFromB_index = -1;
  genpart_KFromKst_index = -1;
  genpart_PiFromKst_index = -1;

  genpart_B_pdg = -1;
  genpart_lep1FromB_pdg = -1;
  genpart_lep2FromB_pdg = -1;
  genpart_KFromB_pdg = -1;
  genpart_KstFromB_pdg = -1;
  genpart_KFromKst_pdg = -1;
  genpart_PiFromKst_pdg = -1;
 
  genpart_B_PtEtaPhiM.clear();
  genpart_lep1_PtEtaPhiM.clear();
  genpart_lep2_PtEtaPhiM.clear();
  genpart_K_PtEtaPhiM.clear();
  genpart_Kst_PtEtaPhiM.clear();
  genpart_Pi_PtEtaPhiM.clear();


  track_pt.clear(); track_eta.clear(); track_phi.clear(); track_dxy.clear();
  track_norm_chi2.clear(); track_charge.clear();  track_edz.clear(); 
  track_dz.clear(); track_MuCleaned.clear(); track_edxy.clear();
  track_vx.clear(); track_vz.clear(); track_vy.clear(); track_mva.clear();
  
  TTrack_PtEtaPhiM.clear(); TTrack_XYZ.clear();  TTrack_cos.clear();
  TTrack_chi_prob.clear(); TTrack_ObjIndex.clear(); 
  TTrack_TrkIndex.clear(); TTrack_KstarIndex.clear();
  TTrack_kid.clear(); TTrack_piid.clear();
  TTrack_mll.clear(); 
  TTrack_mKst.clear();
  TTrack_Lxy.clear();
  TTrack_eLxy.clear(); TTrack_ObjId.clear();

  Kstpair_PtEtaPhiM.clear();
  Kstpair_cos.clear();
  Kstpair_chi_prob.clear();
  Kstpair_Lxy.clear();
  Kstpair_eLxy.clear();

  Epair_PtEtaPhiM.clear(); Epair_XYZ.clear();  Epair_cos.clear();
  Epair_chi_prob.clear(); Epair_ObjIndex.clear(); Epair_TrkIndex.clear();
  Epair_Lxy.clear(); Epair_eLxy.clear(); Epair_ObjId.clear();

  genMu_pt.clear(); genMu_eta.clear(); genMu_phi.clear(); genMu_ch.clear(); 
  genMu_motherId.clear(); genMu_gmotherId.clear();

}


void SkimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  //Get a few collections to apply basic electron ID
  //Get electrons
  edm::Handle<edm::View<reco::GsfElectron> > electrons;
  iEvent.getByToken(electronsToken_, electrons);

  edm::Handle<edm::ValueMap<bool> > ele_veto_id;
  iEvent.getByToken(eleIdMapVetoToken_ ,ele_veto_id);

  edm::Handle<edm::ValueMap<bool> > ele_soft_id;
  iEvent.getByToken(eleIdMapSoftToken_ ,ele_soft_id);

  edm::Handle<edm::ValueMap<bool> > ele_medium_id;
  iEvent.getByToken(eleIdMapMediumToken_ ,ele_medium_id);

  edm::Handle<edm::ValueMap<bool> > ele_tight_id;
  iEvent.getByToken(eleIdMapTightToken_ ,ele_tight_id);

  edm::Handle<edm::ValueMap<int> > ele_mva_id_value;
  iEvent.getByToken( eleIdMapValueToken_ ,ele_mva_id_value);
 
  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonsToken_,muons);

  //Get conversions
  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(conversionsToken_, conversions);    
  
  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  
  //Get vertices 
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  
  //continue if there are no vertices
  if (vertices->empty()) { if(debugCOUT) std::cout << "no vertices"<< std::endl; return;}

  edm::Handle<vector<reco::Track>> tracks;
  iEvent.getByToken(Tracks_, tracks);

  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgresultsToken_, trigResults);

  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);

  edm::Handle<reco::PFClusterCollection> clusters;
  iEvent.getByToken(clusters_,clusters);


  //new lop Pt ele
  std::vector<edm::Handle<edm::ValueMap<float> > > mvaSeeds;
  if(!lookAtpfEle){
    for (const auto& token : mvaSeeds_){ 
      edm::Handle<edm::ValueMap<float> > h;
      iEvent.getByToken(token, h); 
      mvaSeeds.push_back(h);
    }
  }

  edm::Handle<std::vector<reco::GsfTrack> > lowPtGsfTracks;
  if(!lookAtpfEle) iEvent.getByToken(lowPtGsfTracks_, lowPtGsfTracks);



  KalmanVertexFitter theKalmanFitter(false);
  TransientVertex LLvertex;
  
  //run stuff
  run_number = iEvent.id().run();
  ls = iEvent.luminosityBlock();

  event++; 


  for (VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx){
    bool isFake = vtx->isFake();
    if ( isFake) continue;
    vertex_x.push_back(vtx->x());  vertex_y.push_back(vtx->y());  vertex_z.push_back(vtx->z()); 
  }
  if (vertex_x.size() ==0 ) return;
  reco::TrackBase::Point  vertex_point;
  vertex_point.SetCoordinates(vertex_x[0], vertex_y[0], vertex_z[0]);

  beam_x = theBeamSpot->x0(); beam_y = theBeamSpot->y0(); beam_z = theBeamSpot->z0();   


  std::vector<std::vector<float> > genparts;
  std::vector<std::vector<float> > genmu;
  genparts.clear(); genmu.clear();

  
  Init();


  if(!lookAtpfEle && debugCOUT){
    //new low Pt ele
    std::cout << " mvaSeeds.size() = " << mvaSeeds.size() << " gsfSize = " << lowPtGsfTracks->size() << std::endl;
    for (unsigned int iter=0; iter< mvaSeeds.size(); ++iter){
      std::cout << "  mvaSeed:           " 
		<< int( mvaSeeds[iter].isValid() ? mvaSeeds[iter]->size() : -1 ) 
		<< ", ";
      if ( mvaSeeds[iter].isValid() && !mvaSeeds[iter]->empty() && lowPtGsfTracks.isValid() ) {
	reco::GsfTrackRef gsf(lowPtGsfTracks, 0);
	std::cout << "\"" << mvaSeedTags_[iter].instance() << "\"";
	if ( gsf.isNonnull() ) { std::cout << "(example value: " << float( (*mvaSeeds[iter])[gsf] ) << ")"; }
      }
      std::cout << std::endl;
    }
  }


  if(debugCOUT) std::cout << " analyzer => MC gen  " << std::endl;

  if(!IsData){
   genAnalyze(iEvent,iSetup);

   
   genmu = genMuAnalyze(iEvent,iSetup);
   genMu_pt=genmu[0]; genMu_eta=genmu[1]; genMu_phi=genmu[2]; 
   genMu_ch=genmu[3];  genMu_motherId=genmu[4]; genMu_gmotherId=genmu[5]; 

   //gen muon pt decreasing order => can be replaced by SelectTrg_Object for new MC
   SelectedTrgObj_PtEtaPhiCharge = SimulateTrigger(genMu_pt, genMu_eta, genMu_phi, genMu_ch); 
   
  }

  if(debugCOUT) std::cout << " analyzer => HLT paths  " << std::endl;

  if (IsData){    
    std::pair<std::vector<float>, std::vector<std::vector<std::vector<float> > > > trgresult = HLTAnalyze(iEvent, iSetup, HLTPath_, HLTFilter_);
    trigger1=trgresult.first[0]; trigger2=trgresult.first[1]; trigger3=trgresult.first[2];  
    trigger4=trgresult.first[3]; trigger5=trgresult.first[4]; trigger6=trgresult.first[5];   

    if(trigger1+trigger2+trigger3+trigger4+trigger5+trigger6 == 0) { if(debugCOUT)std::cout << " no trigger " << std::endl; return;}
    TrgObj1_PtEtaPhiCharge=trgresult.second[0]; TrgObj2_PtEtaPhiCharge=trgresult.second[1];  
    TrgObj3_PtEtaPhiCharge=trgresult.second[2]; TrgObj4_PtEtaPhiCharge=trgresult.second[3]; 
    TrgObj5_PtEtaPhiCharge=trgresult.second[4]; TrgObj6_PtEtaPhiCharge=trgresult.second[5];      

    SelectedTrgObj_PtEtaPhiCharge = SelectTrg_Object(TrgObj1_PtEtaPhiCharge, TrgObj2_PtEtaPhiCharge, TrgObj3_PtEtaPhiCharge, 
						     TrgObj4_PtEtaPhiCharge, TrgObj5_PtEtaPhiCharge, TrgObj6_PtEtaPhiCharge);
  }



  if(debugCOUT) std::cout << " analyzer => Muons  " << std::endl;
 
  std::vector<std::shared_ptr<reco::Track> > MuTracks; 
  std::vector<unsigned int> object_container;
  std::vector<unsigned int> object_id;
  int muIndex = 0;


  float ZvertexTrg = -10000;
  //  trk_index = -1;
  if(debugCOUT) std::cout << " muon->size() = " << muons->size() << std::endl;
  //select at reco level muons with some requirements
  for (std::vector<reco::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); mu++){
    if (fabs(mu->eta()) > EtaTrack_Cut) continue;
    if (mu->pt() < PtMu_Cut) continue;

    bool tight = false, soft = false;
    if(vertices.isValid()){
      tight = isTightMuonCustom(*mu,(*vertices)[0]);
      soft = muon::isSoftMuon(*mu,(*vertices)[0]);
    }
    if (QualMu_Cut == 1 && !soft) continue;
    if (QualMu_Cut == 2 && !isMediumMuonCustom(*mu)) continue; 
    if (QualMu_Cut == 3 && !tight) continue;
    nmuons++; 

    muon_pt.push_back(mu->pt()); muon_phi.push_back(mu->phi());
    muon_eta.push_back(mu->eta()); muon_charge.push_back(mu->charge());
    
    const Track * mutrack= mu->bestTrack(); 
    muon_trkpt.push_back(mutrack->pt()); muon_trketa.push_back(mutrack->eta());
    muon_trkphi.push_back(mutrack->phi());
    muon_dz.push_back(mutrack->dz(vertex_point));
    muon_dxy.push_back(mutrack->dxy(vertex_point));
    muon_vx.push_back(mu->vx()); muon_vy.push_back(mu->vy());
    muon_vz.push_back(mu->vz()); muon_edz.push_back(mutrack->dzError());
    muon_edxy.push_back(mutrack->dxyError()); muon_d0.push_back(mutrack->d0());
    muon_ed0.push_back(mutrack->d0Error());
    muon_medium.push_back(isMediumMuonCustom(*mu));
    muon_loose.push_back(isLooseMuonCustom(*mu));
    
    muon_tight.push_back(tight); muon_soft.push_back(soft);
    auto muTrack = std::make_shared<reco::Track>(*mutrack);
    MuTracks.push_back(muTrack);
    object_container.push_back(muIndex);
    object_id.push_back(13);

    const MuonPFIsolation& isol = mu->pfIsolationR04();
    double mu_iso = (isol.sumChargedHadronPt+max(0.,isol.sumNeutralHadronEt+isol.sumPhotonEt-0.5*isol.sumPUPt))/mu->pt();
    muon_iso.push_back(mu_iso);

    float mu_trig_dR = DR(mu->eta(), mu->phi(), SelectedTrgObj_PtEtaPhiCharge[1], SelectedTrgObj_PtEtaPhiCharge[2]);
    if ( mu_trig_dR < MuTrgMatchCone && mu_trig_dR < SelectedMu_DR){
      SelectedMu_DR = SelectedMu_DR;
      SelectedMu_index = muIndex;  
      ZvertexTrg = mu->vz();
    }
    ++muIndex;
    //delete mutrack;
  } // reco muon
  
  
  if (SelectedMu_index == -1 && SkipIfNoMuMatch){
    //Init(); 
    t1->Fill();
    if(debugCOUT) std::cout << " SelectedMu_index == -1 " << std::endl;
    return;
  }

  std::vector<std::shared_ptr<reco::Track> > ElTracks;
  int eleIndex = 0;
  if(debugCOUT) std::cout << " electrons->size() = " << electrons->size() << std::endl;
  for(size_t e = 0; e<electrons->size(); e++){
      const auto el = electrons->ptrAt(e);  
      bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, theBeamSpot->position());
      if (!passConvVeto) { if(debugCOUT) std::cout << " !passConvVeto " << std::endl; continue;}
      if (fabs(el->eta()) > EtaTrack_Cut) { if(debugCOUT) std::cout << " > EtaTrack_Cut " << std::endl; continue;}
      if (el->pt() < PtEl_Cut) { if(debugCOUT) std::cout << " < PtEl_Cut " << std::endl; continue;}

      if (SelectedMu_index != -1 ){
	if (fabs(ZvertexTrg - el->vz()) > ElTrgMuDz_Cut ) { if(debugCOUT) std::cout << " vtxZ " << std::endl; continue;}
	if ( DR(el->eta(),el->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < ElTrgExclusionCone) {
	  if(debugCOUT) std::cout << " dR trig " << std::endl;
	  continue;
	}
	//if ( DR(el->eta(),el->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < ElTrgExclusionCone) continue;
      }
      ++nelectron;
      el_pt.push_back(el->pt()); el_eta.push_back(el->eta());
      el_phi.push_back(el->phi()); el_charge.push_back(el->charge());

      const Track * eltrack = el->bestTrack();
      el_trkpt.push_back(eltrack->pt()); el_trketa.push_back(eltrack->eta());
      el_trkphi.push_back(eltrack->phi());
      el_dz.push_back(eltrack->dz(vertex_point));
      el_dxy.push_back(eltrack->dxy(vertex_point));
      el_edxy.push_back(eltrack->dxyError());
      el_edz.push_back(eltrack->dzError());
      el_vx.push_back(el->vx()); el_vy.push_back(el->vy());
      el_vz.push_back(el->vz()); el_mva_out.push_back(el->mva_e_pi());
      el_mva_iso.push_back(el->mva_Isolated());

      double iso = el->pfIsolationVariables().sumChargedHadronPt + max(0.0, el->pfIsolationVariables().sumNeutralHadronEt+el->pfIsolationVariables().sumPhotonEt - 
								       0.5 * el->pfIsolationVariables().sumPUPt) /el->pt();
      el_iso.push_back(iso); el_veto.push_back((*ele_veto_id)[el]);
      el_soft.push_back((*ele_soft_id)[el]);
      el_medium.push_back((*ele_medium_id)[el]);
      el_tight.push_back((*ele_tight_id)[el]);
      el_mva_map_value.push_back((*ele_mva_id_value)[el]);

      auto ElTrack = std::make_shared<reco::Track>(*eltrack);
      ElTracks.push_back(ElTrack); 
      object_container.push_back(eleIndex);
      ++eleIndex; 
      object_id.push_back(11);
      //  delete eltrack;
  }


  std::vector<std::shared_ptr<reco::GsfTrack> > gsfElTracks;
  std::vector<unsigned int> gsf_container;
  std::vector<unsigned int> gsf_id;
  if(!lookAtpfEle){
  int gsfeleIndex = 0;
  if(debugCOUT) std::cout << " lowPtGsfTracks->size() = " << lowPtGsfTracks->size() << std::endl;

  for (std::vector<reco::GsfTrack>::const_iterator gsfT=lowPtGsfTracks->begin(); gsfT!=lowPtGsfTracks->end(); ++gsfT){
    if(debugCOUT) std::cout << " nuova gsf " << int(gsfT - lowPtGsfTracks->begin()) << " pt = " << gsfT->pt() << std::endl;
    //if(!gsfT->isValid()) continue;

    if (fabs(gsfT->eta()) > EtaTrack_Cut) { if(debugCOUT) std::cout << " > EtaTrack_Cut " << std::endl; continue;}
    if (gsfT->pt() < PtEl_Cut) { if(debugCOUT) std::cout << " < PtEl_Cut " << std::endl; continue;}
    
    if (SelectedMu_index != -1 ){
      if (fabs(ZvertexTrg - gsfT->vz()) > ElTrgMuDz_Cut ) { if(debugCOUT) std::cout << " vtxZ " << std::endl; continue;}
      if ( DR(gsfT->eta(),gsfT->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < ElTrgExclusionCone) {
	if(debugCOUT) std::cout << " dR trig " << std::endl;
	continue;
      }
      //if ( DR(el->eta(),el->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < ElTrgExclusionCone) continue;
    }

    if(debugCOUT) std::cout << " look at BDT " << std::endl;
    if(mvaSeeds.size() == 2 && mvaSeeds[0].isValid() && !mvaSeeds[0]->empty() && 
       mvaSeeds[1].isValid() && !mvaSeeds[1]->empty() && lowPtGsfTracks.isValid()){
      reco::GsfTrackRef gsf(lowPtGsfTracks, int(gsfT - lowPtGsfTracks->begin()));

      gsfTrk_seedBDTunb.push_back(float((*mvaSeeds[0])[gsf]));
      gsfTrk_seedBDTbiased.push_back(float((*mvaSeeds[1])[gsf]));
    }
    else continue;
    ++ngsfTracks;
    if(debugCOUT) std::cout << " gsfSalvata " << std::endl;
    gsfTrk_pt.push_back(gsfT->pt()); 
    gsfTrk_eta.push_back(gsfT->eta());
    gsfTrk_phi.push_back(gsfT->phi()); 
    gsfTrk_charge.push_back(gsfT->charge());
    
    auto gsfTrack = std::make_shared<reco::GsfTrack>(*gsfT);
    gsfElTracks.push_back(gsfTrack); 
    gsf_container.push_back(gsfeleIndex);
    ++gsfeleIndex; 
    gsf_id.push_back(11);
    //  delete eltrack;
  }
  if(debugCOUT) std::cout << " move to cleaned tracks " << std::endl;
  }

  std::vector<std::shared_ptr<reco::Track> > cleanedTracks; 
  std::vector<unsigned int> track_container;
  int trk_index = 0;
  if (!UseOnlyBKeeMCForTriplets){
    for (typename vector<reco::Track>::const_iterator trk=tracks->begin(); trk!=tracks->end(); trk++){
      if (!trk->quality(Track::highPurity)) continue;
      if (trk->pt() < PtTrack_Cut) continue;
      if (fabs(trk->eta()) > EtaTrack_Cut) continue;
      if (trk->charge() == 0) continue;
      if (trk->normalizedChi2() > MaxChi2Track_Cut || trk->normalizedChi2() < MinChi2Track_Cut) continue;
      if (fabs(trk->dxy())/trk->dxyError() < TrackSdxy_Cut) continue;

      // exclude tracks overlapping reco muons
      /*
      double minDR = 1000;
      for (typename vector<reco::Muon>::const_iterator mu = muons->begin(); mu!=muons->end(); mu++){
	double tempDR = DR(mu->eta(),mu->phi(),trk->eta(),trk->phi());
	if (tempDR < minDR) minDR = tempDR;
      }
      if (minDR < MuTrkMinDR_Cut) continue;
      */
      
      //parameters from the muon matched in dr to the trigger muon
      if (SelectedMu_index != -1 ){
	if (fabs(ZvertexTrg - trk->vz()) > TrkTrgMuDz_Cut ) continue;
	if ( DR(trk->eta(),trk->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < TrkTrgExclusionCone) continue;
      }

      //assignments
      track_pt.push_back(trk->pt()); track_eta.push_back(trk->eta());
      track_phi.push_back(trk->phi()); track_charge.push_back(trk->charge());
      track_norm_chi2.push_back(trk->normalizedChi2());
      track_dxy.push_back(trk->dxy(vertex_point)); track_dz.push_back(trk->dz(vertex_point));
      track_edxy.push_back(trk->dxyError()); track_edz.push_back(trk->dzError());  

      auto cleanTrack = std::make_shared<reco::Track>(*trk);
      cleanedTracks.push_back(cleanTrack);
      track_vx.push_back(trk->vx()); track_vy.push_back(trk->vy()); track_vz.push_back(trk->vz());

      track_container.push_back(trk_index);
      trk_index++; ntracks++;  

      //bdt   
      /*
      trk_pt =trk->pt(); trk_eta = trk->eta(); trk_phi = trk->phi(); trk_p =trk->p();
      trk_nhits = trk->found(); trk_high_purity = trk->quality(Track::highPurity);
      trk_chi2red = trk->normalizedChi2();
      float features[7]={trk_pt,trk_eta,trk_phi,trk_p,trk_nhits,trk_high_purity,trk_chi2red };
      float bdt_out = bdt.gbr->GetClassifier(features);
      track_mva.push_back(bdt_out);
      */
    }
  }//!use only BToKee MC for triplets 
  


  if (SaveOnlyTracks) {
    t1->Fill(); 
    if(debugCOUT) std::cout << " SaveOnlyTracks " << std::endl;
    return;
  }
  //create mother ee combination


  if(!lookAtpfEle){

    if(debugCOUT) std::cout << " now triplets " << std::endl;
    // fit track pairs  
    std::vector<std::shared_ptr<reco::GsfTrack> > cleanedObjTracks;
    std::vector<std::shared_ptr<reco::GsfTrack> > cleanedPairTracks;
  
    TLorentzVector vel1, vel2;
    std::vector<reco::TransientTrack> tempTracks;
    std::vector<float> tempPtEtaPhiM, tempXYZ;

    std::vector<std::shared_ptr<reco::GsfTrack> > cleanedObjects; 
    for(auto & vec: gsfElTracks) cleanedObjects.push_back(vec);  

    if (cleanedObjects.size() == 0) { if(debugCOUT) std::cout << " cleanedObjects.size() == 0 " << "gsfElTracks.size = " << gsfElTracks.size() << std::endl; return;}

    for(unsigned int iobj=0; iobj<cleanedObjects.size(); ++iobj){
      auto obj = cleanedObjects.at(iobj);

      if(LeptonFinalStateID == 13 && gsf_id.at(iobj) != 13) continue;
      else if (LeptonFinalStateID == 11 && gsf_id.at(iobj) != 11) continue;

      if(gsfTrk_seedBDTunb.at(iobj) < 3.5) continue;


      //the following two for lepton+track+track
      //for(unsigned int itrk2=0; itrk2<cleanedTracks.size(); itrk2++){
      //	auto trk2 = cleanedTracks.at(itrk2);
      //the following two for lepton+lepton+track
      for(unsigned int itrk2=0; itrk2<cleanedObjects.size(); itrk2++){
      	auto trk2 = cleanedObjects.at(itrk2);
	if (obj->charge()*trk2->charge() == 1) continue;
	if (ObjPtLargerThanTrack && obj->pt() < trk2->pt()) continue;

	float dR_l1l2 = DR(obj->eta(), obj->phi(), trk2->eta(), trk2->phi());
	if(dR_l1l2 < TrkObjExclusionCone && (isKll || LeptonFinalStateID == 13)) continue;
	// if(!isKll && LeptonFinalStateID == 11 && 1.*sharedHits((*obj), (*trk2)) / obj->numberOfValidHits() > 0.5){
	// 	std::cout << " >> fraction shared = " << 1.*sharedHits((*obj), (*trk2))/ obj->numberOfValidHits()  << std::endl;
	// }
	
	vel1.SetPtEtaPhiM(obj->pt(), obj->eta(), obj->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	vel2.SetPtEtaPhiM(trk2->pt(), trk2->eta(), trk2->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);


	if(SelectedMu_index != -1 && LeptonFinalStateID == 13){
	  if (gsf_id.at(iobj) == 13 && DR(obj->eta(),obj->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < MuTrgExclusionCone) continue;
	  if (gsf_id.at(iobj)==13 && fabs(ZvertexTrg- obj->vz()) > MuTrgMuDz_Cut ) continue;
      }
      
	//std::cout << " object_id.at(iobj) = " << object_id.at(iobj) << std::endl;
	
	//inv mass on lepton-track pair
	if ((vel1+vel2).M() > MaxMee_Cut || (vel1+vel2).M() < MinMee_Cut ) continue;   
	
	auto tranobj = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*obj,&(*bFieldHandle)));
	auto trantrk2 = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk2,&(*bFieldHandle)));
	tempTracks.clear(); 
	tempTracks.push_back(*tranobj); tempTracks.push_back(*trantrk2);

	LLvertex = theKalmanFitter.vertex(tempTracks);
	if (!LLvertex.isValid()) continue;

	if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < Probee_Cut)  continue;
	if (SelectedMu_index != -1 && fabs(ZvertexTrg-LLvertex.position().z()) > EpairZvtx_Cut ) continue;

	GlobalError err = LLvertex.positionError();
	GlobalPoint Dispbeamspot(-1*( (theBeamSpot->x0() - LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()), 
				 -1*( (theBeamSpot->y0() - LLvertex.position().y())+ (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);
	
	math::XYZVector pperp((vel1+vel2).Px(),(vel1+vel2).Py(),0);
	math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
	float tempCos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
	if (tempCos < Cosee_Cut) continue;
	
	cleanedObjTracks.push_back(obj);
	cleanedPairTracks.push_back(trk2);
	Epair_ObjIndex.push_back(gsf_container.at(iobj));
	//if lep + track+track
	//Epair_TrkIndex.push_back(track_container.at(itrk2));
	//if lep + lep +track
	Epair_TrkIndex.push_back(gsf_container.at(itrk2));
	Epair_ObjId.push_back(gsf_id.at(iobj));   

	tempPtEtaPhiM.clear(); tempXYZ.clear();
	tempPtEtaPhiM.push_back((vel1+vel2).Pt()); tempPtEtaPhiM.push_back((vel1+vel2).Eta());
	tempPtEtaPhiM.push_back((vel1+vel2).Phi()); tempPtEtaPhiM.push_back((vel1+vel2).M());
	tempXYZ.push_back(LLvertex.position().x());  tempXYZ.push_back(LLvertex.position().y()); 
	tempXYZ.push_back(LLvertex.position().z());
	Epair_PtEtaPhiM.push_back(tempPtEtaPhiM); Epair_XYZ.push_back(tempXYZ); 
	Epair_cos.push_back(tempCos);
	Epair_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom())); 	
	Epair_Lxy.push_back(Dispbeamspot.perp());
	Epair_eLxy.push_back(err.rerr(Dispbeamspot));      
      }// cleanedTracks    
    }//cleaned objects
     

    if (SaveOnlyEPairTracks) {
      t1->Fill(); 
      if(debugCOUT) std::cout << " SaveOnlyEPairTracks " << std::endl;
      return;
    }

    // triplet
    TLorentzVector vK; 
    TLorentzVector vPi;
    int kstarIndex = 0;
    for(unsigned int iobj=0; iobj<cleanedObjTracks.size(); iobj++){
      auto objtrk = cleanedObjTracks.at(iobj);
      auto pairtrk = cleanedPairTracks.at(iobj);
 
      for(unsigned int itrk=0; itrk<cleanedTracks.size(); itrk++){
	auto trk = cleanedTracks.at(itrk);

	if(DR(objtrk->eta(), objtrk->phi(), trk->eta(),trk->phi()) < TrkObjExclusionCone) continue;
	if(DR(pairtrk->eta(), pairtrk->phi(), trk->eta(),trk->phi()) < TrkObjExclusionCone) continue;

	if (trk->pt() < PtKTrack_Cut) continue;
         
	//isKll
	//if(isKll || !isKll){
	if(isKll){
	  if (fabs(trk->dxy(vertex_point))/trk->dxyError() < Ksdxy_Cut) continue;

	  //ele ele kaon
	  vel1.SetPtEtaPhiM(objtrk->pt(),objtrk->eta(),objtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	  vel2.SetPtEtaPhiM(pairtrk->pt(),pairtrk->eta(),pairtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	  vK.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(), KaonMass_);
	  
	  if ((vel1+vel2+vK).M() > MaxMB_Cut || (vel1+vel2+vK).M() < MinMB_Cut) continue;
	  if ((vel1+vel2+vK).Pt() < PtB_Cut) continue;
	  
	  auto tranobj = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*objtrk,&(*bFieldHandle)));
	  auto tranpair = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*pairtrk,&(*bFieldHandle)));
	  auto trantrk = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk,&(*bFieldHandle)));
	  tempTracks.clear();
	  tempTracks.push_back(*tranobj); 
	  tempTracks.push_back(*tranpair);
	  tempTracks.push_back(*trantrk);
	  
	  LLvertex = theKalmanFitter.vertex(tempTracks);
	  if (!LLvertex.isValid()) continue;
	  
	  if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < ProbeeK_Cut) continue;
	  GlobalError err = LLvertex.positionError();
	  GlobalPoint Dispbeamspot( -1 * ((theBeamSpot->x0()-LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),
				    -1 * ((theBeamSpot->y0()-LLvertex.position().y()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);
	  
	  math::XYZVector pperp((vel1+vel2+vK).Px(),(vel1+vel2+vK).Py(), 0);
	  math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(), 0.);
	  float tempCos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
	  if (tempCos < CoseeK_Cut) continue;
	  if (SLxy_Cut > Dispbeamspot.perp()/TMath::Sqrt(err.rerr(Dispbeamspot))) continue;
	  
	  //std::cout << " found triplet " << std::endl;
	  
	  tempPtEtaPhiM.clear(); tempXYZ.clear();
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).Pt());
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).Eta()); 
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).Phi());
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).M());       
	  TTrack_PtEtaPhiM.push_back(tempPtEtaPhiM);
	  tempXYZ.push_back(LLvertex.position().x());
	  tempXYZ.push_back(LLvertex.position().y());
	  tempXYZ.push_back(LLvertex.position().z());
	  TTrack_ObjId.push_back(Epair_ObjId.at(iobj));
	  TTrack_XYZ.push_back(tempXYZ); 
	  TTrack_mll.push_back((vel1+vel2).M());
	  TTrack_ObjIndex.push_back(Epair_ObjIndex.at(iobj)); 
	  TTrack_TrkIndex.push_back(Epair_TrkIndex.at(iobj)); 
	  TTrack_kid.push_back(track_container.at(itrk));
	  TTrack_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()));
	  TTrack_cos.push_back(tempCos);
	  TTrack_Lxy.push_back(Dispbeamspot.perp());
	  TTrack_eLxy.push_back(err.rerr(Dispbeamspot));
	  if (EarlyStop) break;
	}//isKll
	else{//K*ll
	//	if(1 == 2){
	  //auto objtrk = cleanedObjTracks.at(iobj);
	  //auto pairtrk = cleanedPairTracks.at(iobj);
	  //auto trk = cleanedTracks.at(itrk);

	  for(unsigned int iPi=0; iPi<cleanedTracks.size(); ++iPi){
	    if(iPi == itrk) continue;
	    auto PiTrk = cleanedTracks.at(iPi); 

	    if (PiTrk->charge() * trk->charge() == 1) continue;

	    if(DR(objtrk->eta(), objtrk->phi(), PiTrk->eta(),PiTrk->phi()) < TrkObjExclusionCone) continue;
	    if(DR(pairtrk->eta(), pairtrk->phi(), PiTrk->eta(),PiTrk->phi()) < TrkObjExclusionCone) continue;
	    if(DR(trk->eta(), trk->phi(), PiTrk->eta(),PiTrk->phi()) < TrkObjExclusionCone) continue;

	    if(debugCOUT) std::cout << " loop 4th track " << std::endl;

	    if (PiTrk->pt() < PtKTrack_Cut) continue;
	    if (std::fabs(PiTrk->eta()) > EtaTrack_Cut) continue;
	    //maybe yes and also for K track but avoid here
	    //if (fabs(PiTrk->dxy(vertex_point))/PiTrk->dxyError() < Ksdxy_Cut) continue;


	    vK.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(), KaonMass_);
	    vPi.SetPtEtaPhiM(PiTrk->pt(),PiTrk->eta(),PiTrk->phi(), PionMass_);

	    if(debugCOUT) std::cout << " K* built " << std::endl;

	    if ((vPi+vK).M() > MaxKst_Cut || (vPi+vK).M() < MinKst_Cut) continue;
	    if ((vPi+vK).Pt() < PtKst_Cut) continue;

	    auto tranobj = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*objtrk,&(*bFieldHandle)));
	    auto tranpair = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*pairtrk,&(*bFieldHandle)));
	    auto trantrk = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk,&(*bFieldHandle)));
	    auto trantrpi = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*PiTrk,&(*bFieldHandle)));
	    tempTracks.clear();
	    // tempTracks.push_back(*tranobj);
	    // tempTracks.push_back(*tranpair);
	    tempTracks.push_back(*trantrk);
	    tempTracks.push_back(*trantrpi);

	    if(debugCOUT) std::cout << " check K* vertex " << std::endl;

	    LLvertex = theKalmanFitter.vertex(tempTracks);
	    if (!LLvertex.isValid()) continue;

	    if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < ProbKst_Cut) continue;
	    GlobalError err = LLvertex.positionError();
	    GlobalPoint Dispbeamspot( -1 * ((theBeamSpot->x0()-LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),
				      -1 * ((theBeamSpot->y0()-LLvertex.position().y()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);

	    math::XYZVector pperp((vPi+vK).Px(),(vPi+vK).Py(), 0);
	    math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(), 0.);
	    float tempCos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
	    if (tempCos < CosKst_Cut) continue;
	    if (SLxyKst_Cut > Dispbeamspot.perp()/TMath::Sqrt(err.rerr(Dispbeamspot))) continue;

	    if(debugCOUT) std::cout << " fill K* " << std::endl;

	    tempPtEtaPhiM.clear();
            tempPtEtaPhiM.push_back((vPi+vK).Pt());
            tempPtEtaPhiM.push_back((vPi+vK).Eta());
            tempPtEtaPhiM.push_back((vPi+vK).Phi());
            tempPtEtaPhiM.push_back((vPi+vK).M());

	    Kstpair_PtEtaPhiM.push_back(tempPtEtaPhiM); 
	    Kstpair_cos.push_back(tempCos);
	    Kstpair_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()));
	    Kstpair_Lxy.push_back(Dispbeamspot.perp());
	    Kstpair_eLxy.push_back(err.rerr(Dispbeamspot));
	    ++kstarIndex;

	    vel1.SetPtEtaPhiM(objtrk->pt(),objtrk->eta(),objtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	    vel2.SetPtEtaPhiM(pairtrk->pt(),pairtrk->eta(),pairtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);

	    if(debugCOUT) std::cout << " building eeK* mass =  " << (vPi+vK+vel1+vel2).M() << " pT = " << (vel1+vel2+vPi+vK).Pt() << std::endl;

	    if ((vel1+vel2+vPi+vK).M() > MaxBeeKst_Cut || (vPi+vK+vel1+vel2).M() < MinBeeKst_Cut) continue;
            if ((vel1+vel2+vPi+vK).Pt() < PtBeeKst_Cut) continue;

	    if(debugCOUT) std::cout << " eeK* post pT and M  " << std::endl;

	    ///if ok 4tracks fit
	    tempTracks.clear();
	    tempTracks.push_back(*tranobj);
	    tempTracks.push_back(*tranpair);
	    tempTracks.push_back(*trantrk);
	    tempTracks.push_back(*trantrpi);

	    if(debugCOUT) std::cout << " before eeK* vertex  " << std::endl;

	    LLvertex = theKalmanFitter.vertex(tempTracks);
	    if (!LLvertex.isValid()) continue;

	    if(debugCOUT) std::cout << " check eeK* vertex " << std::endl;

	    if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < ProbBeeKst_Cut) continue;
	    GlobalError errB = LLvertex.positionError();
	    GlobalPoint DispbeamspotB( -1 * ((theBeamSpot->x0()-LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),
				      -1 * ((theBeamSpot->y0()-LLvertex.position().y()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);

	    math::XYZVector pperpB((vPi+vel1+vel2+vK).Px(),(vel1+vel2+vPi+vK).Py(), 0);
	    math::XYZVector vperpB(DispbeamspotB.x(),DispbeamspotB.y(), 0.);
	    float tempCosB = vperpB.Dot(pperp)/(vperpB.R()*pperp.R());
	    if (tempCosB < CosBeeKst_Cut) continue;
	    if (SLxyBeeKst_Cut > DispbeamspotB.perp()/TMath::Sqrt(errB.rerr(DispbeamspotB))) continue;

	    if(debugCOUT) std::cout << " filling eeK* " << std::endl;

	    tempPtEtaPhiM.clear(); 
	    tempXYZ.clear();
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).Pt());
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).Eta());
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).Phi());
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).M());
	    TTrack_PtEtaPhiM.push_back(tempPtEtaPhiM);

	    tempXYZ.push_back(LLvertex.position().x());
	    tempXYZ.push_back(LLvertex.position().y());
	    tempXYZ.push_back(LLvertex.position().z());
	    TTrack_ObjId.push_back(Epair_ObjId.at(iobj));
	    TTrack_XYZ.push_back(tempXYZ); 
	    TTrack_mll.push_back((vel1+vel2).M());
	    TTrack_mKst.push_back((vPi+vK).M());

	    TTrack_ObjIndex.push_back(Epair_ObjIndex.at(iobj));
	    TTrack_TrkIndex.push_back(Epair_TrkIndex.at(iobj));
	    TTrack_KstarIndex.push_back(kstarIndex-1);
	    TTrack_kid.push_back(track_container.at(itrk));
	    TTrack_piid.push_back(track_container.at(iPi));
	  
	    TTrack_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()));
	    TTrack_cos.push_back(tempCosB);
	    TTrack_Lxy.push_back(DispbeamspotB.perp());
	    TTrack_eLxy.push_back(errB.rerr(DispbeamspotB));

	  if (EarlyStop) break;
	  }//4th track
	}//Kst
      }// tracks 3rd 
  }// objects l1 and l2      
  
  }// lookAtGsfTracks

  
  if(lookAtpfEle){
  // fit track pairs  
  std::vector<std::shared_ptr<reco::Track> > cleanedObjTracks;
  std::vector<std::shared_ptr<reco::Track> > cleanedPairTracks;
  //  cleanedObjTracks.clear(); cleanedPairTracks.clear(); 
  
  TLorentzVector vel1, vel2;
  std::vector<reco::TransientTrack> tempTracks;
  std::vector<float> tempPtEtaPhiM, tempXYZ;

  std::vector<std::shared_ptr<reco::Track> > cleanedObjects; 
  for(auto & vec: MuTracks) cleanedObjects.push_back(vec);
  for(auto & vec: ElTracks) cleanedObjects.push_back(vec);  

  if (cleanedObjects.size() == 0) { if(debugCOUT) std::cout << " cleanedObjects.size() == 0 " << "ElTracks.size = " << ElTracks.size() << std::endl; return;}

  for(unsigned int iobj=0; iobj<cleanedObjects.size(); iobj++){
    auto obj = cleanedObjects.at(iobj);

    if(LeptonFinalStateID == 13 && object_id.at(iobj) != 13) continue;
    else if (LeptonFinalStateID == 11 && object_id.at(iobj) != 11) continue;

    for(unsigned int itrk2=0; itrk2<cleanedTracks.size(); itrk2++){
      auto trk2 = cleanedTracks.at(itrk2);
      if (obj->charge()*trk2->charge() == 1) continue;
      if (ObjPtLargerThanTrack && obj->pt() < trk2->pt()) continue;

      float dR_l1l2 = DR(obj->eta(), obj->phi(), trk2->eta(), trk2->phi());
      if(dR_l1l2 < TrkObjExclusionCone && (isKll || LeptonFinalStateID == 13)) continue;
      // if(!isKll && LeptonFinalStateID == 11 && 1.*sharedHits((*obj), (*trk2)) / obj->numberOfValidHits() > 0.5){
      // 	std::cout << " >> fraction shared = " << 1.*sharedHits((*obj), (*trk2))/ obj->numberOfValidHits()  << std::endl;
      // }

      vel1.SetPtEtaPhiM(obj->pt(), obj->eta(), obj->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
      vel2.SetPtEtaPhiM(trk2->pt(), trk2->eta(), trk2->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);


      if(SelectedMu_index != -1 && LeptonFinalStateID == 13){
	if (object_id.at(iobj) == 13 && DR(obj->eta(),obj->phi(), SelectedTrgObj_PtEtaPhiCharge[1],SelectedTrgObj_PtEtaPhiCharge[2]) < MuTrgExclusionCone) continue;
	if (object_id.at(iobj)==13 && fabs(ZvertexTrg- obj->vz()) > MuTrgMuDz_Cut ) continue;
      }
      
      //std::cout << " object_id.at(iobj) = " << object_id.at(iobj) << std::endl;

      //inv mass on lepton-track pair
      if ((vel1+vel2).M() > MaxMee_Cut || (vel1+vel2).M() < MinMee_Cut ) continue;   

      auto tranobj = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*obj,&(*bFieldHandle)));
      auto trantrk2 = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk2,&(*bFieldHandle)));
      tempTracks.clear(); 
      tempTracks.push_back(*tranobj); tempTracks.push_back(*trantrk2);

      LLvertex = theKalmanFitter.vertex(tempTracks);
      if (!LLvertex.isValid()) continue;

      if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < Probee_Cut)  continue;
      if (SelectedMu_index != -1 && fabs(ZvertexTrg-LLvertex.position().z()) > EpairZvtx_Cut ) continue;

      GlobalError err = LLvertex.positionError();
      GlobalPoint Dispbeamspot(-1*( (theBeamSpot->x0() - LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()), 
			       -1*( (theBeamSpot->y0() - LLvertex.position().y())+ (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);

      math::XYZVector pperp((vel1+vel2).Px(),(vel1+vel2).Py(),0);
      math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(),0.);
      float tempCos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
      if (tempCos < Cosee_Cut) continue;

      cleanedObjTracks.push_back(obj);
      cleanedPairTracks.push_back(trk2);
      Epair_ObjIndex.push_back(object_container.at(iobj));
      Epair_TrkIndex.push_back(track_container.at(itrk2));
      Epair_ObjId.push_back(object_id.at(iobj));   
      tempPtEtaPhiM.clear(); tempXYZ.clear();
      tempPtEtaPhiM.push_back((vel1+vel2).Pt()); tempPtEtaPhiM.push_back((vel1+vel2).Eta());
      tempPtEtaPhiM.push_back((vel1+vel2).Phi()); tempPtEtaPhiM.push_back((vel1+vel2).M());
      tempXYZ.push_back(LLvertex.position().x());  tempXYZ.push_back(LLvertex.position().y()); 
      tempXYZ.push_back(LLvertex.position().z());
      Epair_PtEtaPhiM.push_back(tempPtEtaPhiM); Epair_XYZ.push_back(tempXYZ); 
      Epair_cos.push_back(tempCos);
      Epair_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom())); 	
      Epair_Lxy.push_back(Dispbeamspot.perp());
      Epair_eLxy.push_back(err.rerr(Dispbeamspot));      
    }// cleanedTracks    
  }//cleaned objects
     


  if (SaveOnlyEPairTracks) {
    t1->Fill(); 
    if(debugCOUT) std::cout << " SaveOnlyEPairTracks " << std::endl;
    return;
  }


  // triplet
  TLorentzVector vK; 
  TLorentzVector vPi;
  int kstarIndex = 0;
  for(unsigned int iobj=0; iobj<cleanedObjTracks.size(); iobj++){
    auto objtrk = cleanedObjTracks.at(iobj);
    auto pairtrk = cleanedPairTracks.at(iobj);
 
      for(unsigned int itrk=0; itrk<cleanedTracks.size(); itrk++){
	auto trk = cleanedTracks.at(itrk);

	if(DR(objtrk->eta(), objtrk->phi(), trk->eta(),trk->phi()) < TrkObjExclusionCone) continue;
	if(DR(pairtrk->eta(), pairtrk->phi(), trk->eta(),trk->phi()) < TrkObjExclusionCone) continue;

	if (trk->pt() < PtKTrack_Cut) continue;
         
	//isKll
	//if(isKll || !isKll){
	if(isKll){
	  if (fabs(trk->dxy(vertex_point))/trk->dxyError() < Ksdxy_Cut) continue;

	  //ele ele kaon
	  vel1.SetPtEtaPhiM(objtrk->pt(),objtrk->eta(),objtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	  vel2.SetPtEtaPhiM(pairtrk->pt(),pairtrk->eta(),pairtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	  vK.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(), KaonMass_);
	  
	  if ((vel1+vel2+vK).M() > MaxMB_Cut || (vel1+vel2+vK).M() < MinMB_Cut) continue;
	  if ((vel1+vel2+vK).Pt() < PtB_Cut) continue;
	  
	  auto tranobj = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*objtrk,&(*bFieldHandle)));
	  auto tranpair = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*pairtrk,&(*bFieldHandle)));
	  auto trantrk = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk,&(*bFieldHandle)));
	  tempTracks.clear();
	  tempTracks.push_back(*tranobj); 
	  tempTracks.push_back(*tranpair);
	  tempTracks.push_back(*trantrk);
	  
	  LLvertex = theKalmanFitter.vertex(tempTracks);
	  if (!LLvertex.isValid()) continue;
	  
	  if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < ProbeeK_Cut) continue;
	  GlobalError err = LLvertex.positionError();
	  GlobalPoint Dispbeamspot( -1 * ((theBeamSpot->x0()-LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),
				    -1 * ((theBeamSpot->y0()-LLvertex.position().y()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);
	  
	  math::XYZVector pperp((vel1+vel2+vK).Px(),(vel1+vel2+vK).Py(), 0);
	  math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(), 0.);
	  float tempCos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
	  if (tempCos < CoseeK_Cut) continue;
	  if (SLxy_Cut > Dispbeamspot.perp()/TMath::Sqrt(err.rerr(Dispbeamspot))) continue;
	  
	  //std::cout << " found triplet " << std::endl;
	  
	  tempPtEtaPhiM.clear(); tempXYZ.clear();
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).Pt());
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).Eta()); 
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).Phi());
	  tempPtEtaPhiM.push_back((vel1+vel2+vK).M());       
	  TTrack_PtEtaPhiM.push_back(tempPtEtaPhiM);
	  tempXYZ.push_back(LLvertex.position().x());
	  tempXYZ.push_back(LLvertex.position().y());
	  tempXYZ.push_back(LLvertex.position().z());
	  TTrack_ObjId.push_back(Epair_ObjId.at(iobj));
	  TTrack_XYZ.push_back(tempXYZ); TTrack_mll.push_back((vel1+vel2).M());
	  TTrack_ObjIndex.push_back(Epair_ObjIndex.at(iobj)); 
	  TTrack_TrkIndex.push_back(Epair_TrkIndex.at(iobj)); 
	  TTrack_kid.push_back(track_container.at(itrk));
	  TTrack_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()));
	  TTrack_cos.push_back(tempCos);
	  TTrack_Lxy.push_back(Dispbeamspot.perp());
	  TTrack_eLxy.push_back(err.rerr(Dispbeamspot));
	  if (EarlyStop) break;
	}//isKll
	else{//K*ll
	//	if(1 == 2){
	  //auto objtrk = cleanedObjTracks.at(iobj);
	  //auto pairtrk = cleanedPairTracks.at(iobj);
	  //auto trk = cleanedTracks.at(itrk);

	  for(unsigned int iPi=0; iPi<cleanedTracks.size(); ++iPi){
	    if(iPi == itrk) continue;
	    auto PiTrk = cleanedTracks.at(iPi); 

	    if (PiTrk->charge() * trk->charge() == 1) continue;

	    if(DR(objtrk->eta(), objtrk->phi(), PiTrk->eta(),PiTrk->phi()) < TrkObjExclusionCone) continue;
	    if(DR(pairtrk->eta(), pairtrk->phi(), PiTrk->eta(),PiTrk->phi()) < TrkObjExclusionCone) continue;
	    if(DR(trk->eta(), trk->phi(), PiTrk->eta(),PiTrk->phi()) < TrkObjExclusionCone) continue;

	    if(debugCOUT) std::cout << " loop 4th track " << std::endl;

	    if (PiTrk->pt() < PtKTrack_Cut) continue;
	    if (std::fabs(PiTrk->eta()) > EtaTrack_Cut) continue;
	    //maybe yes and also for K track but avoid here
	    //if (fabs(PiTrk->dxy(vertex_point))/PiTrk->dxyError() < Ksdxy_Cut) continue;


	    vK.SetPtEtaPhiM(trk->pt(),trk->eta(),trk->phi(), KaonMass_);
	    vPi.SetPtEtaPhiM(PiTrk->pt(),PiTrk->eta(),PiTrk->phi(), PionMass_);

	    if(debugCOUT) std::cout << " K* built " << std::endl;

	    if ((vPi+vK).M() > MaxKst_Cut || (vPi+vK).M() < MinKst_Cut) continue;
	    if ((vPi+vK).Pt() < PtKst_Cut) continue;

	    auto tranobj = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*objtrk,&(*bFieldHandle)));
	    auto tranpair = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*pairtrk,&(*bFieldHandle)));
	    auto trantrk = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*trk,&(*bFieldHandle)));
	    auto trantrpi = std::make_shared<reco::TransientTrack>(reco::TransientTrack(*PiTrk,&(*bFieldHandle)));
	    tempTracks.clear();
	    // tempTracks.push_back(*tranobj);
	    // tempTracks.push_back(*tranpair);
	    tempTracks.push_back(*trantrk);
	    tempTracks.push_back(*trantrpi);

	    if(debugCOUT) std::cout << " check K* vertex " << std::endl;

	    LLvertex = theKalmanFitter.vertex(tempTracks);
	    if (!LLvertex.isValid()) continue;

	    if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < ProbKst_Cut) continue;
	    GlobalError err = LLvertex.positionError();
	    GlobalPoint Dispbeamspot( -1 * ((theBeamSpot->x0()-LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),
				      -1 * ((theBeamSpot->y0()-LLvertex.position().y()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);

	    math::XYZVector pperp((vPi+vK).Px(),(vPi+vK).Py(), 0);
	    math::XYZVector vperp(Dispbeamspot.x(),Dispbeamspot.y(), 0.);
	    float tempCos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
	    if (tempCos < CosKst_Cut) continue;
	    if (SLxyKst_Cut > Dispbeamspot.perp()/TMath::Sqrt(err.rerr(Dispbeamspot))) continue;

	    if(debugCOUT) std::cout << " fill K* " << std::endl;

	    tempPtEtaPhiM.clear();
            tempPtEtaPhiM.push_back((vPi+vK).Pt());
            tempPtEtaPhiM.push_back((vPi+vK).Eta());
            tempPtEtaPhiM.push_back((vPi+vK).Phi());
            tempPtEtaPhiM.push_back((vPi+vK).M());

	    Kstpair_PtEtaPhiM.push_back(tempPtEtaPhiM); 
	    Kstpair_cos.push_back(tempCos);
	    Kstpair_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()));
	    Kstpair_Lxy.push_back(Dispbeamspot.perp());
	    Kstpair_eLxy.push_back(err.rerr(Dispbeamspot));
	    ++kstarIndex;

	    vel1.SetPtEtaPhiM(objtrk->pt(),objtrk->eta(),objtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);
	    vel2.SetPtEtaPhiM(pairtrk->pt(),pairtrk->eta(),pairtrk->phi(), (LeptonFinalStateID == 11) ? ElectronMass_ : MuonMass_);

	    if(debugCOUT) std::cout << " building eeK* mass =  " << (vPi+vK+vel1+vel2).M() << " pT = " << (vel1+vel2+vPi+vK).Pt() << std::endl;

	    if ((vel1+vel2+vPi+vK).M() > MaxBeeKst_Cut || (vPi+vK+vel1+vel2).M() < MinBeeKst_Cut) continue;
            if ((vel1+vel2+vPi+vK).Pt() < PtBeeKst_Cut) continue;

	    if(debugCOUT) std::cout << " eeK* post pT and M  " << std::endl;

	    ///if ok 4tracks fit
	    tempTracks.clear();
	    tempTracks.push_back(*tranobj);
	    tempTracks.push_back(*tranpair);
	    tempTracks.push_back(*trantrk);
	    tempTracks.push_back(*trantrpi);

	    if(debugCOUT) std::cout << " before eeK* vertex  " << std::endl;

	    LLvertex = theKalmanFitter.vertex(tempTracks);
	    if (!LLvertex.isValid()) continue;

	    if(debugCOUT) std::cout << " check eeK* vertex " << std::endl;

	    if (ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()) < ProbBeeKst_Cut) continue;
	    GlobalError errB = LLvertex.positionError();
	    GlobalPoint DispbeamspotB( -1 * ((theBeamSpot->x0()-LLvertex.position().x()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dxdz()),
				      -1 * ((theBeamSpot->y0()-LLvertex.position().y()) + (LLvertex.position().z()-theBeamSpot->z0()) * theBeamSpot->dydz()), 0);

	    math::XYZVector pperpB((vPi+vel1+vel2+vK).Px(),(vel1+vel2+vPi+vK).Py(), 0);
	    math::XYZVector vperpB(DispbeamspotB.x(),DispbeamspotB.y(), 0.);
	    float tempCosB = vperpB.Dot(pperp)/(vperpB.R()*pperp.R());
	    if (tempCosB < CosBeeKst_Cut) continue;
	    if (SLxyBeeKst_Cut > DispbeamspotB.perp()/TMath::Sqrt(errB.rerr(DispbeamspotB))) continue;

	    if(debugCOUT) std::cout << " filling eeK* " << std::endl;

	    tempPtEtaPhiM.clear(); 
	    tempXYZ.clear();
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).Pt());
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).Eta());
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).Phi());
	    tempPtEtaPhiM.push_back((vel1+vel2+vPi+vK).M());
	    TTrack_PtEtaPhiM.push_back(tempPtEtaPhiM);

	    tempXYZ.push_back(LLvertex.position().x());
	    tempXYZ.push_back(LLvertex.position().y());
	    tempXYZ.push_back(LLvertex.position().z());
	    TTrack_ObjId.push_back(Epair_ObjId.at(iobj));
	    TTrack_XYZ.push_back(tempXYZ); 
	    TTrack_mll.push_back((vel1+vel2).M());
	    TTrack_mKst.push_back((vPi+vK).M());

	    TTrack_ObjIndex.push_back(Epair_ObjIndex.at(iobj));
	    TTrack_TrkIndex.push_back(Epair_TrkIndex.at(iobj));
	    TTrack_KstarIndex.push_back(kstarIndex-1);
	    TTrack_kid.push_back(track_container.at(itrk));
	    TTrack_piid.push_back(track_container.at(iPi));
	    TTrack_chi_prob.push_back(ChiSquaredProbability(LLvertex.totalChiSquared(),LLvertex.degreesOfFreedom()));
	    TTrack_cos.push_back(tempCosB);
	    TTrack_Lxy.push_back(DispbeamspotB.perp());
	    TTrack_eLxy.push_back(errB.rerr(DispbeamspotB));

	  if (EarlyStop) break;
	  }//4th track
	}//Kst
      }// tracks 3rd 
  }// objects l1 and l2      
  }//lookAtPfEle

  //  if (TTrack_chi_prob.size() == 0) Init();

  if(debugCOUT) std::cout << " filling histo genpart_B_index = " << genpart_B_index << " TTrack_cos.size() = " << TTrack_cos.size() << std::endl;
  t1->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void SkimAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void SkimAnalyzer::endJob() 
{
  // t1->Fill();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SkimAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
//typedef SkimAnalyzer<reco::RecoEcalCandidate> SkimAnalyzerb;
//DEFINE_FWK_MODULE(SkimAnalyzerb);
DEFINE_FWK_MODULE(SkimAnalyzer);

