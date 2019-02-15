#ifndef skim_class_h
#define skim_class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class skim_class {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          run_number;
   UInt_t          ls;
   vector<float>   *vertex_x;
   vector<float>   *vertex_y;
   vector<float>   *vertex_z;
   Float_t         beam_x;
   Float_t         beam_y;
   Float_t         beam_z;
   Int_t           HLT_path1;
   Int_t           HLT_path2;
   Int_t           HLT_path3;
   Int_t           HLT_path4;
   Int_t           HLT_path5;
   Int_t           HLT_path6;
   vector<vector<float> > *TrgObj1_PtEtaPhiCharge;
   vector<vector<float> > *TrgObj2_PtEtaPhiCharge;
   vector<vector<float> > *TrgObj3_PtEtaPhiCharge;
   vector<vector<float> > *TrgObj4_PtEtaPhiCharge;
   vector<vector<float> > *TrgObj5_PtEtaPhiCharge;
   vector<vector<float> > *TrgObj6_PtEtaPhiCharge;
   vector<float>   *SelectedTrgObj_PtEtaPhiCharge;
   Int_t           SelectedMu_index;
   Float_t         SelectedMu_DR;
   UInt_t          nmuon;
   vector<float>   *muon_pt;
   vector<float>   *muon_eta;
   vector<float>   *muon_phi;
   vector<float>   *muon_charge;
   vector<float>   *muon_dxy;
   vector<float>   *muon_dz;
   vector<float>   *muon_edxy;
   vector<float>   *muon_edz;
   vector<float>   *muon_d0;
   vector<float>   *muon_ed0;
   vector<float>   *muon_vx;
   vector<float>   *muon_vy;
   vector<float>   *muon_vz;
   vector<float>   *muon_iso;
   vector<bool>    *muon_soft;
   vector<bool>    *muon_loose;
   vector<bool>    *muon_medium;
   vector<bool>    *muon_tight;
   vector<float>   *muon_trkpt;
   vector<float>   *muon_trketa;
   vector<float>   *muon_trkphi;
   UInt_t          nelectron;
   vector<float>   *el_pt;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<float>   *el_charge;
   vector<float>   *el_dxy;
   vector<float>   *el_dz;
   vector<float>   *el_edxy;
   vector<float>   *el_edz;
   vector<float>   *el_vx;
   vector<float>   *el_vy;
   vector<float>   *el_vz;
   vector<float>   *el_mvaB_wp;
   vector<float>   *el_mvaUnB_wp;
   vector<float>   *el_iso;
   vector<float>   *el_mvaId;
   vector<float>   *el_trkphi;
   vector<float>   *el_trkpt;
   vector<float>   *el_trketa;

   UInt_t ngsfTracks;
   vector<float>   *gsfTrk_pt;
   vector<float>   *gsfTrk_eta;
   vector<float>   *gsfTrk_phi;
   vector<float>   *gsfTrk_charge;
   vector<float>   *gsfTrk_seedBDTunb;
   vector<float>   *gsfTrk_seedBDTbiased;


   /* vector<float>   *genpart_pt; */
   /* vector<float>   *genpart_phi; */
   /* vector<float>   *genpart_eta; */
   /* vector<float>   *genpart_pdgId; */
   /* vector<float>   *genpart_Bindex; */
   /* vector<float>   *genpart_Daughtindex; */
   /* vector<float>   *genpart_charge; */
   /* vector<float>   *genpart_mother_pdgId; */
   /* vector<float>   *genpart_mother_pt; */
   /* vector<float>   *genpart_mother_phi; */
   /* vector<float>   *genpart_mother_eta; */
   /* vector<float>   *genpart_mother_Bindex; */
   /* vector<float>   *genpart_mother_Daughtindex; */
   /* vector<float>   *genpart_grandmother_pdgId; */
   /* vector<float>   *genpart_grandmother_pt; */
   /* vector<float>   *genpart_grandmother_phi; */
   /* vector<float>   *genpart_grandmother_eta; */
   /* vector<float>   *genpart_grandmother_Bindex; */
   /* vector<float>   *genpart_grandmother_x; */
   /* vector<float>   *genpart_grandmother_y; */
   /* vector<float>   *genpart_grandmother_z; */
   /* vector<float>   *genMu_pt; */
   /* vector<float>   *genMu_eta; */
   /* vector<float>   *genMu_phi; */
   /* vector<float>   *genMu_ch; */
   /* vector<float>   *genMu_Id; */
   /* vector<float>   *genMu_motherId; */
   /* vector<float>   *genMu_gmotherId; */

   Int_t genpart_B_index;
   Int_t genpart_lep1FromB_index;
   Int_t genpart_lep2FromB_index;
   Int_t genpart_KFromB_index;
   Int_t genpart_KstFromB_index;
   Int_t genpart_KFromKst_index;
   Int_t genpart_PiFromKst_index;

   Int_t genpart_B_pdg;
   Int_t genpart_lep1FromB_pdg;
   Int_t genpart_lep2FromB_pdg;
   Int_t genpart_KFromB_pdg;
   Int_t genpart_KstFromB_pdg;
   Int_t genpart_KFromKst_pdg;
   Int_t genpart_PiFromKst_pdg;

   std::vector<float> *genpart_B_PtEtaPhiM;
   std::vector<float> *genpart_lep1_PtEtaPhiM;
   std::vector<float> *genpart_lep2_PtEtaPhiM;
   std::vector<float> *genpart_K_PtEtaPhiM;
   std::vector<float> *genpart_Kst_PtEtaPhiM;
   std::vector<float> *genpart_Pi_PtEtaPhiM;

   UInt_t          ntracks;
   vector<float>   *track_pt;
   vector<float>   *track_eta;
   vector<float>   *track_phi;
   vector<float>   *track_norm_chi2;
   vector<float>   *track_charge;
   vector<float>   *track_dxy;
   vector<float>   *track_dz;
   vector<float>   *track_edxy;
   vector<float>   *track_edz;
   vector<float>   *track_vx;
   vector<float>   *track_vy;
   vector<float>   *track_vz;
   vector<float>   *track_mva;
   vector<vector<float> > *Epair_PtEtaPhiM;
   vector<vector<float> > *Epair_XYZ;
   vector<float>   *Epair_cos;
   vector<float>   *Epair_chi_prob;
   vector<float>   *Epair_ObjIndex;
   vector<float>   *Epair_TrkIndex;
   vector<float>   *Epair_Lxy;
   vector<float>   *Epair_eLxy;
   vector<unsigned int> *Epair_ObjId;
   vector<vector<float> > *TTrack_PtEtaPhiM;
   vector<float>   *TTrack_chi_prob;
   vector<vector<float> > *TTrack_XYZ;
   vector<float>   *TTrack_ObjIndex;
   vector<float>   *TTrack_TrkIndex;
   vector<float>   *TTrack_KstarIndex;
   vector<float>   *TTrack_kid;
   vector<float>   *TTrack_piid;
   vector<float>   *TTrack_mll;
   vector<float>   *TTrack_mKst;
   vector<float>   *TTrack_cos;
   vector<float>   *TTrack_Lxy;
   vector<float>   *TTrack_eLxy;
   vector<unsigned int> *TTrack_ObjId;
   vector<vector<float> > *Kstpair_PtEtaPhiM;
   vector<float>   *Kstpair_chi_prob;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run_number;   //!
   TBranch        *b_ls;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_beam_x;   //!
   TBranch        *b_beam_y;   //!
   TBranch        *b_beam_z;   //!
   TBranch        *b_HLT_path1;   //!
   TBranch        *b_HLT_path2;   //!
   TBranch        *b_HLT_path3;   //!
   TBranch        *b_HLT_path4;   //!
   TBranch        *b_HLT_path5;   //!
   TBranch        *b_HLT_path6;   //!
   TBranch        *b_TrgObj1_PtEtaPhiCharge;   //!
   TBranch        *b_TrgObj2_PtEtaPhiCharge;   //!
   TBranch        *b_TrgObj3_PtEtaPhiCharge;   //!
   TBranch        *b_TrgObj4_PtEtaPhiCharge;   //!
   TBranch        *b_TrgObj5_PtEtaPhiCharge;   //!
   TBranch        *b_TrgObj6_PtEtaPhiCharge;   //!
   TBranch        *b_SelectedTrgObj_PtEtaPhiCharge;   //!
   TBranch        *b_SelectedMu_index;   //!
   TBranch        *b_SelectedMu_DR;   //!
   TBranch        *b_nmuon;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_dxy;   //!
   TBranch        *b_muon_dz;   //!
   TBranch        *b_muon_edxy;   //!
   TBranch        *b_muon_edz;   //!
   TBranch        *b_muon_d0;   //!
   TBranch        *b_muon_ed0;   //!
   TBranch        *b_muon_vx;   //!
   TBranch        *b_muon_vy;   //!
   TBranch        *b_muon_vz;   //!
   TBranch        *b_muon_iso;   //!
   TBranch        *b_muon_soft;   //!
   TBranch        *b_muon_loose;   //!
   TBranch        *b_muon_medium;   //!
   TBranch        *b_muon_tight;   //!
   TBranch        *b_muon_trkpt;   //!
   TBranch        *b_muon_trketa;   //!
   TBranch        *b_muon_trkphi;   //!
   TBranch        *b_nelectron;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_dxy;   //!
   TBranch        *b_el_dz;   //!
   TBranch        *b_el_edxy;   //!
   TBranch        *b_el_edz;   //!
   TBranch        *b_el_vx;   //!
   TBranch        *b_el_vy;   //!
   TBranch        *b_el_vz;   //!
   TBranch        *b_el_mvaB_wp;   //!
   TBranch        *b_el_mvaUnB_wp;   //!
   TBranch        *b_el_iso;   //!
   TBranch        *b_el_mvaId;   //!
   TBranch        *b_el_trkphi;   //!
   TBranch        *b_el_trkpt;   //!
   TBranch        *b_el_trketa;   //!

   TBranch        *b_ngsfTracks;
   TBranch        *b_gsfTrk_pt;
   TBranch        *b_gsfTrk_eta;
   TBranch        *b_gsfTrk_phi;
   TBranch        *b_gsfTrk_charge;
   TBranch        *b_gsfTrk_seedBDTunb;
   TBranch        *b_gsfTrk_seedBDTbiased;


   /* TBranch        *b_genpart_pt;   //! */
   /* TBranch        *b_genpart_phi;   //! */
   /* TBranch        *b_genpart_eta;   //! */
   /* TBranch        *b_genpart_pdgId;   //! */
   /* TBranch        *b_genpart_Bindex;   //! */
   /* TBranch        *b_genpart_Daughtindex;   //! */
   /* TBranch        *b_genpart_charge;   //! */
   /* TBranch        *b_genpart_mother_pdgId;   //! */
   /* TBranch        *b_genpart_mother_pt;   //! */
   /* TBranch        *b_genpart_mother_phi;   //! */
   /* TBranch        *b_genpart_mother_eta;   //! */
   /* TBranch        *b_genpart_mother_Bindex;   //! */
   /* TBranch        *b_genpart_mother_Daughtindex;   //! */
   /* TBranch        *b_genpart_grandmother_pdgId;   //! */
   /* TBranch        *b_genpart_grandmother_pt;   //! */
   /* TBranch        *b_genpart_grandmother_phi;   //! */
   /* TBranch        *b_genpart_grandmother_eta;   //! */
   /* TBranch        *b_genpart_grandmother_Bindex;   //! */
   /* TBranch        *b_genpart_grandmother_x;   //! */
   /* TBranch        *b_genpart_grandmother_y;   //! */
   /* TBranch        *b_genpart_grandmother_z;   //! */
   /* TBranch        *b_genMu_pt;   //! */
   /* TBranch        *b_genMu_eta;   //! */
   /* TBranch        *b_genMu_phi;   //! */
   /* TBranch        *b_genMu_ch;   //! */
   /* TBranch        *b_genMu_Id;   //! */
   /* TBranch        *b_genMu_motherId;   //! */
   /* TBranch        *b_genMu_gmotherId;   //! */

   TBranch        *b_genpart_B_index;
   TBranch        *b_genpart_lep1FromB_index;
   TBranch        *b_genpart_lep2FromB_index;
   TBranch        *b_genpart_KFromB_index;
   TBranch        *b_genpart_KstFromB_index;
   TBranch        *b_genpart_KFromKst_index;
   TBranch        *b_genpart_PiFromKst_index;
   TBranch        *b_genpart_B_pdg;
   TBranch        *b_genpart_lep1FromB_pdg;
   TBranch        *b_genpart_lep2FromB_pdg;
   TBranch        *b_genpart_KFromB_pdg;
   TBranch        *b_genpart_KstFromB_pdg;
   TBranch        *b_genpart_KFromKst_pdg;
   TBranch        *b_genpart_PiFromKst_pdg;
   TBranch        *b_genpart_B_PtEtaPhiM;
   TBranch        *b_genpart_lep1_PtEtaPhiM;
   TBranch        *b_genpart_lep2_PtEtaPhiM;
   TBranch        *b_genpart_K_PtEtaPhiM;
   TBranch        *b_genpart_Kst_PtEtaPhiM;
   TBranch        *b_genpart_Pi_PtEtaPhiM;


   TBranch        *b_ntracks;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_norm_chi2;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_dxy;   //!
   TBranch        *b_track_dz;   //!
   TBranch        *b_track_edxy;   //!
   TBranch        *b_track_edz;   //!
   TBranch        *b_track_vx;   //!
   TBranch        *b_track_vy;   //!
   TBranch        *b_track_vz;   //!
   TBranch        *b_track_mva;   //!
   TBranch        *b_Epair_PtEtaPhiM;   //!
   TBranch        *b_Epair_XYZ;   //!
   TBranch        *b_Epair_cos;   //!
   TBranch        *b_Epair_chi_prob;   //!
   TBranch        *b_Epair_ObjIndex;   //!
   TBranch        *b_Epair_TrkIndex;   //!
   TBranch        *b_Epair_Lxy;   //!
   TBranch        *b_Epair_eLxy;   //!
   TBranch        *b_Epair_ObjId;   //!
   TBranch        *b_TTrack_PtEtaPhiM;   //!
   TBranch        *b_TTrack_chi_prob;   //!
   TBranch        *b_TTrack_XYZ;   //!
   TBranch        *b_TTrack_ObjIndex;   //!
   TBranch        *b_TTrack_TrkIndex;   //!
   TBranch        *b_TTrack_KstarIndex;   //!
   TBranch        *b_TTrack_kid;   //!
   TBranch        *b_TTrack_piid;   //!
   TBranch        *b_TTrack_mll;   //!
   TBranch        *b_TTrack_mKst;   //!
   TBranch        *b_TTrack_cos;   //!
   TBranch        *b_TTrack_Lxy;   //!
   TBranch        *b_TTrack_eLxy;   //!
   TBranch        *b_TTrack_ObjId;   //!
   TBranch        *b_Kstpair_PtEtaPhiM;   //!
   TBranch        *b_Kstpair_chi_prob;

   skim_class(TTree *tree=0);
   virtual ~skim_class();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
};

skim_class::skim_class(TTree *tree) : fChain(0) 
{
   Init(tree);
}

skim_class::~skim_class()
{
   if (!fChain) return;
}

Int_t skim_class::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skim_class::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skim_class::Init(TTree *tree)
{

   // Set object pointer
   vertex_x = 0;
   vertex_y = 0;
   vertex_z = 0;
   TrgObj1_PtEtaPhiCharge = 0;
   TrgObj2_PtEtaPhiCharge = 0;
   TrgObj3_PtEtaPhiCharge = 0;
   TrgObj4_PtEtaPhiCharge = 0;
   TrgObj5_PtEtaPhiCharge = 0;
   TrgObj6_PtEtaPhiCharge = 0;
   SelectedTrgObj_PtEtaPhiCharge = 0;
   muon_pt = 0;
   muon_eta = 0;
   muon_phi = 0;
   muon_charge = 0;
   muon_dxy = 0;
   muon_dz = 0;
   muon_edxy = 0;
   muon_edz = 0;
   muon_d0 = 0;
   muon_ed0 = 0;
   muon_vx = 0;
   muon_vy = 0;
   muon_vz = 0;
   muon_iso = 0;
   muon_soft = 0;
   muon_loose = 0;
   muon_medium = 0;
   muon_tight = 0;
   muon_trkpt = 0;
   muon_trketa = 0;
   muon_trkphi = 0;
   el_pt = 0;
   el_eta = 0;
   el_phi = 0;
   el_charge = 0;
   el_dxy = 0;
   el_dz = 0;
   el_edxy = 0;
   el_edz = 0;
   el_vx = 0;
   el_vy = 0;
   el_vz = 0;
   el_mvaB_wp = 0;
   el_mvaUnB_wp = 0;
   el_iso = 0;
   el_mvaId = 0;
   el_trkphi = 0;
   el_trkpt = 0;
   el_trketa = 0;

   ngsfTracks = 0;
   gsfTrk_pt = 0;
   gsfTrk_eta = 0;
   gsfTrk_phi = 0;
   gsfTrk_charge = 0;
   gsfTrk_seedBDTunb = 0;
   gsfTrk_seedBDTbiased = 0;

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

   genpart_B_PtEtaPhiM = 0;
   genpart_lep1_PtEtaPhiM = 0;
   genpart_lep2_PtEtaPhiM = 0;
   genpart_K_PtEtaPhiM = 0;
   genpart_Kst_PtEtaPhiM = 0;
   genpart_Pi_PtEtaPhiM = 0;


   /* genpart_pt = 0; */
   /* genpart_phi = 0; */
   /* genpart_eta = 0; */
   /* genpart_pdgId = 0; */
   /* genpart_Bindex = 0; */
   /* genpart_Daughtindex = 0; */
   /* genpart_charge = 0; */
   /* genpart_mother_pdgId = 0; */
   /* genpart_mother_pt = 0; */
   /* genpart_mother_phi = 0; */
   /* genpart_mother_eta = 0; */
   /* genpart_mother_Bindex = 0; */
   /* genpart_mother_Daughtindex = 0; */
   /* genpart_grandmother_pdgId = 0; */
   /* genpart_grandmother_pt = 0; */
   /* genpart_grandmother_phi = 0; */
   /* genpart_grandmother_eta = 0; */
   /* genpart_grandmother_Bindex = 0; */
   /* genpart_grandmother_x = 0; */
   /* genpart_grandmother_y = 0; */
   /* genpart_grandmother_z = 0; */
   /* genMu_pt = 0; */
   /* genMu_eta = 0; */
   /* genMu_phi = 0; */
   /* genMu_ch = 0; */
   /* genMu_Id = 0; */
   /* genMu_motherId = 0; */
   /* genMu_gmotherId = 0; */

   track_pt = 0;
   track_eta = 0;
   track_phi = 0;
   track_norm_chi2 = 0;
   track_charge = 0;
   track_dxy = 0;
   track_dz = 0;
   track_edxy = 0;
   track_edz = 0;
   track_vx = 0;
   track_vy = 0;
   track_vz = 0;
   track_mva = 0;
   Epair_PtEtaPhiM = 0;
   Epair_XYZ = 0;
   Epair_cos = 0;
   Epair_chi_prob = 0;
   Epair_ObjIndex = 0;
   Epair_TrkIndex = 0;
   Epair_Lxy = 0;
   Epair_eLxy = 0;
   Epair_ObjId = 0;
   TTrack_PtEtaPhiM = 0;
   TTrack_chi_prob = 0;
   TTrack_XYZ = 0;
   TTrack_ObjIndex = 0;
   TTrack_TrkIndex = 0;
   TTrack_KstarIndex = 0;
   TTrack_kid = 0;
   TTrack_piid = 0;
   TTrack_mll = 0;
   TTrack_mKst = 0;
   TTrack_cos = 0;
   TTrack_Lxy = 0;
   TTrack_eLxy = 0;
   TTrack_ObjId = 0;
   Kstpair_PtEtaPhiM = 0;
   Kstpair_chi_prob = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run_number", &run_number, &b_run_number);
   fChain->SetBranchAddress("ls", &ls, &b_ls);
   fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
   fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
   fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
   fChain->SetBranchAddress("beam_x", &beam_x, &b_beam_x);
   fChain->SetBranchAddress("beam_y", &beam_y, &b_beam_y);
   fChain->SetBranchAddress("beam_z", &beam_z, &b_beam_z);
   fChain->SetBranchAddress("HLT_path1", &HLT_path1, &b_HLT_path1);
   fChain->SetBranchAddress("HLT_path2", &HLT_path2, &b_HLT_path2);
   fChain->SetBranchAddress("HLT_path3", &HLT_path3, &b_HLT_path3);
   fChain->SetBranchAddress("HLT_path4", &HLT_path4, &b_HLT_path4);
   fChain->SetBranchAddress("HLT_path5", &HLT_path5, &b_HLT_path5);
   fChain->SetBranchAddress("HLT_path6", &HLT_path6, &b_HLT_path6);
   fChain->SetBranchAddress("TrgObj1_PtEtaPhiCharge", &TrgObj1_PtEtaPhiCharge, &b_TrgObj1_PtEtaPhiCharge);
   fChain->SetBranchAddress("TrgObj2_PtEtaPhiCharge", &TrgObj2_PtEtaPhiCharge, &b_TrgObj2_PtEtaPhiCharge);
   fChain->SetBranchAddress("TrgObj3_PtEtaPhiCharge", &TrgObj3_PtEtaPhiCharge, &b_TrgObj3_PtEtaPhiCharge);
   fChain->SetBranchAddress("TrgObj4_PtEtaPhiCharge", &TrgObj4_PtEtaPhiCharge, &b_TrgObj4_PtEtaPhiCharge);
   fChain->SetBranchAddress("TrgObj5_PtEtaPhiCharge", &TrgObj5_PtEtaPhiCharge, &b_TrgObj5_PtEtaPhiCharge);
   fChain->SetBranchAddress("TrgObj6_PtEtaPhiCharge", &TrgObj6_PtEtaPhiCharge, &b_TrgObj6_PtEtaPhiCharge);
   fChain->SetBranchAddress("SelectedTrgObj_PtEtaPhiCharge", &SelectedTrgObj_PtEtaPhiCharge, &b_SelectedTrgObj_PtEtaPhiCharge);
   fChain->SetBranchAddress("SelectedMu_index", &SelectedMu_index, &b_SelectedMu_index);
   fChain->SetBranchAddress("SelectedMu_DR", &SelectedMu_DR, &b_SelectedMu_DR);
   fChain->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
   fChain->SetBranchAddress("muon_pt", &muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_eta", &muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_dxy", &muon_dxy, &b_muon_dxy);
   fChain->SetBranchAddress("muon_dz", &muon_dz, &b_muon_dz);
   fChain->SetBranchAddress("muon_edxy", &muon_edxy, &b_muon_edxy);
   fChain->SetBranchAddress("muon_edz", &muon_edz, &b_muon_edz);
   fChain->SetBranchAddress("muon_d0", &muon_d0, &b_muon_d0);
   fChain->SetBranchAddress("muon_ed0", &muon_ed0, &b_muon_ed0);
   fChain->SetBranchAddress("muon_vx", &muon_vx, &b_muon_vx);
   fChain->SetBranchAddress("muon_vy", &muon_vy, &b_muon_vy);
   fChain->SetBranchAddress("muon_vz", &muon_vz, &b_muon_vz);
   fChain->SetBranchAddress("muon_iso", &muon_iso, &b_muon_iso);
   fChain->SetBranchAddress("muon_soft", &muon_soft, &b_muon_soft);
   fChain->SetBranchAddress("muon_loose", &muon_loose, &b_muon_loose);
   fChain->SetBranchAddress("muon_medium", &muon_medium, &b_muon_medium);
   fChain->SetBranchAddress("muon_tight", &muon_tight, &b_muon_tight);
   fChain->SetBranchAddress("muon_trkpt", &muon_trkpt, &b_muon_trkpt);
   fChain->SetBranchAddress("muon_trketa", &muon_trketa, &b_muon_trketa);
   fChain->SetBranchAddress("muon_trkphi", &muon_trkphi, &b_muon_trkphi);
   fChain->SetBranchAddress("nelectron", &nelectron, &b_nelectron);
   fChain->SetBranchAddress("el_pt", &el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_dxy", &el_dxy, &b_el_dxy);
   fChain->SetBranchAddress("el_dz", &el_dz, &b_el_dz);
   fChain->SetBranchAddress("el_edxy", &el_edxy, &b_el_edxy);
   fChain->SetBranchAddress("el_edz", &el_edz, &b_el_edz);
   fChain->SetBranchAddress("el_vx", &el_vx, &b_el_vx);
   fChain->SetBranchAddress("el_vy", &el_vy, &b_el_vy);
   fChain->SetBranchAddress("el_vz", &el_vz, &b_el_vz);
   /* fChain->SetBranchAddress("el_mvaB_wp", &el_mvaB_wp, &b_el_mvaB_wp); */
   /* fChain->SetBranchAddress("el_mvaUnB_wp", &el_mvaUnB_wp, &b_el_mvaUnB_wp); */
   /* fChain->SetBranchAddress("el_iso", &el_iso, &b_el_iso); */
   /* fChain->SetBranchAddress("el_mvaId", &el_mvaId, &b_el_mvaId); */
   fChain->SetBranchAddress("el_trkphi", &el_trkphi, &b_el_trkphi);
   fChain->SetBranchAddress("el_trkpt", &el_trkpt, &b_el_trkpt);
   fChain->SetBranchAddress("el_trketa", &el_trketa, &b_el_trketa);

   fChain->SetBranchAddress("ngsfTracks", &ngsfTracks, &b_ngsfTracks);
   fChain->SetBranchAddress("gsfTrk_pt", &gsfTrk_pt, &b_gsfTrk_pt);
   fChain->SetBranchAddress("gsfTrk_eta", &gsfTrk_eta, &b_gsfTrk_eta);
   fChain->SetBranchAddress("gsfTrk_phi", &gsfTrk_phi, &b_gsfTrk_phi);
   fChain->SetBranchAddress("gsfTrk_charge", &gsfTrk_charge, &b_gsfTrk_charge);
   fChain->SetBranchAddress("gsfTrk_seedBDTunb", &gsfTrk_seedBDTunb, &b_gsfTrk_seedBDTunb);
   fChain->SetBranchAddress("gsfTrk_seedBDTbiased", &gsfTrk_seedBDTbiased, &b_gsfTrk_seedBDTbiased);

   /* fChain->SetBranchAddress("genpart_pt", &genpart_pt, &b_genpart_pt); */
   /* fChain->SetBranchAddress("genpart_phi", &genpart_phi, &b_genpart_phi); */
   /* fChain->SetBranchAddress("genpart_eta", &genpart_eta, &b_genpart_eta); */
   /* fChain->SetBranchAddress("genpart_pdgId", &genpart_pdgId, &b_genpart_pdgId); */
   /* fChain->SetBranchAddress("genpart_Bindex", &genpart_Bindex, &b_genpart_Bindex); */
   /* fChain->SetBranchAddress("genpart_Daughtindex", &genpart_Daughtindex, &b_genpart_Daughtindex); */
   /* fChain->SetBranchAddress("genpart_charge", &genpart_charge, &b_genpart_charge); */
   /* fChain->SetBranchAddress("genpart_mother_pdgId", &genpart_mother_pdgId, &b_genpart_mother_pdgId); */
   /* fChain->SetBranchAddress("genpart_mother_pt", &genpart_mother_pt, &b_genpart_mother_pt); */
   /* fChain->SetBranchAddress("genpart_mother_phi", &genpart_mother_phi, &b_genpart_mother_phi); */
   /* fChain->SetBranchAddress("genpart_mother_eta", &genpart_mother_eta, &b_genpart_mother_eta); */
   /* fChain->SetBranchAddress("genpart_mother_Bindex", &genpart_mother_Bindex, &b_genpart_mother_Bindex); */
   /* fChain->SetBranchAddress("genpart_mother_Daughtindex", &genpart_mother_Daughtindex, &b_genpart_mother_Daughtindex); */
   /* fChain->SetBranchAddress("genpart_grandmother_pdgId", &genpart_grandmother_pdgId, &b_genpart_grandmother_pdgId); */
   /* fChain->SetBranchAddress("genpart_grandmother_pt", &genpart_grandmother_pt, &b_genpart_grandmother_pt); */
   /* fChain->SetBranchAddress("genpart_grandmother_phi", &genpart_grandmother_phi, &b_genpart_grandmother_phi); */
   /* fChain->SetBranchAddress("genpart_grandmother_eta", &genpart_grandmother_eta, &b_genpart_grandmother_eta); */
   /* fChain->SetBranchAddress("genpart_grandmother_Bindex", &genpart_grandmother_Bindex, &b_genpart_grandmother_Bindex); */
   /* fChain->SetBranchAddress("genpart_grandmother_x", &genpart_grandmother_x, &b_genpart_grandmother_x); */
   /* fChain->SetBranchAddress("genpart_grandmother_y", &genpart_grandmother_y, &b_genpart_grandmother_y); */
   /* fChain->SetBranchAddress("genpart_grandmother_z", &genpart_grandmother_z, &b_genpart_grandmother_z); */
   /* fChain->SetBranchAddress("genMu_pt", &genMu_pt, &b_genMu_pt); */
   /* fChain->SetBranchAddress("genMu_eta", &genMu_eta, &b_genMu_eta); */
   /* fChain->SetBranchAddress("genMu_phi", &genMu_phi, &b_genMu_phi); */
   /* fChain->SetBranchAddress("genMu_ch", &genMu_ch, &b_genMu_ch); */
   /* fChain->SetBranchAddress("genMu_Id", &genMu_Id, &b_genMu_Id); */
   /* fChain->SetBranchAddress("genMu_motherId", &genMu_motherId, &b_genMu_motherId); */
   /* fChain->SetBranchAddress("genMu_gmotherId", &genMu_gmotherId, &b_genMu_gmotherId); */

   fChain->SetBranchAddress("genpart_B_index", &genpart_B_index, &b_genpart_B_index);
   fChain->SetBranchAddress("genpart_lep1FromB_index", &genpart_lep1FromB_index, &b_genpart_lep1FromB_index);
   fChain->SetBranchAddress("genpart_lep2FromB_index", &genpart_lep2FromB_index, &b_genpart_lep2FromB_index);
   fChain->SetBranchAddress("genpart_KFromB_index", &genpart_KFromB_index, &b_genpart_KFromB_index);
   fChain->SetBranchAddress("genpart_KstFromB_index", &genpart_KstFromB_index, &b_genpart_KstFromB_index);
   fChain->SetBranchAddress("genpart_KFromKst_index", &genpart_KFromKst_index, &b_genpart_KFromKst_index);
   fChain->SetBranchAddress("genpart_PiFromKst_index", &genpart_PiFromKst_index, &b_genpart_PiFromKst_index);

   fChain->SetBranchAddress("genpart_B_pdg", &genpart_B_pdg, &b_genpart_B_pdg);
   fChain->SetBranchAddress("genpart_lep1FromB_pdg", &genpart_lep1FromB_pdg, &b_genpart_lep1FromB_pdg);
   fChain->SetBranchAddress("genpart_lep2FromB_pdg", &genpart_lep2FromB_pdg, &b_genpart_lep2FromB_pdg);
   fChain->SetBranchAddress("genpart_KFromB_pdg", &genpart_KFromB_pdg, &b_genpart_KFromB_pdg);
   fChain->SetBranchAddress("genpart_KstFromB_pdg", &genpart_KstFromB_pdg, &b_genpart_KstFromB_pdg);
   fChain->SetBranchAddress("genpart_KFromKst_pdg", &genpart_KFromKst_pdg, &b_genpart_KFromKst_pdg);
   fChain->SetBranchAddress("genpart_PiFromKst_pdg", &genpart_PiFromKst_pdg, &b_genpart_PiFromKst_pdg);

   fChain->SetBranchAddress("genpart_B_PtEtaPhiM", &genpart_B_PtEtaPhiM, &b_genpart_B_PtEtaPhiM);
   fChain->SetBranchAddress("genpart_lep1_PtEtaPhiM", &genpart_lep1_PtEtaPhiM, &b_genpart_lep1_PtEtaPhiM);
   fChain->SetBranchAddress("genpart_lep2_PtEtaPhiM", &genpart_lep2_PtEtaPhiM, &b_genpart_lep2_PtEtaPhiM);
   fChain->SetBranchAddress("genpart_K_PtEtaPhiM", &genpart_K_PtEtaPhiM, &b_genpart_K_PtEtaPhiM);
   fChain->SetBranchAddress("genpart_Kst_PtEtaPhiM", &genpart_Kst_PtEtaPhiM, &b_genpart_Kst_PtEtaPhiM);
   fChain->SetBranchAddress("genpart_Pi_PtEtaPhiM", &genpart_Pi_PtEtaPhiM, &b_genpart_Pi_PtEtaPhiM);

   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("track_pt", &track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_eta", &track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", &track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_norm_chi2", &track_norm_chi2, &b_track_norm_chi2);
   fChain->SetBranchAddress("track_charge", &track_charge, &b_track_charge);
   fChain->SetBranchAddress("track_dxy", &track_dxy, &b_track_dxy);
   fChain->SetBranchAddress("track_dz", &track_dz, &b_track_dz);
   fChain->SetBranchAddress("track_edxy", &track_edxy, &b_track_edxy);
   fChain->SetBranchAddress("track_edz", &track_edz, &b_track_edz);
   fChain->SetBranchAddress("track_vx", &track_vx, &b_track_vx);
   fChain->SetBranchAddress("track_vy", &track_vy, &b_track_vy);
   fChain->SetBranchAddress("track_vz", &track_vz, &b_track_vz);
   fChain->SetBranchAddress("track_mva", &track_mva, &b_track_mva);
   fChain->SetBranchAddress("Epair_PtEtaPhiM", &Epair_PtEtaPhiM, &b_Epair_PtEtaPhiM);
   fChain->SetBranchAddress("Epair_XYZ", &Epair_XYZ, &b_Epair_XYZ);
   fChain->SetBranchAddress("Epair_cos", &Epair_cos, &b_Epair_cos);
   fChain->SetBranchAddress("Epair_chi_prob", &Epair_chi_prob, &b_Epair_chi_prob);
   fChain->SetBranchAddress("Epair_ObjIndex", &Epair_ObjIndex, &b_Epair_ObjIndex);
   fChain->SetBranchAddress("Epair_TrkIndex", &Epair_TrkIndex, &b_Epair_TrkIndex);
   fChain->SetBranchAddress("Epair_Lxy", &Epair_Lxy, &b_Epair_Lxy);
   fChain->SetBranchAddress("Epair_eLxy", &Epair_eLxy, &b_Epair_eLxy);
   fChain->SetBranchAddress("Epair_ObjId", &Epair_ObjId, &b_Epair_ObjId);
   fChain->SetBranchAddress("TTrack_PtEtaPhiM", &TTrack_PtEtaPhiM, &b_TTrack_PtEtaPhiM);
   fChain->SetBranchAddress("TTrack_chi_prob", &TTrack_chi_prob, &b_TTrack_chi_prob);
   fChain->SetBranchAddress("TTrack_XYZ", &TTrack_XYZ, &b_TTrack_XYZ);
   fChain->SetBranchAddress("TTrack_ObjIndex", &TTrack_ObjIndex, &b_TTrack_ObjIndex);
   fChain->SetBranchAddress("TTrack_TrkIndex", &TTrack_TrkIndex, &b_TTrack_TrkIndex);
   fChain->SetBranchAddress("TTrack_KstarIndex", &TTrack_KstarIndex, &b_TTrack_KstarIndex);
   fChain->SetBranchAddress("TTrack_kid", &TTrack_kid, &b_TTrack_kid);
   fChain->SetBranchAddress("TTrack_piid", &TTrack_piid, &b_TTrack_piid);
   fChain->SetBranchAddress("TTrack_mll", &TTrack_mll, &b_TTrack_mll);
   fChain->SetBranchAddress("TTrack_mKst", &TTrack_mKst, &b_TTrack_mKst);
   fChain->SetBranchAddress("TTrack_cos", &TTrack_cos, &b_TTrack_cos);
   fChain->SetBranchAddress("TTrack_Lxy", &TTrack_Lxy, &b_TTrack_Lxy);
   fChain->SetBranchAddress("TTrack_eLxy", &TTrack_eLxy, &b_TTrack_eLxy);
   fChain->SetBranchAddress("TTrack_ObjId", &TTrack_ObjId, &b_TTrack_ObjId);
   fChain->SetBranchAddress("Kstpair_PtEtaPhiM", &Kstpair_PtEtaPhiM, &b_Kstpair_PtEtaPhiM);
   fChain->SetBranchAddress("Kstpair_chi_prob", &Kstpair_chi_prob, &b_Kstpair_chi_prob);
   Notify();
}

Bool_t skim_class::Notify()
{

   return kTRUE;
}

#endif // #ifdef skim_class_cxx
