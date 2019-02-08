#include "skim_class.h"
#include "TChain.h"
#include "TH1.h"


float findmin(float a,float b ,float c){
  if (a<b && a<c) return a;
  else if (b<a && b<c) return b;
  else return c;
}

float findmax(float a,float b ,float c){
  if (a>b && a>c) return a;
  else if (b>a && b>c) return b;
  else return c;
}

float Dphi(float phi1, float phi2){
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float DR(float eta1, float phi1, float eta2, float phi2){
  return TMath::Sqrt((eta1-eta2)*(eta1-eta2)+Dphi(phi1,phi2)*Dphi(phi1,phi2));
}



int studyEfficiencyRate(){
  TChain *ccmc = new TChain("demo/mytree");
  TChain *ccdata = new TChain("demo/mytree");

  ccdata->Add("../test/outputResHPt_DATA.root");
  //ccmc->Add("../test/outputResHPt_MC10k_newHLT.root");
  ccmc->Add("../test/outputResHPt_MC1k.root");
  //ccmc->Add("../test/outputResHPt.root");
  skim_class mc;
  mc.Init(ccmc);

  skim_class data;
  data.Init(ccdata);

  //match gen
  TH1F* hBptGen = new TH1F("hBptGen"," ",60,0,60);
  TH1F* hE1ptGen = new TH1F("hE1ptGen"," ",60,0,60);
  TH1F* hE2ptGen = new TH1F("hE2ptGen"," ",60,0,60);
  TH1F* hKptGen = new TH1F("hKptGen"," ",60,0,60);

  TH1F* hDReltrk = new TH1F("hDReltrk"," ",100,0,0.3);
  TH1F* hDRK = new TH1F("hDRK"," ",100,0,0.3);
  TH1F* hDRrecoE = new TH1F("hDRrecoE"," ",100,0,0.3);

  TH1F* hBRecoPt = new TH1F("hBRecoPt"," ",60,0,60);
  TH1F* hE1RecoPt = new TH1F("hE1RecoPt"," ",60,0,60);
  TH1F* hE2RecoPt = new TH1F("hE2RecoPt"," ",60,0,60);
  TH1F* hKRecoPt = new TH1F("hKRecoPt"," ",60,0,60);


  float Cos_cuts[10] = {-1.1, 0, 0.5,0.7, 0.8, 0.9, 0.95,0.99,0.995,0.999};
  float PtB_cuts[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  float Lxy_cuts[10] = {0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6 };
  float M_cuts[10] = {3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.1, 5.27 };
  float Trk_cuts[10] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8 };
  float SigEffCos[10] = {0},BkgEffCos[10]={0};
  float SigEffPtB[10] = {0},BkgEffPtB[10]={0}; 
  float SigEffLxy[10] = {0},BkgEffLxy[10]={0};
  float SigEffM[10] = {0},BkgEffM[10]={0};
  float SigEffTrk[10] = {0},BkgEffTrk[10]={0};

  //signal
  float nEventsMC = ccmc->GetEntries();
  std::cout << " mc events = " << nEventsMC << std::endl;
  

  float ngen = 0, ngenAcceptance = 0;
  float nPt2 = 0; // gen lepton pT > 2
  float reco_obj = 0;
  float reco_vtx = 0,reco_dz = 0,reco_mll = 0,reco_llprob = 0,reco_kpt = 0, reco_Mtri = 0;
  float reco_prob = 0, reco_cos = 0,reco_ptB = 0, reco_Lxy = 0,  nObj=0, nTrk=0, nK=0, nTrk2=0;

  for (int ev=0; ev<nEventsMC; ev++){
    mc.GetEntry(ev);

    //take reco candidate closest to gen
    TLorentzVector vel1, vel2, vK;

    if(mc.genpart_B_index == -1) continue;

    //Kee
    vel1.SetPtEtaPhiM(mc.genpart_lep1_PtEtaPhiM->at(0), mc.genpart_lep1_PtEtaPhiM->at(1), mc.genpart_lep1_PtEtaPhiM->at(2), mc.genpart_lep1_PtEtaPhiM->at(3));
    vel2.SetPtEtaPhiM(mc.genpart_lep2_PtEtaPhiM->at(0), mc.genpart_lep2_PtEtaPhiM->at(1), mc.genpart_lep2_PtEtaPhiM->at(2), mc.genpart_lep2_PtEtaPhiM->at(3));
    vK.SetPtEtaPhiM(mc.genpart_K_PtEtaPhiM->at(0), mc.genpart_K_PtEtaPhiM->at(1), mc.genpart_K_PtEtaPhiM->at(2), mc.genpart_K_PtEtaPhiM->at(3));
    bool BKee = true;
    float Kcharge = (mc.genpart_KFromB_pdg > 0) ? 1 : -1;
    float Ele1charge = (mc.genpart_lep1FromB_pdg > 0) ? 1 : -1;
    float Ele2charge = (mc.genpart_lep2FromB_pdg > 0) ? 1 : -1;
    if (!BKee) continue;
    

    if(std::abs(vel1.Eta()) > 2.4 || std::abs(vel2.Eta()) > 2.4 || std::abs(vK.Eta()) > 2.4 ||
       vel1.Pt() < 0.5 || vel2.Eta() < 0.5 || vK.Eta() < 0.5) continue;
    ++ngenAcceptance;

    if (fabs(vel1.Eta())>2.4 || fabs(vel2.Eta())>2.4 || fabs(vK.Eta())>2.4) continue;
    ngen++;

     hBptGen->Fill((vel1+vel2+vK).Pt()); 
     hE1ptGen->Fill(vel1.Pt());
     hE2ptGen->Fill(vel2.Pt());
     hKptGen->Fill(vK.Pt()); 
     //     std::cout << " >>> (vel1+vel2+vK).M() = " << (vel1+vel2+vK).M() << std::endl;
     //     std::cout << " mc triplets = " << mc.TTrack_cos->size() << std::endl;

     float minDRK = 1000;
     float minDRrecoE = 1000;
     float minDRtrkE = 1000;
     
     int recoE = -1;
     int trkEl = -1;
     int trkK = -1;
     int Bindex = -1;
     
     bool found = false; int selVetx = -1;
     for (unsigned int iB=0; iB<mc.TTrack_cos->size(); iB++)  {
       if (mc.TTrack_ObjId->at(iB) != 11) continue;
       //std::cout << " >> mc.TTrack_ObjId->at(iB) = " << mc.TTrack_ObjId->at(iB) << std::endl;
       
       int recoEd = mc.TTrack_ObjIndex->at(iB);
       int trkEld = mc.TTrack_TrkIndex->at(iB);
       int trkKd = mc.TTrack_kid->at(iB);
       
       float dR_ele1 = DR(vel1.Eta(),vel1.Phi(),mc.el_eta->at(recoEd),mc.el_phi->at(recoEd));
       float dR_ele2 = DR(vel2.Eta(),vel2.Phi(),mc.track_eta->at(trkEld),mc.track_phi->at(trkEld));
       float dR_trk = DR(vK.Eta(),vK.Phi(),mc.track_eta->at(trkKd),mc.track_phi->at(trkKd));
       if ((dR_ele1 + dR_ele2 + dR_trk) < (minDRrecoE+minDRrecoE+minDRtrkE) && //mc.el_charge->at(recoEd) == Ele1charge && 
	   //mc.track_charge->at(trkEld) == Ele2charge &&
	   mc.track_charge->at(trkKd) == Kcharge &&
	   mc.el_charge->at(recoEd) * mc.track_charge->at(trkEld) < 0){                                                                              
	 minDRrecoE = dR_ele1;
	 recoE = recoEd;
	 //recoEle2charge = mc.track_charge->at(itrk);
	 
	 minDRtrkE = dR_ele2;
	 trkEl = trkEld;
	 //recoEle2charge = mc.track_charge->at(itrk);
	 
	 minDRK = dR_trk;
	 trkK = trkKd;
	 //recoEle2charge = mc.track_charge->at(itrk);
	 Bindex = iB;
	 //std::cout << " >>> found triplet " << std::endl;
       }
     }//loop over triplets
     
     if(trkK != -1 && trkEl != -1 && recoE != -1){
       found = true; 
       selVetx = Bindex;
     }
     if(!found) continue;
     
     hDRrecoE->Fill(recoE);
     hDReltrk->Fill(minDRtrkE);
     hDRK->Fill(minDRK);

    ++reco_obj;
    TLorentzVector vrel, vreltrk, vrK;
    vrel.SetPtEtaPhiM(mc.el_pt->at(recoE), mc.el_eta->at(recoE), mc.el_phi->at(recoE), 0.0005);
    vreltrk.SetPtEtaPhiM(mc.track_pt->at(trkEl), mc.track_eta->at(trkEl),mc.track_phi->at(trkEl),0.0005);   
    vrK.SetPtEtaPhiM(mc.track_pt->at(trkK),mc.track_eta->at(trkK),mc.track_phi->at(trkK),0.493);

    hBRecoPt->Fill((vrel+vreltrk+vrK).Pt()); 
    hE1RecoPt->Fill(vrel.Pt()); 
    hE2RecoPt->Fill(vreltrk.Pt()); 
    hKRecoPt->Fill(vrK.Pt()); 

    bool pairProb=false;
    for (unsigned int iB=0; iB<mc.Epair_cos->size(); iB++)  {
      if (mc.Epair_ObjId->at(iB) != 11) continue;
      if (mc.Epair_ObjIndex->at(iB) == recoE && mc.Epair_TrkIndex->at(iB) == trkEl){
	if(mc.Epair_chi_prob->at(iB) > 1.e-34) 
	  pairProb = true;
      }
    }

    if(mc.SelectedMu_index == -1) continue;

    bool setCuts1 = false;
    if (found){ reco_vtx++; }


    if (found && fabs(mc.muon_vz->at(mc.SelectedMu_index) - mc.el_vz->at(recoE) ) < 0.3) { 
      reco_dz++;
      //      std::cout << " >>> (vrel+vreltrk).M() = " << (vrel+vreltrk).M() << std::endl;
    }
    if (found && fabs(mc.muon_vz->at(mc.SelectedMu_index) - mc.el_vz->at(recoE) ) <0.3 && (vrel+vreltrk).M() < 5) reco_mll++;
    if (found && pairProb && fabs(mc.muon_vz->at(mc.SelectedMu_index)-mc.el_vz->at(recoE))<0.3 && (vrel+vreltrk).M()<5) reco_llprob++;
    if (found && pairProb && fabs(mc.muon_vz->at(mc.SelectedMu_index)-mc.el_vz->at(recoE))<0.3 && (vrel+vreltrk).M()<5 && vrK.Pt()>0.4){
      reco_kpt++;
      //      std::cout<< " >>> (vrel+vreltrk+vrK).M() = " << (vrel+vreltrk+vrK).M() << std::endl;
    }
    if (found && pairProb && fabs(mc.muon_vz->at(mc.SelectedMu_index)-mc.el_vz->at(recoE))<0.3 && (vrel+vreltrk).M()<5 && vrK.Pt()>0.4 && 
	3.5 < (vrel+vreltrk+vrK).M() && (vrel+vreltrk+vrK).M()<6 ) {
      reco_Mtri++;
      setCuts1 = true;
    }


    if (setCuts1 && mc.TTrack_chi_prob->at(selVetx) > 0.001){ 
      reco_prob++;

      for (int i=0; i<10; i++){
       if (mc.TTrack_cos->at(selVetx) > Cos_cuts[i]) SigEffCos[i]++;  
      }
    }

    if (setCuts1 && mc.TTrack_chi_prob->at(selVetx)>0.001 && mc.TTrack_cos->at(selVetx)>0.9 ){ 
      reco_cos++;
      for (int i=0; i<10; i++){
	if (mc.TTrack_PtEtaPhiM->at(selVetx).at(0) > PtB_cuts[i]) SigEffPtB[i]++;  
      }
    }

    if (setCuts1 && mc.TTrack_chi_prob->at(selVetx)>0.001 && mc.TTrack_cos->at(selVetx)>0.9 && mc.TTrack_PtEtaPhiM->at(selVetx).at(0) > 5 ){ 
      reco_ptB++;
      for (int i=0; i<10; i++){
	if (mc.TTrack_Lxy->at(selVetx)/TMath::Sqrt(mc.TTrack_eLxy->at(selVetx))>Lxy_cuts[i]) SigEffLxy[i]++;  
      }
    }

    bool setCuts2 = false;
    if (setCuts1 && mc.TTrack_chi_prob->at(selVetx)>0.001 && mc.TTrack_cos->at(selVetx)>0.9 && mc.TTrack_PtEtaPhiM->at(selVetx).at(0)>5 && 
	mc.TTrack_Lxy->at(selVetx)/TMath::Sqrt(mc.TTrack_eLxy->at(selVetx))>4 ) {
      reco_Lxy++;
      setCuts2 = true;
    }


    if (setCuts2){
      for (int i=0; i<10; i++){
	if ( mc.TTrack_PtEtaPhiM->at(selVetx).at(3)>M_cuts[i]) SigEffM[i]++;
      }
    }
    if (setCuts2 && mc.TTrack_PtEtaPhiM->at(selVetx).at(3)>4.1){
      for (int i=0; i<10; i++){
	if ( vreltrk.Pt()>Trk_cuts[i] && vrK.Pt()>Trk_cuts[i]) SigEffTrk[i]++;
      }
    }
     
  }
  //  hSigQ->SetBinContent(1,reco_nosoft);  hSigQ->SetBinContent(2,reco_soft);
 


  std::cout << " \n\n Ngen = " << ngen 
	    << " \n Ngen in acceptance = " << ngenAcceptance 
	    << " \n Nreco triplets genMatched = " << reco_obj 
	    << " \n NrecoTriplets and tagMu = " << reco_vtx 
	    << " \n Ele1-Trg dz = " << reco_dz 
	    << " \n mll < 5 = " << reco_mll 
    //<< " \n ll_vtxCL > 1.e-34 = " << reco_llprob
	    << " \n kpt > 0.4 = " << reco_kpt 
	    << " \n Mtrip in [3.5,6] = " << reco_Mtri 
	    << " \n B_vtxCL > 1.e-3 = " << reco_prob
	    << " \n Bcos > 0.9 = " << reco_cos 
	    << " \n ptB >5  = " << reco_ptB 
	    << " \n Lxy > 4 = " << reco_Lxy << std::endl;

  //  std::cout<<"Pt>2 "<<nPt2/ngen<<" nobj "<<nObj/ngen<<"  ntrk "<<nTrk/nTrk2<<" nK "<<nK/ngen<<endl;


  //bakg
  float nEventsData = ccdata->GetEntries();
  std::cout << "data " << nEventsData << std::endl;

  float nbkg=0,nDz=0,nMll=0,nKpt=0,nMtri=0,nProb=0,nCos=0,nSoft=0,nNoSoft=0,nPtB=0,nLxy=0,nSoftEl=0,nTotal=0;
  float Ntotal=0; int NMutag = 0;
  float Ndz=0; float Nmll=0; float Nkpt=0; float Nmtri=0;
  float Nprob=0; float Nvtx=0; //float Ncos=0; float NptB=0; float Nlxy=0;
  int entries= 1000;//ccdata->GetEntries();
  for (int ev=0; ev<nEventsData; ev++){
    data.GetEntry(ev);
    if (ev%100==0) cout<<ev<<endl;
    //cuts
    int tntrk=0,tnel=0;


    float maxcos=-100,maxPtB=-100,maxLxy=-100,Mtri=100,minTrk=-100;

    bool found=false; 
    bool CutMutag = false;
    bool CutDz=false; bool CutMll=false; bool CutKpt=false;
    bool CutMtri=false; bool CutProb=false; bool CutCos=false; bool CutPtB=false;
    bool CutLxy=false; bool CutSoft=false; bool CutMtri2=false;
    bool Triplet=false;


    if( data.TTrack_cos->size() > 0) {
      Triplet = true;
      Ntotal = data.TTrack_cos->size();
    }

    if (Triplet)nTotal++;
    for (int ivtx=0; ivtx<data.TTrack_cos->size(); ivtx++){
      if (data.TTrack_ObjId->at(ivtx) != 11){ 
	//  std::cout << " data.TTrack_ObjId->at(ivtx) = " << data.TTrack_ObjId->at(ivtx) << std::endl; 
	continue;
      }
      found = true;
      Nvtx++;
      int eid = data.TTrack_ObjIndex->at(ivtx);

      if (data.SelectedMu_index < 0) continue;
      CutMutag = true;

      float ZmuVtx = data.muon_vz->at(data.SelectedMu_index);

      if (fabs(data.el_vz->at(eid) - ZmuVtx)>0.3) continue;
      CutDz=true;
      Ndz++;
      if (data.TTrack_mll->at(ivtx)>5) continue;
      CutMll=true;
      Nmll++;
      int kid = data.TTrack_kid->at(ivtx);
      if (data.track_pt->at(kid) < 0.4) continue;
      CutKpt=true;
      Nkpt++;
      if (data.TTrack_PtEtaPhiM->at(ivtx).at(3)<3.5 || data.TTrack_PtEtaPhiM->at(ivtx).at(3)>6 ) continue;
      CutMtri=true;
      Nmtri++;
      if (data.TTrack_chi_prob->at(ivtx)<0.001)  continue;
      CutProb=true; 
      Nprob++;
      if ( maxcos<data.TTrack_cos->at(ivtx)) maxcos=data.TTrack_cos->at(ivtx);
      if (data.TTrack_cos->at(ivtx)<0.9)  continue;
      CutCos=true;
      if ( maxPtB<data.TTrack_PtEtaPhiM->at(ivtx).at(0)) maxPtB=data.TTrack_PtEtaPhiM->at(ivtx).at(0);
      if (data.TTrack_PtEtaPhiM->at(ivtx).at(0)<5)  continue;
      CutPtB=true;
      if (maxLxy<data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx)))
	maxLxy=data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx));
      if (data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx))<4) 
	continue;
      CutLxy=true;
      
      if (fabs(Mtri-5.27)>fabs(5.27-data.TTrack_PtEtaPhiM->at(ivtx).at(3))) Mtri=data.TTrack_PtEtaPhiM->at(ivtx).at(3); 
      if (data.TTrack_PtEtaPhiM->at(ivtx).at(3)<4.1) continue;
      CutMtri2=true;
      float minpt=data.track_pt->at(kid);
      if (data.track_pt->at(kid) >data.track_pt->at(eid))  minpt=data.track_pt->at(eid);
      if (minTrk<minpt) minTrk=minpt;
    }//loop over triplets
    
    if (found) nbkg++;
    if(found && CutMutag) ++NMutag;
    else continue;

    if (found && CutDz) nDz++;
    if (found && CutDz && CutMll) nMll++;
    if (found && CutDz && CutMll && CutKpt) nKpt++;
    if (found && CutDz && CutMll && CutKpt && CutMtri) nMtri++;
    // if (found && CutDz && CutMll && CutKpt && CutMtri ) nSoftEl++;
    if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb){ nProb++;
      for (int i=0; i<10; i++){
	if (maxcos>Cos_cuts[i]) BkgEffCos[i]++;}}
    
    if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb && CutCos){ nCos++;
      for (int i=0; i<10; i++){
	if (maxPtB>PtB_cuts[i]) BkgEffPtB[i]++;}}
    
    if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb && CutCos && CutPtB){ nPtB++;
      for (int i=0; i<10; i++){
	if (maxLxy>Lxy_cuts[i]) BkgEffLxy[i]++;}}
    
    if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb && CutCos && CutPtB && CutLxy) nLxy++; 
    if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb && CutCos && CutPtB && CutLxy && CutSoft) nSoftEl++;   
    if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb && CutCos && CutPtB && CutLxy){
      for (int i=0; i<10; i++){
	if (Mtri>M_cuts[i]) BkgEffM[i]++;}}
  if (found && CutDz && CutMll && CutKpt && CutMtri && CutProb && CutCos && CutPtB && CutLxy && CutMtri2){
    for (int i=0; i<10; i++){
      if (minTrk>Trk_cuts[i]) BkgEffTrk[i]++;}} 
  }

  
  //  hBkgQ->SetBinContent(1,nNoSoft);  hBkgQ->SetBinContent(2,nSoft); 


 std::cout << " Triplet " << nTotal << " Ntriplets1ele " << nbkg << " ele1Trg dz " << nDz << " mll < 5 " << nMll
	   << " kpt > 0.4 " << nKpt << " Mtrip in [3.5,6] " << nMtri << " VtxCl Prob > 1.e-3 " << nProb 
	   << " Bcos > 0.9 " << nCos << " ptb > 5 " << nPtB << " Lxy > 4 " << nLxy << std::endl;


 std::cout << " \n\n Total = " << nTotal
   //<< " \n Nreco triplets = " << nbkg
           << " \n NrecoTriplets and tagMu = " << NMutag
           << " \n Ele1-Trg dz = " << nDz
           << " \n mll < 5 = " << nMll
   //<< " \n ll_vtxCL > 1.e-34 = NA"
           << " \n kpt > 0.4 = " << nKpt
           << " \n Mtrip in [3.5,6] = " << nMtri
           << " \n B_vtxCL > 1.e-3 = " << nProb
           << " \n Bcos > 0.9 = " << nCos
           << " \n ptB >5  = " << nPtB
           << " \n Lxy > 4 = " << nLxy << std::endl;



 for(int i=0; i<10; i++){
   BkgEffCos[i]= BkgEffCos[i]/ccdata->GetEntries();
   SigEffCos[i]=SigEffCos[i]/ngenAcceptance;
   BkgEffPtB[i]= BkgEffPtB[i]/ccdata->GetEntries();
   SigEffPtB[i]=SigEffPtB[i]/ngenAcceptance;
   BkgEffLxy[i]= BkgEffLxy[i]/ccdata->GetEntries();
   SigEffLxy[i]=SigEffLxy[i]/ngenAcceptance;
   BkgEffM[i]= BkgEffM[i]/ccdata->GetEntries();
   SigEffM[i]=SigEffM[i]/ngenAcceptance;
   BkgEffTrk[i]= BkgEffTrk[i]/ccdata->GetEntries();
   SigEffTrk[i]=SigEffTrk[i]/ngenAcceptance;
 }

 
 TGraph * rocCos = new TGraph(10,BkgEffCos,SigEffCos);
 TCanvas *crocCos = new TCanvas();
 crocCos->cd();
 rocCos->SetMarkerStyle(7);
 rocCos->GetXaxis()->SetTitle("bkg eff. (vtxCL)");
 rocCos->GetYaxis()->SetTitle("signal eff. (vtxCL)");
 rocCos->Draw("ap");

 TGraph * rocPtB = new TGraph(10,BkgEffPtB,SigEffPtB);
 TCanvas *crocPtB = new TCanvas();
 crocPtB->cd();
 rocPtB->SetMarkerStyle(7);
 rocPtB->GetXaxis()->SetTitle("bkg eff. (PtB)");
 rocPtB->GetYaxis()->SetTitle("signal eff. (PtB)");
 rocPtB->Draw("ap");

 TGraph * rocLxy=new TGraph(10,BkgEffLxy,SigEffLxy);
 TCanvas *crocLxy = new TCanvas();
 crocLxy->cd();
 rocLxy->SetMarkerStyle(7);
 rocLxy->GetXaxis()->SetTitle("bkg eff. (Lxy)");
 rocLxy->GetYaxis()->SetTitle("signal eff. (Lxy)");
 rocLxy->Draw("ap");

 TGraph * rocM=new TGraph(10,BkgEffM,SigEffM);
 TCanvas *crocM = new TCanvas();
 crocM->cd();
 rocM->SetMarkerStyle(7);
 rocM->GetXaxis()->SetTitle("bkg eff. (M)");
 rocM->GetYaxis()->SetTitle("signal eff. (M)");
 rocM->Draw("ap");

 TGraph * rocTrk=new TGraph(10,BkgEffTrk,SigEffTrk);
 TCanvas *crocTrk = new TCanvas();
 crocTrk->cd();
 rocTrk->SetMarkerStyle(7);
 rocTrk->GetXaxis()->SetTitle("bkg eff. (Trk)");
 rocTrk->GetYaxis()->SetTitle("signal eff. (Trk)");
 rocTrk->Draw("ap");

TFile * fout= new TFile("FastPlot.root","RECREATE");
 hBptGen->Write();
 hE1ptGen->Write(); 
 hE2ptGen->Write(); 
 hKptGen->Write();
 hDReltrk->Write();
 hDRK->Write();
 hDRrecoE->Write();
 hBRecoPt->Write();
 hE1RecoPt->Write();
 hE2RecoPt->Write();
 hKRecoPt->Write();


return 0;
} 
