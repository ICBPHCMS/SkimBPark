//root studyEfficiencyRate_lowPtEle.C  <= for efficiency and rate on Kee
//root studyEfficiencyRate_lowPtEle.C'(1)'  <= for efficiency and rate on K*ee
//update the names of the input files


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



int studyEfficiencyRate_lowPtEle(int isKstar = 0){
  TChain *ccmc = new TChain("demo/mytree");
  TChain *ccdata = new TChain("demo/mytree");

  if(!isKstar){
    ccdata->Add("../plugins/outputBSkim_dataGeorge_test_lowPtele.root");
    ccmc->Add("../plugins/outputBSkim_KeeVince_test_lowPtele.root");
  }
  else{
    ccmc->Add("../test/outputBSkim_KsteeVince_test_lowPtele.root");
    ccdata->Add("../test/outputBSkim_KsteeVince_test_lowPtele.root");
  }
  skim_class mc;
  mc.Init(ccmc);

  skim_class data;
  data.Init(ccdata);

  //match gen
  TH1F* hBptGen = new TH1F("hBptGen"," ",60,0,60);
  TH1F* hE1ptGen = new TH1F("hE1ptGen"," ",60,0,60);
  TH1F* hE2ptGen = new TH1F("hE2ptGen"," ",60,0,60);
  TH1F* hKptGen = new TH1F("hKptGen"," ",60,0,60);
  TH1F* hPptGen = new TH1F("hPptGen"," ",60,0,60);

  TH1F* hKstarMassGen = new TH1F("hKstarMassGen", "", 100, 0., 2.);
  TH1F* hKstarMassReco = new TH1F("hKstarMassReco", "", 100, 0., 2.);
  TH1F* hllMassGen = new TH1F("hllMassGen", "", 100, 0., 6.);
  TH1F* hllMassReco = new TH1F("hllMassReco", "", 100, 0., 6.);

  TH1F* hDRP = new TH1F("hDRP"," ",100,0,0.3);
  TH1F* hDRK = new TH1F("hDRK"," ",100,0,0.3);
  TH1F* hDReltrk = new TH1F("hDReltrk"," ",100,0,0.3);
  TH1F* hDRrecoE = new TH1F("hDRrecoE"," ",100,0,0.3);
  TH1F* hDRGlobal = new TH1F("hDRGlobal", "", 100,0,0.3);

  TH1F* hBRecoPt = new TH1F("hBRecoPt"," ",60,0,60);
  TH1F* hE1RecoPt = new TH1F("hE1RecoPt"," ",60,0,60);
  TH1F* hE2RecoPt = new TH1F("hE2RecoPt"," ",60,0,60);
  TH1F* hKRecoPt = new TH1F("hKRecoPt"," ",60,0,60);
  TH1F* hPRecoPt = new TH1F("hPRecoPt"," ",60,0,60);


  float Cos_cuts[10] = {-1.1, 0, 0.5,0.7, 0.8, 0.9, 0.95,0.99,0.995,0.999};
  float PtB_cuts[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  float Lxy_cuts[10] = {0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6 };
  float M_cuts[10] = {3.5, 3.7, 3.9, 4.1, 4.3, 4.5, 4.7, 4.9, 5.1, 5.27 };
  float Trk_cuts[10] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8 };
  float SigEffCos[10] = {0}, BkgEffCos[10] = {0};
  float SigEffPtB[10] = {0}, BkgEffPtB[10] = {0}; 
  float SigEffLxy[10] = {0}, BkgEffLxy[10] = {0};
  float SigEffM[10] = {0}, BkgEffM[10] = {0};
  float SigEffTrk[10] = {0}, BkgEffTrk[10] = {0};

  float BDT_cuts[100];

  float SigBDT[100] = {0.};
  float SigBDTbiased[100] = {0.};
  float SigBDTnoPE[100] = {0.};
  float SigBDTbiasednoPE[100] = {0.};
  float SigBDT_DR[100] = {0.};
  float SigBDTbiased_DR[100] = {0.};
  float SigBDTnoPE_DR[100] = {0.};
  float SigBDTbiasednoPE_DR[100] = {0.};

  float BkgBDT[100] = {0.};
  float BkgBDTbiased[100] = {0.};
  float BkgBDTnoPE[100] = {0.};
  float BkgBDTbiasednoPE[100] = {0.};

  for(int ij=0; ij<100; ++ij){
    BDT_cuts[ij] = 0.1*ij;
  }

  // float PDG_BCharged_mass = 5.279;
  // float PDG_B0_mass = 5.279;
  float PDG_B_mass = 5.279;
  float PDG_K0st_mass = 0.89581;
  float KaonMass = 0.493677;
  float PionMass = 0.139570;
  float MuonMass = 0.10565837;
  float ElectronMass = 0.5109989e-3;

  float dRMCmatching_val = 0.06;
  float dzProbeTagMuon_val = 0.3;
  float bdt_val = 6.5;
  float mllCut_val = 100.;
  float k0starMassdiff_minval = 0.; //20MeV in PDG error is 0.2MeV
  float k0starMassdiff_maxval = 100; //20MeV in PDG error is 0.2MeV
  float pTK_val = 0.;
  float mBmax_val = 100.;
  float mBmin_val = 1.;
  float vtxCLB_val = 1.e-3;
  float cosAlpha_val = 0.9;
  float pTB_val = 5.;
  float IPSB_val = 4.;

  //signal
  float nEventsMC = ccmc->GetEntries();
  std::cout << " \n  mc events = " << nEventsMC << std::endl;
  
  float ngen[4] = {0., 0., 0., 0.};
  float ngenAcceptance[4] = {0., 0., 0., 0.};
  float nTrip_matched[4] = {0., 0., 0., 0.};
  float nTrip_matched_dR[4] = {0., 0., 0., 0.};
  float nTrip_matched_MuTag[4] = {0., 0., 0., 0.};
  float ndZ_cut[4] = {0., 0., 0., 0.};
  float nMll_cut[4] = {0., 0., 0., 0.};
  float nK0StarMass_cut[4] = {0., 0., 0., 0.};
  float nK0StarVtxCL_cut[4] = {0., 0., 0., 0.};
  float npTK_cut[4] = {0., 0., 0., 0.};
  float nMB_cut[4] = {0., 0., 0., 0.};
  float nVtxCLB_cut[4] = {0., 0., 0., 0.};
  float nCosAlpha_cut[4] = {0., 0., 0., 0.};
  float npTB_cut[4] = {0., 0., 0., 0.};
  float nIPSB_cut[4] = {0., 0., 0., 0.};


  for (int ev=0; ev<nEventsMC; ev++){
    mc.GetEntry(ev);

    //take reco candidate closest to gen
    TLorentzVector vel1, vel2, vK, vP;

    if(mc.genpart_B_index == -1) continue;

    bool BKee = false;
    float Kcharge = 0.;
    float Picharge = 0.;
    float Ele1charge = 0.;
    float Ele2charge = 0.;


    bool CutGenMatched[2] = {false, false};

    if(!isKstar){
      //Kee
      vel1.SetPtEtaPhiM(mc.genpart_lep1_PtEtaPhiM->at(0), mc.genpart_lep1_PtEtaPhiM->at(1), mc.genpart_lep1_PtEtaPhiM->at(2), mc.genpart_lep1_PtEtaPhiM->at(3));
      vel2.SetPtEtaPhiM(mc.genpart_lep2_PtEtaPhiM->at(0), mc.genpart_lep2_PtEtaPhiM->at(1), mc.genpart_lep2_PtEtaPhiM->at(2), mc.genpart_lep2_PtEtaPhiM->at(3));
      vK.SetPtEtaPhiM(mc.genpart_K_PtEtaPhiM->at(0), mc.genpart_K_PtEtaPhiM->at(1), mc.genpart_K_PtEtaPhiM->at(2), mc.genpart_K_PtEtaPhiM->at(3));
      
      BKee = true;
      Kcharge = (mc.genpart_KFromB_pdg > 0) ? 1 : -1;
      Ele1charge = (mc.genpart_lep1FromB_pdg > 0) ? 1 : -1;
      Ele2charge = (mc.genpart_lep2FromB_pdg > 0) ? 1 : -1;
    }
    else{
      //Kstee
      vel1.SetPtEtaPhiM(mc.genpart_lep1_PtEtaPhiM->at(0), mc.genpart_lep1_PtEtaPhiM->at(1), mc.genpart_lep1_PtEtaPhiM->at(2), mc.genpart_lep1_PtEtaPhiM->at(3));
      vel2.SetPtEtaPhiM(mc.genpart_lep2_PtEtaPhiM->at(0), mc.genpart_lep2_PtEtaPhiM->at(1), mc.genpart_lep2_PtEtaPhiM->at(2), mc.genpart_lep2_PtEtaPhiM->at(3));
      vK.SetPtEtaPhiM(mc.genpart_K_PtEtaPhiM->at(0), mc.genpart_K_PtEtaPhiM->at(1), mc.genpart_K_PtEtaPhiM->at(2), mc.genpart_K_PtEtaPhiM->at(3));
      vP.SetPtEtaPhiM(mc.genpart_Pi_PtEtaPhiM->at(0), mc.genpart_Pi_PtEtaPhiM->at(1), mc.genpart_Pi_PtEtaPhiM->at(2), mc.genpart_Pi_PtEtaPhiM->at(3));

      hKstarMassGen->Fill((vK+vP).M());
      hllMassGen->Fill((vel1+vel2).M());

      BKee = true;
      Kcharge = (mc.genpart_KFromKst_pdg > 0) ? 1 : -1;
      Picharge = (mc.genpart_PiFromKst_pdg > 0) ? 1 : -1;
      Ele1charge = (mc.genpart_lep1FromB_pdg > 0) ? 1 : -1;
      Ele2charge = (mc.genpart_lep2FromB_pdg > 0) ? 1 : -1;
    }

    if (!BKee) continue;
    ++ngen[0];

    if( DR(vel1.Eta(), vel1.Phi(), vel2.Eta(), vel2.Phi()) > 0.02) ++ngen[2];

    if(std::abs(vel1.Eta()) > 2.5 || std::abs(vel2.Eta()) > 2.5 || std::abs(vK.Eta()) > 2.5 ||
       vel1.Pt() < 0.2 || vel2.Pt() < 0.2 || vK.Pt() < 0.2) continue;
    if(isKstar && (vP.Pt() < 0.2 || std::abs(vP.Eta()) > 2.5 )) continue;
    ++ngenAcceptance[0];
    CutGenMatched[0] = true;

    if( DR(vel1.Eta(), vel1.Phi(), vel2.Eta(), vel2.Phi()) > 0.02){
      ++ngenAcceptance[2];
      CutGenMatched[1] = true;
    }


     hBptGen->Fill((vel1+vel2+vK).Pt()); 
     hE1ptGen->Fill(vel1.Pt());
     hE2ptGen->Fill(vel2.Pt());
     hKptGen->Fill(vK.Pt()); 
     if(isKstar) hPptGen->Fill(vP.Pt()); 
     //     std::cout << " >>> (vel1+vel2+vK).M() = " << (vel1+vel2+vK).M() << std::endl;
     //     std::cout << " mc triplets = " << mc.TTrack_cos->size() << std::endl;


     float minDRP[2] = {1000., 1000.};
     float minDRK[2] = {1000., 1000.};
     float minDRrecoE[2] = {1000., 1000.};
     float minDRtrkE[2] = {1000., 1000.};

     int recoE[2] = {-1, -1};
     int trkEl[2] = {-1, -1};
     int trkK[2] = {-1, -1};
     int trkP[2] = {-1, -1};
     int Bindex[2] = {-1, -1};
     int selVetx[2] = {-1, -1};
 
     TLorentzVector genB = isKstar ? vP+vK+vel1+vel2 : vK+vel1+vel2;
     TLorentzVector vrel, vreltrk, vrK, vrP, vrB;
     TLorentzVector vrelPE, vreltrkPE, vrKPE, vrPPE, vrBPE;

     float globaldR[2] = {1000., 1000.};
     float bestPtRatio[2] = {1000., 1000.};

     bool found = false; 
     bool foundPE = false; 

     float dR_BTD[100];
     float dR_BTDbiased[100];
     float dR_BTDnoPE[100];
     float dR_BTDbiasednoPE[100];
     for(int ij=0; ij<100; ++ij){
       dR_BTD[ij] = -1;
       dR_BTDbiased[ij] = -1;
       dR_BTDnoPE[ij] = -1;
       dR_BTDbiasednoPE[ij] = -1;
     }

     for (unsigned int iB=0; iB<mc.TTrack_cos->size(); iB++)  {
       if (mc.TTrack_ObjId->at(iB) != 11) continue;
       //std::cout << " >> mc.TTrack_ObjId->at(iB) = " << mc.TTrack_ObjId->at(iB) << std::endl;
       
       int recoEd = mc.TTrack_ObjIndex->at(iB);
       int trkEld = mc.TTrack_TrkIndex->at(iB);
       int trkKd = mc.TTrack_kid->at(iB);
       int trkPd = isKstar ? mc.TTrack_piid->at(iB) : -1;

       float bdt = mc.gsfTrk_seedBDTunb->at(recoEd);
       float bdtbiased = mc.gsfTrk_seedBDTbiased->at(recoEd);

       float dR_ele1 = DR(vel1.Eta(), vel1.Phi(), mc.gsfTrk_eta->at(recoEd), mc.gsfTrk_phi->at(recoEd));
       float dR_ele2 = DR(vel2.Eta(), vel2.Phi(), mc.gsfTrk_eta->at(trkEld), mc.gsfTrk_phi->at(trkEld));
       float dR_trk = DR(vK.Eta(), vK.Phi(), mc.track_eta->at(trkKd), mc.track_phi->at(trkKd));
       float dR_trp = isKstar ? DR(vP.Eta(), vP.Phi(), mc.track_eta->at(trkPd), mc.track_phi->at(trkPd)) : 0;

       float withBdR = dR_ele1+dR_ele2+dR_trk+dR_trp;

       if(CutGenMatched[0] &&  withBdR < globaldR[0] && bdt > bdt_val){       
       	 bestPtRatio[0] = vrB.Pt()/genB.Pt();
	 globaldR[0] = withBdR;

	 minDRrecoE[0] = dR_ele1;
	 recoE[0] = recoEd;
	 
	 minDRtrkE[0] = dR_ele2;
	 trkEl[0] = trkEld;
	 
	 minDRK[0] = dR_trk;
	 trkK[0] = trkKd;

         minDRP[0] = dR_trp;
         trkP[0] = trkPd;

         Bindex[0] = iB;
	 //std::cout << " >>> found triplet " << std::endl;
       }
       float dRrecoEleEle = DR(mc.gsfTrk_eta->at(recoEd), mc.gsfTrk_phi->at(recoEd), mc.gsfTrk_eta->at(trkEld), mc.gsfTrk_phi->at(trkEld));
       if(CutGenMatched[1] && dRrecoEleEle > 0.02 && withBdR < globaldR[1] && bdt > bdt_val){
	 bestPtRatio[1] = vrB.Pt()/genB.Pt();
	 globaldR[1] = withBdR;

	 minDRrecoE[1] = dR_ele1;
	 recoE[1] = recoEd;
	 
	 minDRtrkE[1] = dR_ele2;
	 trkEl[1] = trkEld;
	 
	 minDRK[1] = dR_trk;
	 trkK[1] = trkKd;

         minDRP[1] = dR_trp;
         trkP[1] = trkPd;

         Bindex[1] = iB;
       }

       for(int ij=0; ij<100; ++ij){
	 if(CutGenMatched[0]){
	   if(bdt > BDT_cuts[ij]){
	     if(withBdR < dR_BTD[ij] || dR_BTD[ij] == -1) dR_BTD[ij] = withBdR;
	   }
	   if(bdtbiased > BDT_cuts[ij]){
	     if(withBdR < dR_BTDbiased[ij] || dR_BTDbiased[ij] == -1) dR_BTDbiased[ij] = withBdR;
	   }
	 }

	 if(CutGenMatched[1] && dRrecoEleEle > 0.02){
	   if(bdt > BDT_cuts[ij]){
	     if(withBdR < dR_BTDnoPE[ij] || dR_BTDnoPE[ij] == -1) dR_BTDnoPE[ij] = withBdR;
	   }
	   if(bdtbiased > BDT_cuts[ij]){
	     if(withBdR < dR_BTDbiasednoPE[ij] || dR_BTDbiasednoPE[ij] == -1) dR_BTDbiasednoPE[ij] = withBdR;
	   }	   
	 }
       }


     }//loop over triplets
     
     if(trkK[0] != -1 && trkEl[0] != -1 && recoE[0] != -1 && (!isKstar || trkP[0] != -1 )){
       found = true; 
       selVetx[0] = Bindex[0];
     }
     if(trkK[1] != -1 && trkEl[1] != -1 && recoE[1] != -1 && (!isKstar || trkP[1] != -1 )){
       foundPE = true; 
       selVetx[1] = Bindex[1];
     }

     if(!found) continue; 
    
     if(mc.SelectedMu_index != -1){
       for(int ij=0; ij<100; ++ij){
	 
	 if(CutGenMatched[0] && dR_BTD[ij] != -1) ++SigBDT[ij];
	 if(CutGenMatched[0] && dR_BTDbiased[ij] != -1) ++SigBDTbiased[ij];
	 if(CutGenMatched[1] && dR_BTDnoPE[ij] != -1) ++SigBDTnoPE[ij];
	 if(CutGenMatched[1] && dR_BTDbiasednoPE[ij] != -1) ++SigBDTbiasednoPE[ij];
	 
	 if(CutGenMatched[0] && ( (isKstar && dR_BTD[ij] < 4.*dRMCmatching_val) ||
				  (!isKstar && dR_BTD[ij] < 3.*dRMCmatching_val)) )  ++SigBDT_DR[ij];
	 if(CutGenMatched[0] && ( (isKstar && dR_BTDbiased[ij] < 4.*dRMCmatching_val) ||
				  (!isKstar && dR_BTDbiased[ij] < 3.*dRMCmatching_val)) ) ++SigBDTbiased_DR[ij];
	 
	 if(CutGenMatched[1] && ( (isKstar && dR_BTDnoPE[ij] < 4.*dRMCmatching_val) ||
				  (!isKstar && dR_BTDnoPE[ij] < 3.*dRMCmatching_val)) )  ++SigBDTnoPE_DR[ij];
	 if(CutGenMatched[1] && ( (isKstar && dR_BTDbiasednoPE[ij] < 4.*dRMCmatching_val) ||
				  (!isKstar && dR_BTDbiasednoPE[ij] < 3.*dRMCmatching_val)) )  ++SigBDTbiasednoPE_DR[ij];
       }
     }


     if(found) ++nTrip_matched[0];     
     if(foundPE) ++nTrip_matched[2];     

     bool CutMatch[2] = {false, false};
     bool CutMuTag[2] = {false, false};
     bool CutDz[2] = {false, false};
     bool CutMll[2] = {false, false};
     bool CutKstar[2] = {false, false};
     bool CutKstarVtxCL[2] = {false, false};
     bool CutKpt[2] = {false, false};
     bool CutMB[2] = {false, false};
     bool CutProb[2] = {false, false};
     bool CutCos[2] = {false, false};
     bool CutPtB[2] = {false, false};
     bool CutLxy[2] = {false, false};

     if(found){
     vrel.SetPtEtaPhiM(mc.gsfTrk_pt->at(recoE[0]), mc.gsfTrk_eta->at(recoE[0]), mc.gsfTrk_phi->at(recoE[0]), ElectronMass);
     vreltrk.SetPtEtaPhiM(mc.gsfTrk_pt->at(trkEl[0]), mc.gsfTrk_eta->at(trkEl[0]),mc.gsfTrk_phi->at(trkEl[0]), ElectronMass);   
     vrK.SetPtEtaPhiM(mc.track_pt->at(trkK[0]),mc.track_eta->at(trkK[0]),mc.track_phi->at(trkK[0]), KaonMass);
     if(isKstar) vrP.SetPtEtaPhiM(mc.track_pt->at(trkP[0]),mc.track_eta->at(trkP[0]),mc.track_phi->at(trkP[0]), PionMass);     
     vrB = isKstar ? vrP+vrK+vrel+vreltrk : vrK+vrel+vreltrk;
     }

     if(foundPE){
     vrelPE.SetPtEtaPhiM(mc.gsfTrk_pt->at(recoE[1]), mc.gsfTrk_eta->at(recoE[1]), mc.gsfTrk_phi->at(recoE[1]), ElectronMass);
     vreltrkPE.SetPtEtaPhiM(mc.gsfTrk_pt->at(trkEl[1]), mc.gsfTrk_eta->at(trkEl[1]),mc.gsfTrk_phi->at(trkEl[1]), ElectronMass);   
     vrKPE.SetPtEtaPhiM(mc.track_pt->at(trkK[1]),mc.track_eta->at(trkK[1]),mc.track_phi->at(trkK[1]), KaonMass);
     if(isKstar) vrPPE.SetPtEtaPhiM(mc.track_pt->at(trkP[1]),mc.track_eta->at(trkP[1]),mc.track_phi->at(trkP[1]), PionMass);     
     vrBPE = isKstar ? vrPPE+vrKPE+vrelPE+vreltrkPE : vrKPE+vrelPE+vreltrkPE;
     }

     hDRrecoE->Fill(minDRrecoE[0]);
     hDReltrk->Fill(minDRtrkE[0]);
     hDRK->Fill(minDRK[0]);
     if(isKstar) hDRP->Fill(minDRP[0]);
     hDRGlobal->Fill(globaldR[0]);

     if((isKstar && globaldR[0] < 4.*dRMCmatching_val) || (globaldR[0] < 3.*dRMCmatching_val)){ CutMatch[0] = true;  ++nTrip_matched_dR[0]; }
     if((isKstar && globaldR[1] < 4.*dRMCmatching_val) || (globaldR[1] < 3.*dRMCmatching_val)){ CutMatch[1] = true;  ++nTrip_matched_dR[2]; }

     hBRecoPt->Fill((vrel+vreltrk+vrK).Pt()); 
     hE1RecoPt->Fill(vrel.Pt()); 
     hE2RecoPt->Fill(vreltrk.Pt()); 
     hKRecoPt->Fill(vrK.Pt()); 
     if(isKstar) {
       hPRecoPt->Fill(vrP.Pt());
       hKstarMassReco->Fill((vrK+vrP).M());
       hllMassReco->Fill((vrel+vreltrk).M());
     }
     
    if(mc.SelectedMu_index == -1) continue;
    if(CutMatch[0]) { CutMuTag[0] = true;   ++nTrip_matched_MuTag[0];}
    if(CutMatch[1]) { CutMuTag[1] = true;  ++nTrip_matched_MuTag[2];}
  
    /*
    if(CutMuTag[0] && fabs(mc.muon_vz->at(mc.SelectedMu_index) - mc.gsfTrk_vz->at(recoE[0]) ) <  dzProbeTagMuon_val) {
      CutDz[0] = true;
      ++ndZ_cut[0];
    }
    if(CutMuTag[1] && fabs(mc.muon_vz->at(mc.SelectedMu_index) - mc.gsfTrk_vz->at(recoE[1]) ) <  dzProbeTagMuon_val){
      CutDz[1] = true;
      ++ndZ_cut[2];
    }
    */

    if( CutMuTag[0] && (vrel+vreltrk).M() < mllCut_val) {
      CutMll[0] = true;
      ++nMll_cut[0];
    }
    if( CutMuTag[1] && (vrelPE+vreltrkPE).M() < mllCut_val) {
      CutMll[1] = true;
      ++nMll_cut[2];
    }


    if(isKstar){
      if( CutMll[0] && (vrK+vrP).M() > k0starMassdiff_minval && (vrK+vrP).M() < k0starMassdiff_maxval){
	CutKstar[0] = true;
	++nK0StarMass_cut[0];
      }
      if( CutMll[1] && (vrKPE+vrPPE).M() > k0starMassdiff_minval && (vrKPE+vrPPE).M() < k0starMassdiff_maxval){
	CutKstar[1] = true;
	++nK0StarMass_cut[2];
      }
    }

    TLorentzVector vrKfromB = isKstar ? vrK+vrP : vrK;
    TLorentzVector vrKfromBPE = isKstar ? vrKPE+vrPPE : vrKPE;

    if(((isKstar && CutKstar[0]) || (!isKstar && CutMll[0]) ) && vrKfromB.Pt() > pTK_val){
      CutKpt[0] = true;
      ++npTK_cut[0];
    }
    if(((isKstar && CutKstar[1]) || (!isKstar && CutMll[1]) ) && vrKfromBPE.Pt() > pTK_val){
      CutKpt[1] = true;
      ++npTK_cut[2];
    }

    if(CutKpt[0])
      for (int i=0; i<10; i++){
	if ( mc.TTrack_PtEtaPhiM->at(selVetx[0]).at(3) > M_cuts[i]) ++SigEffM[i];
      }
    if( CutKpt[0] && ( vrB.M() > mBmin_val && vrB.M() < mBmax_val) ) {
      CutMB[0] = true;
      ++nMB_cut[0];
    }
    if( CutKpt[1] && ( vrBPE.M() > mBmin_val && vrBPE.M() < mBmax_val) ) {
      CutMB[1] = true;
      ++nMB_cut[2];
    }

    if(CutMB[0] && mc.TTrack_chi_prob->at(selVetx[0]) > vtxCLB_val){
      CutProb[0] = true;
      ++nVtxCLB_cut[0];
    }
    if(CutMB[1] && mc.TTrack_chi_prob->at(selVetx[1]) > vtxCLB_val){
      CutProb[1] = true;
      ++nVtxCLB_cut[2];
    }

    if(CutProb[0])
      for (int i=0; i<10; ++i){
	if (mc.TTrack_cos->at(selVetx[0]) > Cos_cuts[i]) ++SigEffCos[i];
      }

    if(CutProb[0] && mc.TTrack_cos->at(selVetx[0]) > cosAlpha_val){
      CutCos[0] = true;
      ++nCosAlpha_cut[0];
    }
    if(CutProb[1] && mc.TTrack_cos->at(selVetx[1]) > cosAlpha_val){
      CutCos[1] = true;
      ++nCosAlpha_cut[2];
    }

    if(CutCos[0])
      for (int i=0; i<10; i++){
	if (mc.TTrack_PtEtaPhiM->at(selVetx[0]).at(0) > PtB_cuts[i]) ++SigEffPtB[i];  
      }


    if (CutCos[0] && mc.TTrack_PtEtaPhiM->at(selVetx[0]).at(0) > pTB_val ){
      CutPtB[0] = true;
      ++npTB_cut[0];
    }
    if (CutCos[1] && mc.TTrack_PtEtaPhiM->at(selVetx[1]).at(0) > pTB_val ){
      CutPtB[1] = true;
      ++npTB_cut[2];
    }

    if(CutPtB[0])
      for (int i=0; i<10; i++){
	if (mc.TTrack_Lxy->at(selVetx[0])/TMath::Sqrt(mc.TTrack_eLxy->at(selVetx[0])) > Lxy_cuts[i]) ++SigEffLxy[i];  
      }
    
    if (CutPtB[0] && mc.TTrack_Lxy->at(selVetx[0])/TMath::Sqrt(mc.TTrack_eLxy->at(selVetx[0])) > IPSB_val ){
      CutLxy[0] = true;
      ++nIPSB_cut[0];
    }
    if (CutPtB[1] && mc.TTrack_Lxy->at(selVetx[1])/TMath::Sqrt(mc.TTrack_eLxy->at(selVetx[1])) > IPSB_val ){
      CutLxy[1] = true;
      ++nIPSB_cut[2];
    }

    if(CutLxy[0])
      for (int i=0; i<10; i++){
	if (vreltrk.Pt() > Trk_cuts[i] && vrK.Pt() > Trk_cuts[i] && (!isKstar || vrP.Pt() > Trk_cuts[i])) ++SigEffTrk[i];
    }
  }//MC events


  for(int ij=0; ij<3; ++ij){
    if(ij == 0) std::cout << "\n  K* with charge exchange " << std::endl;
    else if(ij == 2) std::cout << "\n  K* withOUT charge exchange " << std::endl;
    else continue;

    std::cout << " Ngen = " << ngen[ij] 
	      << " \n Ngen in acceptance = " << ngenAcceptance[ij] 
	      << " \n Nreco triplets genMatched = " << nTrip_matched[ij]
	      << " \n Nreco triplets genMatched in dR = " << nTrip_matched_dR[ij]
	      << " \n NrecoTriplets and tagMu = " << nTrip_matched_MuTag[ij]
	      << " \n Ele1-Trg dz = " << ndZ_cut[ij] 
	      << " \n mll < 5 = " << nMll_cut[ij] << std::endl;
    if(isKstar) std::cout << " d(mK*) < 0.02 = " << nK0StarMass_cut[ij] << std::endl;
    std::cout << " kpt > 0.4 = " << npTK_cut[ij]
	      << " \n Mtrip in [4.1,6] = " << nMB_cut[ij]
	      << " \n B_vtxCL > 1.e-3 = " << nVtxCLB_cut[ij]
	      << " \n Bcos > 0.9 = " << nCosAlpha_cut[ij] 
	      << " \n ptB > 5  = " << npTB_cut[ij]
	      << " \n Lxy > 4 = " << nIPSB_cut[ij] << std::endl;
  }

  //bakg
  float nEventsData = ccdata->GetEntries();
  std::cout << " \n data events " << nEventsData << std::endl;

  for (int ev=0; ev<nEventsData; ev++){
    //for (int ev=0; ev<100; ev++){
    data.GetEntry(ev);
    //if (ev%100==0) std::cout << "reading evt " << ev << std::endl;

    bool Triplet[2] = {false, false};
    bool CutMutag[2] = {false, false};
    bool CutDz[2] = {false, false}; 
    bool CutMll[2] = {false, false}; 
    bool CutKstar[2] = {false, false}; 
    bool CutKstarVtxCL[2] = {false, false};
    bool CutKpt[2] = {false, false};
    bool CutMB[2] = {false, false}; 
    bool CutProb[2] = {false, false}; 
    bool CutCos[2] = {false, false}; 
    bool CutPtB[2] = {false, false};
    bool CutLxy[2] = {false, false}; 


    //cuts
    //int tntrk=0,tnel=0;
    float maxcos[2] = {-100, -100};
    float maxPtB[2] = {-100, -100};
    float maxLxy[2] = {-100, -100};
    float Mtri[2] = {-100, -100};
    float minTrk[2] = {-100, -100};

    bool BDT_found[100];
    bool BDTbiased_found[100];
    bool BDTnoPE_found[100];
    bool BDTbiasednoPE_found[100];
    for(int ij=0; ij<100; ++ij){
      BDT_found[ij] = false;
      BDTbiased_found[ij] = false;
      BDTnoPE_found[ij] = false;
      BDTbiasednoPE_found[ij] = false;
    }

    if( data.TTrack_cos->size() <= 0) continue; 
    Triplet[0] = true;
    
    for (int ivtx=0; ivtx<data.TTrack_cos->size(); ivtx++){
      if (data.TTrack_ObjId->at(ivtx) != 11){ 
	//  std::cout << " data.TTrack_ObjId->at(ivtx) = " << data.TTrack_ObjId->at(ivtx) << std::endl; 
	continue;
      }


      int recoEd = data.TTrack_ObjIndex->at(ivtx);
      int trkEld = data.TTrack_TrkIndex->at(ivtx);

      float bdt = data.gsfTrk_seedBDTunb->at(recoEd);
      float bdtbiased = data.gsfTrk_seedBDTbiased->at(recoEd);
      bool found[2] = {true, true}; 

      if(bdt < bdt_val){
	found[0] = false;
	found[1] = false;
      }

      float dR_ele1ele2 = DR( data.gsfTrk_eta->at(recoEd), data.gsfTrk_phi->at(recoEd), data.gsfTrk_eta->at(trkEld), data.gsfTrk_phi->at(trkEld));
      if(dR_ele1ele2 > 0.02) { /*std::cout << " dR_ele1ele2 = " << dR_ele1ele2 << std::endl; */ Triplet[1] = true;}
      else found[1] = false;

      if (data.SelectedMu_index < 0) continue;
      if(found[0]) CutMutag[0] = true;
      if(found[1]) CutMutag[1] = true;


      for(int ij=0; ij<100; ++ij){
	if(found[0] && bdt > BDT_cuts[ij]) BDT_found[ij] = true;
	if(found[0] && bdtbiased > BDT_cuts[ij]) BDTbiased_found[ij] = true;
	if(found[1] && bdt > BDT_cuts[ij]) BDTnoPE_found[ij] = true;
	if(found[1] && bdtbiased > BDT_cuts[ij]) BDTbiasednoPE_found[ij] = true;
      }



      
      float ZmuVtx = data.muon_vz->at(data.SelectedMu_index);

      /*
      if (found[0] && fabs(data.gsfTrk_vz->at(recoEd) - ZmuVtx) < dzProbeTagMuon_val) CutDz[0] = true;
      else found[0] = false;
      if (found[1] && fabs(data.gsfTrk_vz->at(recoEd) - ZmuVtx) < dzProbeTagMuon_val) CutDz[1] = true;
      else found[1] = false;
      */

      if (found[0] && data.TTrack_mll->at(ivtx) < mllCut_val) CutMll[0] = true;
      else found[0] = false;
      if (found[1] && data.TTrack_mll->at(ivtx) < mllCut_val) CutMll[1] = true;
      else found[1] = false;

      int KstarIndex = isKstar? data.TTrack_KstarIndex->at(ivtx) : -1;
      if(isKstar){
	if (found[0] && data.Kstpair_PtEtaPhiM->at(KstarIndex).at(3) > k0starMassdiff_minval && data.Kstpair_PtEtaPhiM->at(KstarIndex).at(3) < k0starMassdiff_maxval) CutKstar[0] = true;
	else found[0] = false;
	if (found[1] && data.Kstpair_PtEtaPhiM->at(KstarIndex).at(3) > k0starMassdiff_minval && data.Kstpair_PtEtaPhiM->at(KstarIndex).at(3) < k0starMassdiff_maxval) CutKstar[1] = true;
	else found[1] = false;
	// if(data.TTrack_mKst->at(ivtx) != data.Kstpair_PtEtaPhiM->at(KstarIndex).at(3)){
	//   std::cout << " >>> KstarIndex = " << KstarIndex << " data.Kstpair_PtEtaPhiM.size() = " << data.Kstpair_PtEtaPhiM->size() << std::endl;
	//   std::cout << " >>> problem index K* " << std::endl;
	// }

	// if(data.Kstpair_chi_prob->at(KstarIndex) < 1.e-3) continue;
	// CutKstarVtxCL = true;
      }


      int kid = data.TTrack_kid->at(ivtx);
      int pid = isKstar ? data.TTrack_piid->at(ivtx) : -1;
      if(found[0] && ((data.track_pt->at(kid) > pTK_val && !isKstar) || (isKstar && data.Kstpair_PtEtaPhiM->at(KstarIndex).at(0) > pTK_val) ) ) CutKpt[0] = true;
      else found[0] = false;
      if(found[1] && ((data.track_pt->at(kid) > pTK_val && !isKstar) || (isKstar && data.Kstpair_PtEtaPhiM->at(KstarIndex).at(0) > pTK_val) ) ) CutKpt[1] = true;
      else found[1] = false;

      if(found[0] && fabs(Mtri[0] - PDG_B_mass) > fabs( PDG_B_mass - data.TTrack_PtEtaPhiM->at(ivtx).at(3)) ) Mtri[0] = data.TTrack_PtEtaPhiM->at(ivtx).at(3);
      if(found[0] && (data.TTrack_PtEtaPhiM->at(ivtx).at(3) > mBmin_val && data.TTrack_PtEtaPhiM->at(ivtx).at(3) < mBmax_val) ) CutMB[0] = true;
      else found[0] = false;
      if(found[1] && (data.TTrack_PtEtaPhiM->at(ivtx).at(3) > mBmin_val && data.TTrack_PtEtaPhiM->at(ivtx).at(3) < mBmax_val) ) CutMB[1] = true;
      else found[1] = false;

      if(found[0] && data.TTrack_chi_prob->at(ivtx) > vtxCLB_val) CutProb[0] = true; 
      else found[0] = false;
      if(found[1] && data.TTrack_chi_prob->at(ivtx) > vtxCLB_val) CutProb[1] = true; 
      else found[1] = false;

      if (found[0] &&  maxcos[0] < data.TTrack_cos->at(ivtx)) maxcos[0] = data.TTrack_cos->at(ivtx);
      if (found[0] && data.TTrack_cos->at(ivtx) > cosAlpha_val)  CutCos[0] = true;
      else found[0] = false;
      if (found[1] && data.TTrack_cos->at(ivtx) > cosAlpha_val)  CutCos[1] = true;
      else found[1] = false;

      if (found[0] && maxPtB[0] < data.TTrack_PtEtaPhiM->at(ivtx).at(0)) maxPtB[0] = data.TTrack_PtEtaPhiM->at(ivtx).at(0);
      if (found[1] && data.TTrack_PtEtaPhiM->at(ivtx).at(0) > pTB_val) CutPtB[0] = true;
      else found[0] = false;
      if (found[1] && data.TTrack_PtEtaPhiM->at(ivtx).at(0) > pTB_val) CutPtB[1] = true;
      else found[1] = false;

      if (found[0] && maxLxy[0] < data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx)))
	maxLxy[0] = data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx));
      if (found[0] && data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx)) > IPSB_val) CutLxy[0] = true;
      else found[0] = false;
      if (found[1] && data.TTrack_Lxy->at(ivtx)/TMath::Sqrt(data.TTrack_eLxy->at(ivtx)) > IPSB_val) CutLxy[1] = true;
      else found[1] = false;
      
      float minpt = data.track_pt->at(kid);
      if (found[0] && data.track_pt->at(kid) > data.gsfTrk_pt->at(recoEd)) minpt = data.gsfTrk_pt->at(recoEd);
      if (found[0] && minTrk[0] < minpt) minTrk[0] = minpt;
      if(found[0] && isKstar && data.track_pt->at(pid) < minTrk[0]) minTrk[0] = data.track_pt->at(pid);

    }//loop over triplets
    
    for(int ij=0; ij<100; ++ij){
      if(BDT_found[ij]) ++BkgBDT[ij];
      if(BDTbiased_found[ij]) ++BkgBDTbiased[ij];  
      if(BDTnoPE_found[ij]) ++BkgBDTnoPE[ij];
      if(BDTbiasednoPE_found[ij]) ++BkgBDTbiasednoPE[ij];
    }

    for(int ij = 0; ij < 2; ++ij){
      if (Triplet[ij]) ++nTrip_matched[2*ij+1];
      if (CutMutag[ij]) ++nTrip_matched_MuTag[2*ij+1]; 
      if (CutDz[ij]) ++ndZ_cut[2*ij+1];
      if (CutMll[ij])  ++nMll_cut[2*ij+1];
      if (isKstar && CutKstar[ij])  ++nK0StarMass_cut[2*ij+1];
      //if (isKstar && CutKstarVtxCL[ij]) ++nK0StarVtxCL_cut[1];
      if (CutKpt[ij]) ++npTK_cut[2*ij+1];
      if (CutMB[ij]) ++nMB_cut[2*ij+1];
      if (CutProb[ij]) ++nVtxCLB_cut[2*ij+1];
      if (CutCos[ij]) ++nCosAlpha_cut[2*ij+1];
      if (CutPtB[ij]) ++npTB_cut[2*ij+1];
      if (CutLxy[ij]) ++nIPSB_cut[2*ij+1];
    }

    for (int i=0; i<10; i++){
      if (maxcos[0] > Cos_cuts[i]) ++BkgEffCos[i];
      if (maxPtB[0] > PtB_cuts[i]) ++BkgEffPtB[i];
      if (maxLxy[0] > Lxy_cuts[i]) ++BkgEffLxy[i];
      if (Mtri[0] > M_cuts[i]) ++BkgEffM[i];
      if (minTrk[0] > Trk_cuts[i]) ++BkgEffTrk[i];
    }
  }//loop over events



  for(int ij=0; ij<4; ++ij){
    if(ij == 1) std::cout << "\n  data K* with charge exchange " << std::endl;
    //else if(ij == 3) std::cout << "\n  data K* withOUT charge exchange " << std::endl;
    else continue;
    
    std::cout << " Total = " << nTrip_matched[ij]
              << " \n NrecoTriplets and tagMu = " << nTrip_matched_MuTag[ij]
              << " \n Ele1-Trg dz = " << ndZ_cut[ij]
              << " \n mll < 5 = " << nMll_cut[ij] << std::endl;
    if(isKstar) std::cout << " d(mK*) < 0.02 = " << nK0StarMass_cut[ij] << std::endl;
    std::cout << " kpt > 0.4 = " << npTK_cut[ij]
              << " \n Mtrip in [4.1,6] = " << nMB_cut[ij]
              << " \n B_vtxCL > 1.e-3 = " << nVtxCLB_cut[ij]
              << " \n Bcos > 0.9 = " << nCosAlpha_cut[ij]
              << " \n ptB > 5  = " << npTB_cut[ij]
              << " \n Lxy > 4 = " << nIPSB_cut[ij] << std::endl;
  }

  
  for(int i=0; i<10; i++){
    BkgEffCos[i] = BkgEffCos[i] / nEventsData;
    SigEffCos[i] = SigEffCos[i] / ngenAcceptance[0];
    
    BkgEffPtB[i] = BkgEffPtB[i] / nEventsData;
    SigEffPtB[i] = SigEffPtB[i] / ngenAcceptance[0];
    
    BkgEffLxy[i] = BkgEffLxy[i] / nEventsData;
    SigEffLxy[i] = SigEffLxy[i] / ngenAcceptance[0];
    
    BkgEffM[i] = BkgEffM[i] / nEventsData;
    SigEffM[i] = SigEffM[i] / ngenAcceptance[0];
    
    BkgEffTrk[i] = BkgEffTrk[i] / nEventsData;
    SigEffTrk[i] = SigEffTrk[i] / ngenAcceptance[0];

  }

  for(int ij=0; ij<100; ++ij){
    SigBDT_DR[ij] = SigBDT_DR[ij]/ngenAcceptance[0];
    SigBDTbiased_DR[ij] = SigBDTbiased_DR[ij]/ngenAcceptance[0];
    SigBDTnoPE_DR[ij] = SigBDTnoPE_DR[ij]/ngenAcceptance[1];
    SigBDTbiasednoPE_DR[ij] = SigBDTbiasednoPE_DR[ij]/ngenAcceptance[1];

    BkgBDT[ij] = BkgBDT[ij]/nEventsData;
    BkgBDTbiased[ij] = BkgBDTbiased[ij]/nEventsData;
    BkgBDTnoPE[ij] = BkgBDTnoPE[ij]/nEventsData;
    BkgBDTbiasednoPE[ij] = BkgBDTbiasednoPE[ij]/nEventsData;
  }

  /*
  TGraph * rocBDT = new TGraph(100, BkgBDT, SigBDT_DR);
  TCanvas *crocBDT = new TCanvas();
  crocBDT->cd();
  rocBDT->SetMarkerStyle(7);
  rocBDT->GetXaxis()->SetTitle("bkgRate (leading ele BDT)");
  rocBDT->GetYaxis()->SetTitle("sigEff gen-matched (leading ele BDT)");
  rocBDT->Draw("ap");


  TGraph * rocBDTbiased = new TGraph(100, BkgBDTbiased, SigBDTbiased_DR);
  TCanvas *crocBDTbiased = new TCanvas();
  crocBDTbiased->cd();
  rocBDTbiased->SetMarkerStyle(7);
  rocBDTbiased->GetXaxis()->SetTitle("bkgRate  (leading ele BDTbiased)");
  rocBDTbiased->GetYaxis()->SetTitle("sigEff gen-matched (leading ele BDTbiased)");
  rocBDTbiased->Draw("ap");


  TGraph * rocBDTnoPE = new TGraph(100, BkgBDTnoPE, SigBDTnoPE_DR);
  TCanvas *crocBDTnoPE = new TCanvas();
  crocBDTnoPE->cd();
  rocBDTnoPE->SetMarkerStyle(7);
  rocBDTnoPE->GetXaxis()->SetTitle("bkgRate noPE (leading ele BDT)");
  rocBDTnoPE->GetYaxis()->SetTitle("sigEff gen-matched noPE (leading ele BDT)");
  rocBDTnoPE->Draw("ap");


  TGraph * rocBDTbiasednoPE = new TGraph(100, BkgBDTbiasednoPE, SigBDTbiasednoPE_DR);
  TCanvas *crocBDTbiasednoPE = new TCanvas();
  crocBDTbiasednoPE->cd();
  rocBDTbiasednoPE->SetMarkerStyle(7);
  rocBDTbiasednoPE->GetXaxis()->SetTitle("bkgRate noPE (leading ele BDTbiased)");
  rocBDTbiasednoPE->GetYaxis()->SetTitle("sigEff gen-matched noPE (leading ele BDTbiased)");
  rocBDTbiasednoPE->Draw("ap");
  */

  /* not reviewed
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
  */

  TFile * fout= new TFile("FastPlot.root","RECREATE");
  hBptGen->Write();
  hE1ptGen->Write(); 
  hE2ptGen->Write(); 
  hKptGen->Write();
  hPptGen->Write();

  hDReltrk->Write();
  hDRrecoE->Write();
  hDRK->Write();
  hDRP->Write();
  hDRGlobal->Write();

  hBRecoPt->Write();
  hE1RecoPt->Write();
  hE2RecoPt->Write();
  hKRecoPt->Write();
  hPRecoPt->Write();
  
  hKstarMassReco->Write();
  hKstarMassGen->Write();

  hllMassReco->Write();
  hllMassGen->Write();

  return 0;
} 
