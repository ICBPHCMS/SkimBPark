#include "skim_class.h"
#include "TChain.h"
#include "TH1.h"
#include "match_jpsi.h"
#include "plot_jpsi.h"
#include "skimhelper.h"
#include "skim_newe_class.h"

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

int ObjSingleE_Eff(bool write=false,TString name="efftest.root"){
TChain *ccmc=new TChain("demo/mytree");
//10X
// ccmc->Add("output10XMC_oldGsf.root");
// 94X new gsf
ccmc->Add("/eos/cms/store/group/cmst3/user/gkaratha/BKee_SkimmingNtuples_lowptE_10X_ID_94XMC/CRAB_UserFiles/crab_BKee_SkimmingNtuples_lowptE_10X_ID_94XMC/190127_121737/0000/*.root");
//94X old gsf
//ccmc->Add("/eos/cms/store/group/cmst3/user/gkaratha/BKee_SkimmingNtuples_lowptE_10X_ID_94XMC_normalGsfTracks/CRAB_UserFiles/crab_BKee_SkimmingNtuples_lowptE_10X_ID_94XMC_normalGsfTracks/190127_154944/0000/*.root");

 int nentries=ccmc->GetEntries();

 nentries=35000;
 bool useMVA=true;
 bool useBiased=false;
  bool probe=true;
 float probeCone=0.4;
 float probeDz=0.3;
 float PtCut=0;
skim_class mc;
mc.Init(ccmc,"0");

//unbiased
 float rangeBDT[3]={1.03,1.75,2.61};
//match gen
TH1F * hDenB=new TH1F("hDenB"," ",30,0,10);
TH1F * hSigBL=new TH1F("hSigBL"," ",30,0,10);
TH1F * hSigBM=new TH1F("hSigBM"," ",30,0,10);
TH1F * hSigBT=new TH1F("hSigBT"," ",30,0,10);


/*
TH1F * hBkgB=new TH1F("hBkggB"," ",50,-4,6);
TH1F * hBkgUnB=new TH1F("hBkgUnB"," ",50,-4,6);*/
TH1F * hDR=new TH1F("hDR"," ",100,0,0.1);


//signal
 cout<<"total "<<ccmc->GetEntries()<<" run on "<<nentries<<endl;
 float matchedL=0,matchedM=0,matchedT=0,total=0;
 

//biased
 if(useBiased){rangeBDT[0]=-0.48; rangeBDT[1]=0.76; rangeBDT[2]=1.83;}
 float effBDT[3]={0};
 for (int ev=0; ev<nentries; ev++){
  mc.GetEntry(ev);
  if (ev%500==0) cout<<ev<<endl;
  int MuIndex=mc.SelectedMu_index;
  if (probe && MuIndex==-1) continue;
  for(int iel=0; iel<mc.genMu_Id->size(); iel++){
    if (fabs(mc.genMu_Id->at(iel))==13) continue;
    if (mc.genMu_pt->at(iel)<PtCut) continue;
    if (fabs(mc.genMu_eta->at(iel))>2.4) continue;
    if (probe){
          if (DR(mc.muon_eta->at(MuIndex),mc.muon_phi->at(MuIndex),mc.genMu_eta->at(iel),mc.genMu_phi->at(iel))<probeCone) continue;
       }
    hDenB->Fill(mc.genMu_pt->at(iel)); total++;
    float minDR=1000;  int indexE=0;
    for (int irel=0; irel<mc.nelectron; irel++){
       if (mc.el_charge->at(irel)!=mc.genMu_ch->at(iel)) continue;
       if (minDR<DR(mc.genMu_eta->at(iel),mc.genMu_phi->at(iel),mc.el_eta->at(irel),mc.el_phi->at(irel))) continue;
        minDR=DR(mc.genMu_eta->at(iel),mc.genMu_phi->at(iel),mc.el_eta->at(irel),mc.el_phi->at(irel));
  }
    hDR->Fill(minDR);
   if (minDR<0.03){
       if (useMVA){
          if ( useBiased && mc.el_mvaB_wp->at(indexE)>-0.48) {
                 hSigBL->Fill(mc.genMu_pt->at(iel)); matchedL++; }
          if ( useBiased && mc.el_mvaB_wp->at(indexE)>0.76) {
                 hSigBM->Fill(mc.genMu_pt->at(iel)); matchedM++; }
          if ( useBiased && mc.el_mvaB_wp->at(indexE)>1.83) {
                 hSigBT->Fill(mc.genMu_pt->at(iel)); matchedT++; }
          if ( !useBiased && mc.el_mvaUnB_wp->at(indexE)>1.03) {
                 hSigBL->Fill(mc.genMu_pt->at(iel)); matchedL++; }
          if ( !useBiased && mc.el_mvaUnB_wp->at(indexE)>1.75) {
                 hSigBM->Fill(mc.genMu_pt->at(iel)); matchedM++; }
          if ( !useBiased && mc.el_mvaUnB_wp->at(indexE)>2.61) {
                 hSigBT->Fill(mc.genMu_pt->at(iel)); matchedT++; }
       }
      else{ 
       hSigBL->Fill(mc.genMu_pt->at(iel)); matchedL++;}
          
  }
}
}
if (!useMVA) 
   cout<<"total "<<total<<" eff "<<matchedL/total<<endl;
else
  cout<<"total  "<<total<<" eff "<<matchedL/total<<"  "<<matchedM/total<<"  "<<matchedT/total<<endl;
// cout<<"total  "<<total<<" eff "<<matchedL<<"  "<<matchedM<<"  "<<matchedT<<endl;
TCanvas * cDR=canvas_1plot(hDR,"bDR",true,"#Delta R","");


hSigBL->Divide(hDenB); hSigBM->Divide(hDenB); hSigBT->Divide(hDenB);
if (!useMVA) 
   TCanvas * cEff=canvas_1plot(hSigBL,"bDR",false,"p_{T}","");
else
   TCanvas * cBDT=canvas_3plot(hSigBL,hSigBM,hSigBT,"bBDT",false,false,-1,-1,"p_{T} (gen e)","","Loose","Medium","Tight");


 if (write){
   TFile * fout=new TFile(name,"RECREATE");
   hDR->Write();
   hSigBL->Write();
   if (useMVA){
   hSigBM->Write(); hSigBT->Write(); }
 }
return 0;
} 
