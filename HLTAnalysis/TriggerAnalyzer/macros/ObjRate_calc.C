#include "skim_class.h"
#include "TChain.h"
#include "TH1.h"
#include "match_jpsi.h"
#include "plot_jpsi.h"

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
/*
float signal[5](float x[5], std::vector<float> * var, bool BestMin){
  float y[5]={0};
  float N=0;
  for(int j =0; j<var->size(); j++){ 
    N++;
    for(int i=0; i<5; i++){
      if (fabs(var->at(j))>x[i] && !BestMin) y[i]=y[i]+1;
      if (fabs(var->at(j))<x[i] && BestMin) y[i]=y[i]+1;
     }
  }
  for(int i=0; i<5; i++)y[i]= 

*/
int ObjRate_calc(){

TChain *ccdata=new TChain("demo/mytree");


 ccdata->Add("outputdata_test.root");

skim_class data;
// skimObj_class data;
data.Init(ccdata);

cout<<"data "<<ccdata->GetEntries()<<endl;

 int nbkg=0,npass=0;
 for (int ev=0; ev<ccdata->GetEntries(); ev++){
 data.GetEntry(ev);
 nbkg++;
 if (ev%1000==0) cout<<ev<<endl;
 //cuts
 int itrg=data.SelectedMu_index;
 if (data.SelectedMu_index<0) continue;
 float ZmuVtx=data.muon_vz->at(data.SelectedMu_index);

 int temp_pass=0;
 for(int ivtx=0; ivtx<data.TTrack_cos->size(); ivtx++){
   //   if (data.TTrack_ObjId->at(ivtx)!=11) continue;
       temp_pass++;
           }
 if (temp_pass>0) npass++; 
 }

 cout<<npass<<"  "<<nbkg<<endl;


return 0;
} 
