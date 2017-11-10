#include "XECCylFit_Func.h"

Double_t thcoSUS=16*pow(10,-6);
Double_t dT=130;

void ShrinkVector(Double_t *shrvec,Double_t *RawPos, Double_t *SOPos);

void Faroshrink(){
  TFile* fin = new TFile("FAROMPPC_ip_mesh.root","read");
  TTree* tin = (TTree*)fin->Get("tip");
  TFile* fout = new TFile("FAROMPPC_ip_mesh_shrink.root","recreate");
  TTree* tout = new TTree("tip","tip");
  Double_t RawPos[3];

  tin->SetBranchAddress("XPos",&RawPos[0]);
  tin->SetBranchAddress("YPos",&RawPos[1]);
  tin->SetBranchAddress("ZPos",&RawPos[2]);


  Int_t channel;
  Double_t ShrinkPos[3];
  //Bool_t DataQual;
  tout->Branch("XPos",&ShrinkPos[0]);
  tout->Branch("YPos",&ShrinkPos[1]);
  tout->Branch("ZPos",&ShrinkPos[2]);
  tout->Branch("channel",&channel);
  //tout->Branch("DataQual",&DataQual);

  Double_t ShrinkOrigin[3]={0,-1100,0};

  for (size_t iCh = 0; iCh < NMPPC; iCh++) {
    tin->GetEntry(iCh);
    //std::cout<<"Raw X:"<<RawPos[0]<<std::endl;
    Double_t shrvec[3];
    ShrinkVector(shrvec,RawPos,ShrinkOrigin);
    //std::cout<<"shrink x:"<<shrvec[0]<<"shrink y:"<<shrvec[1]<<"shrink z:"<<shrvec[2]<<std::endl;
    for (int i = 0; i < 3; i++) {
      ShrinkPos[i]=RawPos[i]+shrvec[i];
    }
    //DataQual=true;
    channel=iCh;
    tout->Fill();
  }
  tout->Write();
  fout->Close();
}

void ShrinkVector(Double_t *shrvec,Double_t *RawPos, Double_t *SOPos){
  for (int i = 0; i < 3; i++) {
    shrvec[i]=thcoSUS*dT*(SOPos[i]-RawPos[i]);
  }
}


// Double_t Distance(Double_t *RawPos, Double_t *SOPos){
//   Double_t sqdist;
//   for (int i = 0; i < 3; i++) {
//     sqdist+=pow(RawPos[i]-SOPos[i],2);
//   }
//   return TMath::Sqrt(sqdist);
// }
