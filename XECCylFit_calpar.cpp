#include "XECCylFit_Func.h"

/*
* 1 by Shinji
* 2 by Mitsutaka
* 3 by Mitsutaka, center and direction
* 4 by Mitsutaka, # of parameter = 6 (not 7)
* 5 by Mitsutaka, based on 3, z_dir = 1
*/

void XECCylFit_calpar(void){

  TFile* fout=new TFile("FAROMPPC_ip_mesh.root","recreate");
  TTree* tout=new TTree("tip","tip");

  Int_t ChNum;
  Double_t XPos;
  Double_t YPos;
  Double_t ZPos;

  tout->Branch("ChNum",&ChNum);
  tout->Branch("XPos",&XPos);
  tout->Branch("YPos",&YPos);
  tout->Branch("ZPos",&ZPos);

  Double_t *trans = new Double_t[3];
  TCanvas* canvas1=new TCanvas("canvas1","graph",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","ZPhi",600,600);
  TCanvas* canvas3 = new TCanvas("canvas3","Totally interpolated",600,600);
  TGraph2D *gr = new TGraph2D();
  TGraph *grZPhi =new TGraph();
  TGraph *grInterpolation =new TGraph();
  TGraph2D *grTransIP =new TGraph2D();

  trans[0] = 0;
  trans[1] = 0;
  trans[2] = 0;

  Int_t id;
  Double_t a, b, c;
  vector<Int_t> id_v;
  vector<Double_t> a_v, b_v, c_v;
  //ifstream data("cyl_data_20170311_1.dat");
  ifstream data("XECFAROSiPMOnCOBRA.dat");
  //ifstream data("mppc_cobra.dat");
  while(data>>id>>a>>b>>c){
    id_v.push_back(id);
    a_v.push_back(a);
    b_v.push_back(b);
    c_v.push_back(c);
  }
  nwf= id_v.size();
  data.close();
  for (Int_t i = 0; i < nwf; i++){
    Int_t tmpid=id_v.at(i);
    Int_t tmprow=tmpid/NColumn;
    Int_t tmpcolumn=tmpid%NColumn;

    WFUsedMPPC[tmpid] = true;
    WFMPPCX[tmpid] = a_v.at(i);
    WFMPPCY[tmpid] = b_v.at(i);
    WFMPPCZ[tmpid] = c_v.at(i);
    //gr->SetPoint(gr->GetN(), WFMPPCX[tmpid], WFMPPCY[tmpid], WFMPPCZ[tmpid]);
  }

  //canvas1->cd();
  //gr->Draw("p0");


  for (mode = 0; mode < 8; mode++) {
    TMinuit* gMinuit= new TMinuit(7);
    gMinuit->SetFCN(fcn);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0]=1;

    Double_t initial[7]={r_set, trans[0], trans[1], trans[2], 0, 0, 1};
    Double_t vstart[7]={r_set, trans[0], trans[1], trans[2], 0, 0, 1};
    Double_t step[7]={1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

    gMinuit->mnparm(0, "R", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "centerX", vstart[1], step[1], 0, 0, ierflg);
    gMinuit->mnparm(2, "centerY", vstart[2], step[2], 0, 0, ierflg);
    gMinuit->mnparm(3, "centerZ", vstart[3], step[3], 0, 0, ierflg);
    gMinuit->mnparm(4, "directionX", vstart[4], step[4], -1, 1, ierflg);
    gMinuit->mnparm(5, "directionY", vstart[5], step[5], -1, 1, ierflg);
    gMinuit->mnparm(6, "directionZ", vstart[6], step[6], -1, 1, ierflg);

    gMinuit->FixParameter(3);
    arglist[0]=500;
    arglist[1]=1.;

    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    Double_t *Par1st = new Double_t[7];
    Double_t *Err1st = new Double_t[7];
    Double_t temp_par, temp_err, chi2, ndf;

    for (Int_t i = 0; i < 7; i++){
      gMinuit->GetParameter(i, temp_par, temp_err);
      Par1st[i] = temp_par;
      Err1st[i] = temp_err;
      cout << initial[i] << " " << temp_par << "+-" << temp_err << ":" << temp_par - initial[i] << endl;
    }

    Double_t fedm,errdef,fmin;
    Int_t npari,nparx,istat;
    //gMinuit -> mnstat(fmin,fedm,errdef,npari,nparx,istat);

    Double_t R=Par1st[0];
    Double_t Center[3]={Par1st[1],Par1st[2],Par1st[3]};
    TVector3 ZAxis(Par1st[4],Par1st[5],Par1st[6]);
    //ZAxis.Print();
    for (Int_t i = 0; i < nwf; i++){
      Int_t tmpid=id_v.at(i);
      if(criteria(tmpid,mode)==true){
        if (WFUsedMPPC[tmpid]==true) {
          Double_t coord[3]={WFMPPCX[tmpid],WFMPPCY[tmpid],WFMPPCZ[tmpid]};
          Double_t Cyl[3];
          Cartes2Cyl(coord,Cyl,Center,ZAxis,R);
          Double_t PhiCyl=TMath::RadToDeg()*Cyl[1];
          Double_t ZCyl=Cyl[2];
          ZAllMPPC[tmpid]=ZCyl;
          PhiAllMPPC[tmpid]=PhiCyl;
          //std::cout<<"Theta: "<<PhiCyl<<" Z: "<<ZCyl<<std::endl;
          //grZPhi->SetPoint(grZPhi->GetN(),ZCyl,PhiCyl);
        }
      }
    }

    TMinuit *gMinuit_Zmesh =new TMinuit(2);
    gMinuit_Zmesh->SetFCN(Zmesh);
    Int_t ierflgZ=0;
    Double_t arglistZ[10];
    arglistZ[0]=500;
    arglistZ[1]=1;
    gMinuit_Zmesh->mnparm(0, "ZOffset",-300, 0.1, 0, 0, ierflg);
    gMinuit_Zmesh->mnparm(1, "ZSpace",15.1, 0.1, 0, 0, ierflg);
    gMinuit_Zmesh->mnexcm("MIGRAD", arglistZ, 2, ierflgZ);

    Double_t ZOffset,ZOffsetErr,ZSpace,ZSpaceErr;
    gMinuit_Zmesh->GetParameter(0, ZOffset, ZOffsetErr);
    gMinuit_Zmesh->GetParameter(1, ZSpace, ZSpaceErr);

    TMinuit *gMinuit_Phimesh =new TMinuit(2);
    gMinuit_Phimesh->SetFCN(Phimesh);
    Int_t ierflgPhi=0;
    Double_t arglistPhi[10];
    arglistPhi[0]=500;
    arglistPhi[1]=1;
    gMinuit_Phimesh->mnparm(0, "PhiOffset",120 , 0.1, 0, 0, ierflgPhi);
    gMinuit_Phimesh->mnparm(1, "PhiSpace",1 , 0.1, 0, 0, ierflgPhi);
    gMinuit_Phimesh->mnexcm("MIGRAD", arglistPhi, 2, ierflgPhi);

    Double_t PhiOffset,PhiOffsetErr,PhiSpace,PhiSpaceErr;
    gMinuit_Phimesh->GetParameter(0, PhiOffset, PhiOffsetErr);
    gMinuit_Phimesh->GetParameter(1, PhiSpace, PhiSpaceErr);

    for (Int_t i = 0; i < NMPPC; i++) {
      Int_t tmprow=i/NColumn;
      Int_t tmpcolumn=i%NColumn;
      Double_t CylZPos=tmpcolumn*ZSpace+ZOffset;
      Double_t CylPhiPos=tmprow*PhiSpace+PhiOffset;
      if(criteria(i,mode)==true){
        //grInterpolation->SetPoint(grInterpolation->GetN(),CylZPos,CylPhiPos);
        Double_t tmpcyl[3]={R,TMath::DegToRad()*CylPhiPos,CylZPos};
        Double_t tmpcartes[3]={};
        Cyl2Cartes(tmpcartes,tmpcyl,Center,ZAxis,R);
        ChNum=i;
        XPos=tmpcartes[0];
        YPos=tmpcartes[1];
        ZPos=tmpcartes[2];
        tout->Fill();
        //grTransIP->SetPoint(grTransIP->GetN(),XPos,YPos,ZPos);
      }
    }
  }


  //canvas2->cd();
  grZPhi->SetMarkerStyle(20);
  grZPhi->SetMarkerSize(2);
  grZPhi->SetMarkerColor(kRed);
  //grZPhi->Draw("ap");
  grInterpolation->SetMarkerStyle(20);
  grInterpolation->SetMarkerColor(kBlue);
  grInterpolation->Draw("same p");

  //canvas3->cd();
  grTransIP->SetMarkerStyle(20);
  grTransIP->SetMarkerColor(kRed);
  grTransIP->Draw("p0");
  //gr->Draw("same p0");

  tout->Write();
  fout->Close();

  return;
}
