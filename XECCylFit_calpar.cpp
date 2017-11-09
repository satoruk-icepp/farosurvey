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
  Bool_t DataQual;

  tout->Branch("channel",&ChNum);
  tout->Branch("XPos",&XPos);
  tout->Branch("YPos",&YPos);
  tout->Branch("ZPos",&ZPos);
  tout->Branch("DataQual",&DataQual);

  Double_t *trans = new Double_t[3];
  TCanvas* canvas1=new TCanvas("canvas1","graph",600,600);
  TCanvas* canvas2=new TCanvas("canvas2","ZPhi",2400,1200);
  TCanvas* canvas3 = new TCanvas("canvas3","Totally interpolated",600,600);
  TGraph2D *grWF = new TGraph2D();
  TGraph2D *grTransIP =new TGraph2D();
  TGraph* grZPhi[NPart];
  TGraph *grInterpolation[NPart];

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
    grWF->SetPoint(grWF->GetN(), WFMPPCX[tmpid], WFMPPCY[tmpid], WFMPPCZ[tmpid]);
  }
  canvas2->Divide(4,2);
  for (mode = 0; mode < NPart; mode++) {
    grZPhi[mode]=new TGraph();
    grInterpolation[mode]=new TGraph();

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
          grZPhi[mode]->SetPoint(grZPhi[mode]->GetN(),ZCyl,PhiCyl);
        }
      }
    }

    TMinuit *gMinuit_Zmesh =new TMinuit(3);
    gMinuit_Zmesh->SetFCN(Zmesh);
    Int_t ierflgZ=0;
    Double_t arglistZ[10];
    arglistZ[0]=500;
    arglistZ[1]=1;
    gMinuit_Zmesh->mnparm(0, "ZOffset",-300, 0.1, 0, 0, ierflg);
    gMinuit_Zmesh->mnparm(1, "ZSpace",15.1, 0.1, 0, 0, ierflg);
    gMinuit_Zmesh->mnparm(2, "ZTilt",0, 0.1, 0, 0, ierflg);
    gMinuit_Zmesh->mnexcm("MIGRAD", arglistZ, 2, ierflgZ);

    Double_t ZOffset,ZOffsetErr,ZSpace,ZSpaceErr,ZTilt,ZTiltErr;
    gMinuit_Zmesh->GetParameter(0, ZOffset, ZOffsetErr);
    gMinuit_Zmesh->GetParameter(1, ZSpace, ZSpaceErr);
    gMinuit_Zmesh->GetParameter(2, ZTilt, ZTiltErr);

    TMinuit *gMinuit_Phimesh =new TMinuit(3);
    gMinuit_Phimesh->SetFCN(Phimesh);
    Int_t ierflgPhi=0;
    Double_t arglistPhi[10];
    arglistPhi[0]=500;
    arglistPhi[1]=1;
    gMinuit_Phimesh->mnparm(0, "PhiOffset",120 , 0.1, 0, 0, ierflgPhi);
    gMinuit_Phimesh->mnparm(1, "PhiSpace",1 , 0.1, 0, 0, ierflgPhi);
    gMinuit_Phimesh->mnparm(2, "PhiTilt",0, 0.1, 0, 0, ierflg);
    gMinuit_Phimesh->mnexcm("MIGRAD", arglistPhi, 2, ierflgPhi);

    Double_t PhiOffset,PhiOffsetErr,PhiSpace,PhiSpaceErr,PhiTilt,PhiTiltErr;
    gMinuit_Phimesh->GetParameter(0, PhiOffset, PhiOffsetErr);
    gMinuit_Phimesh->GetParameter(1, PhiSpace, PhiSpaceErr);
    gMinuit_Phimesh->GetParameter(2, PhiTilt, PhiTiltErr);

    for (Int_t i = 0; i < NMPPC; i++) {
      Int_t tmprow=i/NColumn;
      Int_t tmpcolumn=i%NColumn;
      Double_t CylZPos=tmpcolumn*ZSpace+ZOffset+tmprow*ZTilt;
      Double_t CylPhiPos=tmprow*PhiSpace+PhiOffset+tmpcolumn*PhiTilt;
      if(criteria(i,mode)==true){
        grInterpolation[mode]->SetPoint(grInterpolation[mode]->GetN(),CylZPos,CylPhiPos);
        Double_t tmpcyl[3]={R,TMath::DegToRad()*CylPhiPos,CylZPos};
        Double_t tmpcartes[3]={};
        Cyl2Cartes(tmpcartes,tmpcyl,Center,ZAxis,R);
        IPMPPCX[i]=tmpcartes[0];
        IPMPPCY[i]=tmpcartes[1];
        IPMPPCZ[i]=tmpcartes[2];
        // grTransIP->SetPoint(grTransIP->GetN(),XPos,YPos,ZPos);
      }
    }

    canvas2->cd(mode+1);
    TString strside;
    if(mode%2==0){
      strside="US";
    }else{
      strside="DS";
    }
    Int_t CFRP=TMath::Floor(mode/2);
    TString strCFRP;
    switch (CFRP) {
      case 0:
      strCFRP="CFRP: A";
      break;
      case 1:
      strCFRP="CFRP: B";
      break;
      case 2:
      strCFRP="CFRP: C";
      break;
      case 3:
      strCFRP="CFRP: D";
      break;
      default:
      strCFRP="unknown";
      break;
    }

    TString MeshTitle="Mesh Fitting";
    TString RawTitle = "Well-Fitted MPPCs";
    //TString TopTitle=MeshTitle+" "+strCFRP+" "+strside+";";
    TString TopTitle=RawTitle+" "+strCFRP+" "+strside+";";
    TString AxisTitle="Z [mm];#phi [deg]";
    TString TotalTitle=TopTitle+AxisTitle;
    grZPhi[mode]->SetTitle(TotalTitle);
    grZPhi[mode]->SetMarkerStyle(20);
    grZPhi[mode]->SetMarkerSize(1);
    grZPhi[mode]->SetMarkerColor(kRed);
    grZPhi[mode]->Draw("ap");
    grInterpolation[mode]->SetMarkerStyle(20);
    //grZPhi[mode]->SetMarkerSize(0.1);
    grInterpolation[mode]->SetMarkerColor(kBlue);
    grInterpolation[mode]->Draw("same p");

  }


  for (int i = 0; i < NMPPC; i++) {
    ChNum=i;
    XPos=IPMPPCX[i];
    YPos=IPMPPCY[i];
    ZPos=IPMPPCZ[i];
    DataQual=true;
    tout->Fill();
  }
  canvas1->cd();
  grWF->SetMinimum(-500);
  grWF->SetMaximum(500);
  grWF->GetXaxis()->SetLimits(-800,0);
  grWF->GetYaxis()->SetLimits(-800,800);
  grWF->SetMarkerStyle(20);
  grWF->SetMarkerColor(kRed);
  grWF->SetMarkerSize(0.5);
  grWF->Draw("p0");
  grTransIP->SetMarkerStyle(20);
  grTransIP->SetMarkerColor(kBlue);
  grTransIP->SetMarkerSize(0.3);
  grTransIP->Draw("same p0 ");
  canvas1->Print("./images/interpolation_with_mesh.pdf");

  canvas3->cd();
  grZPhi[4]->Draw("ap");
  //grInterpolation[4]->Draw("same p");

  tout->Write();
  fout->Close();

  return;
}
