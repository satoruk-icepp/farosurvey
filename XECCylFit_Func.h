#define NColumn 44
#define NRow 93
#define NMPPC 4092
#define NCFRP 4
#define NPart 8
Bool_t MC = false;
Int_t nmax = 600;
Double_t r_set = 650;
Double_t z_set = 0;
Double_t err = 0.15;
Double_t ZAllMPPC[NMPPC];
Double_t PhiAllMPPC[NMPPC];
Int_t CFRPOrigin[NCFRP+1]={0,24,47,70,93};
Int_t mode;
Int_t nwf;

Double_t WFMPPCX[NMPPC];
Double_t WFMPPCY[NMPPC];
Double_t WFMPPCZ[NMPPC];
Bool_t WFUsedMPPC[NMPPC];


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
//void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
Double_t fitfunc(Double_t x, Double_t y, Double_t z, Double_t *par, Int_t iter);
Bool_t criteria(Int_t iCh,Int_t i);
Double_t Zdevmesh(Int_t channel, Double_t *par, Int_t iter);
Double_t Phidevmesh(Int_t channel,Double_t *par, Int_t iter);
void Zmesh(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void Phimesh(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void Cartes2Cyl(Double_t *Cartes,Double_t *Cyl,Double_t *Center,TVector3 ZAxis,Double_t R);
void Cyl2Cartes(Double_t *Cartes,Double_t *Cyl,Double_t *Center,TVector3 ZAxis,Double_t R);

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  // Int_t npnt=12;
  //Int_t npnt = nwf;
  Int_t iter = 1;
  Double_t chisq=0;
  Double_t delta;
  for(int i=0; i<NMPPC;i++){
    if(criteria(i,mode)==true){
      if(WFUsedMPPC[i]==true){
        delta=fitfunc(WFMPPCX[i], WFMPPCY[i], WFMPPCZ[i], par, iter);
        chisq += delta* delta;
      }
    }
  }
  f = chisq;
}
/*
void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
// Int_t npnt=12;
Int_t npnt = nwf;
Int_t iter = 2;
Double_t chisq=0;
Double_t delta;
for(int i=0; i<npnt;i++){
if(WFUsedMPPC[i]==true){
delta=fitfunc(WFMPPCX[i], WFMPPCY[i], WFMPPCZ[i], par, iter);
chisq += delta* delta;
}
}
f = chisq;
}
*/

Double_t fitfunc(Double_t x, Double_t y, Double_t z, Double_t *par, Int_t iter){

  Double_t AxisL = TMath::Sqrt(par[4] * par[4] + par[5] * par[5] + par[6] * par[6]);
  Double_t AxisDir[3] = {par[4]/AxisL, par[5]/AxisL, par[6]/AxisL};
  Double_t Center[3] = {par[1], par[2], par[3]};
  Double_t A[3] = {x - Center[0], y - Center[1], z - Center[2]};

  Double_t k = A[0]*AxisDir[0] + A[1]*AxisDir[1] +A[2]*AxisDir[2];
  // Double_t Foot[3] = {Center[0] + k * AxisDir[0], Center[1] + k * AxisDir[1], Center[2] + k * AxidDir[2]};
  Double_t l = TMath::Sqrt((A[0] - k*AxisDir[0]) *(A[0] - k*AxisDir[0]) + (A[1] - k*AxisDir[1]) *(A[1] - k*AxisDir[1]) + (A[2] - k*AxisDir[2]) *(A[2] - k*AxisDir[2]));
  // cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] << " " << par[5] << " " << par[6] << " " <<TMath::Sqrt((l-par[0])*(l-par[0])) << endl;
  Double_t value;
  switch (iter){
    case 1:
    // return TMath::Sqrt((l-par[0])*(l-par[0]))/0.2;
    value= (l-par[0])*(l-par[0])/(err * err);
    break;
    case 2:
    value= k*k +(l-par[0])*(l-par[0])/(err * err);
    //return TMath::Sqrt(k*k +(l-par[0])*(l-par[0]))/0.2;
    //return TMath::Sqrt(k*k);
    break;
  }
  return value;
}

Bool_t criteria(Int_t iCh,Int_t part){
  Int_t Row=iCh/NColumn;
  Int_t Column=iCh%NColumn;
  Bool_t cri_col=false;
  Bool_t cri_row=false;
  Bool_t cri_total=false;
  if(part%2==0){
    if(Column<22){
      cri_col=true;
    }
  }else{
    if(Column>=22){
      cri_col=true;
    }
  }
  Int_t CFRP=part/2;
  if(Row>=CFRPOrigin[CFRP]&&Row<CFRPOrigin[CFRP+1]){
    cri_row=true;
  }

  if(cri_col==true&&cri_row==true){
    cri_total=true;
    //std::cout<<"selected: "<<iCh<<std::endl;
  }

  return cri_total;
}

Double_t Zdevmesh(Int_t channel, Double_t *par, Int_t iter){
  Double_t Column=(double)(channel%NColumn);
  Double_t ZOffset=par[0];
  Double_t ZSpace=par[1];
  Double_t ZMesh=Column*ZSpace+ZOffset;
  Double_t deviation=ZMesh-ZAllMPPC[channel];
  //std::cout<<"channel: "<<channel<<" Column: "<<Column<<" Faro: "<<ZAllMPPC[channel]<<" mesh: "<<ZMesh<<std::endl;
  return deviation;

}

Double_t Phidevmesh(Int_t channel,Double_t *par, Int_t iter){
  Int_t Row = channel/NColumn;
  Double_t PhiOffset=par[0];
  Double_t PhiSpace=par[1];
  Double_t PhiMesh=(double)Row*PhiSpace+PhiOffset;
  Double_t deviation=PhiMesh-PhiAllMPPC[channel];
  return deviation;
}

void Zmesh(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  // Int_t npnt=12;
  Int_t iter = 1;
  Double_t chisq=0;
  Double_t delta;
  for(int i=0; i<NMPPC;i++){
    if(criteria(i,mode)==true){
      if(WFUsedMPPC[i]==true){
        delta=Zdevmesh(i, par, iter);
        chisq += delta* delta;
      }
    }
  }
  f = chisq;
}

void Phimesh(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  // Int_t npnt=12;
  Int_t iter = 1;
  Double_t chisq=0;
  Double_t delta;

  for(int i=0; i<NMPPC;i++){
    if(criteria(i,mode)==true){
      if(WFUsedMPPC[i]==true){
        //std::cout<<"channel: "<<MPPCIndex[i]<<std::endl;
        delta=Phidevmesh(i, par, iter);
        chisq += delta* delta;
      }
    }
  }
  f = chisq;
}

void Cartes2Cyl(Double_t *Cartes,Double_t *Cyl,Double_t *Center,TVector3 ZAxis,Double_t R){
  //par[0]:centerx
  //par[1]:centerY
  //ZAxis.X():vectorX
  //ZAxis.Y():vectorY
  //ZAxis.Z():vectorZ

  TVector3 XAxis(ZAxis.Z(),0,-ZAxis.X());
  TVector3 YAxis(-ZAxis.X(),(ZAxis.X()*ZAxis.X()+ZAxis.Z()*ZAxis.Z())/ZAxis.Y(),-ZAxis.Z());
  TVector3 uZAxis,uXAxis, uYAxis;
  uZAxis=ZAxis.Unit();
  uXAxis=XAxis.Unit();
  uYAxis=YAxis.Unit();


  TVectorD *RelativeCartes=new TVectorD(3);
  for (int i = 0; i < 3; i++) {
    (*RelativeCartes)[i]=(Double_t) (Cartes[i]-Center[i]);
  }
  //ZAxis.Print();
  TMatrixD M(3,3);
  M[0][0]=uXAxis.X();
  M[0][1]=uYAxis.X();
  M[0][2]=uZAxis.X();
  M[1][0]=uXAxis.Y();
  M[1][1]=uYAxis.Y();
  M[1][2]=uZAxis.Y();
  M[2][0]=uXAxis.Z();
  M[2][1]=uYAxis.Z();
  M[2][2]=uZAxis.Z();

  TDecompSVD* SVD=new TDecompSVD(M);
  TMatrixD U_T =SVD->GetU();
  U_T.T();
  TMatrixD V =SVD->GetV();
  TVectorD S =SVD->GetSig();
  TMatrixD S_inv(3,3);
  for(int i = 0; i < 3; i++) {
    S_inv[i][i]=1/S[i];
  }
  Int_t Layers=3;
  //multiply V S足1
  TMatrixD VS_inv(Layers,Layers);
  VS_inv.Mult(V,S_inv);
  //calculate M足1
  TMatrixD M_inv(Layers,Layers);
  M_inv.Mult(VS_inv,U_T);
  //M_inv.Print();
  //M.Print();
  //RelativeCartes->Print();
  *RelativeCartes *= M_inv;
  Cyl[0] = R;
  Cyl[1] = TMath::ATan2((*RelativeCartes)[1],(*RelativeCartes)[0]);
  Cyl[2] = (*RelativeCartes)[2];
  //RelativeCartes->Print();
  //std::cout<<"it works!!"<<std::endl;
}

void Cyl2Cartes(Double_t *Cartes,Double_t *Cyl,Double_t *Center,TVector3 ZAxis,Double_t R){
  TVector3 XAxis(ZAxis.Z(),0,-ZAxis.X());
  TVector3 YAxis(-ZAxis.X(),(ZAxis.X()*ZAxis.X()+ZAxis.Z()*ZAxis.Z())/ZAxis.Y(),-ZAxis.Z());
  TVector3 uZAxis,uXAxis, uYAxis;
  uZAxis=ZAxis.Unit();
  uXAxis=XAxis.Unit();
  uYAxis=YAxis.Unit();

  //ZAxis.Print();
  TMatrixD M(3,3);
  M[0][0]=uXAxis.X();
  M[0][1]=uYAxis.X();
  M[0][2]=uZAxis.X();
  M[1][0]=uXAxis.Y();
  M[1][1]=uYAxis.Y();
  M[1][2]=uZAxis.Y();
  M[2][0]=uXAxis.Z();
  M[2][1]=uYAxis.Z();
  M[2][2]=uZAxis.Z();

  TVectorD *RelativeCyl=new TVectorD(3);
  (*RelativeCyl)[0]=R*TMath::Cos(Cyl[1]);
  (*RelativeCyl)[1]=R*TMath::Sin(Cyl[1]);
  (*RelativeCyl)[2]=Cyl[2];

  *RelativeCyl *= M;
  for (int i = 0; i < 3; i++) {
    Cartes[i]=(*RelativeCyl)[i]+Center[i];
  }
}
