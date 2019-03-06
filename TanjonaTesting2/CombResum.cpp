#include <numeric>
#include "CombResum.h"



long double ConstResum::MUF=125.09;
long double ConstResum::MUR=125.09;
long double ConstResum::as=0.118;
long double ConstResum::sigma0Higgs=1.;
long double ConstResum::Q=125.09;
long double ConstResum::Gf = 0.0000116637;// Fermi constant in GeV^-2
long double ConstResum::GeVtopb = 389379660.;// GeV^-2 to pb conversion factor == (hc)^2 

std::complex<long double> Tmatch(std::complex<long double> N, long double xp){
  return (std::pow(N,3)*std::pow(xp,2)/(1.+std::pow(N,3)*std::pow(xp,2)));
}

std::complex<long double> T1(std::complex<long double> N, long double xp){
  return (1.);
}

std::complex<long double> T0(std::complex<long double> N, long double xp){
  return (0.);
}


CombResum::CombResum(int ordres, int ordmatch, int channel, string PDFNAME, bool Wilson,
				long double mH, long double Nc, long double Nf, long double Mur, long double Muf)
:_Lumi(LHAPDF::mkPDF(PDFNAME),Muf,Nf){
   _Nc=Nc;_Nf=Nf;
 //Initialise constants
  _Ca=_Nc; 
  _Cf=(_Nc*_Nc-1.)/(2.*_Nc);
  ConstResum::MUF=Muf;
  ConstResum::MUR=Mur;
  ConstResum::Q=mH;
 
   
   //REMEMBER!!! ALL in unity as
  beta_0=(11.*_Ca-2.*_Nf)/(12.*M_PIl);
  beta_1=((17.*_Ca*_Ca-5.*_Ca*_Nf-3.*_Cf*_Nf)*2./3.)/(16.*M_PIl*M_PIl);
  beta_2=((2857./54.*_Ca*_Ca*_Ca+(_Cf*_Cf-205./18.*_Cf*_Ca-1415./54.*_Ca*_Ca)*_Nf
  +(11./9.*_Cf+79./54.*_Ca)*_Nf*_Nf))/std::pow(4.*M_PIl,3);
  ConstResum::as=_Lumi.get_alphaS(Mur);
  cout << "as = " << ConstResum::as << endl;
  _ordres=ordres;_ordmatch=ordmatch;_channel=channel;_Wilson=Wilson;_Nc=Nc;_Nf=Nf;
  _Fix=new FixedptResum(ordres, ordmatch,channel, Wilson, Nc, Nf);
  _Joint=new Joint(ordres, ordmatch, channel, Wilson, Nc, Nf);
   ConstResum::sigma0Higgs=std::sqrt(2.)*ConstResum::Gf*ConstResum::as*ConstResum::as/(576.*M_PIl);
  
  
  
  
  const long double zeta2=gsl_sf_zeta_int(2);
  const long double zeta3=gsl_sf_zeta_int(3);
  const long double zeta4=gsl_sf_zeta_int(4);
  
  
  
  
  
  _MatchFun=Tmatch;
  
 
  
  
  
}

CombResum::~CombResum(){
  delete _Fix;
  delete _Joint;
}


void CombResum::SetMatchingFunction(std::complex<long double> (Func)(std::complex<long double>, long double)){
  _MatchFun=Func;
}

void CombResum::SetMUF(long double Muf){
  ConstResum::MUF=Muf;
  _Lumi.Cheb_Lum(ConstResum::MUF);
}

void CombResum::SetMUR(long double Mur){
  ConstResum::MUR=Mur;
  ConstResum::as=_Lumi.get_alphaS(ConstResum::MUR);
  ConstResum::sigma0Higgs=std::sqrt(2.)*ConstResum::Gf*ConstResum::as*ConstResum::as/(576.*M_PIl);
}

void CombResum::SetChannel(int channel){
  _channel=channel;
  _Fix->SetChannel(channel);
  _Joint->SetChannel(channel);
}

void CombResum::SetOrdRes(int ordres){
  _ordres=ordres;
  _Fix->SetOrdRes(ordres);
  _Joint->SetOrdRes(ordres);
}

void CombResum::SetOrdMatch(int ordmatch){
  _ordmatch=ordmatch;
  _Fix->SetOrdMatch(ordmatch);
  _Joint->SetOrdMatch(ordmatch);
}


long double CombResum::alpha_s_muR(long double mu)
{
  long double X=1.+alpha_Mz*beta_0*std::log(std::pow(mu,2)/std::pow(Mz,2));
  long double ris=0.0;
  ris= alpha_Mz/X-beta_1/beta_0*std::pow(alpha_Mz/X,2.)*std::log(X)+std::pow(alpha_Mz/X,3)*(beta_2/beta_0*(1.-X)
  +std::pow(beta_1/beta_0,2)*(std::pow(std::log(X),2.)-std::log(X)-1.+X));
  return ris;
}

struct ptstruct {
  CombResum *CombR; 
  long double xp;
};

std::complex<long double> ResCrossSecFixpt(std::complex<long double> N, void *p){
  ptstruct par=*(ptstruct *)p;
  return(par.CombR->FixptPartRes(N,par.xp));
}

std::complex<long double> ResCrossSecJoint(std::complex<long double> N,std::complex<long double> lchi, void *p){
  ptstruct par=*(ptstruct *)p;
  return(par.CombR->JointPartRes(N,lchi,par.xp));
}

std::complex<long double> MatchCrossSecFixpt(std::complex<long double> N, void *p){
  ptstruct par=*(ptstruct *)p;
  return(par.CombR->FixptPartMatch(N,par.xp));
}

std::complex<long double> MatchCrossSecJoint(std::complex<long double> N, void *p){
  ptstruct par=*(ptstruct *)p;
  return(par.CombR->JointPartMatch(N,par.xp));
}

long double CombResum::ResummedCrossSection(long double CMS, long double xp, int matching, long double *err){
  //Matching Function
  //matching==0 only Joint
  //matching==1 only FixedptResum
  //matching!=0 || matching !=1 used matching function defined in CombResum
  std::function<std::complex<long double>(std::complex<long double>, long double)> _MatchFun2=_MatchFun;
  if (matching==0) _MatchFun=T0;
  if (matching==1) _MatchFun=T1;
  //Fixed pt Part
  long double FixptPart=0.,errFixptPart=0.;
  long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/CMS,2.);
  ptstruct Fixptpar;
  Fixptpar.xp=xp;
  Fixptpar.CombR=this;
  FixptPart=integration::InverseMellin_path(3,ResCrossSecFixpt,tauprime,2.,1.5,&Fixptpar, &errFixptPart);
  FixptPart*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  //Joint Part
  long double Joint=0.,errJoint=0.;
  if (xp < std::pow(3.5/ConstResum::Q,2.)){
    cout << " low pt " << endl;
    Borel::BorelJointC(3,tauprime,xp,ConstResum::as*beta_0,ResCrossSecJoint,&Fixptpar,&Joint,&errJoint,2.,3.,1.0); //Good convergence at small pt not at high pt
  }
  else{
    Borel::BorelJointCp(3,tauprime,xp,ConstResum::as*beta_0,ResCrossSecJoint,&Fixptpar,&Joint,&errJoint,2.,3.,1.0);//Quick and good convergence for  pt>4... C=1,2,3... no difference in this region
  }
  Joint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  errJoint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  //cout << "Joint= " << Joint << " +- " << errJoint << endl;
  _MatchFun=_MatchFun2;
  *err=errJoint+errFixptPart;
  return (FixptPart+Joint);
  
  
}

std::vector<long double> CombResum::ResummedCrossSection(long double CMS, std::vector<long double> xp,int matching, std::vector<long double> *err){
  //Matching Function
  //matching==0 only Joint
  //matching==1 only FixedptResum
  //matching!=0 || matching !=1 used matching function defined in CombResum
  std::function<std::complex<long double>(std::complex<long double>, long double)> _MatchFun2=_MatchFun; 
  if (matching==0) _MatchFun=T0;
  if (matching==1) _MatchFun=T1;
  
  std::vector<long double> ris,error;
  int i=0;
  for (auto pt : xp){
    std::cout << i << endl ;
    //Fixed pt Part
    long double FixptPart=0., errFixptPart=0.;
    long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+pt)+std::sqrt(pt))/CMS,2.);
    ptstruct Fixptpar;
    Fixptpar.xp=pt;
    Fixptpar.CombR=this;
    FixptPart=integration::InverseMellin_path(3,ResCrossSecFixpt,tauprime,2.,1.5,&Fixptpar,&errFixptPart);
    long double Joint=0.,errJoint=0.;
    if (pt < std::pow(3.5/ConstResum::Q,2.)){
      cout << " low pt " << endl;
      Borel::BorelJointC(3,tauprime,pt,ConstResum::as*beta_0,ResCrossSecJoint,&Fixptpar,&Joint,&errJoint,2.,3.,1.0); //Good convergence at small pt not at high pt
    }
    else{
      Borel::BorelJointCp(3,tauprime,pt,ConstResum::as*beta_0,ResCrossSecJoint,&Fixptpar,&Joint,&errJoint,2.,3.,1.0);//Quick and good convergence for  pt>4... C=1,2,3... no difference in this region
    }
    ris.push_back((FixptPart+Joint)*tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb);
    error.push_back((errFixptPart+errJoint)*tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb);
    
    i++;
  }
  _MatchFun=_MatchFun2;
  *err=error;
  return ris;
}
std::vector<std::vector<long double>> CombResum::ResummedCrossSection(std::vector<long double> CMS, std::vector<long double> xp,int matching,
  std::vector<std::vector<long double>> *err){
  //Matching Function
  //matching==0 only Joint
  //matching==1 only FixedptResum
  //matching!=0 || matching !=1 used matching function defined in CombResum
  std::function<std::complex<long double> (std::complex<long double>, long double)> _MatchFun2=_MatchFun;
  if (matching==0) _MatchFun=T0;
  if (matching==1) _MatchFun=T1;
  
  std::vector<std::vector<long double>> ris, error;
  int i=0;
  for (auto cms : CMS ){
    cout << i << " " ;
    std::vector<long double> ris2,error2;
    for (auto pt : xp){
      //Fixed pt Part
      long double FixptPart=0.,errFixptPart=0.;
      long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+pt)+std::sqrt(pt))/cms,2.);
      ptstruct Fixptpar;
      Fixptpar.xp=pt;
      Fixptpar.CombR=this;
      FixptPart=integration::InverseMellin_path(3,ResCrossSecFixpt,tauprime,2.,1.5,&Fixptpar, &errFixptPart);
      //Joint Part
      long double Joint=0.,errJoint=0.;
      if (pt < std::pow(3.5/ConstResum::Q,2.)){
	cout << " low pt " << endl;
	Borel::BorelJointC(3,tauprime,pt,ConstResum::as*beta_0,ResCrossSecJoint,&Fixptpar,&Joint,&errJoint,2.,3.,1.0); //Good convergence at small pt not at high pt
      }
      else{
	Borel::BorelJointCp(3,tauprime,pt,ConstResum::as*beta_0,ResCrossSecJoint,&Fixptpar,&Joint,&errJoint,2.,3.,1.0);//Quick and good convergence for  pt>4... C=1,2,3... no difference in this region
      }
      ris2.push_back((FixptPart+Joint)*tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb);
      error2.push_back((errFixptPart+errJoint)*tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb);
    }
    ris.push_back(ris2);
    error.push_back(error2);
    i++;
  }
  _MatchFun=_MatchFun2;
  *err=error;
  return ris;
}



long double CombResum::MatchingCrossSection(long double CMS, long double xp, int matching, long double *err){
  //Matching Function
  //matching==0 only Joint
  //matching==1 only FixedptResum
  //matching!=0 || matching !=1 used matching function defined in CombResum
  std::function<std::complex<long double>(std::complex<long double>, long double)> _MatchFun2=_MatchFun;
  if (matching==0) _MatchFun=T0;
  if (matching==1) _MatchFun=T1;
  //Fixed pt Part
  long double FixptPart=0.,errFixptPart=0.;
  long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/CMS,2.);
  ptstruct Fixptpar;
  Fixptpar.xp=xp;
  Fixptpar.CombR=this;
  FixptPart=integration::InverseMellin_path(3,MatchCrossSecFixpt,tauprime,2.,1.5,&Fixptpar, &errFixptPart);
  FixptPart*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  //Joint Part
  long double Joint=0.,errJoint=0.;
  Joint =integration::InverseMellin_path(3,MatchCrossSecJoint,tauprime,2.,1.5,&Fixptpar, &errJoint);
  Joint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  errJoint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  cout << "Joint= " << Joint << " +- " << errJoint << endl;
  _MatchFun=_MatchFun2;
  *err=errJoint+errFixptPart;
  return (FixptPart+Joint);
  
  
}

std::vector<long double> CombResum::MatchingCrossSection(long double CMS, std::vector<long double> xp,int matching, std::vector<long double> *err){
  //Matching Function
  //matching==0 only Joint
  //matching==1 only FixedptResum
  //matching!=0 || matching !=1 used matching function defined in CombResum
  std::function<std::complex<long double>(std::complex<long double>, long double)> _MatchFun2=_MatchFun; 
  if (matching==0) _MatchFun=T0;
  if (matching==1) _MatchFun=T1;
  
  std::vector<long double> ris,error;
  int i=0;
  for (auto pt : xp){
    std::cout << i << endl ;
    //Fixed pt Part
    long double FixptPart=0., errFixptPart=0.;
    long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+pt)+std::sqrt(pt))/CMS,2.);
    ptstruct Fixptpar;
    Fixptpar.xp=pt;
    Fixptpar.CombR=this;
    //cout << "CIAO" << endl;
    FixptPart=integration::InverseMellin_path(3,MatchCrossSecFixpt,tauprime,2.,1.5,&Fixptpar, &errFixptPart);
    //cout << FixptPart << " CIAO2" << endl;
    FixptPart*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
    //Joint Part
    long double Joint=0.,errJoint=0.;
    Joint =integration::InverseMellin_path(3,MatchCrossSecJoint,tauprime,2.,1.5,&Fixptpar, &errJoint);
    Joint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
    errJoint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
    ris.push_back((FixptPart+Joint));
    error.push_back((errFixptPart+errJoint));
    
    i++;
  }
  _MatchFun=_MatchFun2;
  *err=error;
  return ris;
}
std::vector<std::vector<long double>> CombResum::MatchingCrossSection(std::vector<long double> CMS, std::vector<long double> xp,int matching,
  std::vector<std::vector<long double>> *err){
  //Matching Function
  //matching==0 only Joint
  //matching==1 only FixedptResum
  //matching!=0 || matching !=1 used matching function defined in CombResum
  std::function<std::complex<long double> (std::complex<long double>, long double)> _MatchFun2=_MatchFun;
  if (matching==0) _MatchFun=T0;
  if (matching==1) _MatchFun=T1;
  
  std::vector<std::vector<long double>> ris, error;
  int i=0;
  for (auto cms : CMS ){
    cout << i << " " ;
    std::vector<long double> ris2,error2;
    for (auto pt : xp){
      //Fixed pt Part
      long double FixptPart=0.,errFixptPart=0.;
      long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+pt)+std::sqrt(pt))/cms,2.);
      ptstruct Fixptpar;
      Fixptpar.xp=pt;
      Fixptpar.CombR=this;
      FixptPart=integration::InverseMellin_path(3,MatchCrossSecFixpt,tauprime,2.,1.5,&Fixptpar, &errFixptPart);
      FixptPart*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
      //Joint Part
      long double Joint=0.,errJoint=0.;
      Joint =integration::InverseMellin_path(3,MatchCrossSecJoint,tauprime,2.,1.5,&Fixptpar, &errJoint);
      Joint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
      errJoint*=tauprime*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
      ris2.push_back((FixptPart+Joint));
      error2.push_back((errFixptPart+errJoint));
    }
    ris.push_back(ris2);
    error.push_back(error2);
    i++;
  }
  _MatchFun=_MatchFun2;
  *err=error;
  return ris;
}





//Core Functions

std::complex<long double> CombResum::FixptPartRes(std::complex<long double> N, long double xp){
  std::vector<std::complex<long double>> Fixchannel,Lumi;std::complex<long double> zero(0.,0.);
  Fixchannel=_Fix-> ComputeFixedptResum(N,xp);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return (_MatchFun(N,xp)*std::inner_product(Fixchannel.begin(),Fixchannel.end(),Lumi.begin(),zero));
}

std::complex<long double> CombResum::JointPartRes(std::complex<long double> N, std::complex<long double> lchi, long double xp){
  std::vector<std::complex<long double>>Jointchannel,Lumi;std::complex<long double> zero(0.,0.);
  Jointchannel=_Joint->ComputeJointRes(N,lchi,xp);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return ((1.-_MatchFun(N,xp))*std::inner_product(Jointchannel.begin(),Jointchannel.end(),Lumi.begin(),zero));
}

std::complex<long double> CombResum::FixptPartMatch(std::complex<long double> N, long double xp){
  std::vector<std::complex<long double>> Fixchannel,Lumi;std::complex<long double>zero(0.,0.);
  Fixchannel=_Fix->ComputeMatching(N,xp);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return (_MatchFun(N,xp)*std::inner_product(Fixchannel.begin(),Fixchannel.end(),Lumi.begin(),zero));
}

std::complex<long double> CombResum::JointPartMatch(std::complex<long double> N, long double xp){
  std::vector<std::complex<long double>>Jointchannel,Lumi;std::complex<long double> zero(0.,0.);
  Jointchannel=_Joint->ComputeMatching(N,xp);
  //cout << Jointchannel[0] << " " << Jointchannel[1] << " " << Jointchannel[2] << endl;
  Lumi=_Lumi.Higgs_Lum_N(N);
  //cout << (1.-_MatchFun(N,xp))*std::inner_product(Jointchannel.begin(),Jointchannel.end(),Lumi.begin(),zero) << endl;
  return ((1.-_MatchFun(N,xp))*(Jointchannel[0]*Lumi[0]+Jointchannel[1]*Lumi[1]+Jointchannel[2]*Lumi[2]));
}


//pt resummation Standard Borel Approach own implementation used for cross check

std::complex<long double> ResCrossSecpt(std::complex<long double> N,std::complex<long double> lchi, void *p){
  ptstruct par=*(ptstruct *)p;
  return(par.CombR->ptPartRes(N,lchi,par.xp));
}

std::complex<long double> CombResum::ptPartRes(std::complex<long double> N, std::complex<long double> lchi, long double xp){
  std::vector<std::complex<long double>>Jointchannel,Lumi;std::complex<long double> zero(0.,0.);
  Jointchannel=_Joint->ComputeptRes(N,lchi,xp);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return (std::inner_product(Jointchannel.begin(),Jointchannel.end(),Lumi.begin(),zero));
}


long double CombResum::ResummedCrossSectionpt(long double CMS, long double xp, long double *err){
  long double ptRes=0.,errptRes=0.;
  ptstruct Fixptpar;
  Fixptpar.xp=xp;
  Fixptpar.CombR=this;
  long double tauprime=std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/CMS,2.);
  Borel::BorelJointpt(3,tauprime,xp,ConstResum::as*beta_0,ResCrossSecpt,&Fixptpar,&ptRes,&errptRes,2.,3.,1.0);
  ptRes*=std::pow(ConstResum::Q/CMS,2)*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  errptRes*=std::pow(ConstResum::Q/CMS,2)*ConstResum::sigma0Higgs*ConstResum::GeVtopb;
  *err=errptRes;
  return (ptRes);
}


