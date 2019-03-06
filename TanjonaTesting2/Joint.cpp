#include "Joint.h"


const long double P0[12]={0.,2.833327652357064,-1.8330909624101532,0.7181879191200942,
      -0.0280403220046588,-0.181869786537805,0.532318519269331,-1.07281686995035,1.38194913357518,
      -1.11100841298484,0.506649587198046,-0.100672390783659};
   
const long double P1[4]={11./6.,-3.,3./2.,-1./3.};

double _Li2(double x, void *p) {
  return (std::pow(1.-x,3.)*gsl_sf_dilog(x));
}


Joint::Joint(const int ordres, const int ordmatch, int channel, bool Wilson, long double Nc, long double Nf): AP(Nf,Nc){
  _ordres=ordres; _ordmatch=ordmatch; _channel=channel;_Wilson=Wilson;_Nc=Nc;_Nf=Nf;
  
  //Initialise constants
  _Ca=_Nc; 
  _Cf=(_Nc*_Nc-1.)/(2.*_Nc);
  
  const long double zeta2=gsl_sf_zeta_int(2);
  const long double zeta3=gsl_sf_zeta_int(3);
  const long double zeta4=gsl_sf_zeta_int(4);
  
  
  //REMEMBER!!! ALL in unity as
  beta_0=(11.*_Ca-2.*_Nf)/(12.*M_PIl);
  beta_1=((17.*_Ca*_Ca-5.*_Ca*_Nf-3.*_Cf*_Nf)*2./3.)/(16.*M_PIl*M_PIl);
  beta_2=((2857./54.*_Ca*_Ca*_Ca+(_Cf*_Cf-205./18.*_Cf*_Ca-1415./54.*_Ca*_Ca)*_Nf
  +(11./9.*_Cf+79./54.*_Ca)*_Nf*_Nf))/std::pow(4.*M_PIl,3);
  
  Apt1g=_Ca/M_PIl;
  Apt2g=(_Ca/2.*(_Ca*(67./18.-zeta2)-5./9.*_Nf))/std::pow(M_PIl,2);
  Apt3g=(_Ca/4.*(_Ca*_Ca*(15503./648.-67./9.*zeta2-11.*zeta3+11./2.*zeta4)
  +_Cf*_Nf*(-55./24.+2.*zeta3)+_Ca*_Nf*(-2051./324.+10./9.*zeta2)-25./81.*_Nf*_Nf))/std::pow(M_PIl,3);
  
  Bpt1g=-2.*beta_0;
  Bpt2g=(-2.*(-8./3.*_Ca*_Nf+(32./3.+12.*zeta3)*_Ca*_Ca-2.*_Cf*_Nf)/(std::pow(4.,2)))/std::pow(M_PIl,2)+beta_0*_Ca*zeta2/M_PIl;
  
  Dpt2g=(_Ca*_Ca*(-101./27.+7./2.*zeta3)+14./27.*_Ca*_Nf)/std::pow(M_PIl,2);
  
  //Hpt1g=(3.*_Ca*zeta2+11./2.)/M_PIl;
  Hpt1g=(3.*_Ca*zeta2)/M_PIl;
  Hpt2g=(_Ca*_Ca*(93./16.+67./12.*zeta2-55./18.*zeta3+65./8.*zeta4)+_Ca*_Nf*(-5./3.-5./6.*zeta2-4./9.*zeta3))/std::pow(M_PIl,2);
  
  Apt1q=_Cf/M_PIl;
  
  if (Wilson){
    W1=0;
    W2=0;
  }
  else{
    W1=0;
    W2=0;
  }
  
  ULL=new std::complex<long double>*[2];
  ULL[0]=new std::complex<long double>[2];
  ULL[1]=new std::complex<long double>[2];
  
  UNLL=new std::complex<long double>*[2];
  UNLL[0]=new std::complex<long double>[2];
  UNLL[1]=new std::complex<long double>[2];
  
  UNNLL=new std::complex<long double>*[2];
  UNNLL[0]=new std::complex<long double>[2];
  UNNLL[1]=new std::complex<long double>[2];
  
  V1=new std::complex<long double>*[2];
  V1[0]=new std::complex<long double>[2];
  V1[1]=new std::complex<long double>[2];
  
  if (_ordres>2){
    std::cout << "Resummation is implemented only up to NNLL! Set resummation order to NNLL" << std::endl;
    _ordres=2;
  }
  if (_ordmatch>2){
    std::cout << "Matching is implemented only up to NNLO! Set matching order to NNLO" << std::endl;
    _ordmatch=2;
  }
  
  integration::ComputeFiniteDifferenceCoefficients(&Coeff,10,6);
  
  Li2_coeff.push_back(0.);
  double il = 0.;
  for (int i = 1; i < (orderexp + 1); i++) {
    il = static_cast<double> (i);
    Li2_coeff.push_back(1./(il*il));
  }
}

Joint::~Joint(){
}


void Joint::SetOrdRes(int ordres){
  _ordres=ordres;
}
void Joint::SetOrdMatch(int ordmatch){
  _ordmatch=ordmatch;
}
int Joint::SetChannel(int channel){
  _channel=channel;
}

//Sudakov Higgs
std::complex<long double> Joint::Sudakov_g(std::complex<long double> N, std::complex<long double> lchi, int ord){
  const long double LR=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUR,2.));
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lchimb=lchi-2.*ConstResum::as*beta_0*EulerGamma;
  
  const std::complex<long double> LC= lchi/ConstResum::as/beta_0;
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  std::complex<long double> Sud_g1_g(0.,0.),Sud_g2_g(0.,0.),Sud_g3_g(0.,0.);
  if (ord>=0) {
    Sud_g1_g=Apt1g/beta_0/beta_0*(lchi+std::log(1.-lchi));
  }
  if (ord>=1) {
    Sud_g2_g=Apt1g*beta_1/std::pow(beta_0,3)*((lchi+std::log(1.-lchi))/(1.-lchi)+0.5*std::pow(std::log(1.-lchi),2))
					    -Apt2g/beta_0/beta_0*(std::log(1.-lchi)+lchi/(1.-lchi))+Bpt1g/beta_0*std::log(1.-lchi)
					    +Apt1g/beta_0*(std::log(1.-lchi)+lchi/(1.-lchi))*LR;
  }
  if (ord>=2) {
    gsl_sf_result redilog,imdilog;
    long double r=std::abs(std::exp((lN-lchi)/(ConstResum::as*beta_0)));
    long double theta=std::arg(std::exp((lN-lchi)/(ConstResum::as*beta_0)));
    gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
    //std::complex<long double> risdilog(redilog.val,imdilog.val);
    std::complex<long double> risdilog(0.,0);
    if(_dilog==true){
      risdilog=redilog.val+II*imdilog.val;
    }
    Sud_g3_g=Apt1g*beta_1*beta_1/(2.*std::pow(beta_0,4))*((lchi+std::log(1.-lchi))/(std::pow(1.-lchi,2))*(lchi+(1.-2.*lchi)*std::log(1.-lchi)))
	+Apt1g*beta_2/std::pow(beta_0,3)*(((2.-3.*lchi)*lchi)/(2.*std::pow(1.-lchi,2))+std::log(1.-lchi))
	-Apt2g*beta_1/std::pow(beta_0,3)*(((2.-3.*lchi)*lchi)/(2.*std::pow(1.-lchi,2))+(1.-2.*lchi)*std::log(1.-lchi)/std::pow(1.-lchi,2))
	+Bpt1g*beta_1/beta_0*(lchi+std::log(1.-lchi))/(1.-lchi)-Apt3g*lchi*lchi/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
	-Bpt2g/beta_0*lchi/(1.-lchi)+Apt1g*lN/(1.-lN)*risdilog+Bpt1g*lchi/(1.-lchi)*LR+Apt2g/beta_0*lchi*lchi/std::pow(1.-lchi,2)*LR
	+Apt1g*beta_1/std::pow(beta_0,2)*((lchi*(1.-lchi)+(1.-2.*lchi)*std::log(1.-lchi))/std::pow(1.-lchi,2))*LR-Apt1g/2.*std::pow(lchi/(1.-lchi)*LR,2);    
  }
  return std::exp(1./ConstResum::as*Sud_g1_g+Sud_g2_g+ConstResum::as*Sud_g3_g); 
}

//Hard part Higgs
std::complex<long double> Joint::Hgggq1(std::complex<long double> N) {
  return (_Cf/2./(N+1.))/M_PIl;
}

std::complex<long double> Joint::Hggggreg2(std::complex<long double> N) {
  MellinFunc Func;
  std::complex<long double> hggggreg2;
  std::complex<long double> Nm1=N-1.;
  std::complex<long double> N1=N+1.;
  std::complex<long double> N2=N+2.;
  std::complex<long double> N3=N+3.;
  std::complex<long double> N4=N+4.;
  
  hggggreg2=1./12.*_Ca*(_Ca-_Nf)*Func.Logminus(N1)+1./216.*(18.*_Cf*_Nf*(Func.Logz3(N)+Func.Logz3(N1))
	   -36.*_Ca*_Ca*(Func.Logz3minusplus(N)+2.*Func.Logz3minusplus(N1)-Func.Logz3minusplus(N2)
	   -2.*Func.Logz3minusplus(N3)+Func.Logz3minusplus(N4)))+_Ca*_Ca*gsl_sf_zeta_int(2)*
	   (Func.Logplusplus(Nm1)+2.*Func.Logplusplus(N)+3.*Func.Logplusplus(N1)+2.*Func.Logplusplus(N2)+Func.Logplusplus(N3))
	   -_Ca*_Ca/3.*(Func.Logplus3plus(Nm1)+2.*Func.Logplus3plus(N)+3.*Func.Logplus3plus(N1)+2.*Func.Logplus3plus(N2)
	   +Func.Logplus3plus(N3))+1./24.*(2.*_Ca*_Nf*(Func.Logz2(N)+Func.Logz2(N1))+3.*_Cf*_Nf*
	   (3.*Func.Logz2(N)+Func.Logz2(N1))+_Ca*_Ca*(25.*Func.Logz2(N)-11.*Func.Logz2(N1)+44.*Func.Logz2(N2)))
	   +_Ca*_Ca/2.*(Func.Logz2Logminusminus(Nm1)-2.*Func.Logz2Logminusminus(N)+3.*Func.Logz2Logminusminus(N1)
	   -2.*Func.Logz2Logminusminus(N2)+Func.Logz2Logminusminus(N3))-_Ca*_Ca/2.*(Func.Logz2Logplusplus(Nm1)
	   +2.*Func.Logz2Logplusplus(N)+3.*Func.Logz2Logplusplus(N1)+2.*Func.Logz2Logplusplus(N2)+Func.Logz2Logplusplus(N3))
	   +_Ca*_Ca/2.*(Func.LogzLogminus2minus(Nm1)-2.*Func.LogzLogminus2minus(N)+3.*Func.LogzLogminus2minus(N1)
	   -2.*Func.LogzLogminus2minus(N2)+Func.LogzLogminus2minus(N3))
	   +_Ca*_Ca*(Func.LogzLogplus2plus(Nm1)+2.*Func.LogzLogplus2plus(N)+3.*Func.LogzLogplus2plus(N1)
	   +2.*Func.LogzLogplus2plus(N2)+Func.LogzLogplus2plus(N3))+1./216.*((78.*_Ca*_Nf+324.*_Cf*_Nf-2319.*_Ca*_Ca)*Func.Logz(N)
	   -216.*_Ca*_Ca*Func.Logz(Nm1)+(60.*_Ca*_Nf+324.*_Cf*_Nf-447.*_Ca*_Ca)*Func.Logz(N1)-1608.*_Ca*_Ca*Func.Logz(N2)
	   +216.*_Ca*_Ca*(Func.Li2mzLogzplus(Nm1)+2.*Func.Li2mzLogzplus(N)+3.*Func.Li2mzLogzplus(N1)+2.*Func.Li2mzLogzplus(N2)
	   + Func.Li2mzLogzplus(N3)))+_Ca*_Ca*(3.*Func.Li2zLogzplusminus(Nm1)-Func.Li2zLogzplusminus(N)+3.*Func.Li2zLogzplusminus(N1)
	   +Func.Li2zLogzplusminus(N2)-3.*Func.Li2zLogzplusminus(N3)+Func.Li2zLogzplusminus(N4))-_Ca*_Ca*(5.*Func.Li3zregminusplus(Nm1)
	   -Func.Li3zregminusplus(N)+5.*Func.Li3zregminusplus(N1)+Func.Li3zregminusplus(N2)-5.*Func.Li3zregminusplus(N3)+Func.Li3zregminusplus(N4))
	   +1./216.*((_Ca*_Ca*(-8272.+4.*N*(-2145.+N*(-1386.+211.*N)))+108.*_Cf*(-16.+7.*N*(1.+N))*_Nf+2.*_Ca*(332.+N*(327.+(141.-74.*N)*N))*_Nf)
	   /((Nm1)*N*N1*N2)+72.*_Ca*_Ca*(11.*Func.Li2minus(Nm1)-12.*Func.Li2minus(N)+12.*Func.Li2minus(N1)-11.*Func.Li2minus(N2))
	   -216.*_Ca*_Ca*(Func.Li3mzplus(Nm1)+2.*Func.Li3mzplus(N)+3.*Func.Li3mzplus(N1)+2.*Func.Li3mzplus(N2)+Func.Li3mzplus(N3))
	   +432.*_Ca*_Ca*(Func.Li3zoverplusplus(Nm1)+2.*Func.Li3zoverplusplus(N)+3.*Func.Li3zoverplusplus(N1)+2.*Func.Li3zoverplusplus(N2)+Func.Li3zoverplusplus(N3))
	   -108.*_Ca*_Ca*gsl_sf_zeta_int(3)*(-2.*Func.plus(Nm1)+17.*Func.plus(N)+22.*Func.plus(N1)+10.*Func.plus(N2)+12.*Func.plus(N3)));
	   return hggggreg2/std::pow(M_PIl,2);
}

std::complex<long double> Joint::Hgggq2(std::complex<long double> N){
  MellinFunc Func;
  std::complex<long double> hgggq2;
  std::complex<long double> Nm1=N-1.;
  std::complex<long double> N1=N+1.;
  std::complex<long double> N2=N+2.;
  std::complex<long double> N3=N+3.;
  std::complex<long double> N4=N+4.;
  hgggq2=_Cf/72.*(_Ca*(-152.*Func.Logminus(Nm1)+152.*Func.Logminus(N)-43.*Func.Logminus(N1))
	+4.*_Nf*(5.*Func.Logminus(Nm1)-5.*Func.Logminus(N)+Func.Logminus(N1))
	+9.*_Cf*(16.*Func.Logminus(Nm1)-16.*Func.Logminus(N)+5.*Func.Logminus(N1)))
	+_Cf/48.*(_Ca*(-22.*Func.Logminus2(Nm1)+22.*Func.Logminus2(N)-5.*Func.Logminus2(N1))
	+3.*_Cf*(6.*Func.Logminus2(Nm1)-6.*Func.Logminus2(N)+Func.Logminus2(N1))
	+2.*_Nf*(2.*Func.Logminus2(Nm1)-2.*Func.Logminus2(N)+Func.Logminus2(N1)))
	-(_Ca-_Cf)*(_Cf*(2.*Func.Logminus3(Nm1)-2.*Func.Logminus3(N)+Func.Logminus3(N1)))/24.
	-1./48.*_Cf*(_Cf*(-2.*Func.Logz3(N)+Func.Logz3(N1))+2.*_Ca*(2.*Func.Logz3(N)+Func.Logz3(N1)))
	-_Ca*_Cf*gsl_sf_zeta_int(2)/4.*(2.*Func.Logplus(Nm1)+2.*Func.Logplus(N)+Func.Logplus(N1))
	+_Ca*_Cf/12.*(2.*Func.Logplus3(Nm1)+2.*Func.Logplus3(N)+Func.Logplus3(N1))
	-1./96.*_Cf*(3.*_Cf*(4.*Func.Logz2(N)+3.*Func.Logz2(N1))-2.*_Ca*(36.*Func.Logz2(N)+9.*Func.Logz2(N1)
	+8.*Func.Logz2(N2)))+_Ca*_Cf/8.*(2.*Func.Logz2Logminus(Nm1)-2.*Func.Logz2Logminus(N)
	+Func.Logz2Logminus(N1))-_Ca*_Cf/8.*(2.*Func.Logz2Logplus(Nm1)+2.*Func.Logz2Logplus(N)
	+Func.Logz2Logplus(N1))+_Ca*_Cf/12.*(22.*Func.Li2z(Nm1)+24.*Func.Li2z(N)-9.*Func.Li2z(N1)+4.*Func.Li2z(N2))
	+1./8.*_Ca*_Cf*Func.Li2z2(N1)+_Cf/144.*(45.*_Cf*(-3.*Func.Logz(N)+Func.Logz(N1))-2.*_Ca*(72.*Func.Logz(Nm1)
	+321.*Func.Logz(N)-6.*Func.Logz(N1)+88.*Func.Logz(N2)))+_Ca*_Cf/12*(-22.*Func.LogzLogminus(Nm1)
	+24.*Func.LogzLogminus(N)-9.*Func.LogzLogminus(N1)+4.*Func.LogzLogminus(N2))
	+_Ca*_Cf/8.*(2.*Func.LogzLogminus2(Nm1)-2.*Func.LogzLogminus2(N)+Func.LogzLogminus2(N1))+1./4.*_Ca*_Cf*Func.LogzLogplus(N1)
	+_Ca*_Cf/2.*(2.*Func.Li2zLogz(Nm1)-2.*Func.Li2zLogz(N)+Func.Li2zLogz(N1))+_Ca*_Cf/8.*(2.*Func.Li2z2Logz(Nm1)
	+2.*Func.Li2z2Logz(N)+Func.Li2z2Logz(N1))-_Ca*_Cf/2.*(2.*Func.Li3z(Nm1)-4.*Func.Li3z(N)+Func.Li3z(N1))
	-3.*_Ca*_Cf/16.*(2.*Func.Li3z2(Nm1)+2.*Func.Li3z2(N)+Func.Li3z2(N1))-_Ca*_Cf/2.*(2.*Func.Li3overplus(Nm1)
	+2.*Func.Li3overplus(N)+Func.Li3overplus(N1))+(_Cf*(-N2*(27.*_Cf*Nm1*(-10.+3.*N)-4.*(56.+N*(43.+13.*N))*_Nf)
	+4.*_Ca*(-2014+72.*6.*gsl_sf_zeta_int(2)+540.*gsl_sf_zeta_int(3)+N*(-2905.+39.*6.*gsl_sf_zeta_int(2)+702.*gsl_sf_zeta_int(3)
	+N*(-1137.+20.*N+63.*6.*gsl_sf_zeta_int(2)+24.*6.*gsl_sf_zeta_int(2)*N+54.*(18.+7.*N)*gsl_sf_zeta_int(3))))))/(432*Nm1*N*N1*N2);
  return hgggq2/std::pow(M_PIl,2);
  
}

std::complex<long double> Joint::Hggqq2(std::complex<long double> N){
  return (_Cf*_Cf*(4.+N*(1.+N)*(8.+(-3.+N)*N))/(4.*N*N*std::pow(1.-N*N,2)))/std::pow(M_PIl,2);
}


void Joint::ComputeEvolution(std::complex<long double> N, std::complex<long double> lchi){
  AP.ComputeGamma(N,1);
  const std::complex<long double> rad=std::sqrt(std::pow(AP.gg0-AP.SS0,2)+4.*AP.gS0*AP.Sg0);
  std::complex<long double> lnEPLL(0.,0.),lnEMLL(0.,0.),lnEPNLL(0.,0.),lnEMNLL(0.,0.),lnEPNNLL(0.,0.),lnEMNNLL(0.,0.);
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lchimb=lchi-2.*ConstResum::as*beta_0*EulerGamma;
  const std::complex<long double> lNmb=2.*ConstResum::as*beta_0*std::log(N);
  
  
  const std::complex<long double> LC= lchi/ConstResum::as/beta_0;
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  const long double LR=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUR,2.));
  const long double LF=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUF,2.));
  if (_ordres >=0){
    lnEPLL=-Apt1g/2./beta_0/beta_0*lNmb*std::log(1.-lchi);
    lnEMLL=-Apt1q/2./beta_0/beta_0*lNmb*std::log(1.-lchi);

    ULL[0][0]=(std::exp(1./ConstResum::as*lnEMLL));
    ULL[0][1]=0.+II*0.;
    ULL[1][0]=0.+II*0.;
    ULL[1][1]=(std::exp(1./ConstResum::as*lnEPLL));
  }
  if (_ordres >=1){
    lnEPNLL=AP.plus0/beta_0*(std::log(1.-lchi)+ConstResum::as*beta_0*LF)-Apt1g/2./beta_0/beta_0*lNmb*
    (beta_1/beta_0*(std::log(1.-lchi))/(1.-lchi)+beta_0*lchi/(1.-lchi)*LR);
    lnEMNLL=AP.minus0/beta_0*(std::log(1.-lchi)+ConstResum::as*beta_0*LF)-Apt1q/2./beta_0/beta_0*lNmb*
    (beta_1/beta_0*(std::log(1.-lchi))/(1.-lchi)+beta_0*lchi/(1.-lchi)*LR);

    UNLL[0][0]=(std::exp(lnEMNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEPNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    UNLL[0][1]=(std::exp(lnEPNLL)-std::exp(lnEMNLL))*AP.Sg0/rad;
    UNLL[1][0]=(std::exp(lnEPNLL)-std::exp(lnEMNLL))*AP.gS0/rad;
    UNLL[1][1]=(std::exp(lnEPNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEMNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    V1[0][0]=(-beta_1*(4.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.SS0+beta_0*beta_0*(AP.gS1*AP.Sg0-AP.gS0*AP.Sg1+beta_1*AP.SS0)
	    -std::pow(beta_0,3)*AP.SS1+beta_0*(2.*AP.gg1*AP.gS0*AP.Sg0-(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)*(AP.gg0-AP.SS0)
	    +(2.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.SS1))/(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[0][1]=(-std::pow(beta_0,3)*AP.Sg1-beta_1*AP.Sg0*(4.*AP.Sg0*AP.gS0+std::pow(AP.gg0-AP.SS0,2))+beta_0*AP.Sg0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)
	    +(AP.gg0-AP.SS0)*(AP.gg1-AP.SS1))+beta_0*beta_0*(AP.Sg1*(-AP.gg0+AP.SS0)+AP.Sg0*(beta_1+AP.gg1-AP.SS1)))
	    /(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[1][0]=(-std::pow(beta_0,3)*AP.gS1-beta_1*AP.gS0*(4.*AP.Sg0*AP.gS0+std::pow(AP.gg0-AP.SS0,2))+beta_0*AP.gS0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)
	    +(AP.SS0-AP.gg0)*(AP.SS1-AP.gg1))+beta_0*beta_0*(AP.gS1*(+AP.gg0-AP.SS0)+AP.gS0*(beta_1-AP.gg1+AP.SS1)))
	    /(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[1][1]=(-beta_1*(4.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.gg0+beta_0*beta_0*(-AP.gS1*AP.Sg0+AP.gS0*AP.Sg1+beta_1*AP.gg0)
	    -std::pow(beta_0,3)*AP.gg1+beta_0*(2.*AP.SS1*AP.Sg0*AP.gS0-(AP.Sg1*AP.gS0+AP.Sg0*AP.gS1)*(AP.SS0-AP.gg0)
	    +(2.*AP.Sg0*AP.gS0+std::pow(AP.SS0-AP.gg0,2))*AP.gg1))/(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
  }
  if (_ordres >=2){
    lnEPNNLL=AP.plus0/beta_0*(std::log(1.-lchi)+ConstResum::as*(beta_1/beta_0*(std::log(1.-lchi)/(1.-lchi))+beta_0*lchi/(1.-lchi)*LR+beta_0*LF)
    +ConstResum::as*ConstResum::as*(-0.5*beta_0*beta_0*LR*LR+0.5*beta_0*beta_0*std::pow(-LF+LR,2)+LF*beta_1))
    -ConstResum::as*Apt1g/2./beta_0/beta_0*lNmb*((2.*(beta_1*beta_1-beta_0*beta_2)*lchi+2.*beta_1*beta_1*std::log(1.-lchi)
    -beta_1*beta_1*std::pow(std::log(1.-lchi),2))/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
    +LR*(beta_1*(2.-lchi)*lchi-beta_1*std::log(1.-lchi))/std::pow(1.-lchi,2)+beta_0*beta_0*LR*LR*(-2.+lchi)*lchi/(2.*std::pow(1.-lchi,2))
    );
    lnEMNNLL=AP.minus0/beta_0*(std::log(1.-lchi)+ConstResum::as*(beta_1/beta_0*(std::log(1.-lchi)/(1.-lchi))+beta_0*lchi/(1.-lchi)*LR+beta_0*LF)
      +ConstResum::as*ConstResum::as*(-0.5*beta_0*beta_0*LR*LR+0.5*beta_0*beta_0*std::pow(-LF+LR,2)+LF*beta_1))
    -ConstResum::as*Apt1q/2./beta_0/beta_0*lNmb*((2.*(beta_1*beta_1-beta_0*beta_2)*lchi+2.*beta_1*beta_1*std::log(1.-lchi)
    -beta_1*beta_1*std::pow(std::log(1.-lchi),2))/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
    +LR*(beta_1*(2.-lchi)*lchi-beta_1*std::log(1.-lchi))/std::pow(1.-lchi,2)+beta_0*beta_0*LR*LR*(-2.+lchi)*lchi/(2.*std::pow(1.-lchi,2))
    );
    
    
    UNNLL[0][0]=(std::exp(lnEMNNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEPNNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    UNNLL[0][1]=(std::exp(lnEPNNLL)-std::exp(lnEMNNLL))*AP.Sg0/rad;
    UNNLL[1][0]=(std::exp(lnEPNNLL)-std::exp(lnEMNNLL))*AP.gS0/rad;
    UNNLL[1][1]=(std::exp(lnEPNNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEMNNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
  }
  return;
}

std::vector<std::complex<long double> > Joint::ComputeJointRes(std::complex<long double> N, std::complex<long double> lchi, long double xp){
  /* channel legends
   * 0 all channels
   * 1 gg only
   * 2 gq only
   * 3 qq only (means all flavours since Higgs resummation is flavour blind)
   */
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> LC= lchi/ConstResum::as/beta_0;
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  
  const std::complex<long double> lchimb=lchi-2.*ConstResum::as*beta_0*EulerGamma;
  const std::complex<long double> lNmb=2.*ConstResum::as*beta_0*std::log(N);
  
  const long double LR=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUR,2.));
  const long double LF=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUF,2.));
  std::complex<long double> RggNLL,RggNNLL,RgqNLL;
  std::vector<std::complex<long double> > ris;
  ComputeEvolution(N,lchi);
  gsl_sf_result redilog,imdilog;
  long double r=std::abs(std::exp((lN-lchi)/(ConstResum::as*beta_0)));
  long double theta=std::arg(std::exp((lN-lchi)/(ConstResum::as*beta_0)));
  gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
  std::complex<long double> risdilog(0.,0);
  //std::complex<long double> risdilog(redilog.val,imdilog.val);
  _dilog=false;
  if (xp > std::pow(3.5/ConstResum::Q,2.)){
    _dilog=true;
    risdilog=redilog.val+II*imdilog.val;
  }
  //LL resummation
  switch(_ordres){
    case(0):{
      if ((_channel==0)||(_channel==1)){
	ris.push_back(std::pow(ULL[1][1],2)*Sudakov_g(N,lchi,0));
      }
      else ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      break;
    }
    case(1):{
      RggNLL=std::exp(1./2./beta_0/beta_0*lNmb*(Apt2g-beta_1/beta_0*Apt1g)*lchi/(1.-lchi));

      
      RgqNLL=std::exp(-std::log(1.-lchi));
      if((_channel==0)||(_channel==1)){
	ris.push_back((std::pow(RggNLL*UNLL[1][1],2)
	+ConstResum::as*((Hpt1g+Apt1g*risdilog-2.*beta_0*LR)*std::pow(ULL[1][1],2)))*Sudakov_g(N,lchi,1));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==2)){
	ris.push_back((std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]
	+ConstResum::as*(Hgggq1(N)*ULL[1][1]*ULL[0][0]-std::pow(ULL[1][1],2)*V1[1][0]
	+ULL[1][1]*ULL[0][0]*V1[1][0]))*Sudakov_g(N,lchi,1));

      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==3)){
	ris.push_back((std::pow(RggNLL*UNLL[1][0],2))*Sudakov_g(N,lchi,1));
      }
      else ris.push_back((0.,0.));
      break;
    }
    case(2):{
      RggNLL=std::exp(1./2./beta_0/beta_0*lNmb*(Apt2g-beta_1/beta_0*Apt1g)*lchi/(1.-lchi));
           
      RggNNLL=std::exp(ConstResum::as*(V1[1][1]*(lchi/(1.-lchi)-ConstResum::as*beta_0*LF)
      +1./std::pow(beta_0,3)*(beta_1*beta_1/2./beta_0*Apt1g-beta_2*Apt1g/2.-beta_1*Apt2g/2.+beta_0*Apt3g/2.)*lNmb
      *(2.-lchi)*lchi/(2.*std::pow(1.-lchi,2))+1./2./beta_0*(Apt2g-beta_1/beta_0*Apt1g)*lNmb*lchi*(2.-lchi)/std::pow(1.-lchi,2)*LR
      -beta_1*(1./2./beta_0/beta_0*lNmb*(Apt2g-beta_1/beta_0*Apt1g))/beta_0*std::log(1.-lchi)/std::pow(1.-lchi,2)
      -1./4.*Dpt2g/beta_0*lN));
      
      RgqNLL=std::exp(-std::log(1.-lchi));
      if((_channel==0)||(_channel==1)){
	ris.push_back((std::pow(RggNNLL*UNNLL[1][1],2)
	+ConstResum::as*(2.*Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][1]+(Hpt1g+Apt1g*risdilog)*std::pow(RggNLL*UNLL[1][1],2)
	+2.*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][1]*V1[1][0]-2.*RggNLL*RggNLL*UNLL[1][1]*UNLL[1][0]*V1[0][1]
	-2.*std::pow(RggNLL*UNLL[1][1],2)*beta_0*LR)+ConstResum::as*ConstResum::as*((Hggggreg2(N)+Hpt2g-3.*beta_0*(Hpt1g+Apt1g*risdilog)*LR
	-2.*beta_1*LR+3.*beta_0*beta_0*LR*LR)*std::pow(ULL[1][1],2)))*Sudakov_g(N,lchi,2));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==2)){
	ris.push_back((std::pow(RggNNLL,2)*UNNLL[1][1]*UNNLL[1][0]+ConstResum::as*(Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][0]
	+std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[1][1]+(Hpt1g+risdilog*Apt1g)*std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]
	  -std::pow(RggNLL*UNLL[1][1],2)*V1[1][0]+RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][0]*V1[1][0]-std::pow(RggNLL*UNLL[1][0],2)*V1[0][1]
	  -std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[0][0]-2.*beta_0*std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*LR)
	+ConstResum::as*ConstResum::as*(Hgggq2(N)*ULL[1][1]*ULL[0][0]
	-3.*Hgggq1(N)*beta_0*ULL[1][1]*ULL[0][0]*LR-ULL[1][1]*ULL[0][0]*V1[1][0]*beta_0*LF))*Sudakov_g(N,lchi,2));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==3)){
	ris.push_back((std::pow(RggNNLL*UNNLL[1][0],2)+2.*ConstResum::as*(Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][0]*UNLL[0][0]
	+(Hpt1g+Apt1g*risdilog)*std::pow(RggNLL*UNLL[1][0],2)+std::pow(RggNLL*UNLL[1][0],2)*V1[1][1]
	  -std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[1][0]+RggNLL*RgqNLL*UNLL[1][0]*UNLL[0][0]*V1[1][0]
	  -std::pow(RggNLL*UNLL[1][0],2)*V1[0][0]-2.*beta_0*LR*std::pow(RggNLL*UNLL[1][0],2))
	+ConstResum::as*ConstResum::as*(Hggqq2(N)*std::pow(ULL[0][0],2)))*Sudakov_g(N,lchi,2));
      }
      else ris.push_back((0.,0.));
      break;
    }
  }
  return ris; 
}

//Standard pt Resummation implemented as cross check

void Joint::ComputeEvolutionpt(std::complex<long double> N, std::complex<long double> lchi){
  AP.ComputeGamma(N,1);
  const std::complex<long double> rad=std::sqrt(std::pow(AP.gg0-AP.SS0,2)+4.*AP.gS0*AP.Sg0);
  std::complex<long double> lnEPLL(0.,0.),lnEMLL(0.,0.),lnEPNLL(0.,0.),lnEMNLL(0.,0.),lnEPNNLL(0.,0.),lnEMNNLL(0.,0.);
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lchimb=lchi-2.*ConstResum::as*beta_0*EulerGamma;
  const std::complex<long double> lNmb=2.*ConstResum::as*beta_0*std::log(N);
  
  
  const std::complex<long double> LC= lchi/ConstResum::as/beta_0;
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  const long double LR=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUR,2.));
  const long double LF=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUF,2.));
  if (_ordres >=0){
    lnEPLL=0.;
    lnEMLL=0.;
    ULL[0][0]=(std::exp(1./ConstResum::as*lnEMLL));
    ULL[0][1]=0.+II*0.;
    ULL[1][0]=0.+II*0.;
    ULL[1][1]=(std::exp(1./ConstResum::as*lnEPLL));
  }
  if (_ordres >=1){
    lnEPNLL=AP.plus0/beta_0*(std::log(1.-lchi));
    lnEMNLL=AP.minus0/beta_0*(std::log(1.-lchi));
    UNLL[0][0]=(std::exp(lnEMNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEPNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    UNLL[0][1]=(std::exp(lnEPNLL)-std::exp(lnEMNLL))*AP.Sg0/rad;
    UNLL[1][0]=(std::exp(lnEPNLL)-std::exp(lnEMNLL))*AP.gS0/rad;
    UNLL[1][1]=(std::exp(lnEPNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEMNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    V1[0][0]=(-beta_1*(4.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.SS0+beta_0*beta_0*(AP.gS1*AP.Sg0-AP.gS0*AP.Sg1+beta_1*AP.SS0)
	    -std::pow(beta_0,3)*AP.SS1+beta_0*(2.*AP.gg1*AP.gS0*AP.Sg0-(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)*(AP.gg0-AP.SS0)
	    +(2.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.SS1))/(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[0][1]=(-std::pow(beta_0,3)*AP.Sg1-beta_1*AP.Sg0*(4.*AP.Sg0*AP.gS0+std::pow(AP.gg0-AP.SS0,2))+beta_0*AP.Sg0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)
	    +(AP.gg0-AP.SS0)*(AP.gg1-AP.SS1))+beta_0*beta_0*(AP.Sg1*(-AP.gg0+AP.SS0)+AP.Sg0*(beta_1+AP.gg1-AP.SS1)))
	    /(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[1][0]=(-std::pow(beta_0,3)*AP.gS1-beta_1*AP.gS0*(4.*AP.Sg0*AP.gS0+std::pow(AP.gg0-AP.SS0,2))+beta_0*AP.gS0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)
	    +(AP.SS0-AP.gg0)*(AP.SS1-AP.gg1))+beta_0*beta_0*(AP.gS1*(+AP.gg0-AP.SS0)+AP.gS0*(beta_1-AP.gg1+AP.SS1)))
	    /(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[1][1]=(-beta_1*(4.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.gg0+beta_0*beta_0*(-AP.gS1*AP.Sg0+AP.gS0*AP.Sg1+beta_1*AP.gg0)
	    -std::pow(beta_0,3)*AP.gg1+beta_0*(2.*AP.SS1*AP.Sg0*AP.gS0-(AP.Sg1*AP.gS0+AP.Sg0*AP.gS1)*(AP.SS0-AP.gg0)
	    +(2.*AP.Sg0*AP.gS0+std::pow(AP.SS0-AP.gg0,2))*AP.gg1))/(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
  }
  if (_ordres >=2){
    
    lnEPNNLL=AP.plus0/beta_0*(std::log(1.-lchi)+ConstResum::as*(beta_1/beta_0*(std::log(1.-lchi)/(1.-lchi))+beta_0*lchi/(1.-lchi)*LR+beta_0*LF));
    lnEMNNLL=AP.minus0/beta_0*(std::log(1.-lchi)+ConstResum::as*(beta_1/beta_0*(std::log(1.-lchi)/(1.-lchi))+beta_0*lchi/(1.-lchi)*LR+beta_0*LF));
    
    UNNLL[0][0]=(std::exp(lnEMNNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEPNNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    UNNLL[0][1]=(std::exp(lnEPNNLL)-std::exp(lnEMNNLL))*AP.Sg0/rad;
    UNNLL[1][0]=(std::exp(lnEPNNLL)-std::exp(lnEMNNLL))*AP.gS0/rad;
    UNNLL[1][1]=(std::exp(lnEPNNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEMNNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
  }
  return;
}


std::vector<std::complex<long double> > Joint::ComputeptRes(std::complex<long double> N, std::complex<long double> lchi, long double xp){
  /* channel legends
   * 0 all channels
   * 1 gg only
   * 2 gq only
   * 3 qq only (means all flavours since Higgs resummation is flavour blind)
   */
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> LC= lchi/ConstResum::as/beta_0;
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  
  const std::complex<long double> lchimb=lchi-2.*ConstResum::as*beta_0*EulerGamma;
  const std::complex<long double> lNmb=2.*ConstResum::as*beta_0*std::log(N);
  
  const long double LR=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUR,2.));
  const long double LF=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUF,2.));
  std::complex<long double> RggNLL,RggNNLL,RgqNLL;
  std::vector<std::complex<long double> > ris;
  ComputeEvolutionpt(N,lchi);
  gsl_sf_result redilog,imdilog;
  long double r=std::abs(std::exp((lN-lchi)/(ConstResum::as*beta_0)));
  long double theta=std::arg(std::exp((lN-lchi)/(ConstResum::as*beta_0)));
  gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
  std::complex<long double> risdilog(0.,0);
  _dilog=false;
  //LL resummation
  switch(_ordres){
    case(0):{
      if ((_channel==0)||(_channel==1)){
	ris.push_back(std::pow(ULL[1][1],2)*Sudakov_g(N,lchi,0));
      }
      else ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      break;
    }
    case(1):{
      RggNLL=std::exp(0.);
      if((_channel==0)||(_channel==1)){
	ris.push_back((std::pow(RggNLL*UNLL[1][1],2)
	+ConstResum::as*((Hpt1g-2.*beta_0*LR)*std::pow(ULL[1][1],2)))*Sudakov_g(N,lchi,1));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==2)){
	ris.push_back((std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]
	+ConstResum::as*(Hgggq1(N)*ULL[1][1]*ULL[0][0]-std::pow(ULL[1][1],2)*V1[1][0]
	+ULL[1][1]*ULL[0][0]*V1[1][0]))*Sudakov_g(N,lchi,1));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==3)){
	ris.push_back((std::pow(RggNLL*UNLL[1][0],2))*Sudakov_g(N,lchi,1));
      }
      else ris.push_back((0.,0.));
      break;
    }
    case(2):{
      
      RggNLL=std::exp(0.);
      
      RggNNLL=std::exp(ConstResum::as*(V1[1][1]*(lchi/(1.-lchi))-1./4.*Dpt2g/beta_0*lN));
      
      RgqNLL=std::exp(-std::log(1.-lchi));
      if((_channel==0)||(_channel==1)){
	ris.push_back((std::pow(RggNNLL*UNNLL[1][1],2)
	+ConstResum::as*(2.*Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][1]+(Hpt1g+Apt1g*risdilog)*std::pow(RggNLL*UNLL[1][1],2)
	+2.*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][1]*V1[1][0]-2.*RggNLL*RggNLL*UNLL[1][1]*UNLL[1][0]*V1[0][1]
	-2.*std::pow(RggNLL*UNLL[1][1],2)*beta_0*LR)+ConstResum::as*ConstResum::as*((Hggggreg2(N)+Hpt2g-3.*beta_0*(Hpt1g+Apt1g*risdilog)*LR
	-2.*beta_1*LR+3.*beta_0*beta_0*LR*LR)*std::pow(ULL[1][1],2)))*Sudakov_g(N,lchi,2));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==2)){
	ris.push_back((std::pow(RggNNLL,2)*UNNLL[1][1]*UNNLL[1][0]+ConstResum::as*(Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][0]
	+std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[1][1]+(Hpt1g+risdilog*Apt1g)*std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]
	  -std::pow(RggNLL*UNLL[1][1],2)*V1[1][0]+RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][0]*V1[1][0]-std::pow(RggNLL*UNLL[1][0],2)*V1[0][1]
	  -std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[0][0]-2.*beta_0*std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*LR)
	+ConstResum::as*ConstResum::as*(Hgggq2(N)*ULL[1][1]*ULL[0][0]
	-3.*Hgggq1(N)*beta_0*ULL[1][1]*ULL[0][0]*LR-ULL[1][1]*ULL[0][0]*V1[1][0]*beta_0*LF))*Sudakov_g(N,lchi,2));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==3)){
	ris.push_back((std::pow(RggNNLL*UNNLL[1][0],2)+2.*ConstResum::as*(Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][0]*UNLL[0][0]
	+(Hpt1g+Apt1g*risdilog)*std::pow(RggNLL*UNLL[1][0],2)+std::pow(RggNLL*UNLL[1][0],2)*V1[1][1]
	  -std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[1][0]+RggNLL*RgqNLL*UNLL[1][0]*UNLL[0][0]*V1[1][0]
	  -std::pow(RggNLL*UNLL[1][0],2)*V1[0][0]-2.*beta_0*LR*std::pow(RggNLL*UNLL[1][0],2))
	+ConstResum::as*ConstResum::as*(Hggqq2(N)*std::pow(ULL[0][0],2)))*Sudakov_g(N,lchi,2));
      }
      else ris.push_back((0.,0.));
      break;
    }
  }
  return ris; 
}

//Matching

//Useful Core Function For Matching
std::complex<long double> GenFunLCHI(std::complex<long double> N,long double xi, long double xp){
  const long double EulerGamma=0.57721566;
  return(2.*std::pow(xp,-0.5-0.5*xi)*std::pow(N,1.+xi)*LBesselK(1.+xi,2.*N*std::sqrt(xp))*Gamma_inv(-xi)*std::exp(2.*EulerGamma*xi));
}

std::complex<long double> GenFunLCHI2(std::complex<long double> N,long double xi, long double xp, long double j){
  const long double EulerGamma=0.57721566;
  return(2.*std::pow(xp,-0.5-0.5*xi+0.5*j)*std::pow(N,1.+xi+j)*LBesselK(1.+xi-j,2.*N*std::sqrt(xp))*Gamma_inv(j-xi)*std::exp(2.*EulerGamma*xi));
}


std::complex<long double> Joint::LCHI(int power, std::complex<long double> N, long double xp){
  long double hh=0.01;
  int order=6;
  
  if (power == 0)
    return GenFunLCHI(N,0,xp);

  if(power>DerMAX){
    std::cout<< "ERROR: NO COEFFICIENT EVALUATED AT THIS ORDER OF DERIVATIVE" << std::endl;
    return (0.,0.);
  }
  std::vector<long double> al;
  al.push_back(0.0);
  int NN= order+power;
  for (int i=1;i<(NN/2+1);i++){
    al.push_back((static_cast<long double> (i))*hh);
    al.push_back(-(static_cast<long double> (i))*hh);
  }
  std::complex<long double> ris;
  for (int i=0;i<(NN);i++){
    if (Coeff[power][NN-1][i]!=0)
      ris+=Coeff[power][NN-1][i]*GenFunLCHI(N,al[i],xp);
  }
  ris/=std::pow(hh,(long double) (power));
  return ris;
}

std::complex<long double> Joint::GenFunLi2(std::complex<long double> N,long double xi, long double xp){
  std::complex<long double> sum;
  long double il;
  for (int i = 0; i < (orderexp+1); i++) {
    il = static_cast<long double> (i);
    sum += Li2_coeff[i]*GenFunLCHI2(N,xi,xp,il);
  }
  return sum;
}

std::complex<long double> Joint::Li2Integral(int power, std::complex<long double> N, long double xp){
  long double hh=0.01;
  int order=6;
  
  if (power == 0)
    return GenFunLi2(N,0,xp);

  if(power>DerMAX){
    std::cout<< "ERROR: NO COEFFICIENT EVALUATED AT THIS ORDER OF DERIVATIVE" << std::endl;
    return (0.,0.);
  }
  std::vector<long double> al;
  al.push_back(0.0);
  int NN= order+power;
  for (int i=1;i<(NN/2+1);i++){
    al.push_back((static_cast<long double> (i))*hh);
    al.push_back(-(static_cast<long double> (i))*hh);
  }
  std::complex<long double> ris;
  for (int i=0;i<(NN);i++){
    if (Coeff[power][NN-1][i]!=0)
      ris+=Coeff[power][NN-1][i]*GenFunLi2(N,al[i],xp);
  }
  ris/=std::pow(hh,(long double) (power));
  return ris;
}




std::vector<std::complex<long double> > Joint::ComputeMatching(std::complex<long double> N, long double xp){
  /* channel legends
   * 0 all channels
   * 1 gg only
   * 2 gq only
   * 3 qq only (means all flavours since Higgs resummation is flavour blind)
   */
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> LN=2.*std::log(Nbar);
  const long double LR=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUR,2.));
  const long double LF=std::log(ConstResum::Q*ConstResum::Q/std::pow(ConstResum::MUF,2.));
  std::complex<long double> risdilog(0.,0), risdiloglogchi(0.,0.), risdiloglogchi2(0.,0.);
  _dilog=false;
  if (xp > std::pow(3.5/ConstResum::Q,2.)){
    _dilog=true;
  }
  if (_dilog) {
	    risdilog =Li2Integral(0,N,xp);
	    risdiloglogchi = Li2Integral(1,N,xp);
	    risdiloglogchi2 = Li2Integral(2,N,xp);
  }
  AP.ComputeGamma(N,1);
  
  std::vector<std::complex<long double> > match;
  switch (_ordres){
    case (0):{
      switch(_ordmatch){
	case(0):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back((0.,0.));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
	case(1):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back(ConstResum::as*(-Apt1g*LCHI(2,N,xp)/2.+Apt1g*LCHI(1,N,xp)*LN));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
	case(2):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back(ConstResum::as*(-Apt1g*LCHI(2,N,xp)/2.+Apt1g*LCHI(1,N,xp)*LN)
	  +ConstResum::as*ConstResum::as*(1./8.*Apt1g*Apt1g*LCHI(4,N,xp)-0.5*Apt1g*Apt1g*LCHI(3,N,xp)*LN
	  +0.5*Apt1g*Apt1g*LCHI(2,N,xp)*LN*LN-1./3.*Apt1g*LCHI(3,N,xp)*beta_0+0.5*beta_0*Apt1g*LCHI(2,N,xp)*LN));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
     }
  }
    case(1):{
      switch(_ordmatch){
	case(0):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back((0.,0.));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
	case(1):{
	  if ((_channel==0)||(_channel==1))
	  //match.push_back(ConstResum::as*((Apt1g*risdilog)-Apt1g*LCHI(2,N,xp)/2.-2.*LR*beta_0+LCHI(1,N,xp)*(-Bpt1g-2.*AP.gg0)+2.*LF*AP.gg0));
	  match.push_back(ConstResum::as*(-Apt1g*LCHI(2,N,xp)/2.-2.*LR*beta_0+LCHI(1,N,xp)*(-Bpt1g-2.*AP.gg0)+2.*LF*AP.gg0));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2))
	  match.push_back(ConstResum::as*(Hgggq1(N)-LCHI(1,N,xp)*AP.gS0+LF*AP.gS0));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
	case(2):{
	  std::complex<long double> rad=std::sqrt((AP.gg0-AP.SS0)*(AP.gg0-AP.SS0)+4.*AP.gS0*AP.Sg0);
	  if ((_channel==0)||(_channel==1))
	  match.push_back(ConstResum::as*((Hpt1g+Apt1g*risdilog)-Apt1g*LCHI(2,N,xp)/2.-2.*LR*beta_0+LCHI(1,N,xp)*(-Bpt1g-2.*AP.gg0)+2.*LF*AP.gg0)
	  +ConstResum::as*ConstResum::as*(Apt1g*Apt1g*LCHI(4,N,xp)/8.+1./6.*Apt1g*LCHI(3,N,xp)*(3.*Bpt1g-2.*beta_0+6.*AP.gg0)+LF*LF*(2.*AP.gg0*AP.gg0+AP.gS0*AP.Sg0)
	  +0.5*LCHI(2,N,xp)*(-Apt2g+Bpt1g*Bpt1g-Apt1g*Hpt1g-Bpt1g*beta_0+3.*Apt1g*LR*beta_0+4.*Bpt1g*AP.gg0-2.*Apt1g*LF*AP.gg0-2.*beta_0*AP.gg0+4.*AP.gg0*AP.gg0
	  +2.*AP.gS0*AP.Sg0)-0.5*Apt1g*Apt1g*risdiloglogchi2+LCHI(1,N,xp)*((-2.*beta_0*(-Apt2g*LN+4.*LF*AP.gg0*AP.gg0+Bpt1g*(Hpt1g-2*LR*beta_0+2.*LF*AP.gg0)+2.*LF*AP.gS0*AP.Sg0)
	  *rad-Apt1q*LN*(LR*beta_0*beta_0-beta_1)*(-AP.gg0+AP.SS0+rad)-Apt1g*LN*(-2.*Hpt1g*beta_0*rad+LR*beta_0*beta_0*(AP.gg0-AP.SS0+5.*rad)
	  +beta_1*(-AP.gg0+AP.SS0+rad)))/(2.*beta_0*rad))+Apt1g*(-Bpt1g+Apt1g*LN)*risdiloglogchi));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2))
	  match.push_back(ConstResum::as*(Hgggq1(N)-LCHI(1,N,xp)*AP.gS0+LF*AP.gS0)+ConstResum::as*ConstResum::as*(0.5*Apt1g*LCHI(3,N,xp)*AP.gS0+0.5*LF*LF*AP.gS0*(3.*AP.gg0+AP.SS0)
	  +0.5*LCHI(2,N,xp)*(-Apt1g*(Hgggq1(N)+LF*AP.gS0)+AP.gS0*(2.*Bpt1g-beta_0+3.*AP.gg0+AP.SS0))+0.5*LCHI(1,N,xp)*(-6.*LF*AP.gg0*AP.gS0-2.*Bpt1g*(Hgggq1(N)+LF*AP.gS0)
	  -(Apt1g*LN*LR*beta_0*AP.gS0-Apt1q*LN*LR*beta_0*AP.gS0-Apt1g*LN*beta_1/beta_0*AP.gS0+Apt1q*LN*beta_1/beta_0*AP.gS0)/rad
	  +1./beta_0/beta_0*((Apt1g+Apt1q)*Hgggq1(N)*LN*beta_0*beta_0+(-Apt1g+Apt1q)*LN*beta_1*AP.gS0-2.*LF*beta_0*beta_0*AP.gS0*AP.SS0)
	  +((Apt1g-Apt1q)*LN*(beta_0*beta_0*AP.gS1-AP.gS0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)+(AP.gg0-AP.SS0)*(AP.gg1-AP.SS1))
	  +beta_0*(AP.gg1*AP.gS0-AP.gg0*AP.gS1+AP.gS1*AP.SS0-AP.gS0*AP.SS1)))/(beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-(AP.gg0-AP.SS0)*(AP.gg0-AP.SS0))))));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3))
	  match.push_back(ConstResum::as*ConstResum::as*(LCHI(2,N,xp)*AP.gS0*AP.gS0-2.*LCHI(1,N,xp)*LF*AP.gS0*AP.gS0+LF*LF*AP.gS0*AP.gS0));
	  else match.push_back((0.,0.));
	}
      }
    }
    case(2):{
      switch(_ordmatch){
	case(0):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back((0.,0.));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
	case(1):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back(ConstResum::as*(Hpt1g-Bpt1g*LCHI(1,N,xp)-Apt1g*LCHI(2,N,xp)/2.-2.*LR*beta_0-2.*LCHI(1,N,xp)*AP.gg0+2.*LF*AP.gg0+Apt1g*risdilog));
	  if ((_channel==0)||(_channel==2))
	  match.push_back(ConstResum::as*(Hgggq1(N)+(-LCHI(1,N,xp)+LF)*AP.gS0));
	  else match.push_back((0.,0.));
	  match.push_back((0.,0.));
	}
	case(2):{
	  if ((_channel==0)||(_channel==1))
	  match.push_back(ConstResum::as*(Hpt1g-Bpt1g*LCHI(1,N,xp)-Apt1g*LCHI(2,N,xp)/2.-2.*LR*beta_0-2.*LCHI(1,N,xp)*AP.gg0+2.*LF*AP.gg0+Apt1g*risdilog)
	  +ConstResum::as*ConstResum::as*((Hggggreg2(N)+Hpt2g)-Bpt2g*LCHI(1,N,xp)-Apt2g*LCHI(2,N,xp)/2.+Bpt1g*Bpt1g*LCHI(2,N,xp)/2.-0.5*Apt1g*Hpt1g*LCHI(2,N,xp)
	  +Apt1g*Apt1g/8.*LCHI(4,N,xp)-Dpt2g*LN/2.-1./3.*Apt1g*beta_0*LCHI(3,N,xp)-3.*Hpt1g*LR*beta_0+3./2.*Apt1g*LCHI(2,N,xp)*LR*beta_0
	  +3.*LR*LR*beta_0*beta_0-2.*LR*beta_1-2.*Hpt1g*LCHI(1,N,xp)*AP.gg0+Apt1g*LCHI(3,N,xp)*AP.gg0+2.*Hpt1g*LF*AP.gg0
	  -Apt1g*LCHI(2,N,xp)*LF*AP.gg0-LCHI(2,N,xp)*beta_0*AP.gg0+std::pow(-LF+LR,2)*beta_0*AP.gg0+6.*LCHI(1,N,xp)*LR*beta_0*AP.gg0-4.*LF*LR*beta_0*AP.gg0
	  -LR*LR*beta_0*AP.gg0+2.*LCHI(2,N,xp)*AP.gg0*AP.gg0-4.*LCHI(1,N,xp)*LF*AP.gg0*AP.gg0+2.*LF*LF*AP.gg0*AP.gg0+0.5*Bpt1g*
	  (-2.*Hpt1g*LCHI(1,N,xp)+Apt1g*LCHI(3,N,xp)-LCHI(2,N,xp)*beta_0+6.*LR*beta_0*LCHI(1,N,xp)+4.*LCHI(2,N,xp)*AP.gg0-4.*LF*AP.gg0*LCHI(1,N,xp))
	  -2.*LCHI(1,N,xp)*AP.gg1+2.*LF*AP.gg1+(-LF)*(-2.*Hgggq1(N)+(LCHI(1,N,xp)-LF)*AP.gS0)*AP.Sg0+(-2.*Hgggq1(N)*LCHI(1,N,xp)+(LCHI(2,N,xp)-LF*LCHI(1,N,xp))*AP.gS0)*AP.Sg0
	  -0.5*Apt1g*(2.*Bpt1g*risdiloglogchi+Apt1g*risdiloglogchi2-2.*LN*beta_0*risdilog+6.*LR*beta_0*risdilog+4.*risdiloglogchi*AP.gg0-4.*LF*AP.gg0*risdilog)));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2))
	  match.push_back(ConstResum::as*(Hgggq1(N)+(LF-LCHI(1,N,xp))*AP.gS0)+ConstResum::as*ConstResum::as*0.5*(2.*Hgggq2(N)+2.*Hgggq1(N)*LCHI(1,N,xp)*beta_0-6.*Hgggq1(N)*LR*beta_0
	  -2.*Hgggq1(N)*LCHI(1,N,xp)*AP.gg0
	  +2.*Hgggq1(N)*LF*AP.gg0-2.*Hpt1g*LCHI(1,N,xp)*AP.gS0+2.*Hpt1g*LF*AP.gS0-LCHI(2,N,xp)*beta_0*AP.gS0+std::pow(-LF+LR,2)*beta_0*AP.gS0+6.*LCHI(1,N,xp)*LR*beta_0*AP.gS0
	  -4.*LF*LR*beta_0*AP.gS0-LR*LR*beta_0*AP.gS0+3.*LCHI(2,N,xp)*AP.gg0*AP.gS0-6.*LCHI(1,N,xp)*LF*AP.gg0*AP.gS0+3.*LF*LF*AP.gg0*AP.gS0+2.*Bpt1g*(-Hgggq1(N)*LCHI(1,N,xp)
	  +(LCHI(2,N,xp)-LF*LCHI(1,N,xp))*AP.gS0)+Apt1g*(-Hgggq1(N)*LCHI(2,N,xp)+(LCHI(3,N,xp)-LF*LCHI(1,N,xp))*AP.gS0)-2.*LCHI(1,N,xp)*AP.gS1+2.*LF*AP.gS1
	  +(-LF)*(-2.*Hgggq1(N)+(LCHI(1,N,xp)-LF)*AP.gS0)*AP.SS0+(-2.*Hgggq1(N)*LCHI(1,N,xp)+(LCHI(2,N,xp)-LF*LCHI(1,N,xp))*AP.gS0)*AP.SS0
	  +2.*Apt1g*(LF*risdilog-risdiloglogchi)*AP.gS0));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3))
	  match.push_back(ConstResum::as*ConstResum::as*(Hggqq2(N)+(-LF)*AP.gS0*(-2.*Hgggq1(N)+(LCHI(1,N,xp)-LF)*AP.gS0)+AP.gS0*(-2.*Hgggq1(N)*LCHI(1,N,xp)+(LCHI(2,N,xp)-LF*LCHI(1,N,xp))*AP.gS0)));
	  else match.push_back((0.,0.));
	}
      }
    }
  }
  //cout << match[0] << endl;
  return match;
}










































































