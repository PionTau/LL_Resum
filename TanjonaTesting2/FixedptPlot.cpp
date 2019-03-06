#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <LHAPDF/LHAPDF.h>
#include <vector>
#include <string>
#include "CombResum.h"


using namespace std;
using namespace LHAPDF;

int main(){
  long double CMS=13000.;
  long double ptstart=30.;
  long double ptend=250.;
  long double mH=125.09;
  int NUM=100;
  long double rate=(ptend-ptstart)/((long double) NUM);
  stringstream name;
  name << "graph/FixedptPlot_" << CMS/1000. << "_prova.dat";
  ofstream OUT((name.str()).c_str());
  std::vector<long double> xp;
  for (int i=0;i<NUM;i++){
    xp.push_back(std::pow((ptstart+i*rate)/mH,2));
  }
  std::vector<long double> sigma,sigma2;
  CombResum Final(2.,2.,0.,"PDF4LHC15_nnlo_100",true,mH,3.,5.,mH/2.,mH/2.);
  sigma=Final.ResummedCrossSection(CMS,xp,1);
  sigma2=Final.ResummedCrossSection(CMS,xp,2);
  for (int i=0;i<xp.size();i++){
    OUT << CMS << "\t"<< ptstart+i*rate << "\t" << sigma[i] << "\t" << sigma2[i] << endl;
  }
  OUT.close();
  
  
}