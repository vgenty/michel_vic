#include "ChiFit.h"


Double_t ChiFit::fit_function(float x,Double_t *par) { return par[0]+par[1]*x; }

void ChiFit::calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  double chisq = 0;
  for (int i=0;i<_iNum; i++) {
    double delta  = (_y[i]-fit_function(_x[i],par))/_errory[i];
    chisq += delta*delta;
  }
  f = chisq;
  return;
}
