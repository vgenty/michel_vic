#ifndef CHIFIT_H
#define CHIFIT_H

#include "TROOT.h"
class ChiFit
{
public:
  static ChiFit& getInstance()
  {
    static ChiFit    instance; // Guaranteed to be destroyed.
    // Instantiated on first use.
    return instance;
  }

  
  Double_t fit_function(float x,Double_t *par);
  void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
    
  Int_t _iNum;
  Double_t *_x;
  Double_t *_y;
  Double_t *_errory;
  
  
private:
  ChiFit() {};                   // Constructor? (the {} brackets) are needed here.

  // C++ 11
  // =======
  // We can use the better technique of deleting the methods
  // we don't want.
  ChiFit(ChiFit const&)          = delete;
  void operator=(ChiFit const&)  = delete;



};

#endif
