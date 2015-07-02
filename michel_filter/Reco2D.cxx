#ifndef RECO2D_CXX
#define RECO2D_CXX

#include "Reco2D.h"

unsigned int Reco2D::nCk( unsigned int n, unsigned int k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  
  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}


Double_t Reco2D::coeff(Double_t k, Double_t N) {
  auto m = (N - 3.0)/2.0;
  return 1.0/pow(2,2*m+1) * (nCk(2*m,m-k+1) - nCk(2*m,m-k-1));
}


std::vector<Double_t> Reco2D::windowed_means(int window_size;, int p_above, int p_below,
					     const std::vector<ahit>    & data,
					     const std::vector<HitIdx_t>& order) {

  auto w = window_size + 2;
  w = (w-1)/2;

  auto num = data.size();

  
  

}


#endif
