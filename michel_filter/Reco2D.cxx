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


Double_t Reco2D::smooth_derive(const std::vector<Double_t> f,
			       const std::vector<Double_t> x,
			       Int_t N) {
 
  // N should def be odd.
  
  auto M   = int{(N - 1)/2};
  auto tot = double{0.0};
  
  for(int k = 0; k < M; ++k)
    tot += coeff(k+1,N) * (f[k+M] - f[M - 1 - k])/(x[k + M] - x[M - 1 - k]) * 2 * (k+1);
  
  return tot;
  
  
}

//This should definitely become a templated function but for now lets get moving who cares
std::vector<Double_t> Reco2D::windowed_means(int window_size, Double_t p_above, Double_t p_below,
					     const std::vector<ahit>    & data,
					     const std::vector<HitIdx_t>& order) {

  // auto printvec = [](const std::vector<double>& d)
  //   { 
  //     std::cout << "{ ";
  //     for(const auto& e : d)
  //     	std::cout << e << " ";
  //     std::cout <<" }\n";
  //   };

  auto w = window_size + 2;
  w = (unsigned int)((w - 1)/2);
  
  auto num = order.size();
  std::vector<Double_t> mean_window;
  std::vector<Double_t> means;
  
  auto charge = [](const ahit& h){ return h.hit.Integral(); };
     
  for(int i = 1; i <= num; ++i) {
    if(i < w) {
      //std::cout << "    a    \n";
      for(int j = 0; j < 2 * (i%w) - 1; ++j)
	means.push_back(charge(data[order[j]]));
    }else if (i > num - w + 1){
      //std::cout << "    b    \n";
      for(int j = num - 2*((num - i)%w)-1 ; j < num; ++j)
	means.push_back(charge(data[order[j]]));
    }    else{
      //std::cout << "    c    \n";
      for(int j = i - w; j < i + w - 1; ++j)
	means.push_back(charge(data[order[j]]));
      
    }
    // std::cout << "\ndata before...:";
    // printvec(means);
    
    if(means.size() > 3) cut(means,p_above,1);
    
    // std::cout << "\ndata after...:";
    // printvec(means);
    
    mean_window.push_back(calc_mean(means));
    means.clear();
  }
  
  //printvec(mean_window);
  //std::cout << "we made it!!!\n";
  return mean_window;
  
}

inline Double_t Reco2D::calc_mean(std::vector<Double_t> &data) {
  auto sum = double{0.0};
  for(const auto& d : data) sum += d;
  return sum/data.size();
}

inline void Reco2D::cut(std::vector<Double_t>& data,
			double frac, bool above) {
  
  auto size = data.size();
  int to_stay = floor(frac*size);
  
  //sort the array based on charge
  std::sort(data.begin(),data.end(),
	    [](const Double_t& a, const Double_t& b) -> bool
	    {
	      return a < b;	      
	    });
  
  //std::cout << "\tfrac : " << frac << " to_stay : " << to_stay << " above : " << above << " size : " << size << "\n";

  if(above) {
    data.erase(data.begin() + to_stay, data.end());
  }
  else {
    //to_stay = 1 - to_stay;
    data.erase(data.begin(),data.begin()+to_stay);
  }
  
}

std::pair<size_t,size_t> Reco2D::DetEVtx(const std::vector<Double_t>& q,
					 const std::vector<Double_t>& dqds) {

  
  auto the_max = 0.0;
  size_t loc = 9999;
  size_t candidate_loc = 9999;
  
  for(size_t i = 0; i < q.size(); ++i) {
    if(q[i] > the_max) {
      the_max = q[i]; candidate_loc = i;
    }
  }
  

  //get lowest dqds
  auto the_min = 0.0;
  size_t dqdsloc = 9999;
  size_t dqdscandidate_loc = 9999;
  
  for(size_t i = 0; i < dqds.size(); ++i) {
    if(dqds[i] < the_min) {
      the_min = dqds[i]; dqdscandidate_loc = i;
    }
  }
  

  //if the dqdscandidate_loc is within 30 of the candidate_loc fine...
  if(abs(dqdscandidate_loc - candidate_loc) < 30)
    return std::make_pair(candidate_loc,dqdscandidate_loc);  
  else
    return std::make_pair(999,999);
  

}


size_t Reco2D::REALDetEVtx(std::vector<ahit> h,
			   std::vector<HitIdx_t> o,
			   size_t mean_michel_vtx) {
   

  auto window_size = 20;
  auto left  = 0;
  auto right = 0;
  
  
  right = o.size() - mean_michel_vtx;
  left  = mean_michel_vtx;
  
  auto smallest = [](const size_t& r, const size_t& l)
    { if (r < l) return r;
      if (l < r) return l;
    };

  
  if(smallest(right,left) < window_size)
    window_size = smallest(right,left);

  
  auto k   = 0.0;
  auto idx = 0;
  // std::cout << "a\n";
  // std::cout << "window start" << mean_michel_vtx - window_size;
  // std::cout << "   window size " << window_size << "\n";

  for(int window = mean_michel_vtx - window_size;
      window < mean_michel_vtx + window_size; ++window){
    
    auto c = h[o[window]].hit.Integral();
    if(c > k) {
      k = c; idx = window;
    }
    
  }
    
  return idx;

}

#endif
