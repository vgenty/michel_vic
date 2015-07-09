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
//I just realized we could use TSpectrum::SmoothMarkov as well...
//you can also use TSpectrum to estimate the the pointwise background, which may be useful
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

  
  auto candidate_loc     = find_max(q);
  auto dqdscandidate_loc = find_min(dqds); 
  
  std::cout << "candidate_loc " << candidate_loc << std::endl;
  //if the dqdscandidate_loc is within 20 of the candidate_loc fine...
  if(abs(dqdscandidate_loc - candidate_loc) < 20)
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

  if(window_size == 0)
    window_size++; //this only works because I pop...
  
  std::cout << "window size is... " <<  window_size << std::endl;
  std::cout << "mean_michel_vtx is.. " << mean_michel_vtx << "\n";
 
  
  for(int window = mean_michel_vtx - window_size;
      window < mean_michel_vtx + window_size; ++window){
    auto c = h[o[window]].hit.Integral();
    if(c > k) {
      k = c; idx = window;
    }
    
  }
  
  return idx;

}


//these following two methods need to be templated
size_t Reco2D::find_max(const std::vector<Double_t>& data) {

  auto the_max = 0.0;
  size_t candidate_loc = 9999;
  
  for(size_t i = 0; i < data.size(); ++i) {
    if(data[i] > the_max) {
      the_max = data[i]; candidate_loc = i;
    }
  }
   
  return candidate_loc;
}

size_t Reco2D::find_min(const std::vector<Double_t>& data) {
  
  //get lowest dqds
  auto the_min = 0.0;
  size_t dqdscandidate_loc = 9999;
  
  for(size_t i = 0; i < data.size(); ++i) {
    if(data[i] < the_min) {
      the_min = data[i]; dqdscandidate_loc = i;
    }
  }

  return dqdscandidate_loc;

}

void Reco2D::tag_michel(ClusterYPlane*& c, //for now this DOES have 1 michel b/c of filter
			size_t idx,      // of chosen in michel in orderd_pts
			bool forward,    //higher/lower in orderedpts
			const larlite::event_hit *evt_hits) { //all the hits
  auto E = double{0.0};
  auto L = double{0.0};
  
  //determine the radius of the circle here
  Double_t radius = 0.0;
  
  //for(const auto& i : c->_ordered_pts) {
  
  std::vector<size_t>       michel_idxs;
  std::vector<larlite::hit> michel_hits;

  for(size_t i = 0; i < c->_ordered_pts.size(); ++i)
    //c->_ahits[c->_ordered_pts[i]].hit.Integral();
    
    if(forward)  {
      if(i >= idx) {
	radius += c->_ds[i];
	michel_idxs.push_back(c->_ordered_pts[i]);
      }
    }
    else {
      if( i <= idx ) {
	radius += c->_ds[i];
	michel_idxs.push_back(c->_ordered_pts[i]);
      }
    }
  
  bool there = false;
  std::cout << "michel idxs . size() : " << michel_idxs.size() << "\n";
  std::cout << "_ahits . size() : "      << c->_ahits.size() << "\n";
  std::cout << "_ordered . size() : "    << c->_ordered_pts.size() << "\n";
  
  int k = 0;
  //  std::vector<larlite::hit> cluster_hits;
  std::vector<size_t> cluster_hits;
  std::cout << "michelidxs {";
  for(auto o : michel_idxs)
    std::cout << o << ",";
  std::cout << "}\n";
  std::cout << "pushing back: {";

  for(size_t i = 0; i < c->_ahits.size(); ++i) {
    for(size_t j = 0; j < michel_idxs.size(); ++j) {
      if( i == michel_idxs[j]) {
	std::cout << " *  ";
	k = 1;
	there = true;
	break;
      }
    }
    //std::cout << " " << k << " ";
    if(k  == 0) {
      std::cout << i << ",";
    
      cluster_hits.push_back(i);
    }
    //there = false;
    k = 0;
  }
  std::cout << "}";
  std::cout << "cluster_hits . size() : " << cluster_hits.size() << "\n";
  std::cout << "radius... " << radius << "\n";
  //michel_idxs.clear();
  
  //erase the michels out of here
  // if(forward)
  //   ahits.erase(ahits.begin() + idx, ahits.end());
  // else
  //   ahits.erase(ahits.begin(), ahits.begin() + idx + 1);
  std::cout << "num hits in this event... " << evt_hits->size() << "\n";
  
  //there = false;
  k = 0;
  std::cout << "idx.. " << idx << "\n";
  
  for(const auto& m : michel_idxs)
    std::cout << "the distance between me and michels is "
	      << distance(c->_ahits[m].hit,c->_ahits[c->_ordered_pts[idx]].hit) << "\n";

  // for(const auto& m : cluster_hits)
  //   std::cout << "the distance between cluster_hits and michels is "
  // 	      << distance(m,c->_ahits[c->_ordered_pts[idx]].hit) << "\n";
  
  
  int counter = 0;
  int total_view_hits = 0;
  for(const auto& ehit : *evt_hits) {   // loop over all the hits
    if(ehit.View() == 2){ // look at Y plane only
      total_view_hits++;
      //for(const auto& ahitz : cluster_hits) { //loop over all cluster hits that are not michels;
      for(const auto& a : cluster_hits) { //loop over all cluster hits that are not michels;
	auto ahitz = c->_ahits[a].hit;
	//if ((ehit < ahitz) == false)  { //whoa does this fail???
	if ((ehit == ahitz))  { //whoa does this fail???
	// if(ehit.StartTick() == ahitz.StartTick() &&
	//    ehit.EndTick()   == ahitz.EndTick()) {
	  //there = true;
	  k = 1;
	  std::cout << "matched " << a << ", ";
	  counter++;
	  // std::cout << "kaka the distance is... " 
	  // 	    << distance(ehit,c->_ahits[c->_ordered_pts[idx]].hit) << "\n";
	  
	  if(distance(ehit,c->_ahits[c->_ordered_pts[idx]].hit) == 0)
	    std::cout << "culprit " << a << std::endl;
	  
	  std::cout << std::endl;
	  break;
	}
      }
      
      if( k == 0 )
	std::cout << " the distance is... " 
		  << distance(ehit,c->_ahits[c->_ordered_pts[idx]].hit) << "\n";
      if( k == 0 && (distance(ehit,c->_ahits[c->_ordered_pts[idx]].hit) <= radius) )
  	michel_hits.push_back(ehit);
      
    }//end Y plane
    //there = false;
    k =0;
  } //end all hits
  std::cout << "counter : " << counter << "\n";
  std::cout << "total_view_hits : " << total_view_hits << "\n";
  std::cout << "b" << std::endl;
  L = radius;
  //as a check lets loop through the michel hits and check to see that the idx is in there...
  
  
  for(const auto& h : michel_hits) {
    E += h.SummedADC();
  }
  
  
  // for(const auto& id : michel_idxs) 
  //   E += c->_ahits[id].hit.Integral();
  
  //c->_michel = new Michel(E,L,c->_ahits[c->_ordered_pts[idx]].vec,michel_idxs.size());

  std::cout << "E " << E << std::endl;
  c->_michel = new Michel(E,L,c->_ahits[c->_ordered_pts[idx]].vec,michel_hits.size());
  
  c->_michel->dump();

}


inline Double_t Reco2D::distance(const larlite::hit& a, const larlite::hit& b) {
  
  auto x1 = a.WireID().Wire * 0.3;
  auto x2 = a.PeakTime()    * 0.0802814; 

  auto y1 = b.WireID().Wire * 0.3;
  auto y2 = b.PeakTime()    * 0.0802814; 
  
  return sqrt((x1-y1)*(x1-y1) + (x2-y2)*(x2-y2));
}

#endif
