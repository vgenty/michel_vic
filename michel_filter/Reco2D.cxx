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
//you can also use TSpectrum to estimate the the pointwise background, which may be useful (as a reminder~)
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
    
    if(means.size() > 3) cut(means,p_above,1);
    
    mean_window.push_back(calc_mean(means));
    means.clear();
  }
  
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
  
  if(above) 
    data.erase(data.begin() + to_stay, data.end());
  else 
    data.erase(data.begin(),data.begin()+to_stay);
  
  
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
  
  // int w_left  = window_size;
  // int w_right = window_size;
  
  auto k   = 0.0;
  auto idx = 0;
  // std::cout << "a\n";
  // std::cout << "window start" << mean_michel_vtx - window_size;
  // std::cout << "   window size " << window_size << "\n";

  if(window_size == 0)
    window_size++; //this only works because I pop...
  
  //this loop is most contentious!!
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


void Reco2D::tag_muon(ClusterYPlane*& c,  
		      size_t idx,        
		      bool forward,      
		      const larlite::event_hit *evt_hits) {
  //kas771
  std::vector<ahit> muon_hits;
  std::vector<double> distance;
   
   if (forward == true){
     for (int i = 0; i < c->_ordered_pts.size();i++){
       if (i < idx){
	 muon_hits.push_back(c->_ahits[c->_ordered_pts[i]]);
	 distance.push_back(c->_s[i]);
       }
     }
     
   }

   else{
     for (int i = 0; i < c->_ordered_pts.size();i++){
       if (i > idx){
	 muon_hits.push_back(c->_ahits[c->_ordered_pts[i]]);
	 distance.push_back(c->_s[i]);
       }
     }
   }
   
   c-> _muon= new Muon(muon_hits, distance);
}



  
void Reco2D::tag_michel(ClusterYPlane*& c, //for now this DOES have 1 michel b/c of filter
			size_t idx,        // of chosen in michel in orderd_pts
			bool forward,      //higher/lower in orderedpts
			const larlite::event_hit *evt_hits) { //all the hits
  std::cout << "idx : " << idx;
  std::cout << "forward: " << forward;
  
  auto E = double{0.0};
  
  //determine the radius of the circle here
  Double_t radius = 0.0;
  
  std::vector<size_t>       michel_idxs;
  std::vector<larlite::hit> michel_hits;
  std::cout << "ordered_pts . size() " << c->_ordered_pts.size() << std::endl;
  std::cout << "michel... idx... {";
  for(size_t i = 0; i < c->_ordered_pts.size(); ++i) {
    if(forward)  {
      if(i >= idx) {
	michel_idxs.push_back(c->_ordered_pts[i]);

	if(i < c->_ordered_pts.size() - 1) radius += c->_ds[i];
	else radius += 0.3;
	std::cout << " i " << c->_ordered_pts[i] << " " ;
      }
    }
    else {
      if( i <= idx ) {
	michel_idxs.push_back(c->_ordered_pts[i]);
	
	if(i < c->_ordered_pts.size() - 1) radius += c->_ds[i];
	else radius += 0.3;
	std::cout << " i " << c->_ordered_pts[i] << " " ;
      }
    }
    //std::cout << "radius : " << radius << " c->_ds[i] : " << c->_ds[i] << "\n";
  }
  std::cout << " }\n";
  if(radius < 0.3) radius = 0.3;
  
  // //remove duplicate hits...............
  // int w = 0;
  // //auto _ahits_copy = c->ahits();
  // std::vector<ahit> _ahits_copy;
  // //std::vector<size_t> duplicates;
  // bool ggg = true;
  // std::vector<ahit>::iterator a1,a2;
  // int total = 0;
  // int called = 0;
  // while(1) {
  //   //for(const auto& ahit1 : c->_ahits) {
  //   for(a1 = c->_ahits.begin(); a1 != c->_ahits.end(); ++a1) {
  //     for(a2 = c->_ahits.begin(); a2 != c->_ahits.end(); ++a2) {
  // 	//for(const auto& ahit2 : c->_ahits) {
  // 	if((*a1).hit == (*a2).hit) {
  // 	  w++;
  // 	}
  //     }
  //     if(w > 1) {
  // 	c->_ahits.erase(a1);
  // 	called++;
  // 	w = 0;
  // 	std::cout << "found a dup\n";
  // 	break;
  //     }
  //     w = 0;
  //   }
  //   if(called == 0)
  //     break;
  //   called = 0;
  // }
  
  bool there = false;
  
  std::vector<larlite::hit> cluster_hits;
  
  for(size_t i = 0; i < c->_ahits.size(); ++i) {
    for(size_t j = 0; j < michel_idxs.size(); ++j) {
      if( i == michel_idxs[j] ) {
	there = true;
	break;
      }						
    }
    if(!there)  {
      cluster_hits.push_back(c->_ahits[i].hit);
    }
    there = false;
  }
  
  // int f = 0;
  // int fdupe = 0;
  
  // for(const auto& ahit1 : c->_ahits) {
  //   for(const auto& ahit2 : c->_ahits) {
  //     if(ahit1.hit == ahit2.hit) {
  // 	f++;
  //     }
  //   }
  //   if( f > 1 )
  //     fdupe++;
  //   f = 0;
  // }
  // std::cout << "fdupes: " << fdupe << "\n";
  
  //std::cout << "}\n";
  
  //remove duplicate hits...............
  //int w = 0;
  // int called = 0;
  //   std::vector<size_t>::iterator a1,a2;
  
  // while(1) {
  //   //for(const auto& ahit1 : cluster_hits) {
  //   for(a1 = cluster_hits.begin(); a1 != cluster_hits.end(); ++a1) {
  //     for(a2 = cluster_hits.begin(); a2 != cluster_hits.end(); ++a2) {
  //   	if((*a1) == (*a2)) {
  // 	  w++;
  // 	}
  //     }
  //     if(w > 1) {
  // 	cluster_hits.erase(a1);
  // 	called++;
  // 	w = 0;
  // 	std::cout << "found a dup\n";
  // 	break;
  //     }
  //     w = 0;
  //   }
  //   if(called == 0)
  //     break;
  //   called = 0;
  // }
  


  std::cout << "cluster hits . size() " << cluster_hits.size() << "\n";
  std::cout << "michel_idx   . size() " << michel_idxs.size() << "\n";
  
  //if(radius < 0.3) radius = 0.3;
  
  //there = false;
  
  // if(michel_idxs.size() == 1)
  //   std::cout << " before any farther ~ dist to idx is ~ " << distance(c->_ahits[michel_idxs[0]].hit,
  // 								       c->_ahits[c->_ordered_pts[idx]].hit)
  // 	      << "\n";
  // std::cout << "cluster hits.sioze() " << cluster_hits.size() << "\n";

  // std::cout << "checking if michell in evt_hits...\n";
  // int counter =0;
  // for(const auto& ehit : *evt_hits) {
  //   if(ehit.View() == 2) {
  //     for(int i = 0 ; i < michel_idxs.size(); ++i) { //loop over all cluster hits that are not michels;
  // 	if(ehit == c->_ahits[michel_idxs[i]].hit) {
  // 	  counter++;
  // 	  std::cout << "found counter: " << counter << " michels at i "<< michel_idxs[i] << std::endl;
  // 	}
  //     }
  //   }
  // }


  ////////check for dupes between cluster_hits and michel
  ////////check for dupes between cluster_hits and itself...


  std::vector<larlite::hit>::iterator a1;
  std::cout << "...checking if michel in cluster_hits..\n";
  
  int  w      = 0;
  bool dup    = true;
  int  ndupes = 0;
  
  while(1) {
    for(a1 = cluster_hits.begin(); a1 != cluster_hits.end(); ++a1) {
      for(int i = 0 ; i < michel_idxs.size(); ++i) { //loop over all cluster hits that are not michels;
	if(*a1 == c->_ahits[michel_idxs[i]].hit) {
	  w++;
	  break;  
	}
      }
      if(w > 0)
	break;
      w = 0;
    }
    if(w > 0) {
      cluster_hits.erase(a1);
      w = 0;
      dup = true;
      ndupes++;
    }
    w = 0;
     
    if(!dup)
      break;
    
    dup = false;
  }
  
  std::cout << "encountered : " << ndupes << "\n";


  std::vector<larlite::hit>::iterator a2, a3;
  std::cout << "...checking if cluster_hits in cluster_hits..\n";
  w = 0;
  dup = false;
  int ndupes2 = 0;
  
  while(1) {
    for(a2 = cluster_hits.begin(); a2 != cluster_hits.end(); ++a2) {
      for(a3 = cluster_hits.begin(); a3 != cluster_hits.end(); ++a3) {
	if(*a2 == *a3) {
	  w++;
	}
	if(w > 1)
	  break;  
      }
      if(w > 1)
	break;
      w = 0;
    }
    if(w > 1) {
      cluster_hits.erase(a2);
      w = 0;
      dup = true;
      ndupes2++;
    }
    w = 0;
    
    if(!dup)
      break;
    
    dup = false;
  }
  
  std::cout << "encountered : " << ndupes + ndupes2<< "\n";
  
  for(const auto& ehit : *evt_hits) {   // loop over all the hits
    if(ehit.View() == 2){ // look at Y plane only
      for(const auto& ahitz : cluster_hits) { //loop over all cluster hits that are not michels;
	//auto ahitz = c->_ahits[eh].hit;
	if (ehit == ahitz)  { //custom hit.h
	  // if (ehit.Channel()   == ahitz.Channel() &&
	  //     ehit.StartTick() == ahitz.StartTick() &&
	  //     ehit.EndTick()   == ahitz.EndTick()  )  { //custom hit.h
	  there = true;
	  break;
	}
      }
      // if(!there)
      // 	std::cout << "not there but distance is..." << distance(ehit,c->_ahits[c->_ordered_pts[idx]].hit) << std::endl;
      if( !there && (distance(ehit,c->_ahits[c->_ordered_pts[idx]].hit) <= radius) )
  	michel_hits.push_back(ehit);
      
    } //end Y plane
    there = false;
  } //end all hits
  
  for(const auto& h : michel_hits)
    E += h.Integral();

  std::cout << "michel_hits.size() : " << michel_hits.size() << std::endl;
  
  
  c->_michel = new Michel(E,radius,c->_ahits[c->_ordered_pts[idx]].vec,
			  michel_hits,michel_hits.size());
  
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
