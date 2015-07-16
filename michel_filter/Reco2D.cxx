#ifndef RECO2D_CXX
#define RECO2D_CXX

#include "Reco2D.h"

// Double_t Reco2D::fit_function(float x,Double_t *par) { return par[0]+par[1]*x; }

// //void Reco2D::calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
// Double_t Reco2D::chi_square();
// {
//   double chisq = 0;
//   for (int i=0;i<_iNum; i++) {
//     double delta  = (_y[i]-fit_function(_x[i],par))/_errory[i];
//     chisq += delta*delta;
//   }
//   f = chisq;
//   return;
// }

void Reco2D::do_chi(ClusterYPlane*& c, //probably void right just so c is updated
		    Int_t window_size)
{
 
  std::vector<Double_t> chi;
  
  
  //prepare MINUIT
  // TMinuit *ptMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
  // ptMinuit->SetPrintLevel(-1);

  // Double_t arglist[10];
  // arglist[0] = -1;
  // Int_t ierflg = 0;

  //go here to find the command you want to use...
  //http://cern-tex.web.cern.ch/cern-tex/minuit/node18.html

  // ptMinuit->Command("SET PRINT", arglist, 1);
  // static Double_t vstart[2] = {1, 1};
  // static Double_t step[2]   = {0.1 , 0.1};
  // ptMinuit->mnparm(0, "a1", vstart[0], step[0],0,0,ierflg);
  // ptMinuit->mnparm(1, "a2", vstart[1], step[1],0,0,ierflg);


  //pointer to member function fucked me here :-(
  //  ptMinuit->SetFCN(calc_chi_square);
  
  // void (Reco2D::*func)(Int_t&, Double_t*, Double_t&, Double_t*, Int_t);
  // func = &Reco2D::calc_chi_square;
  
  // ptMinuit->SetFCN(func);
  
  //https://root.cern.ch/root/roottalk/roottalk04/1099.html
  //ptMinuit->SetFCN(&_chifit->calc_chi_square);

  
  // auto func1 = std::bind(&Reco2D::calc_chi_square, this);
  // ptMinuit->SetFCN(func1);
  
  // auto fff = std::mem_fn(&Reco2D::calc_chi_square);
  // ptMinuit->SetFCN(fff);
  
  
  

  // arglist[0] = 500;
  // arglist[1] = 1.;

  //get windows of ordered_pts

  
  
  // TFitter minuit(2);
  // minuit.SetFCN(calc_chi_square);
  // static Double_t vstart[2] = {1, 1};
  // static Double_t step[2]   = {0.1 , 0.1};
  // minuit.SetParameter(0, "a1", vstart[0], step[0],0,0,ierflg);
  // minuit.SetParameter(1, "a2", vstart[1], step[1],0,0,ierflg);

  // ROOT::Math::Minimizer* min = 
  //   ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");  
  // min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  // min->SetMaxIterations(10000);  // for GSL
  // min->SetTolerance(0.001);
  // min->SetPrintLevel(-1);  

  //ROOT::Fit::Fitter::MinuitFCN_t fcn
  //ROOT::Math::Functor f(&,2);

  // static Double_t vstart[2] = {1, 1};
  // static Double_t step[2]   = {0.1 , 0.1};
  // min->SetVariable(0, "a1", vstart[0], step[0],0,0);
  // min->SetVariable(1, "a2", vstart[1], step[1],0,0);

  
  //fuck me I give up
  TGraphErrors *graph;
  TF1 *tf;

  auto chi_data = get_windows(c->_ordered_pts,window_size);
  
  for(const auto& cd : chi_data) {
    
    const unsigned int SIZE = cd.size();

    if(SIZE == 1) {
      chi.push_back(0.0);
      continue;
    }
    
    Double_t x[SIZE];
    Double_t y[SIZE];
    Double_t xerr[SIZE];
    Double_t yerr[SIZE];
    
    
    // std::cout << "{";
    // for(const auto& b2: cd) {
    //   std::cout << b2 << ",";
    // }
    // std::cout << "}\n";
    
    for(unsigned int i = 0; i < SIZE; ++i) {
      x[i]    = c->_ahits[cd[i]].vec->X();
      y[i]    = c->_ahits[cd[i]].vec->Y();
      xerr[i] = 0.3;
      yerr[i] = 1.0;
    }
    
    // _iNum = SIZE; //should be reference?
    // _x = x; 
    // _y = y; 
    // _errory = yerr;

    
    graph = new TGraphErrors(SIZE,x,y,xerr,yerr);
    tf    = new TF1("aho","[0] + [1]*x");
    tf->SetParameter(0,1);
    tf->SetParameter(1,1);
    graph->Fit(tf,"F 0 N Q");
    // Now ready for minimization step, then get chi
    //ptMinuit->mnexcm("MIGRAD", arglist,2,ierflg);
    //minuit.ExecuteCommand("MIGRAD",0,0)

    // Double_t amin,edm,errdef;
    // Int_t nvpar,nparx,icstat;    
    // ptMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    //minuit.GetStats(amin,edm,errdef,nvpar,nparx,icstat);

    
    double amin = tf->GetChisquare()/(SIZE - 1);
    chi.push_back(amin);
    //std::cout << "SIZE: " << SIZE << " chi: " << amin/(SIZE - 1) << "\n";
    
    delete graph;
    delete tf;
    //ptMinuit->mnrset(1);
    //minuit.GetMinuit->mnrset(1);
  }

  c->_chi2 = chi;
  //delete ptMinuit;
} 

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

template<typename T>
std::vector<std::vector<T> > Reco2D::get_windows(const std::vector<T>& the_thing, int window_size)
{
  
  std::vector<std::vector<T> > data;
  std::vector<T> inner;
  
  auto w = window_size + 2;
  w = (unsigned int)((w - 1)/2);
  auto num = the_thing.size();
  
  for(int i = 1; i <= num; ++i) {
    if(i < w) {
      for(int j = 0; j < 2 * (i%w) - 1; ++j)
	inner.push_back(the_thing[j]);
    }else if (i > num - w + 1){
      for(int j = num - 2*((num - i)%w)-1 ; j < num; ++j)
	inner.push_back(the_thing[j]);
    }else{
      for(int j = i - w; j < i + w - 1; ++j)
	inner.push_back(the_thing[j]);
    }
    data.push_back(inner);
    inner.clear();
  }

  return data;
  
}

//This should definitely become a templated function but for now lets get moving who caresb
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

int Reco2D::find_max(const std::vector<Double_t>& data, const std::vector<int> ref) {

  auto the_max = 0.0;
  size_t candidate_loc = 9999;

  bool found = false;
  int f = 0;
  for(size_t i = 0; i < data.size(); ++i) {
    if( f < ref.size()) {
      for (const auto& idx: ref){
	if (i == idx){
	  found = true;
	  f++;
	  break;
	}
      }
    }
    
    if(!found && data[i] > the_max) {
      the_max = data[i]; candidate_loc = i;
    }
    found = false;
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
			const larlite::event_hit *evt_hits,//all the hits
			 Double_t _min_rad){ 
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
  if(radius < _min_rad) radius = _min_rad;
  
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
  
  //for(size_t i = 0; i < c->_ahits.size(); ++i) { //vic 07142015
  for(size_t i = 0; i < c->_ordered_pts.size(); ++i) {
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
  
  // std::cout << "encountered : " << ndupes << "\n";


  // std::vector<larlite::hit>::iterator a2, a3;
  // std::cout << "...checking if cluster_hits in cluster_hits..\n";
  // w = 0;
  // dup = false;
  // int ndupes2 = 0;
  
  // while(1) {
  //   for(a2 = cluster_hits.begin(); a2 != cluster_hits.end(); ++a2) {
  //     for(a3 = cluster_hits.begin(); a3 != cluster_hits.end(); ++a3) {
  // 	if(*a2 == *a3) {
  // 	  w++;
  // 	}
  // 	if(w > 1)
  // 	  break;  
  //     }
  //     if(w > 1)
  // 	break;
  //     w = 0;
  //   }
  //   if(w > 1) {
  //     cluster_hits.erase(a2);
  //     w = 0;
  //     dup = true;
  //     ndupes2++;
  //   }
  //   w = 0;
    
  //   if(!dup)
  //     break;
    
  //   dup = false;
  // }
  
  // std::cout << "encountered : " << ndupes + ndupes2<< "\n";
  
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

std::vector<int> Reco2D::chi_max_pos(const ClusterYPlane *c,const int num_maxs){
  std::vector<int> the_maxs;
  
  int count = 0;
  while(count< num_maxs){
    the_maxs.push_back(find_max(c-> _chi2, the_maxs));
    count++;
  }
  
  return the_maxs;
}

//this one is to find all the local maximums above a pedestal value
//takes the cluster, whether forward or back, a window size for calculation, and an rms cutoff value

std::vector<int> Reco2D::chi_max_pos(ClusterYPlane *c, bool forward, int window, 
				     float cutoff, float rise_edge, float fall_edge, float threshold){ 
  auto peaks = Reconstruct(c->_chi2, forward, window, cutoff, rise_edge, fall_edge, threshold);
  return peaks;
}
// const std::vector<size_t> Reco2D::chi_max_pos(ClusterYPlane *c, bool forward, int window, 
// 					      float cutoff, float rise_edge, float fall_edge, float threshold) const{
  
//   auto peaks = Reconstruct(c->_chi2, forward, window, cutoff, rise_edge, fall_edge, threshold);
//   return peaks;
// }

int Reco2D::find_peak(const std::vector<Double_t>& data, int istart, int iend)
  {
    auto the_max = double{0.0};
    int cl = 4096;
    
    for(int i = istart; i < iend; ++i) {
      if(data[i] > the_max) { the_max = data[i]; cl = i; }
    }
    
    return cl;
  }

std::vector<int> Reco2D::Reconstruct( const std::vector<Double_t>& chi2, bool forward, 
				      int window, float cutoff, float rise_edge, float fall_edge, float threshold)
{
    std::vector<int> result;
    std::cout << "e1\n";
    auto ped_info = PedEstimate(chi2,forward, window, cutoff);
    std::cout << "e2\n";
    bool found_pulse = false; 
    size_t t = 0;
    
    while (  t < chi2.size() ) {
      if(chi2[t] > (ped_info.first  + rise_edge * ped_info.second)  && 
	chi2[t] > threshold + ped_info.first &&
	 !found_pulse)
	found_pulse = true;
      
      if(found_pulse) {
        size_t  t_end = t;

	while(1) {
	  if(t_end == chi2.size() - 1)
	    break;
	  
	  if(chi2[t_end] <= (ped_info.first + fall_edge * ped_info.second) &&
	     (chi2[t_end] <= threshold + ped_info.first))
	    break;
	  else 
	    ++t_end;
	}
	
	result.push_back(find_peak(chi2, t, t_end));	
	
	while(t < t_end) ++t; //secretly increases t...
	
      }
      ++t;
      found_pulse = false;
    }
    
    return result;
  }



std::pair<float,float> Reco2D::PedEstimate(const std::vector<Double_t>& chi2, bool start, int window, float cutoff) {
  float mean = 0;
  float rms = 0;
  int n = chi2.size();
  int k = 0;
  std::cout << "e4\n";
  //number of points to consider in calculation;
  //int window = 10;
  
  bool below = false;
  std::cout << "e4.1 here we go start is :" << start << "\n";
  //need minimum number to calculate
  if (n >= window) {
    std::cout << "e4.2\n";
    if (start == true){
      std::cout << "e4.3\n";
      while (below == false && k+window < n){
	std::cout << "e4.4\n";
	auto mean_rms = getrms(chi2, k, k+window, window);
	mean = mean_rms.first;
	rms  = mean_rms.second;
	
	if (rms < cutoff){
	  below = true;
	}
	  k++;
      }
    }
    else{
      std::cout << "e4.21\n";
      
      std::cout << "e4.22\n";
      k = n-1;
      while (below == false && k - window > 0){
	std::cout << "e4.23\n";
	auto mean_rms =  getrms(chi2, k-window, k, window);
	mean = mean_rms.first;
	rms  = mean_rms.second;
	std::cout << "e4.24 aho\n";
	if (rms < cutoff){
	  below = true;
	}
	  k--;
      }
    }
  }
  
  std::cout << "e6\n";
  //returns <mean, rms>, or 0,0 if nothing in vector or bad index or not below cutoff
  std::cout << "mean: " << mean << " rms: " << rms <<"\n";
  return std::pair<float,float>(mean,rms);
}

std::pair<float,float> Reco2D::getrms (const std::vector<Double_t>& chi2, int k, int m, int window) {
  float mean = 0;
  float rms = 0;
  for (int i  = k; i < m; i++) mean += chi2.at(i);
  
  mean = mean/window;

  for (int i= k; i < m; i++){
    float diff = chi2.at(i) - mean;
    rms += diff*diff;
  }
  
  rms = sqrt(rms/window);
  return  std::pair<float,float>(mean,rms);
}


#endif
