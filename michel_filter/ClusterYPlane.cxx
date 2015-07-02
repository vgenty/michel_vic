#ifndef CLUSTERYPLANE_CXX
#define CLUSTERYPLANE_CXX

#include "ClusterYPlane.h"

ClusterYPlane::ClusterYPlane(const ClusterYPlane& other)
{


  _ahits    = other._ahits;
  _clusters = other._clusters;
 
  // _start = other._start;
  // _end  = other._end;
  
  // _ordered_pts =other._ordered_pts;
  // _ds = other._ds;
  // _s = other._s;
  //quickly sort the hits based on x location
  sort_hits();
  
  //set start and end point
  set_start_end();
  
  //order the points
  order_points();

}


ClusterYPlane::ClusterYPlane(std::vector<larlite::hit>     in_hits,
			     std::vector<larlite::cluster> in_clusters)
{

  _clusters = in_clusters;
  _ahits.resize(in_hits.size());
  
  auto xy = [](const larlite::hit& h,Double_t *a) 
    { a[0] = h.WireID().Wire * 0.3;
      a[1] = h.PeakTime() * 0.0802814; };
  
  for(unsigned int i = 0; i < in_hits.size(); ++i) {
    Double_t a[2]; xy(in_hits[i],a);
    _ahits[i].hit = in_hits[i];
    _ahits[i].vec = new TVector2(a);
  }
  
  
  //quickly sort the hits based on x location
  sort_hits();
  
  //set start and end point
  set_start_end();
  
  //order the points
  order_points();
  
}

ClusterYPlane::ClusterYPlane(std::vector<larlite::hit>     in_hits,
			     larlite::cluster              in_cluster)
{ 
  
  _clusters.push_back(in_cluster);
  _ahits.resize(in_hits.size());
  
  //Lambda encapsulates xy conversion
  
    // auto xy = [](const larlite::hit& h){ 
    //   Double_t c[2] = { h.WireID().Wire * 0.3,
    // 		       h.PeakTime() * 0.0802814 }; return c; };
  auto xy = [](const larlite::hit& h,Double_t *a) 
    { a[0] = h.WireID().Wire * 0.3;
      a[1] = h.PeakTime() * 0.0802814; };
  
  for(unsigned int i = 0; i < in_hits.size(); ++i) {
    Double_t a[2]; xy(in_hits[i],a);
    _ahits[i].hit = in_hits[i];
    _ahits[i].vec = new TVector2(a);
  }
  
  
    //quickly sort the hits based on x location
  sort_hits();
  
  //set start and end point
  set_start_end();
  
  //order the points
  order_points();
  
}


void ClusterYPlane::order_points() {
  
  std::vector<HitIdx_t> all_pts(_ahits.size() - 1,0);
  for(size_t b = 0 ; b < all_pts.size(); ++b)
    all_pts[b] = b+1;
  
  
  _ordered_pts.push_back(0);
  
  bool aho = true;
  Double_t zz = 0.0;
  size_t j = 0;
  size_t cnt = 0;
  std::vector<HitIdx_t>::iterator idxholder;
  Double_t closest = 9999.9;
  Double_t stot = 0.0;
  
   while(aho) {
  
     for(std::vector<HitIdx_t>::iterator itr = all_pts.begin();
	   itr != all_pts.end(); ++itr) {
       zz = distance(_ahits[_ordered_pts[cnt]].vec,_ahits[*itr].vec);
       // std::cout << zz << "  between  " << _ordered_pts[cnt]
       // 		 << " and " << *itr << "\n";
       if(zz < closest && zz < 0.3*6) {
	 idxholder = itr;
	 closest   = zz;
	 j = 1;
	 //std::cout << "!!best!! " << std::endl; 
       }
     }
     
     
     if(j) {
       
       auto o = *idxholder;
       _ordered_pts.push_back(o);
       all_pts.erase(idxholder); //this is scary, erase does idxholder++ !
       
       _ds.push_back(closest);
       
       stot += closest;
       _s.push_back(stot);
       
       
       std::cout << "Found closest point to "
		 << _ordered_pts[cnt] << " is " 
		  << o << " at closest = "
		  << closest << "\n";
       
       closest = 9999.9;
       zz      = 0.0;
       cnt    += 1;
     }
     
     
     if(all_pts.size() == 0 || j == 0)
       aho = false;
       
       
     j = 0;
   }
  
   // for(const auto& id : _ordered_pts)
   //   std::cout <<  id << " , ";
       
   // std::cout << "\n";
   
   // std::cout << "_ordered_pts.size() : " << _ordered_pts.size() << " "
   // 	     << "_ahits.size()       : " << _ahits.size() << "\n";
   
  
}


Double_t ClusterYPlane::distance(const TVector2* a,
				 const TVector2* b) {
  
  return ((*a - *b).Mod());
}

ClusterYPlane ClusterYPlane::operator+(const ClusterYPlane* other) {
  
  // (this->_ahits).insert   ( this->_ahits.end(), 
  // 			    other->_ahits.begin(), other->_ahits.end() );
  // (this->_clusters).insert( this->_clusters.end(), 
  // 			    other->_clusters.begin(), other->_clusters.end() );
  //scary
  
  
  std::vector<larlite::hit>     send_hits;
  std::vector<larlite::cluster> send_clusters;
  
  for(auto const& h : this->_ahits)
    send_hits.push_back(h.hit);
  for(auto const& h : other->_ahits)
    send_hits.push_back(h.hit);

  for(auto const& c : this->_clusters)
    send_clusters.push_back(c);
  for(auto const& c : other->_clusters)
    send_clusters.push_back(c);
  
  ClusterYPlane baby(send_hits,send_clusters);
  // std::cout << "found a baby ... \n";
  // baby.dump();
  return baby;
  
}

bool ClusterYPlane::touching(const ClusterYPlane* other) {
  // std::cout << "in touching...\n";

  // std::cout << "(" << this->_start.vec->X() 
  // 	    << "," << this->_start.vec->Y()
  // 	    << ")\n";

  // std::cout << "(" << other->_start.vec->X() 
  // 	    << "," << other->_start.vec->Y()
  // 	    << ")\n";
  
  if(near(this->_start.vec,other->_end.vec) ||
     near(this->_end.vec  ,other->_start.vec))
    return true;

  
  return false;
}

bool ClusterYPlane::near(const TVector2* a, const TVector2* b) {
  //std::cout << "in near...\n";
  if(abs(a->X() - b->X()) < 1.0 &&
     abs(a->Y() - b->Y()) < 3.75) //hard cutoff, user should set not me :(
    return true;
  return false;
  
}

void ClusterYPlane::dump() {
  // std::vector<ahit>             _ahits;
  // std::vector<larlite::cluster> _clusters;
  
  // ahit _start;
  // ahit _end;
  
  // std::vector<HitIdx_t>  _ordered_pts;
  // std::vector<Double_t>  _ds;
  // std::vector<Double_t>  _s;
  
  std::cout << "\n==start dump==\n";
  std::cout << "A cluster with " << _ahits.size() << " hits in " << _clusters.size() << "\n";
  std::cout << "The start point is at (" << _start.vec->X() << "," << _start.vec->Y() << ")\n";
  std::cout << "The end point is at ("   << _end.vec->X() << ","   << _end.vec->Y()   << ")\n";
  std::cout << _ordered_pts.size() << " of the hits are ordered and nearby ";
  std::cout << "\n==sned dump==\n";
  
  

}
#endif
