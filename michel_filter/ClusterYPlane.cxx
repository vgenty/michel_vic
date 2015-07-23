#ifndef CLUSTERYPLANE_CXX
#define CLUSTERYPLANE_CXX

#include "ClusterYPlane.h"

ClusterYPlane::ClusterYPlane(const ClusterYPlane& other)
{

  _ahits    = other._ahits;
  _clusters = other._clusters;
  
  _nX       = other._nX; //set local nears wire
  _nX       = other._nY; //set local nears time
  _d_cutoff = other._d_cutoff;

  if(_ahits.empty()) throw std::exception();

  //quickly sort the hits based on x location
  sort_hits();
  
  //order the points
  order_points();
  
  //set start and end point
  set_start_end();
  
  //calculate distances
  calculate_distances();
 
}


ClusterYPlane::ClusterYPlane(const std::vector<larlite::hit>&     in_hits,
			     const std::vector<larlite::cluster>& in_clusters,
			     const Double_t near_X, const Double_t near_Y,
			     const Double_t d_cut)
{
  
  //necessary params
  _nX = near_X;
  _nY = near_Y;
  _d_cutoff = d_cut;

  _clusters = in_clusters;
  _ahits.resize(in_hits.size());

  if(_ahits.empty()) throw std::exception();
  
  auto xy = [](const larlite::hit& h,Double_t *a) 
    { a[0] = h.WireID().Wire * 0.3;
      a[1] = h.PeakTime()    * 0.0802814; }; // currently missing trigger offset
  
  for(unsigned int i = 0; i < in_hits.size(); ++i) {
    Double_t a[2]; xy(in_hits[i],a);
    _ahits[i].hit = in_hits[i];
    _ahits[i].vec = TVector2(a);
  }

  //quickly sort the hits based on x location
  sort_hits();
  
  //order the points
  order_points();

  //set start and end point
  set_start_end();
  
  //calculate distances
  calculate_distances();
  
   
}

ClusterYPlane::ClusterYPlane(const std::vector<larlite::hit>&     in_hits,
			     const larlite::cluster&              in_cluster,
			     const Double_t near_X,
			     const Double_t near_Y,
			     const Double_t d_cut)
{ 
  //necessary params
  _nX = near_X;
  _nY = near_Y;
  _d_cutoff = d_cut;
  
  _clusters.push_back(in_cluster);
  _ahits.resize(in_hits.size());
  if(_ahits.empty()) throw std::exception();
  
  //Lambda encapsulates xy conversion

  auto xy = [](const larlite::hit& h,Double_t *a) 
    { a[0] = h.WireID().Wire * 0.3;
      a[1] = h.PeakTime() * 0.0802814; }; // no trigger offset
  
  for(unsigned int i = 0; i < in_hits.size(); ++i) {
    Double_t a[2]; xy(in_hits[i],a);
    _ahits[i].hit = in_hits[i];
    _ahits[i].vec = TVector2(a);
  }
  
  
    //quickly sort the hits based on x location
  sort_hits();
  
  //order the points
  order_points();

  //set start and end point
  set_start_end();
  
  //calculate distances
  calculate_distances();
  
}

void ClusterYPlane::calculate_distances() {
  
  auto tot_dist = double{0.0};
  
  _s .push_back(0.0);
  _ds.push_back(0.0);

  for(size_t u = 0; u < _ordered_pts.size() - 1; ++u) {
    auto zz = distance(_ahits[_ordered_pts[ u ]].vec,
		       _ahits[_ordered_pts[u+1]].vec);
    _ds.push_back(zz);

    tot_dist += zz;

    _s.push_back(tot_dist);
  }
  
  
}

void ClusterYPlane::set_start_end() {
  if(_ordered_pts.empty()) {
    std::cout<<"SHIT: "<<_ordered_pts.size()<<"/"<<_ahits.size()<<std::endl;
    throw std::exception();
    return ;
  }
  _start = _ahits[_ordered_pts.front()];
  _end   = _ahits[_ordered_pts.back() ];
}


void ClusterYPlane::order_points() {
  // try to order the points starting from left to right (lower wire -> higher)
  // then try from right to left (higher -> lower)
  // take whichever cluster is larger...
  
  //feed do_ordering the start index inside of ahits

  if(_ahits.empty()) throw std::exception();
  auto right_to_left = do_ordering(_ahits.size() - 1); // right to left
  auto left_to_right = do_ordering(0);                 // left  to right
  
  
  if(right_to_left.size() > left_to_right.size())
    _ordered_pts = right_to_left;
  else
    _ordered_pts = left_to_right;

  std::cout<<_ordered_pts.size()<<"/"<<_ahits.size()<<std::endl;
}


std::vector<HitIdx_t> ClusterYPlane::do_ordering(const size_t start_idx) {
  
  std::vector<HitIdx_t> all_pts(_ahits.size() - 1,0);
  for(size_t b = 0 ; b < all_pts.size(); ++b)
    all_pts[b] = b+1;

  std::vector<HitIdx_t> the_order;
  the_order.reserve(_ahits.size());
  the_order.push_back(start_idx);

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
       zz = distance(_ahits[the_order[cnt]].vec,_ahits[*itr].vec);
       
       //       if(zz < closest && zz < 0.3*20) { //hard cutoff here to avoid delta ray
       if(zz < closest && zz < _d_cutoff) { //hard cutoff here to avoid delta ray
	 idxholder = itr;
	 closest   = zz;
	 j = 1;
       }
     }
     
     if(j) {
       auto o = *idxholder;
       the_order.push_back(o);
       all_pts.erase(idxholder); //this is scary, erase does idxholder++ !
       closest = 9999.9;
       zz      = 0.0;
       cnt    += 1;
     }
     
     if(all_pts.size() == 0 || j == 0)
       aho = false;
     
     j = 0;
   }
   
   return the_order;
}

Double_t ClusterYPlane::distance(const TVector2& a,
				 const TVector2& b) {
  
  return ((a - b).Mod());
}

ClusterYPlane ClusterYPlane::operator+(const ClusterYPlane& other) {
  
  // (this->_ahits).insert   ( this->_ahits.end(), 
  // 			    other._ahits.begin(), other._ahits.end() );
  // (this->_clusters).insert( this->_clusters.end(), 
  // 			    other._clusters.begin(), other._clusters.end() );
  //scary
  
  std::vector<larlite::hit>     send_hits;
  std::vector<larlite::cluster> send_clusters;

  send_hits.reserve(this->_ahits.size() + other._ahits.size());
  send_clusters.reserve(this->_clusters.size() + other._clusters.size());
  
  for(auto const& h : this->_ahits)
    send_hits.push_back(h.hit);
  for(auto const& h : other._ahits)
    send_hits.push_back(h.hit);

  for(auto const& c : this->_clusters)
    send_clusters.push_back(c);
  for(auto const& c : other._clusters)
    send_clusters.push_back(c);
  
  ClusterYPlane merged(send_hits,send_clusters,_nX,_nY,_d_cutoff);
  // std::cout << "found a baby ... \n";
  // baby.dump();
  return merged;
  
}

bool ClusterYPlane::touching(const ClusterYPlane& other) {
  
  if(near(this->_start.vec,other._end.vec) ||
     near(this->_end.vec  ,other._start.vec))
    return true;
  
  
  return false;
}

bool ClusterYPlane::near(const TVector2& a, const TVector2& b) {
  
  if(abs(a.X() - b.X()) < _nX &&
     abs(a.Y() - b.Y()) < _nY) //hard cutoff, user should set not me :(
    return true;
  return false;
  
}

size_t ClusterYPlane::find_closest_hit(const TVector2& point) {

  auto dist = 999.0;
  auto idx  = size_t{0};
  
  for(size_t i = 0; i < _ordered_pts.size(); ++i) {
    auto d = distance(_ahits[_ordered_pts[i]].vec,point);
    if( d < dist) { dist = d; idx = i; }
  } 
  
  _michel_location = idx;
  _michel_dist     = dist;
  return idx;
}

int ClusterYPlane::match(const TVector2& michel_loc) {
  
  ///in cosmics
  if(_has_michel)
    return 2; //2 I already have a michel...
  
  auto dist = 999.0;
  auto idx  = size_t{0};
  
  for(size_t i = 0; i < _ordered_pts.size(); ++i) {
    auto d = distance(_ahits[_ordered_pts[i]].vec,michel_loc);
    if( d < dist) { dist = d; idx = i; }
  } 
  
  //distance to this michel is dist... , what is the cut off for this??
  if(dist < 4 ) { //? is this reasonable...
    _michel_location = idx;
    _michel_dist     = dist;
    _has_michel      = true;
  }
  else { return 0; } //michel isn't right for me :(
  
  return 1; //1 is matched with accepted parameters
}

void ClusterYPlane::dump() {
  std::cout << "\t\n==start dump==\n";
  std::cout << "A cluster with " << _ahits.size() << " hits in " << _clusters.size() << "\n";
  std::cout << "The start point is at (" << _start.vec.X() << "," << _start.vec.Y() << ")\n";
  std::cout << "The end point is at ("   << _end.vec.X() << ","   << _end.vec.Y()   << ")\n";
  std::cout << _ordered_pts.size() << " of the hits are ordered and nearby ";
  std::cout << "I may or may not have a michel, do i? " << _has_michel << "\n";
  if(_has_michel) 
    std::cout << "ok I do, it is at a distance : " << _michel_dist << " from one of the hits in my cham\n";
  std::cout << "\t\n==end dump==\n";
}
#endif

