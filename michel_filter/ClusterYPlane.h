#ifndef CLUSTERYPLANE_H
#define CLUSTERYPLANE_H

//C++
#include <iostream>

//larLite
#include "DataFormat/cluster.h"

//Particles
#include "Michel.h"
#include "Muon.h"

typedef size_t HitIdx_t;

class ClusterYPlane{

public:

  ClusterYPlane(){}
  
  /// ctors
  
  ClusterYPlane(const std::vector<larlite::hit>&     in_hits,
		const std::vector<larlite::cluster>& in_clusters,
		const Double_t near_X, const Double_t near_Y,
		const Double_t d_cut);
  
  ClusterYPlane(const std::vector<larlite::hit>&     in_hits,
		const larlite::cluster&              in_cluster,
		const Double_t near_X, const Double_t near_Y,
		const Double_t d_cut);
  
  /// Default destructor
  virtual ~ClusterYPlane(){}
  
  //copy ctor THIS must Go under destructor for some reason.
  ClusterYPlane(const ClusterYPlane& other);
  
  //everything is public save from getters.
  //private: 
  
  /////////RECO VARIABLES////////////
  bool _has_michel = false;
  
  size_t   _michel_location;
  Double_t _michel_dist;
  
  std::vector<ahit>             _ahits;
  std::vector<larlite::cluster> _clusters;
  
  ahit _start;
  ahit _end;
  
  Double_t _nX;
  Double_t _nY;
  Double_t _d_cutoff;

  //all of these "ordered" and will have same size at _ordered_pts
  std::vector<HitIdx_t>  _ordered_pts;
  std::vector<Double_t>  _ds;
  std::vector<Double_t>  _s;
  std::vector<Double_t>  _chi2;
  std::vector<Double_t>  _t_means;
  std::vector<Double_t>  _t_dqds;

  //Particles we will create based on algo
  Michel _michel;
  Muon   _muon;
  
  
  //Setters caleld by constructor only
  void set_start_end();
  void order_points();
  void calculate_distances();
  
  
  Double_t distance(const TVector2& a,
		    const TVector2& b); //Get distance between two TVector2

  std::vector<HitIdx_t> do_ordering(const size_t start_idx); //from the given index, order the points while stepping around.

  bool near    (const TVector2& a, const TVector2& b); //there two points are "near" based on given criteria
  bool touching(const ClusterYPlane& other);           //used in merging, tell me if this cluster is touching another

  int match(const TVector2& michel_loc); //Used in the cosmic ana. Does't do anything yet.

  size_t find_closest_hit(const TVector2& point); //Return closest his in ordered_pts to point
  
  //Inline Methods
  inline void sort_hits() {
    std::sort(_ahits.begin(), _ahits.end(), 
	      [](const ahit & a, const ahit & b) -> bool
	      { 
		return a.vec.X() < b.vec.X();
	      });
  };
  
  void dump(); //Print out some class variables
  

  
  //Operator methods
  ClusterYPlane operator+(const ClusterYPlane& other);
  
};

#endif
/** @} */ // end of doxygen group 

