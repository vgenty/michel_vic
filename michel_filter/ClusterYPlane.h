/**
 * \file ClusterYPlane.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class ClusterYPlane
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/
#ifndef CLUSTERYPLANE_H
#define CLUSTERYPLANE_H

//C++
#include <iostream>

//larLite
#include "DataFormat/cluster.h"
#include "DataFormat/hit.h"

//ROOT
#include "TVector2.h"

//vic
#include "Michel.h"

/**
   \class ClusterYPlane
   User defined class ClusterYPlane ... these comments are used to generate
   doxygen documentation!
 */


typedef size_t HitIdx_t;

struct ahit {
  larlite::hit hit;
  TVector2 *vec;
};


class ClusterYPlane{

public:

  /// ctors
  ClusterYPlane(std::vector<larlite::hit>     in_hits,
		larlite::cluster              in_cluster);  
  ClusterYPlane(std::vector<larlite::hit>     in_hits,
		std::vector<larlite::cluster> in_clusters);  

  
  /// Default destructor
  virtual ~ClusterYPlane(){}
  
  //copy ctor? THIS must Go under destructor wtf??
  ClusterYPlane(const ClusterYPlane& other);
  

  //everything is public : - )
  //private:
  
  
  //class variables
  bool _has_michel = false;
  
  size_t   _michel_location;
  Double_t _michel_dist;
  
  std::vector<ahit>             _ahits;
  std::vector<larlite::cluster> _clusters;
  
  ahit _start;
  ahit _end;
  
  std::vector<HitIdx_t>  _ordered_pts;
  std::vector<Double_t>  _ds;
  std::vector<Double_t>  _s;

  
  //Michel
  Michel *_michel;
  
  
  //Operator methods
  ClusterYPlane operator+(const ClusterYPlane* other);
  
  //Usual methods
  void order_points();
  void dump();
  void calculate_distances();

  Double_t distance(const TVector2* a,
		    const TVector2* b);

  std::vector<HitIdx_t> do_ordering(const size_t start_idx);

  bool near(const TVector2* a, const TVector2* b);
  bool touching(const ClusterYPlane* other);

  int match(const TVector2* michel_loc);
  size_t find_closest_hit(const TVector2* point);
  
  //Inline Methods
  inline void sort_hits() {
    std::sort(_ahits.begin(), _ahits.end(), 
	      [](const ahit & a, const ahit & b) -> bool
	      { 
		return a.vec->X() < b.vec->X();
	      });
  };
  
  
  inline void set_start_end() {
    _start = _ahits[_ordered_pts.front()];
    _end   = _ahits[_ordered_pts.back() ];
  }

};

#endif
/** @} */ // end of doxygen group 
