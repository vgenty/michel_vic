/**
 * \file Reco2D.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class Reco2D
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/
#ifndef RECO2D_H
#define RECO2D_H

//C++
#include <iostream>
#include <cmath>

//Vic
#include "ClusterYPlane.h"

class Reco2D{

public:

  Reco2D(){}

  ~Reco2D(){}
  
  //everything is public really would just like to hold reco methods here...
  
  Double_t coeff(Double_t k, Double_t N);
  unsigned int nCk( unsigned int n, unsigned int k );
  
  std::vector<Double_t> windowed_means(int window_size;, int p_above, int p_below,
				       const std::vector<ahit>    & data,
				       const std::vector<HitIdx_t>& order);
  
  
};

#endif
/** @} */ // end of doxygen group 

