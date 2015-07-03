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
  
  std::vector<Double_t> windowed_means(int window_size, Double_t p_above, Double_t p_below,
				       const std::vector<ahit>    & data,
				       const std::vector<HitIdx_t>& order);

  
  Double_t smooth_derive(const std::vector<Double_t> f,
			 const std::vector<Double_t> x,
			 Int_t N);
			 
			 
  std::pair<size_t,size_t> DetEVtx(const std::vector<Double_t>& q,
				   const std::vector<Double_t>& dqds);
  size_t REALDetEVtx(std::vector<ahit> h,
		     std::vector<HitIdx_t> o,
		     size_t mean_michel_vtx);
  
  //inline methods
  inline void cut(std::vector<Double_t>& data,
		  double frac, bool above);
  
  inline Double_t calc_mean(std::vector<Double_t> &data);
};

#endif
/** @} */ // end of doxygen group 

