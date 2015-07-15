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
#include <cmath>

//Vic
#include "ClusterYPlane.h"

#include "ChiFit.h"

#include "TMinuit.h"
#include "TFitter.h"
#include "TF1.h"
#include "TGraphErrors.h"

// void mywrapper(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//   Reco2D::myfittingfunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
// }

class Reco2D{
  
public:


  Reco2D(){ _chifit = &ChiFit::getInstance();  }

  ~Reco2D(){}


  ChiFit *_chifit;
  Int_t _iNum;
  Double_t *_x;
  Double_t *_y;
  Double_t *_errory;

  Double_t fit_function(float x,Double_t *par);
  
  void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);



  //everything is public really would just like to hold reco methods here...
  
  
  //for fitting TMinuit requires globally scoped BS
  //end TMINUIT


  void do_chi(ClusterYPlane*& c,
	      Int_t window_size);
  
  Double_t coeff(Double_t k, Double_t N);
  unsigned int nCk( unsigned int n, unsigned int k );
  
  //why is template dying?
  //template<typename T>
  std::vector<std::vector<size_t> > get_windows(const std::vector<size_t>& the_thing, int window_size);

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

  size_t find_max(const std::vector<Double_t>& data);
  size_t find_min(const std::vector<Double_t>& data);
  void tag_michel(ClusterYPlane*& c, //for now this DOES have 1 michel b/c of filter
		  size_t idx,      // of chosen in michel in orderd_pts
		  bool forward,    //higher/lower in orderedpts
		  const larlite::event_hit *evt_hits,  //all the hits
		  Double_t _min_rad);

  void tag_muon(ClusterYPlane*& c, //for now this DOES have 1 michel b/c of filter
		  size_t idx,      // of chosen in michel in orderd_pts
		  bool forward,    //higher/lower in orderedpts
		  const larlite::event_hit *evt_hits); //all the hits

  //inline methods
  inline void cut(std::vector<Double_t>& data,
		  double frac, bool above);
  
  inline Double_t calc_mean(std::vector<Double_t> &data);
  inline Double_t distance(const larlite::hit& a, const larlite::hit& b);


  };

#endif
/** @} */ // end of doxygen group 

