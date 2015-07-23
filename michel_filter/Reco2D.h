#ifndef RECO2D_H
#define RECO2D_H

//C++
#include <cmath>

//Vic
#include "ClusterYPlane.h"

//ROOT
#include "TF1.h"
#include "TGraphErrors.h"

class Reco2D{
  
public:


  Reco2D(){} // _chifit = &ChiFit::getInstance(); 

  ~Reco2D(){}
  
  
  //Chi variables
  void do_chi(ClusterYPlane& c, const Int_t window_size, size_t ref_index=0);
  std::vector<int> chi_max_pos(const ClusterYPlane& c,const int num_maxs);
  

  //Various window calculation utilities
  Double_t     coeff(Double_t k, Double_t N);
  unsigned int nCk  (unsigned int n, unsigned int k );
  
  template<typename T>
  std::vector<std::vector<T> > get_windows   (const std::vector<T>& the_thing, const int window_size);
  std::vector<Double_t>        windowed_means(int window_size, Double_t p_above, Double_t p_below,
					      const std::vector<ahit>    & data,
					      const std::vector<HitIdx_t>& order);
  Double_t smooth_derive(const std::vector<Double_t> f, 
			 const std::vector<Double_t> x, 
			 Int_t N);
  
  
  //Calculation of michel start point functions
  std::pair<size_t,size_t> DetEVtx(const std::vector<Double_t>& q,
				   const std::vector<Double_t>& dqds);
  
  size_t REALDetEVtx(const std::vector<ahit>& h,
		     const std::vector<HitIdx_t>& o,
		     size_t mean_michel_vtx);


  
  size_t find_max(const std::vector<Double_t>& data);
  int find_max(const std::vector<Double_t>& data,const std::vector<int>& ref);



  size_t find_min(const std::vector<Double_t>& data);

  Double_t tag_michel(ClusterYPlane& c,                 
		      size_t idx,                        
		      bool forward,                       
		      const larlite::event_hit *evt_hits,  
		      Double_t _min_rad);
  
  void tag_muon(ClusterYPlane& c,                    
		size_t idx,                          
		bool forward,                        
		const larlite::event_hit *evt_hits); 

  
  //For finding peaks (max/min) in data vector
  std::vector<int> find_max_pos     (const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);
  std::vector<int> find_min_pos     (const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);
  std::vector<int> Reconstruct_Maxes(const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);
  std::vector<int> Reconstruct_Mins (const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);

  std::vector<int> scan_forward  (const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);
  std::vector<int> scan_backward (const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);
  std::vector<int> scan_forward2 (const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);
  std::vector<int> scan_backward2(const std::vector<Double_t>& data, float rise_edge, float fall_edge, float threshold, float ped_mean, float ped_rms);


  std::pair<float,float> PedEstimate(const std::vector<Double_t>& data, bool start, int window, float cutoff);
  std::pair<float,float> getrms     (const std::vector<Double_t>& data, int k,      int m,      int window);
  
  int find_max_peak(const std::vector<Double_t>& data, int istart, int iend);
  int find_min_peak(const std::vector<Double_t>& data, int istart, int iend);
  
  
  //Inline methods
  inline void cut(std::vector<Double_t>& data,
		  double frac, bool above);
  
  inline Double_t calc_mean(std::vector<Double_t> &data);
  inline Double_t distance(const larlite::hit& a, const larlite::hit& b);
  


  };

#endif

// #include "ChiFit.h"

// #include "TMinuit.h"
// #include "TFitter.h"
// #include "TF1.h"
// #include "TGraphErrors.h"

// void mywrapper(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
//   Reco2D::myfittingfunction(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
// }



// ChiFit *_chifit;
// Int_t _iNum;
// Double_t *_x;
// Double_t *_y;
// Double_t *_errory;

// Double_t fit_function(float x,Double_t *par);
  
// void calc_chi_square(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  
//everything is public really would just like to hold reco methods here...
  
//for fitting TMinuit requires globally scoped bullshit why even bother
//end TMINUIT
  

