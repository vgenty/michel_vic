/**
 * \file Michel2DAna.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class Michel2DAna
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/

#ifndef LARLITE_MICHEL2DANA_H
#define LARLITE_MICHEL2DANA_H

//C++
#include <string>
#include <algorithm>

//larlite 
#include "Analysis/ana_base.h"
#include "DataFormat/event_ass.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/DetectorProperties.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/simch.h"
#include "DataFormat/shower.h" // not sure if needed

#include "MCComp/MCMatchAlg.h"

//ROOT
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include <TStopwatch.h>

//Vic's
#include "Reco2D.h"

namespace larlite {
  /**
     \class Michel2DAna
     User custom analysis class made by me
   */
  class Michel2DAna : public ana_base{
  
  public:

    /// Default constructor
    
    Michel2DAna(std::string cp) : 
      _output_tree     (nullptr),
      _cluster_producer(cp)
    { _name="Michel2DAna"; _fout=0;}


    /// Default destructor
    virtual ~Michel2DAna(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    // Setter methods for various reco parameters as per ATLAS coding guidelines
    void set_min_merge_cluster_size(const int i)    { _min_cluster_size = i; }
    void set_min_proto_cluster_size(const Int_t i)   { _min_proto_cluster = i;}
    void set_n_window_size(const int i)       { _n_window_size    = i; }
    void set_window_cutoff(const Double_t i)  { _window_cutoff    = i; }
    void set_truncated_shave(const int i)     { _truncated_shave  = i; }
    void set_min_rad(const Double_t i)        { _min_michel_rad   = i; }

    void set_near_X(const Double_t i)         { _nearX = i; }
    void set_near_Y(const Double_t i)         { _nearY = i; }
    void set_d_cutoff(const Double_t i)       { _d_cutoff = i;}

    void set_chi2_rise(const Double_t i)     { _chi2_rise = i;}
    void set_chi2_fall(const Double_t i)     { _chi2_fall = i;}
    void set_chi2_threshold(const Double_t i){  _chi2_thresh = i;}
    
    void set_tmean_rise(const Double_t i){_tmean_rise = i;}
    void set_tmean_fall(const Double_t i ){_tmean_fall = i;}
    void set_tmean_threshold(const Double_t i){ _tmean_thresh = i;}

    void set_tdqds_rise(const Double_t i){_tdqds_rise = i;}
    void set_tdqds_fall(const Double_t i){_tdqds_fall = i;}
    void set_tdqds_threshold(const Double_t i){ _tdqds_thresh = i;}
    
    

  private:
   
    
    Int_t win = 0;
    ::btutil::MCMatchAlg fBTAlg;
    
    std::vector<ClusterYPlane> _clusters; 
    std::string _cluster_producer;
    
    
    Double_t _nearX = 0;
    Double_t _nearY = 0;
    Double_t _d_cutoff = 0;
    Int_t    _min_proto_cluster = 4;

    Double_t  _chi2_rise = 5;
    Double_t  _chi2_fall = 5;
    Double_t  _chi2_thresh = 0;

    Double_t  _tmean_rise = 5;
    Double_t  _tmean_fall = 5;
    Double_t  _tmean_thresh = 0;

    Double_t  _tdqds_rise = 5;
    Double_t  _tdqds_fall = 5;
    Double_t  _tdqds_thresh = 0;

    //Variables going into tree...
    TTree *_output_tree;
    
    Double_t _tX;
    Double_t _tY;
    Double_t _rX;
    Double_t _rY;
    
    Bool_t _has_michel;
    Bool_t _tru_id;
    Bool_t _mis_id;

    
    std::vector<Double_t> _ahits_X_copy;
    std::vector<Double_t> _ahits_Y_copy;
    std::vector<Double_t> _charges_copy;
    std::vector<size_t>   _ordered_pts_copy;
    std::vector<Double_t> _mean_charges_copy;
    std::vector<Double_t> _dqds_copy;
    std::vector<Double_t> _s_copy;
    std::vector<Double_t> _chi2_copy;

    Int_t _michel_hits;
    
    //0715
    Double_t _num_hits;
    Double_t _num_wires;
    Double_t _num_hits_p_wire;
    
    Double_t _startX;
    Double_t _startY;
    Double_t _endX;
    Double_t _endY;
    
    Double_t _michel_E;
    Double_t _michel_L;
    Double_t _true_michel_E;
    Double_t _simch_shower_michel_E;
    Double_t _simch_notshower_michel_E;

    Double_t _d_m_h;
    
    Double_t _reco_Q_o_mc_Q;
    Double_t _mcQ_frac;
    Double_t _MeV_scale;
    Double_t _true_michel_Det;

    Double_t _Q_tot_p2;
    Double_t _Q_u_p2;

    std::vector<Double_t> _large_frac_shower_hits_X;
    std::vector<Double_t> _large_frac_shower_hits_Y;
    std::vector<Double_t> _ALL_hits_p2_X;
    std::vector<Double_t> _ALL_hits_p2_Y;

    std::vector<int> _the_chi_max_peak;
    int _num_chi_max_peaks;

    std::vector<int> _the_tmean_max_peak;
    int _num_tmean_max_peaks;
    
    std::vector<int> _the_tdqds_min_peak;
    int _num_tdqds_min_peaks;


    Double_t _lifetime_correction;

    Double_t _matched_max_s;
    Double_t _matched_min_s;

    float _tmean_ped_mean ;
    float _tmean_ped_rms ;
    
    float _tdqds_ped_mean ;
    float _tdqds_ped_rms ;
    
    
    //simchannel
    
    Double_t _simch_michel_true_shower_E    ;
    Double_t _simch_michel_false_shower_E   ;
    
    Double_t _simch_plane_true_shower_E     ;
    Double_t _simch_plane_false_shower_E    ;
    
    Double_t _simch_ordered_true_shower_E   ;
    Double_t _simch_ordered_false_shower_E  ;
    
    Double_t _simch_cluster_true_shower_E   ;
    Double_t _simch_cluster_false_shower_E  ;

    Double_t  _michel_L_true;
    
    int  _min_hits_to_edge;

    Int_t _evt = 0;

    //THE Reco object
    Reco2D r2d;
    
    //Geo variables
    Double_t _time2cm;
    Double_t _wire2cm;
    Double_t _ne2ADC;
    
    //Python setters with various variables
    int      _min_cluster_size = 25;
    int      _n_window_size    = 25;
    Double_t _window_cutoff    = 0.25;
    int      _truncated_shave  = 2;

    Double_t _min_michel_rad = 0.3;
    
    bool determine_forward(bool& ddiirr, 
			   size_t mean_michel_vtx,
			   size_t real_michel_vtx,
			   const ClusterYPlane& c);

     std::pair<int, int> find_match_peaks(const ClusterYPlane& c, std::vector<int>& the_tmean_max_peaks,
					  std::vector<int>& the_tdqds_min_peaks, int range);
   
    
    
    
    //General Methods

    void printvec(std::vector<int> v);
    
    void diagnostic(int min, int max, std::vector<int> v);
      
    bool convert_2d(const event_hit     *evt_hits,
		    const event_cluster *evt_clusters,
		    const event_ass     *evt_ass_data);
    
    bool find_projected_start(TVector2& p, 
			      const event_mcshower* evt_mcshower);
    
    void check_cluster_boundaries();
    
    
    std::pair<Double_t,Double_t> get_summed_mcshower_other(const ::btutil::MCBTAlg& aho,
							   const std::vector<larlite::hit>& hits,
							   bool plane2hits);
    
    void clear_all();

    int  N_to_edge(const ClusterYPlane& c, int tmean_max_ind);

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
