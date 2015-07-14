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
    void set_min_cluster_size(const int i)    { _min_cluster_size = i; }
    void set_n_window_size(const int i)       { _n_window_size    = i; }
    void set_window_cutoff(const Double_t i)  { _window_cutoff    = i; }
    void set_truncated_shave(const int i)     { _truncated_shave  = i; }
    void set_min_rad(const Double_t i) {_min_michel_rad = i;}
    
  private:
   
    
    ::btutil::MCMatchAlg fBTAlg;
    
    std::vector<ClusterYPlane*> _clusters; 
    std::string _cluster_producer;
    
    //Variables going into tree...
    TTree *_output_tree;
    
    Double_t _tX;
    Double_t _tY;
    Double_t _rX;
    Double_t _rY;

    std::vector<Double_t> _ahits_X_copy;
    std::vector<Double_t> _ahits_Y_copy;
    std::vector<Double_t> _charges_copy;
    std::vector<size_t>   _ordered_pts_copy;
    std::vector<Double_t> _mean_charges_copy;
    std::vector<Double_t> _dqds_copy;
    std::vector<Double_t> _s_copy;

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
    
    //simchannel
    
    Double_t _simch_michel_true_shower_E     ;
    Double_t _simch_michel_false_shower_E   ;
    
    Double_t _simch_plane_true_shower_E     ;
    Double_t _simch_plane_false_shower_E    ;
    
    Double_t _simch_ordered_true_shower_E   ;
    Double_t _simch_ordered_false_shower_E  ;
    
    Double_t _simch_cluster_true_shower_E   ;
    Double_t _simch_cluster_false_shower_E  ;
    
    
    

    Int_t _evt = 0;

    //THE Reco object
    Reco2D* r2d;
    
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

   
    
    
    
    //General Methods
    bool convert_2d(const event_hit     *evt_hits,
		    const event_cluster *evt_clusters,
		    const event_ass     *evt_ass_data);
    
    bool find_projected_start(TVector2*& p, 
			      const event_mcshower* evt_mcshower);
    
    void check_cluster_boundaries();
    
    
    std::pair<Double_t,Double_t> get_summed_mcshower_other(const ::btutil::MCBTAlg& aho,
							   const std::vector<larlite::hit>& hits,
							   bool plane2hits);
    
    void clear_all();
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
