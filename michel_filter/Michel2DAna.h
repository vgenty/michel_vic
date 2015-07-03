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
#include "DataFormat/mcshower.h"

//ROOT
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"

//Vic's
#include "ClusterYPlane.h"
#include "Reco2D.h"

namespace larlite {
  /**
     \class Michel2DAna
     User custom analysis class made by SHELL_USER_NAME
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

    
  private:
    
    TTree *_output_tree;
    TBranch *_aho;
    
    std::vector<ClusterYPlane*> _clusters; 
    
    std::string _cluster_producer;
    
    Double_t _tX;
    Double_t _tY;
    Double_t _rX;
    Double_t _rY;


    Reco2D* r2d;
    
    //Temporary stuff...
    // TGraph *tgMeans;
    // TGraph *tgPoints;
    // TGraph *tgdqdx;
    
    //geo
    Double_t _time2cm;
    Double_t _wire2cm;
    
    //crap
    Double_t _sss;
    
    //reco2d
    
    
    //Methods
    bool convert_2d(const event_hit     *evt_hits,
		    const event_cluster *evt_clusters,
		    const event_ass     *evt_ass_data);
    
    
    void check_cluster_boundaries();
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
