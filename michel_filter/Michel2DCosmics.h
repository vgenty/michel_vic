/**
 * \file Michel2DCosmics.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class Michel2DCosmics
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/

#ifndef LARLITE_MICHEL2DCOSMICS_H
#define LARLITE_MICHEL2DCOSMICS_H

#include "Michel2DAna.h"

namespace larlite {
  /**
     \class Michel2DCosmics
     User custom analysis class made by SHELL_USER_NAME
   */
  class Michel2DCosmics : public ana_base{
  
  public:
    
    /// Default constructor
    Michel2DCosmics(std::string cp) :
      _cluster_producer(cp)
    { _name="Michel2DCosmics"; _fout=0;}
    //m2Da = new Michel2DAna(cp);}
    
    /// Default destructor
    virtual ~Michel2DCosmics(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  private:
    TTree* _output_tree;
    int _evt = 0;
    
    std::vector<Double_t> _tX;
    std::vector<Double_t> _tY;
    
    std::vector<std::vector<Double_t> > _ahits_X_copy;
    std::vector<std::vector<Double_t> > _ahits_Y_copy;

    std::vector<std::vector<size_t> >  _ordered_pts_copy;
    
    //geo
    Double_t _time2cm;
    Double_t _wire2cm;
    
    std::string _cluster_producer;
    std::vector<ClusterYPlane*> _clusters; 

    bool convert_2d(const event_hit     *evt_hits,
		    const event_cluster *evt_clusters,
		    const event_ass     *evt_ass_data);
    
    
    void check_cluster_boundaries();
    void clear_all();
    bool find_projected_starts(std::vector<TVector2*>& p, //is this correect??
			       const event_mcshower* evt_mcshower);
  

    
    //Michel2DAna *m2Da; //I probably could inherit but why not just make heap instance
    
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
