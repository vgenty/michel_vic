/**
 * \file TriClusterAna.h
 *
 * \ingroup tricluster
 * 
 * \brief Class def header for a class TriClusterAna
 *
 * @author vgenty
 */

/** \addtogroup tricluster

    @{*/

#ifndef LARLITE_TRICLUSTERANA_H
#define LARLITE_TRICLUSTERANA_H

// Larlite includes
#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/Geometry.h"

//Local headers
#include "Artist.h"

namespace larlite {
  /**
     \class TriClusterAna
     User custom analysis class made by SHELL_USER_NAME
   */
  class TriClusterAna : public ana_base{
  
  public:

    /// Default constructor
    TriClusterAna(){ _name="TriClusterAna"; _fout=0;}

    /// Default destructor
    virtual ~TriClusterAna(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    
    inline Artist& get_artist() { return fArtist;}

  protected:

  private:
    Artist fArtist;

    double fWire2cm;
    double fTime2cm;
    
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
