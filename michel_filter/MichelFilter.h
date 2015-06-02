/**
 * \file MichelFilter.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class MichelFilter
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/

#ifndef LARLITE_MICHELFILTER_H
#define LARLITE_MICHELFILTER_H

#include "Analysis/ana_base.h"
#include "DataFormat/mctrack.h"

namespace larlite {
  /**
     \class MichelFilter
     User custom analysis class made by SHELL_USER_NAME
   */
  class MichelFilter : public ana_base{
  
  public:

    /// Default constructor
    MichelFilter(){ _name="MichelFilter"; _fout=0;}

    /// Default destructor
    virtual ~MichelFilter(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

  protected:
    
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
