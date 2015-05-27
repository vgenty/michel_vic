/**
 * \file Triangles.h
 *
 * \ingroup triangles
 * 
 * \brief Class def header for a class Triangles
 *
 * @author vgenty
 */

/** \addtogroup triangles

    @{*/

#ifndef LARLITE_TRIANGLES_H
#define LARLITE_TRIANGLES_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class Triangles
     User custom analysis class made by SHELL_USER_NAME
   */
  class Triangles : public ana_base{
  
  public:

    /// Default constructor
    Triangles(){ _name="Triangles"; _fout=0;}

    /// Default destructor
    virtual ~Triangles(){}

    //Standard routine
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
