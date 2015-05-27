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

//C++ includes
// #include <map>
// #include <vector>

// Larlite includes
#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/Geometry.h"

// ROOT includes
#include "TGraphErrors.h"
#include "TF1.h"

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
    

  private:

    TF1 *fTf;
    std::vector<TGraphErrors*> fTg;
    
    std::map<UChar_t,std::vector<std::pair<double,double> > > fHits_xy;
    std::map<UChar_t,std::vector<std::pair<double,double> > > fHits_xy_err;

    
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
