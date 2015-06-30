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

//larlite 
#include "Analysis/ana_base.h"

//ROOT
#include "TTree.h"


namespace larlite {
  /**
     \class Michel2DAna
     User custom analysis class made by SHELL_USER_NAME
   */
  class Michel2DAna : public ana_base{
  
  public:

    /// Default constructor
    Michel2DAna() : 
      _output_tree(nullptr)
    { _name="Michel2DAna"; _fout=0;}

    /// Default destructor
    virtual ~Michel2DAna(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    
  private:
    
    TTree *_output_tree;
    
    
    
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
