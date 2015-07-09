/**
 * \file Michel.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class Michel
 *
 * @author vgenty
 */

/** \addtogroup michel_filter

    @{*/
#ifndef MICHEL_H
#define MICHEL_H

#include <iostream>

#include "TVector2.h"

/**
   \class Michel
   User defined class Michel ... these comments are used to generate
   doxygen documentation!
 */
class Michel{

public:

  /// Default constructor
  Michel(Double_t c, Double_t l, TVector2* s, Int_t nh) :
    _charge(c),
    _length(l),
    _start (*s),
    _num_hits(nh)
  {}
  
  /// Default destructor
  ~Michel(){}
  
  
  Double_t _charge;
  Double_t _length;
  TVector2 _start;
  Int_t    _num_hits;
  
  void dump();
};

#endif
/** @} */ // end of doxygen group 

