/**
 * \file Muon.h
 *
 * \ingroup michel_filter
 * 
 * \brief Class def header for a class Muon
 *
 * @author kathrynsutton
 */

/** \addtogroup michel_filter

    @{*/
#ifndef MUON_H
#define MUON_H

#include <iostream>

#include "TVector2.h"
#include "DataFormat/hit.h"
#include "AHIT.h"

/**
   \class Muon
   User defined class Muon ... these comments are used to generate
   doxygen documentation!
 */
class Muon{

public:

  /// Default constructor
  Muon(const std::vector<ahit> & h, std::vector<double> d ){
    _hits = h;
    _distance = d;

     double sum = 0.0;
     std::for_each(d.begin(), d.end(), [&] (double dist){
	 sum += dist;
       });

     _length = sum;

     double sum2 = 0.0;
     std::for_each(h.begin(), h.end(), [&] (ahit myhit){
	 sum2 += myhit.hit.Integral();
       });

     _charge = sum2;

     if (h[0].hit.Integral() <  h[h.size()-1].hit.Integral()){
       _start  = *h[0].vec;
       _end = *h[h.size()-1].vec;
     }

     else{
        _start  =  *h[h.size()-1].vec;
	_end = *h[0].vec;
     }
   
  }

  /// Default destructor
  ~Muon(){}

  std::vector<ahit> _hits;
  std::vector<Double_t> _distance;

  Double_t _charge;
  Double_t _length;
  
  TVector2 _start;
  TVector2 _end;
  Int_t    _num_hits;
  
  void dump();
  

};

#endif
/** @} */ // end of doxygen group 

