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
#include <sstream>
#include <string>

// Larlite includes
#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"
#include "LArUtil/GeometryUtilities.h"
#include "LArUtil/Geometry.h"

// ROOT includes
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1D.h"

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
    
    TF1 *ONE;
    TF1 *TWO;
    TF1 *THREE;
    
    
    //Functions
    inline double distance(const double a, const double c,
			   const std::pair<double,double>& xy) {
      return fabs(-1.0*a*xy.first + xy.second - 1.0*c)/(sqrt(a*a+1));
    }
    
    inline std::pair<double,double> location(const double a, const double c,
					     const std::pair<double,double>&xy) {
      return std::make_pair(((xy.first + a*xy.second) - a * c) / (a*a + 1),
			    (-1.0*a*(-1.0*xy.first - a*xy.second) + c) / (a*a + 1));
    }
    
    bool opening_direction(double frac,
			   const size_t i,
			   std::vector<size_t>& sxidx);
    
    std::vector<size_t> sort_indexes(const std::vector<double> &v);
    std::vector<size_t> sort_indexes(const std::vector<std::pair<double,double> > &v);
    

    // size_t Triangle::the_index(double percentage,
    // 			       std::vector<double>& dists,
    // 			       std::vector<double>& idx);
    
    //Vars
    TH1D* fTh1d;
    
    std::vector<TF1*> fTf;
    std::vector<TGraphErrors*> fTg;
    
    
    std::map<UChar_t,std::vector<std::pair<double,double> > > fHits_xy;
    std::map<UChar_t,std::vector<double > > fHits_charge;

    std::map<UChar_t,std::vector<std::pair<double,double> > > fHits_xy_err;
    std::map<UChar_t,size_t> fHits_num;
    std::map<int,std::pair<double,double> > fFit_params;
    std::map<int,std::vector<double> > fHit_distances;
    
    std::map<UChar_t,std::vector<double> > fHits_dist;

    std::map<int,std::pair<double,double> > fLeft_right_dist;
    std::map<int,std::pair<double,double> > fLeft_right_dist_rms;
    std::map<int,std::pair<double,double> > fLeft_right_xs;
    
    //Boundaries
    std::map<int,std::pair<double,double> > fFirst_point;
    std::map<int,std::pair<double,double> > fPoint_on_line;
      
    std::map<int,std::pair<double,double> > fPerp_unit_vector;
    std::map<int,std::pair<double,double> > fPoint_above;
    std::map<int,std::pair<double,double> > fPoint_below;
    
    
    std::map<int,std::pair<double,double> > fOne;
    std::map<int,std::pair<double,double> > fTwo;
    std::map<int,std::pair<double,double> > fThree;

    
    
    //Boundary functions
    inline double one  (const int i, const double x) { return fOne[i].first * x   + fOne[i].second; }
    inline double two  (const int i, const double x) { return fTwo[i].first * x   + fTwo[i].second; }
    inline double three(const int i, const double x) { return fThree[i].first * x + fThree[i].second; }
    
    
    bool inside_boundaries(size_t i, const std::pair<double,double>& xy);
    
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
