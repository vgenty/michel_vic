/**
 * \file Artist.h
 *
 * \ingroup tricluster
 * 
 * \brief Class def header for a class Artist
 *
 * @author vgenty
 */

/** \addtogroup tricluster

    @{*/
#ifndef ARTIST_H
#define ARTIST_H

#include <iostream>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>

#include "Triangle.h"

#include "TF1.h"
#include "TGraphErrors.h"


/**
   \class Artist
   User defined class Artist ... these comments are used to generate
   doxygen documentation!
 */
class Artist{

public:

  /// Default constructor
  Artist(){}

  /// Default destructor
  ~Artist(){}

  void set_num_triangles(const int n) { fNumTriangles = n; }
  
  void set_hit_point    (const int p,std::pair<double,double> xy);
  void set_hit_point_err(const int p,std::pair<double,double> xyerr);
  void set_hit_charge   (const int p,const double charge);

  void erase();
  bool opening_direction(double frac,
			 const int p);
  void fit_hits(const int p);
  void create_triangle(const int p);
  void choose_best_triangle();


  inline const std::vector<std::pair<double,double> >& get_hits        (const int p) { return fHits_xy[p];     }
  inline const std::vector<std::pair<double,double> >& get_hits_err    (const int p) { return fHits_xy_err[p]; }
  inline const std::vector<double>                   & get_hits_charge (const int p) { return fHits_charge[p]; }
  inline const Triangle*                               get_triangle    (const int p) { return fTriangles[p];   }
  inline int                                           get_chosen      ()            { return fChosen; }
private:
  
  //
  //Functions
  //
  inline double distance(const int p,const int idx);
  inline std::pair<double,double> location(const int p,
					   const int idx);
  
  
  std::vector<int> sort_indexes(const std::vector<double> & v);
  std::vector<int> sort_indexes(const std::vector<std::pair<double,double> > &v,int coord);

  std::pair<double,double> find_first_hit(const int p);
  std::pair<double,double> find_point_on_line(const int p);
  std::pair<double,double> find_point_above(const int p,
					    const double scale);
  std::pair<double,double> find_point_below(const int p,
					    const double scale);
  
  
  
  //
  //Variables
  //


  //Hits
  std::map<int,std::vector<std::pair<double,double> > > fHits_xy;
  std::map<int,std::vector<std::pair<double,double> > > fHits_xy_err;
  std::map<int,std::vector<double> >                    fHits_charge;
  std::map<int,int>                                     fHits_num;
  
  //Distance from hit to line, averages, RMS

  std::map<int,std::vector<int> > fSorted_proj_x;
  
  std::map<int,std::vector<double> >                    fHits_dist;
  std::map<int,std::vector<std::pair<double,double> > > fHits_proj;
  std::map<int,std::pair<double,double> > fHits_avg_dist;
  std::map<int,std::pair<double,double> > fHits_rms_dist;
  std::map<int,std::pair<double,double> > fHits_avg_x;
  
  //Graph and Fit, only one of each
  TF1* fTf;
  TGraphErrors* fTg;

  std::map<int,double>  fFit_slope;
  std::map<int,double>  fFit_yinter;
  std::map<int,double>  fFit_slope_err;
  std::map<int,double>  fFit_yinter_err;

  std::map<int,std::pair<double,double> > fFit_perp_vec;
  
  //Points
  std::map<int,std::pair<double,double> > fPoint_first;
  std::map<int,std::pair<double,double> > fPoint_on_line;
  std::map<int,std::pair<double,double> > fPoint_above;
  std::map<int,std::pair<double,double> > fPoint_below;
  
  // Triangles
  std::map<int,bool> fRight;
  int fNumTriangles;
  std::map<int,Triangle*> fTriangles;
  
  int fChosen;
  
};

#endif
/** @} */ // end of doxygen group 

