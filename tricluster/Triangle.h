/**
 * \file Triangle.h
 *
 * \ingroup tricluster
 * 
 * \brief Class def header for a class Triangle
 *
 * @author vgenty
 */

/** \addtogroup tricluster

    @{*/
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <iostream>
#include <vector>

class Triangle{

public:

  /// Default constructor
  Triangle(const std::pair<double,double>& First,
	   const std::pair<double,double>& Above,
	   const std::pair<double,double>& Below,
	   const int v);
  
  /// Default destructor
  ~Triangle(){}
  
  void check_boundaries(const std::vector<std::pair<double,double> >& xy );

  //Set chosen
  inline void set_chosen() { fChosen = true; }
  
  inline const bool                       is_chosen()  { return fChosen;     }
  inline const std::vector<int>         & get_fHit()   { return fHits;       }
  inline const std::pair<double,double> & line_one()   { return fLine_one;   }
  inline const std::pair<double,double> & line_two()   { return fLine_two;   }
  inline const std::pair<double,double> & line_three() { return fLine_three; }

  
protected:
  
private:
  
  bool fRight;  //direction of the shower, fRight = true if shower opens to right
  bool fChosen; //this triangle is chosen as the best of the bunch by the artist
  int  fView;
  std::pair<double,double> fLine_one;
  std::pair<double,double> fLine_two;
  std::pair<double,double> fLine_three;
  
  std::vector<int> fHits; //Hit indicies which lie inside the triangle
  
  //Boundary functions
  inline double above (const double& x) { return fLine_one.first   * x   + fLine_one.second; } 
  inline double below (const double& x) { return fLine_two.first   * x   + fLine_two.second; }
  inline double hypot (const double& x) { return fLine_three.first * x   + fLine_three.second; }
  
};

#endif
/** @} */ // end of doxygen group 

