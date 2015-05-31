#ifndef TRIANGLE_CXX
#define TRIANGLE_CXX

#include "Triangle.h"

Triangle::Triangle(const std::pair<double,double>& First,
		   const std::pair<double,double>& Above,
		   const std::pair<double,double>& Below,
		   const int  v,
		   const bool r) {
  
  fLine_one   = std::make_pair(  (Above.second - First.second)
			       / (Above.first - First.first  ),
			         (Above.second - First.second) 
			       / (Above.first - First.first  ) *
			         (-1.0) * First.first + First.second );

  fLine_two   = std::make_pair(  (Below.second - First.second)
			       / (Below.first - First.first  ),
			         (Below.second - First.second) 
			       / (Below.first - First.first  ) *
			         (-1.0) * First.first + First.second );

  fLine_three = std::make_pair(  (Above.second - Below.second)
			       / (Above.first - Below.first  ),
			         (Above.second - Below.second) 
			       / (Above.first - Below.first  ) *
			         (-1.0) * Below.first + Below.second );

  fChosen = false;
  fView   = v;
  fRight  = r;
}

void Triangle::check_boundaries(const std::vector<std::pair<double,double> >& xy ){
  
  for(unsigned int k = 0; k < xy.size(); ++k) {
    if(xy[k].second <= above(xy[k].first)) {
      if(xy[k].second >= below(xy[k].first)) {
	if(fRight) {
	  if(fLine_three.first > 0) { //
	    if(xy[k].second >= hypot(xy[k].first)) {
	      fHits.push_back(k); // write the index to fHits
	    }
	  } else {
	    if(xy[k].second <= hypot(xy[k].first)) {
	      fHits.push_back(k); // write the index to fHits
	    }
	  }
	} else { // if opens to the left
	  if(fLine_three.first < 0) {
	    if(xy[k].second >= hypot(xy[k].first)) {
	      fHits.push_back(k); // write the index to fHits
	    }
	  } else {
	    if(xy[k].second <= hypot(xy[k].first)) {
	      fHits.push_back(k); // write the index to fHits
	    }
	  }
	}
      }
    }
  }



}


#endif
