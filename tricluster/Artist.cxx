#ifndef ARTIST_CXX
#define ARTIST_CXX

#include "Artist.h"

//
// Teach std::pair how to add, subtract, dot product, and multiply by scalar
//

std::pair<double,double> operator+(const std::pair<double,double> & l,const std::pair<double,double> & r) {   
  return {l.first+r.first,l.second+r.second};                                    
} 
std::pair<double,double> operator-(const std::pair<double,double> & l,const std::pair<double,double> & r) {   
  return {l.first-r.first,l.second-r.second};                                    
} 
double operator*(const std::pair<double,double> & l,const std::pair<double,double> & r) {   
  return l.first * r.first + l.second *  r.second;
} 
std::pair<double,double> operator*(const double c,const std::pair<double,double> & r) {   
  return {c * r.first, c *  r.second};                                    
} 
  
//
// Setters for hits
//
void Artist::set_hit_point(const int p,std::pair<double,double> xy) {
  fHits_xy[p].push_back(xy);
}
void Artist::set_hit_point_err(const int p,std::pair<double,double> xyerr){
  fHits_xy_err[p].push_back(xyerr);
}
void Artist::set_hit_charge(const int p,const double charge) {
  fHits_charge[p].push_back(charge);
}

//
// Inline functions for distance to line, and nearest location
//

inline double Artist::distance(const int p,const int idx) {
  return fabs(-1.0*fFit_slope[p]*fHits_xy[p][idx].first 
	      + fHits_xy[p][idx].second 
	      - 1.0*fFit_yinter[p]) /
    (sqrt(fFit_slope[p]*fFit_slope[p]+1));
}

inline std::pair<double,double> Artist::location(const int p,
						 const int idx) {
  return std::make_pair(
			((fHits_xy[p][idx].first 
			  + fFit_slope[p]* fHits_xy[p][idx].second)
			 - fFit_slope[p] * fFit_yinter[p])
			/ (fFit_slope[p] * fFit_slope[p] + 1)
			,
			(-1.0*fFit_slope[p] * (-1.0*fHits_xy[p][idx].first 
					       - fFit_slope[p] 
					       * fHits_xy[p][idx].second)
			 + fFit_yinter[p]) 
			/ (fFit_slope[p]*fFit_slope[p] + 1)
			);
}

//
// Sorting functions
//
std::vector<int> Artist::sort_indexes(const std::vector<double> & v) {
  
  std::vector<int> idx(v.size());
  for (unsigned int i = 0; i != idx.size(); ++i) idx[i] = i;
  
  // sort indexes by comparing values in vector
  std::sort(idx.begin(), idx.end(),
	    [&v](int i1, int i2) {return v[i1] < v[i2];}); //sort low to high
  
  return idx;
}


std::vector<int> Artist::sort_indexes(const std::vector<std::pair<double,double> > &v,int coord) {
  //im not sure how intensive this piece of code is...
  
  // copy the original index locations
  std::vector<int> idx(v.size());
  for (unsigned int i = 0; i != idx.size(); ++i) idx[i] = i;
  
  // sort indexes by comparing values in vector
  if(!coord) //sort by x's
    std::sort(idx.begin(), idx.end(),
	      [&v](int i1, int i2) {return v[i1].first < v[i2].first;}); //sort low to high
  else       //sort by y's
    std::sort(idx.begin(), idx.end(),
	 [&v](int i1, int i2) {return v[i1].second < v[i2].second;}); //sort low to high
  
  return idx;
}  

//
// Determine the first hit location, this will be used 
// as the start of the triangle
//
std::pair<double,double> Artist::find_first_hit(const int p) {
    
    if (fRight[p]) 
      return fHits_xy[p][fSorted_proj_x[p][0]];
    else
      return fHits_xy[p][fSorted_proj_x[p][fHits_num[p] - 1]];
	
}

//
// Determine the point on the line that is the average of the last 
// frac points, hopefully fRight is correct!!!
//
std::pair<double,double> Artist::find_point_on_line(const int p) {
  
  if(fRight[p])
    return std::make_pair(fHits_avg_x[p].second,
			  fFit_slope[p] * fHits_avg_x[p].second + fFit_yinter[p]);
  else
    return std::make_pair(fHits_avg_x[p].first,
			  fFit_slope[p] * fHits_avg_x[p].first + fFit_yinter[p]);
  
}
std::pair<double,double> Artist::find_point_above(const int p,
						  const double scale) {
  if(fRight[p])
    return fPoint_on_line[p] + (fHits_avg_dist[p].second * scale * fFit_perp_vec[p]);
  else
    return fPoint_on_line[p] + (fHits_avg_dist[p].first  * scale * fFit_perp_vec[p]);

}
std::pair<double,double> Artist::find_point_below(const int p,
						  const double scale) {
  
  if(fRight[p])
    return fPoint_on_line[p] - (fHits_avg_dist[p].second * scale * fFit_perp_vec[p]);
  else
    return fPoint_on_line[p] - (fHits_avg_dist[p].first  * scale * fFit_perp_vec[p]);
    
}

bool Artist::opening_direction(double frac,
			       const int p) {
    
  int num = floor(frac*fHits_num[p]);
  double lp_avg = 0.0, rp_avg = 0.0, lp_rms = 0.0, rp_rms = 0.0;
  double lp_avg_xs = 0.0, rp_avg_xs = 0.0;
  
  //calculate the perp distances from the line to all the points
  //calculate the xy projection onto the line (closest point on line)
  
  std::vector<double> distances;
  for(int idx = 0; idx < fHits_num[p]; ++idx) {
    fHits_dist[p].push_back(distance(p,idx));
    fHits_proj[p].push_back(location(p,idx));
  }
  
  
  //Sort projection of hits on line from smallest to largest X
  fSorted_proj_x[p] = sort_indexes(fHits_proj[p],0);

  for(int lp = 0; lp < num; ++lp){ //look left (lp) to right
    lp_avg    += fHits_dist[p][fSorted_proj_x[p][lp]];
    rp_avg    += fHits_dist[p][fSorted_proj_x[p][fHits_num[p] - 1 - lp]];
    lp_avg_xs += fHits_xy[p][fSorted_proj_x[p][lp]].first;
    rp_avg_xs += fHits_xy[p][fSorted_proj_x[p][fHits_num[p] - 1 - lp]].first;
    
  }
  
  lp_avg /= num; rp_avg /= num; lp_avg_xs /= num; rp_avg_xs /= num;
  
  for(int lp = 0; lp < num; ++lp){ // avoid using pow() for no reason
    lp_rms += (fHits_dist[p][fSorted_proj_x[p][lp]] - lp_avg)*(fHits_dist[p][fSorted_proj_x[p][lp]] - lp_avg);
    rp_rms += (fHits_dist[p][fSorted_proj_x[p][fHits_num[p] - 1 - lp]]-rp_avg)*(fHits_dist[p][fSorted_proj_x[p][fHits_num[p] - 1 - lp]]-rp_avg);
    }
    
  
  lp_rms /= num;         rp_rms /= num;
  lp_rms = sqrt(lp_rms); rp_rms = sqrt(rp_rms);
  
  fHits_avg_dist[p] = std::make_pair(lp_avg,rp_avg);
  fHits_rms_dist[p] = std::make_pair(lp_rms,rp_rms);
  fHits_avg_x[p]    = std::make_pair(lp_avg_xs,rp_avg_xs);
  
  if (rp_rms > lp_rms) 
    return true; //opens to the right
  else
    return false;
    
  }

void Artist::fit_hits(const int p) {

  fHits_num[p] = fHits_xy[p].size();
  
  fTg = new TGraphErrors();
  
  fTf = new TF1("","[0]+[1]*x",0,1100);
  fTf->SetParameter(0,400);
  fTf->SetParameter(1,-1.0);
  
  for(int j = 0; j < fHits_num[p]; ++j) {
    fTg->SetPoint(j,fHits_xy[p][j].first,fHits_xy[p][j].second);
    fTg->SetPointError(j,fHits_xy_err[p][j].first,fHits_xy_err[p][j].second);
   }
  
  fTg->Fit(fTf,"QCFN");
  fFit_slope[p]      = fTf->GetParameter(1);
  fFit_yinter[p]     = fTf->GetParameter(0);
  fFit_slope_err[p]  = fTf->GetParError(1);
  fFit_yinter_err[p] = fTf->GetParError(0);
  
  fFit_perp_vec[p] = 1.0/(1.0 + fFit_slope[p]*fFit_slope[p]) * std::make_pair(-1.0,1.0*fFit_slope[p]);

  if(fFit_perp_vec[p].second < 0) 
    fFit_perp_vec[p] = -1.0 * fFit_perp_vec[p];

  
  delete fTg; fTg = 0;
  delete fTf; fTf = 0;
}

void Artist::create_triangle(const int p) {
  
  auto frac = 0.10;
  
  fRight[p] = opening_direction(frac,p);
  
  fPoint_first[p]   = find_first_hit    (p);
  fPoint_on_line[p] = find_point_on_line(p);
  fPoint_above[p]   = find_point_above  (p,5.0);
  fPoint_below[p]   = find_point_below  (p,5.0);
  
  fTriangles[p] = new Triangle(fPoint_first[p],
			       fPoint_above[p],
			       fPoint_below[p],
			       p,
			       fRight[p]);
  
  fTriangles[p]->check_boundaries(fHits_xy[p]);
}

void Artist::choose_best_triangle() {
  
  //Artist will choose the best triangle based
  //on the fit with the lowest slope error...
  
  int chosen = 0;
  double chosen_e = 1000.0;
  
  for(unsigned int i = 0; i < fHits_num.size(); ++i) {
    if(fFit_slope_err[i] < chosen_e) {
      chosen_e = fFit_slope_err[i];
      chosen = i;
    }
  }
  
  fTriangles[chosen]->set_chosen();
  fChosen = chosen;
}

//
// Clear function
//

void Artist::erase() {
  
  //Hits
  fHits_xy.clear();
  fHits_xy_err.clear();
  fHits_charge.clear();
  fHits_num.clear();
  
  //Distance from hit to line, averages, RMS
  fSorted_proj_x.clear();
  
  fHits_dist.clear();
  fHits_proj.clear();
  fHits_avg_dist.clear();
  fHits_rms_dist.clear();
  fHits_avg_x.clear();
  
  fFit_slope.clear();
  fFit_yinter.clear();
  fFit_slope_err.clear();
  fFit_yinter_err.clear();
  
  fFit_perp_vec.clear();
  
  //Points
  fPoint_first.clear();
  fPoint_on_line.clear();
  fPoint_above.clear();
  fPoint_below.clear();
  
  // Triangles
  fRight.clear();
  fTriangles.clear();
  fNumTriangles = 0;

}

#endif
