Tri#ifndef LARLITE_TRIANGLES_CXX
#define LARLITE_TRIANGLES_CXX

#include "Triangles.h"

namespace larlite {
  
  //functions for addition, subtraction, multiplication of two std::pairs
  //Addition of two pair
  //template <typename T,typename U>                                                   
  std::pair<double,double> operator+(const std::pair<double,double> & l,const std::pair<double,double> & r) {   
    return {l.first+r.first,l.second+r.second};                                    
  } 
  std::pair<double,double> operator-(const std::pair<double,double> & l,const std::pair<double,double> & r) {   
    return {l.first-r.first,l.second-r.second};                                    
  } 
  //multiplication
  double operator*(const std::pair<double,double> & l,const std::pair<double,double> & r) {   
    return l.first * r.first + l.second *  r.second;
  } 
  std::pair<double,double> operator*(const double c,const std::pair<double,double> & r) {   
    return {c * r.first, c *  r.second};                                    
  } 

  std::pair<double,double> Triangles::find_first_hit(const bool right,int i) {
    
    //Project all the hits on to the line
    std::vector<double > projection; //first is index distance second is x coordinate on the line

    for(const auto& xy : fHits_xy[i])
      projection.push_back(location(fFit_params[i].second,fFit_params[i].first,xy).first);
    
    //Sort projection from smallest to largest
    auto sorted_proj = sort_indexes(projection);
    
    if (right) 
      return fHits_xy[i][sorted_proj[0]];
    else
      return fHits_xy[i][sorted_proj[fHits_num[i] - 1]];
	
    
  }
  
  void Triangles::reset() {
    //Hit information    
    fHits_xy.clear();
    fHits_charge.clear();
    fHits_xy_err.clear();
    fHits_num.clear();
    
    

    fHit_distances.clear();
    fHits_dist.clear();
    
    fLeft_right_dist.clear();
    fLeft_right_dist_rms.clear();
    fLeft_right_xs.clear();
    
    fFirst_point.clear();
    fPoint_on_line.clear();
    
    fPerp_unit_vector.clear();
    fPoint_above.clear();
    fPoint_below.clear();
    
    
    fOne.clear();
    fTwo.clear();
    fThree.clear();

    for(int m = 0; m < 3; ++m) {
      delete fTg[m];
      fTg[m] = 0;
      delete fTf[m];
      fTf[m] = 0;
    }
      

  }
  
  bool Triangles::inside_boundaries(size_t i, const std::pair<double,double>& xy) {
    //std::cout << "Checking point " << "(" << xy.first << "," << xy.second << ") \n";
    //std::cout << one(i,xy.first) << std::endl;
    //std::cout << two(i,xy.first) << std::endl;

    if(xy.second <= one(i,xy.first)) {
      //std::cout << " is below the first line.." << std::endl;
      if(xy.second >= two(i,xy.first)) {
	//std::cout << " is above the second... " << std::endl;
	if(fThree[i].first > 0) { //
	  if(xy.second >= three(i,xy.first)) {
	    //std::cout << " and is above the hypoteneuse, it's in!!\n";
	    return true;
	  }
	} else {
	  if(xy.second <= three(i,xy.first)) {
	    //std::cout << " and is below the hypoteneuse, it's in!!\n";
	    return true;
	  }
	}
      }
      
    }
    
    return false;
  }
  
  bool Triangles::opening_direction(double frac,
				    const size_t i) {
    
    int num = floor(frac*fHits_num[i]);
    double lp_avg = 0.0, rp_avg = 0.0, lp_rms = 0.0, rp_rms = 0.0;
    double lp_avg_xs = 0.0, rp_avg_xs = 0.0;
    
    //rms = variance here, save a step by not taking the square root
    //average the distances to the line of the first frac points

    //calculate the distances from the line to all the points
    std::vector<double> distances;
    for(auto const& the_hit : fHits_xy[i])
      fHit_distances[i].push_back(distance(fFit_params[i].second,fFit_params[i].first,the_hit));

    
    //calculate the projections of the points on the line, use the farthest right and left projections along with their distances
    std::vector<double > projection; //first is index distance second is x coordinate on the line
    
    for(const auto& xy : fHits_xy[i])
      projection.push_back(location(fFit_params[i].second,fFit_params[i].first,xy).first);
    
    //Sort projection from smallest to largest
    auto sorted_proj = sort_indexes(projection);

    
    
    for(size_t lp = 0; lp < num; ++lp){
      lp_avg += fHit_distances[i][sorted_proj[lp]];
      rp_avg += fHit_distances[i][sorted_proj[fHits_num[i] - 1 - lp]];
      lp_avg_xs += fHits_xy[i][sorted_proj[lp]].first;
      rp_avg_xs += fHits_xy[i][sorted_proj[fHits_num[i] - 1 - lp]].first;

    }
    lp_avg /= num; rp_avg /= num; lp_avg_xs /= num; rp_avg_xs /= num;
    for(size_t lp = 0; lp < num; ++lp){ // avoid using pow() for no reason
      lp_rms += (fHit_distances[i][sorted_proj[lp]] - lp_avg)*(fHit_distances[i][sorted_proj[lp]] - lp_avg);
      rp_rms += (fHit_distances[i][sorted_proj[fHits_num[i] - 1 - lp]]-rp_avg)*(fHit_distances[i][sorted_proj[fHits_num[i] - 1 - lp]]-rp_avg);
    }
    
    
    lp_rms /= num;         rp_rms /= num;
    lp_rms = sqrt(lp_rms); rp_rms = sqrt(rp_rms);
    
    fLeft_right_dist[i]     = std::make_pair(lp_avg,rp_avg);
    fLeft_right_dist_rms[i] = std::make_pair(lp_rms,rp_rms);
    fLeft_right_xs[i]       = std::make_pair(lp_avg_xs,rp_avg_xs);
    
    if (rp_rms > lp_rms) 
      return true; //opens to the right
    else
      return false;
    
  }

  //sorts std::vector<double>
  std::vector<size_t> Triangles::sort_indexes(const std::vector<double> & v) {
    
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    
    // sort indexes by comparing values in vector
	sort(idx.begin(), idx.end(),
	   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];}); //sort low to high
    
    return idx;
  }
  
  //Sorts by pair.first
  std::vector<size_t> Triangles::sort_indexes(const std::vector<std::pair<double,double> > &v,int coord) {
    //im not sure how intensive this piece of code is...
    
    // copy the original index locations
    std::vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
    
    // sort indexes by comparing values in vector
    if(!coord) //sort by x's
      sort(idx.begin(), idx.end(),
	   [&v](size_t i1, size_t i2) {return v[i1].first < v[i2].first;}); //sort low to high
    else       //sort by y's
      sort(idx.begin(), idx.end(),
	   [&v](size_t i1, size_t i2) {return v[i1].second < v[i2].second;}); //sort low to high

    return idx;
  }  
  
  
  bool Triangles::initialize() {
    
    fTg.resize(3);
    fTf.resize(3);
    

    //tTf = new TF1("")
    fTh1d = new TH1D("first",";;",100,0,100);
    
    return true;
  }
  
  bool Triangles::analyze(storage_manager* storage) {
    reset();
    
    size_t chosen;
    double chosen_e = 10000.0;
    
    std::cout << "Event: " << storage->get_entries_read() << "\n";
    
    auto evt_hits = storage->get_data<event_hit>("gaushit"); 
    
    // >> Initialize only once
    const auto wire2cm   = ::larutil::GeometryUtilities::GetME()->WireToCm();
    const auto time2cm   = ::larutil::GeometryUtilities::GetME()->TimeToCm();
    
    
    // no hits just return do nothing
    std::cout << evt_hits->size() << "\n";
    if(evt_hits->size() < 2)
      return true;
    
    for(auto const& h : *evt_hits) {
      // Locations of hits
      double x      = h.WireID().Wire   * wire2cm;
      double y      = h.PeakTime()      * time2cm;
      
      // Arbitrary scaling of x/y errpr bars. Hits with little charge deposition are ``least important"
      double scale = 50.0/h.PeakAmplitude();
      //std::cout << "amp : " << h.PeakAmplitude() << "  sigma: " << h.SigmaPeakAmplitude() << "\n";
      
      double errx   =   scale;  //0.3 * scale;
      double erry   =   scale;  //h.SigmaPeakTime() * time2cm * scale;
      //std::cout << errx << " " << erry << "\n";
      UChar_t plane = ::larutil::Geometry::GetME()->ChannelToPlane(h.Channel());

      fHits_xy[plane].push_back(std::pair<double,double>(x,y));
      fHits_charge[plane].push_back(h.Integral());
      fHits_xy_err[plane].push_back(std::pair<double,double>(errx,erry));
      //fHits_dist[plane].push_back(0.0);
      // print out the plane
      // printf("%u",plane);
      fHits_num[plane] += 1;
    }
    
    
    
    //for each plane
    for(int i = 0; i < 3; ++i) {
      std::stringstream tg, tf;
      tg << "tgraph_plane_" << i;
      tf << "tfit_plane_" << i;
      
      fTg[i] = new TGraphErrors();
      fTg[i]->SetName(tg.str().c_str());
      fTf[i] = new TF1(tf.str().c_str(),"[0]+[1]*x",0,1100); //come back and make sure the ranges are kosher
      //fTf[i] = new TF1(tf.str().c_str(),"1.0++x",0,1100); //special linear fitting symbol ++

      fTf[i]->SetParameter(0,400);
      fTf[i]->SetParameter(1,-1.0);

      for(int j = 0; j < fHits_xy[i].size(); ++j) {
	fTg[i]->SetPoint(j,
			 fHits_xy[i][j].first,
			 fHits_xy[i][j].second);
	fTg[i]->SetPointError(j,
			      fHits_xy_err[i][j].first,
			      fHits_xy_err[i][j].second);
	
      }
      
      fTg[i]->Fit(fTf[i],"QCFN");
      //fTg[i]->Fit(fTf[i],"QFN");
      fFit_params[i] = std::make_pair(fTf[i]->GetParameter(0),fTf[i]->GetParameter(1));
      
      
      //Sort the hits by x location
      //soted_xs returns the index of fHits_xy with the most left wire# (sorted_xs[0])
      auto sorted_xs = sort_indexes(fHits_xy[i],0);
      auto sorted_ys = sort_indexes(fHits_xy[i],1);

      
      //Determine whether the shower is pointing down past pi/4 (currently no up)
      //      bool extreme_up_down = fabs(fFit_params[i].second) > 1.0 ? true : false;
      //std::cout << extreme_up_down << std::endl;

      //Determine which way the shower opens to the left
      //or to the right, bool right is true of it opens to the right
      bool right = opening_direction(0.10,i);
      
      // take the start of the shower (first point on the left or the right
      // depending on how the shower opens, find the point on the line down the way
      // EDIT: taking the first point sometimes fails really really bad. Find the first point starting from
      // the left that is below 1 stdev of the average distance to the line
      
      // if(right) {
      // 	for(size_t q = 0; q < fHits_num[i]; ++q) {
      // 	  if(distance(fFit_params[i].second,fFit_params[i].first,fHits_xy[i][sorted_xs[q]]) < 
      // 	     2*fLeft_right_dist_rms[i].first) {
      // 	    fFirst_point[i] = fHits_xy[i][sorted_xs[q]];
      // 	    break;
      // 	  }  
      // 	}
      // } else {
      // 	for(size_t q = 0; q < fHits_num[i]; ++q) {
      // 	  if(distance(fFit_params[i].second,fFit_params[i].first,fHits_xy[i][sorted_xs[fHits_num[i] - 1 - q]]) < 
      // 	     2*fLeft_right_dist_rms[i].second) {
      // 	    fFirst_point[i] = fHits_xy[i][sorted_xs[fHits_num[i] - 1 - q]];
      // 	    break;
      // 	  }  
      // 	}
      // }
      
      //fFirst_point[i]   = right ? fHits_xy[i][sorted_xs[i]] : fHits_xy[i][sorted_xs[fHits_num[i]-1]];
      fFirst_point[i]   = find_first_hit(right,i);
      
      
      
      fPoint_on_line[i] = right ? std::make_pair(fLeft_right_xs[i].second,
					      fFit_params[i].first + // y = m*x + b
					      fFit_params[i].second *
					      fLeft_right_xs[i].second) : std::make_pair(fLeft_right_xs[i].first,
											 fFit_params[i].first +
											 fFit_params[i].second *
											 fLeft_right_xs[i].first);
      
      fPerp_unit_vector[i] = 5.0/(1.0 + fFit_params[i].second*fFit_params[i].second) * std::make_pair(-1.0,1.0*fFit_params[i].second);

if(fPerp_unit_vector[i].second < 0) fPerp_unit_vector[i] = -1.0*fPerp_unit_vector[i];

      fPoint_above[i]      = right ? fPoint_on_line[i] + (fLeft_right_dist[i].second * fPerp_unit_vector[i]) : fPoint_on_line[i] + (fLeft_right_dist[i].first * fPerp_unit_vector[i]);
      fPoint_below[i]      = right ? fPoint_on_line[i] - (fLeft_right_dist[i].second * fPerp_unit_vector[i]) : fPoint_on_line[i] - (fLeft_right_dist[i].first * fPerp_unit_vector[i]);
      
      
      //Make the lines first is y inter, second is slope
      fOne[i] = std::make_pair((fPoint_above[i].second - fFirst_point[i].second)/
			       (fPoint_above[i].first - fFirst_point[i].first),
			       (fPoint_above[i].second - fFirst_point[i].second)/
			       (fPoint_above[i].first - fFirst_point[i].first)*(-1.0)*fFirst_point[i].first +
			       fFirst_point[i].second);
      fTwo[i] = std::make_pair((fPoint_below[i].second - fFirst_point[i].second)/
			       (fPoint_below[i].first - fFirst_point[i].first),
			       (fPoint_below[i].second - fFirst_point[i].second)/
			       (fPoint_below[i].first - fFirst_point[i].first)*(-1.0)*fFirst_point[i].first +
			       fFirst_point[i].second);
      fThree[i] = std::make_pair((fPoint_above[i].second - fPoint_below[i].second)/
				 (fPoint_above[i].first - fPoint_below[i].first),
				 (fPoint_above[i].second - fPoint_below[i].second)/
				 (fPoint_above[i].first - fPoint_below[i].first)*(-1.0)*fPoint_below[i].first +
				 fPoint_below[i].second);
			       
      
      std::cout << "\n\n\tPlane " << i << " diagnostics \n" << std::endl;
      std::cout << "right " << right << "\n";
      std::cout << " best fit y-inter " << fTf[i]->GetParameter(0) << " slope: " << fTf[i]->GetParameter(1) << "\n";
      std::cout << "average x location of first .1 vals " << fLeft_right_xs[i].first << "\n";
      std::cout << "average x location of last .1 vals "  << fLeft_right_xs[i].second << "\n";
      std::cout << "average distance to line of first .1 vals " << fLeft_right_dist[i].first << "\n";
      std::cout << "average distance to line of last .1 vals "  << fLeft_right_dist[i].second << "\n";
      std::cout << "first point: " << " (" << fFirst_point[i].first << "," << fFirst_point[i].second << ")\n";
      std::cout << "point on line: " << " (" << fPoint_on_line[i].first << "," << fPoint_on_line[i].second << ")\n";
      std::cout << "point above: " << " (" << fPoint_above[i].first << "," << fPoint_above[i].second << ")\n";
      std::cout << "point below: " << " (" << fPoint_below[i].first << "," << fPoint_below[i].second << ")\n";
      std::cout << "perp direction " << "(" << fPerp_unit_vector[i].first << "," << fPerp_unit_vector[i].second << ")\n";
      std::cout << "The CHI for this PLOT is " << fTf[i]->GetChisquare() << "\n";
      std::cout << "The slope err for this PLOT is " << fTf[i]->GetParError(1) << "\n";
      
      if(fTf[i]->GetParError(1) < chosen_e) {
	chosen_e = fTf[i]->GetParError(1);
	chosen = i;
      }
      
      
      std::cout << "\t\n I thought view: " << chosen << " was the best\n";
      //Stop here and check the output this could be horrible we don't know
      if(i == 1) {
      	ONE    = new TF1("ONE",  "[0]+[1]*x",0,1100);
      	TWO    = new TF1("TWO",  "[0]+[1]*x",0,1100);
      	THREE  = new TF1("THREE","[0]+[1]*x",0,1100);

      	ONE->SetParameter(0,fOne[i].second);
      	ONE->SetParameter(1,fOne[i].first);
      	TWO->SetParameter(0,fTwo[i].second);
      	TWO->SetParameter(1,fTwo[i].first);
      	THREE->SetParameter(0,fThree[i].second);
      	THREE->SetParameter(1,fThree[i].first);
	
	
      }
      
      
      //
      //OLD
      //
      // Find the highest 90% index, the rest have distances lower than this value
      // int index = floor(0.1*fHits_dist[i].size());
      // double max_distance = fHits_dist[i][sorted_dists[index]];
      
      // std::cout << max_distance << "\n";
      
      
      // //average the 5% farthest points on the line and return the point on the
      // //line that is closest to the the point with the most average distance
      
      // auto loc = location(fTf[i]->GetParameter(1),
      // 			  fTf[i]->GetParameter(0),
      // 			  fHits_xy[i][the_index(0.05,fHists_dist[i],sorted_dists)]);
      
      // //get unit vector of the tangent line
      // std::pair<double,double> unit_tangent_vec = std::make_pair(-1.0/sqrt((1+fTf[i]->GetParameter(1)*fTf[i]->GetParameter(1))),
      // 								 -1.0*fTf[i]->GetParameter(1)/sqrt(1+fTf[i]->GetParameter(1)*fTf[i]->GetParameter(1)));
      
      // //create the two points above and below the line that we will use to create the triangle
      // std::pair<double,double> above = std::make_pair(loc.first + max_distance*unit_tangent_vec.first,
      // 						      loc.second + max_distance*unit_tangent_vec.second);
      // std::pair<double,double> below = std::make_pair(loc.first - max_distance*unit_tangent_vec.first,
      // 						      loc.second - max_distance*unit_tangent_vec.second);
      
      
      // //find the first and the last hit
      // std::pair<double,double> the_first = fHits_xy[i][sorted_xs[sorted_xs.size()-1]];
      // std::pair<double,double> the_last  = fHits_xy[i][sorted_xs[0]];
      
      
      
      // //find the 
      // //Create the bounding box function which will be used to test of the hit is inside the triangle
      
      
      // // find the third farthest point on the line for no real reason 
      // // and compute the closest XY value on that line
      
      // // auto farthest = location(fTf[i]->GetParameter(1),
      // // 			       fTf[i]->GetParameter(0),
      // // 			       fHits_xy[i][sorted_dists[2]]);
      
      
      
	}
    
  
      //Now I have my chosen view....
    
      //Loop over the hits in this view add up the charge in the triangle

    double tot_charge = 0.0;
    double excess_charge = 0.0;

    for(int k = 0; k < fHits_num[chosen]; ++k)
      if(inside_boundaries(chosen,fHits_xy[chosen][k]))
	tot_charge += fHits_charge[chosen][k];
      else
	excess_charge += fHits_charge[chosen][k];
    
    bool extreme_up_down = fabs(fFit_params[chosen].second) > 1.0 ? true : false;
    
    std::cout << "\tChosen plane:       " << chosen << std::endl;
    std::cout << "\tMy tot_charge:      " << tot_charge << "\n";
    std::cout << "\tAll charge in event " << tot_charge+excess_charge << "\n";
    std::cout << "\tExtreme up down?    " << extreme_up_down << std::endl;
    
    // std::cout << "\n Fit diagonistics: \n";
    // std::cout << fTf[0]->GetParameter(0) << " " << fTf[0]->GetParameter(1) << "\n";
    // std::cout << fTf[1]->GetParameter(0) << " " << fTf[1]->GetParameter(1) << "\n";
    // std::cout << fTf[2]->GetParameter(0) << " " << fTf[2]->GetParameter(1) << "\n";
    
    return true;
  }
  
  bool Triangles::finalize() {
    std::cout << "In finalize... \n";
    for(auto k : fTg)
      k->Write();
    for(auto k : fTf) 
      k->Write();
    
    fTh1d->Write();
    
    ONE->Write();
    TWO->Write();
    THREE->Write();
    
    return true;
  }

}


#endif
