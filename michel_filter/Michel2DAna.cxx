#ifndef LARLITE_MICHEL2DANA_CXX
#define LARLITE_MICHEL2DANA_CXX

#include "Michel2DAna.h"
#include <iomanip>

namespace larlite {
 
  bool Michel2DAna::initialize() {
    
    _wire2cm   = ::larutil::GeometryUtilities::GetME() ->WireToCm();
    _time2cm   = ::larutil::GeometryUtilities::GetME() ->TimeToCm();    
    _ne2ADC    = ::larutil::DetectorProperties::GetME()->ElectronsToADC();
    
    _output_tree = new TTree("out_tree","aho_tree");

    _output_tree->Branch("true_X" , &_tX , "tX/D");
    _output_tree->Branch("true_Y" , &_tY , "tY/D");
    _output_tree->Branch("reco_X" , &_rX , "rX/D");
    _output_tree->Branch("reco_Y" , &_rY , "rY/D"); 
    
    _output_tree->Branch("ahits_X"     , "std::vector<Double_t>" , &_ahits_X_copy);
    _output_tree->Branch("ahits_Y"     , "std::vector<Double_t>" , &_ahits_Y_copy);
    _output_tree->Branch("charges"     , "std::vector<Double_t>" , &_charges_copy);
    _output_tree->Branch("ordered_pts" , "std::vector<size_t>"   , &_ordered_pts_copy);
    _output_tree->Branch("mean_charges", "std::vector<Double_t>" , &_mean_charges_copy);
    _output_tree->Branch("dqds"        , "std::vector<Double_t>" , &_dqds_copy);
    _output_tree->Branch("s"           , "std::vector<Double_t>" , &_s_copy);
    _output_tree->Branch("startX"      , &_startX, "startX/D");
    _output_tree->Branch("startY"      , &_startY, "startY/D");
    _output_tree->Branch("endX"        , &_endX,   "endX/D"  );
    _output_tree->Branch("endY"        , &_endY,   "startX/D");

    _output_tree->Branch("_michel_E"      , &_michel_E      , "_michel_E/D");
    _output_tree->Branch("_michel_L"      , &_michel_L      , "_michel_L/D");
    _output_tree->Branch("d_michel_hit"   , &_d_m_h         , "d_michel_hit/D");
    _output_tree->Branch("_true_michel_E" , &_true_michel_E , "_true_michel_E/D");
    _output_tree->Branch("_reco_Q_o_mc_Q" , &_reco_Q_o_mc_Q , "_reco_Q_o_mc_Q/D");
    
    
    _output_tree->Branch("_mcQ_frac", &_mcQ_frac, "_mcQ_frac/D");

    _output_tree->Branch("_Q_tot_p2" , &_Q_tot_p2 , "_Q_tot_p2/D");
    _output_tree->Branch("_Q_u_p2" , &_Q_u_p2 , "_Q_u_p2/D");
    _output_tree->Branch("_MeV_scale", &_MeV_scale, "_MeV_scale/D");
    _output_tree->Branch("_true_michel_Det", &_true_michel_Det, "_true_michel_Det/D");


    _output_tree->Branch("_simch_michel_true_shower_E",&_simch_michel_true_shower_E,"_simch_michel_true_shower_E/D");
    _output_tree->Branch("_simch_michel_false_shower_E",&_simch_michel_false_shower_E,"_simch_michel_false_shower_E/D");
    _output_tree->Branch("_simch_plane_true_shower_E",&_simch_plane_true_shower_E,"_simch_plane_true_shower_E/D");
    _output_tree->Branch("_simch_plane_false_shower_E",&_simch_plane_false_shower_E,"_simch_plane_false_shower_E/D");
    _output_tree->Branch("_simch_ordered_true_shower_E",&_simch_ordered_true_shower_E,"_simch_ordered_true_shower_E/D");
    _output_tree->Branch("_simch_ordered_false_shower_E",&_simch_ordered_false_shower_E,"_simch_ordered_false_shower_E/D");
    _output_tree->Branch("_simch_cluster_true_shower_E",&_simch_cluster_true_shower_E,"_simch_cluster_true_shower_E/D");
    _output_tree->Branch("_simch_cluster_false_shower_E",&_simch_cluster_false_shower_E,"_simch_cluster_false_shower_E/D");
    
    
    _output_tree->Branch("_large_frac_shower_hits_X", "std::vector<Double_t>" , &_large_frac_shower_hits_X);
    _output_tree->Branch("_large_frac_shower_hits_Y", "std::vector<Double_t>" , &_large_frac_shower_hits_Y);
    _output_tree->Branch("_ALL_hits_p2_X", "std::vector<Double_t>" , &_ALL_hits_p2_X);
    _output_tree->Branch("_ALL_hits_p2_Y", "std::vector<Double_t>" , &_ALL_hits_p2_Y);


    // _simch_shower_michel_E
    return true;
    
  }
  
  bool Michel2DAna::analyze(storage_manager* storage) {
    clear_all(); 
    
    //True
    auto evt_mcshower = storage->get_data<event_mcshower>("mcreco");
    auto evt_simch    = storage->get_data<event_simch>("largeant");
    if(!evt_simch || !(evt_simch->size())) {
      print(msg::kERROR,__FUNCTION__,"SimChannel data product not found!");
      return false;
    }
    


    //Reco
    auto evt_hits     = storage->get_data<event_hit>    ("gaushit");
    auto evt_clusters = storage->get_data<event_cluster>(_cluster_producer);
    auto evt_ass_data = storage->get_data<event_ass>    (_cluster_producer);
	
    
    //Create ClusterYPlane objects, merge joining clusters
    if(!convert_2d(evt_hits,
		   evt_clusters,
		   evt_ass_data)) return false;
    
    
    //lambda returns cluster with most number of ordered_pts
    auto largest = [](const std::vector<ClusterYPlane*>& _cl)
      {
	size_t idx = 0; size_t size = 0; size_t ret = 0;
	for(const auto& c : _cl) {
	  if(c->_ordered_pts.size() > size) {
	    ret = idx; size = c->_ordered_pts.size();
	  }
	  ++idx;
	}
	return ret;
      };
    
    
    std::vector<Double_t> truncated_mean;
    std::vector<Double_t> truncated_dqds;
    
    auto c = _clusters[largest(_clusters)];
    
    if(c->_ordered_pts.size() < _min_cluster_size) //only this cluster size accepted.
      return false;
    
    //do truncated mean
    truncated_mean = r2d->windowed_means(_n_window_size,_window_cutoff,0,
					 c->_ahits, c->_ordered_pts);
    
    //cut off the bullshit on the edges of mean
    truncated_mean.erase(truncated_mean.begin(),
			 truncated_mean.begin() + _truncated_shave);
    
    truncated_mean.erase(truncated_mean.end() - _truncated_shave,
			 truncated_mean.end());
    
    int s = 3; // must be odd, currently has no setter
    
    for(int o = 0; o < s; ++o) truncated_dqds.push_back(0.0);
    
    //do smooth differentiation
    for(int i = s; i < truncated_mean.size() - s + 1; ++i) {
      std::vector<Double_t> f(truncated_mean.begin() + i - s,truncated_mean.begin() + i + s);
      std::vector<Double_t> x(c->_s.begin() + i - s + _truncated_shave,c->_s.begin() + i + s + _truncated_shave);
      truncated_dqds.push_back(r2d->smooth_derive(f,x,2*s+1));
    }
    
    for(int o = 0; o < s; ++o) truncated_dqds.push_back(0.0);
    
    auto mean_michel_vtx = r2d->DetEVtx(truncated_mean,
					truncated_dqds); //should return index of highest charge
    
    if(mean_michel_vtx.first == 999) return false;
    
    auto real_michel_vtx = r2d->REALDetEVtx(c->_ahits,
					    c->_ordered_pts,
					    mean_michel_vtx.first);
    
    auto the_vtx = size_t{0}; //reco vtx

    // bool forward;
    // if(mean_michel_vtx.first < mean_michel_vtx.second) {
    //   the_vtx = real_michel_vtx + 1;
    //   forward = true;
    // }
    // else {
    //   std::cout << " michel was behind me... " << std::endl;
    //   the_vtx = real_michel_vtx - 1;
    //   forward = false;
    // }
    // std::cout << "mean_michel_vtx  " << mean_michel_vtx.first << std::endl;
    // std::cout << "real_michel_vtx  " << real_michel_vtx << std::endl;
    
    //maybe what is better is to count the number of points before or after...    
    bool forward;
    bool ddiirr;

    ddiirr
      = (double)mean_michel_vtx.first > (double)(c->_ordered_pts.size()/2.0);

    Double_t part1 = 0;
    Double_t part2 = 0;
    
    Double_t p1 = 0;
    Double_t p2 = 0;

    for (size_t i = 0; i< c->_ordered_pts.size(); i++){
      if (i < real_michel_vtx)      p1++;
      else if (i > real_michel_vtx) p2++;
    }
    
    for (size_t i = 0; i< c->_ordered_pts.size(); i++){
      if      (i < real_michel_vtx)
	part1 += c->_ahits[c->_ordered_pts[i]].hit.Integral();
      else if (i > real_michel_vtx)
	part2 += c->_ahits[c->_ordered_pts[i]].hit.Integral();    
    }
    
    
    Double_t n_cutoff = 2;
    Double_t c_cutoff = 1.15;
    Int_t    w_cutoff = 10;

    std::cout << "\tp1/p2 : " << std::setprecision(15) << p1/p2 << " \n";
    std::cout << "\tpart1/part2 : " << std::setprecision(15) << part1/part2 << " \n";
    
    ////FIRST////
    if( p1/p2 > n_cutoff || p1/p2 < 1/n_cutoff) {
      
      if (part1 > part2) {
	ddiirr = true;
      } else {
	ddiirr = false; 
      }
      std::cout << "n_cutoff...\n";
    }
    
    ////SECOND/////
    else if (part1/part2 > c_cutoff || part1/part2 < 1/c_cutoff) {
      
      if (part1 > part2) {
	ddiirr = true;
      } else {
	ddiirr = false; 
      }
      std::cout << "c_cutoff...\n";
    }
    ////THIRD/////
    
    
    else if(p1 > w_cutoff && p2 > w_cutoff) {
      
      part1 = 0.0;
      part2 = 0.0;
      
      for (size_t i = real_michel_vtx - w_cutoff; i < real_michel_vtx + w_cutoff; i++){
    	if (i < real_michel_vtx)
    	  part1 += c->_ahits[c->_ordered_pts[i]].hit.Integral();
    	else if (i > real_michel_vtx) 
    	  part2 += c->_ahits[c->_ordered_pts[i]].hit.Integral();
	
      }
      
      if (part1 > part2) {
	ddiirr = true;
      } else {
	ddiirr = false; 
      }
      
      std::cout << "window_cutoff...\n";
      
    }
    
    
    ////FOURTH/////
    else {

      std::cout << "\n~~~~~~~~~~~!!!!!!!!!!!! WE FAIL !!!!!~~~~~~~~~~~~\n";
      ddiirr = false;
      
      
    }
    
    
    // else if(p1 > 10 && p2 > 10) {
      
    //   for (size_t i = real_michel_vtx - 10; i < real_michel_vtx + 10; i++){
    // 	if (i < real_michel_vtx)
    // 	  part1 += c->_ahits[c->_ordered_pts[i]].hit.Integral();
    // 	else if (i > real_michel_vtx) 
    // 	  part2 += c->_ahits[c->_ordered_pts[i]].hit.Integral();
	
    //   }
    //}

    ///////THIRD///////
    // else if (part1 > part2) {
    //   ddiirr = true;
    // } else {
    //   ddiirr = false; 
    // }
    
    std::cout << " \npart1 : " << part1 << " part2: " << part2 << "\n";
    std::cout << " \np1 : " << p1 << " p2: " << p2 << "\n";
    
    if(ddiirr) {
      the_vtx = real_michel_vtx + 1;
      forward = true;
    } else {
      if(real_michel_vtx != 0)
	the_vtx = real_michel_vtx - 1;
      else
	the_vtx = real_michel_vtx;
      forward = false;
    }
    
    // Get the closest reconstructed hit to the start of the mcshower
    TVector2 *proj_start = nullptr;
    if(!find_projected_start(proj_start,evt_mcshower)) //this actually updated proj_start pointer
      return false;

    auto real_michel = c->find_closest_hit(proj_start); //find closest hit to projection
    
    auto thit = c->_ahits[c->_ordered_pts[the_vtx]];    //get hit of the reco vtx point
    
    r2d->tag_michel(c,the_vtx,forward,evt_hits, _min_michel_rad);
    r2d->tag_muon(c,the_vtx,forward,evt_hits);

    _Q_u_p2 = c-> _muon-> _charge;

    //double plane_charge
    
    /////////COMPARE TO MC////////////////
    
    std::cout << "comparing the aho to MC\n";
    std::cout << "creating the g4_tackids...\n";
    //From shower quality 
    auto _mc_energy_min = 0;
    auto _mc_energy_max = 65; // ?? MeV ??
   
    std::vector<std::vector<unsigned int> > g4_trackid_v;
    std::vector<unsigned int> mc_index_v;
    g4_trackid_v.reserve(evt_mcshower->size());
    
    for(size_t mc_index=0; mc_index<evt_mcshower->size(); ++mc_index) {
	auto const& mcs = (*evt_mcshower)[mc_index];
      
      if(mcs.MotherPdgCode() == 13 &&
	 mcs.Process() == "muMinusCaptureAtRest" &&
	 mcs.DetProfile().E()/mcs.Start().E() > 0.95) {
	
	double energy = mcs.DetProfile().E();
      
	std::vector<unsigned int> id_v;
	id_v.reserve(mcs.DaughterTrackID().size());
      
	if( _mc_energy_min < energy && energy < _mc_energy_max ) {
	  for(auto const& id : mcs.DaughterTrackID()) {
	    if(id == mcs.TrackID()) continue;
	    id_v.push_back(id);
	  }
	  id_v.push_back(mcs.TrackID());
	  g4_trackid_v.push_back(id_v);
	  mc_index_v.push_back(mc_index);
	}
      
      }
    }
    if(g4_trackid_v.size() == 0)
      return false;

    std::cout << "doing ass\n";
    
    event_hit* ev_hit = nullptr;
    auto const& ass_hit_v = storage->find_one_ass(evt_clusters->id(),ev_hit,evt_clusters->name());
    
    std::cout << "Building MC map... \n";
    
    if(!fBTAlg.BuildMap(g4_trackid_v, *evt_simch, *ev_hit, ass_hit_v)) {
      print(msg::kERROR,__FUNCTION__,"Failed to build back-tracking map for MC...");
      return false;
    }
    auto aho = fBTAlg.BTAlg();
    

    // for(const auto& ii : g4_trackid_v) {
    //   std::cout << " g4_track_id: { ";
    //   for(const auto& kk : ii) {
    // 	std::cout << " : " << kk << " " ;
    //   }
    //   std::cout << " }";
    // }
    auto reco_michel_hits  = get_summed_mcshower_other(aho,c->_michel->_hits,1);

    double plane_charge = 0.0;
    
    std::vector<larlite::hit> plane2hits;
    for(const auto& h : *evt_hits) {
      if(h.View() == 2) {
	plane2hits.push_back(h);
	plane_charge += h.Integral();
	_ALL_hits_p2_X.push_back(h.WireID().Wire * 0.3);
	_ALL_hits_p2_Y.push_back(h.PeakTime() * 0.0802814);
      }
    }
    
    _Q_tot_p2 = plane_charge;
    
    auto all_hits = get_summed_mcshower_other(aho,plane2hits,1);
    
    std::vector<larlite::hit> ordered_cluster_hits;
    for(const auto& idx : c->_ordered_pts)
      ordered_cluster_hits.push_back(c->_ahits[idx].hit);
    auto ordered_hits = get_summed_mcshower_other(aho,ordered_cluster_hits,0);
    
    std::vector<larlite::hit> all_cluster_hits;
    for(const auto& h : c->_ahits)
      all_cluster_hits.push_back(h.hit);
    auto cluster_hits = get_summed_mcshower_other(aho,all_cluster_hits,0);
    
    _simch_michel_true_shower_E     = reco_michel_hits.first;
    _simch_michel_false_shower_E    = reco_michel_hits.second;
    
    _simch_plane_true_shower_E      = all_hits.first;
    _simch_plane_false_shower_E     = all_hits.second;
    
    _simch_ordered_true_shower_E    = ordered_hits.first;
    _simch_ordered_false_shower_E   = ordered_hits.second;
    
    _simch_cluster_true_shower_E    = cluster_hits.first;
    _simch_cluster_false_shower_E   = cluster_hits.second;

    
    ///////////////////////////////////////////
    ////////////WRITE OUT TO TTREE////////////
    /////////////////////////////////////////
    
    _tX = c->_ahits[c->_ordered_pts[real_michel]].vec->X();
    _tY = c->_ahits[c->_ordered_pts[real_michel]].vec->Y();
    _rX = c->_ahits[c->_ordered_pts[the_vtx]].vec->X();
    _rY = c->_ahits[c->_ordered_pts[the_vtx]].vec->Y();

    for(const auto& hit : c->_ahits) {
      _ahits_X_copy.push_back(hit.vec->X());
      _ahits_Y_copy.push_back(hit.vec->Y());
      _charges_copy.push_back(hit.hit.Integral());
    }
    
    for(const auto& pts : c->_ordered_pts)
      _ordered_pts_copy.push_back(pts);
    
    _mean_charges_copy = truncated_mean;
    _dqds_copy         = truncated_dqds;
    _s_copy            = c->_s;

    _startX = c->_start.vec->X();
    _startY = c->_start.vec->Y();
    _endX = c->_end.vec->X();
    _endY = c->_end.vec->Y();
    
    _d_m_h = c->_michel_dist;
    _michel_E = c->_michel->_charge;
    _michel_L = c->_michel->_length;
    
    _reco_Q_o_mc_Q = 0.0;
    _reco_Q_o_mc_Q = _michel_E / _simch_plane_true_shower_E;


    
    _mcQ_frac = _simch_michel_true_shower_E/( _simch_plane_true_shower_E);
    _MeV_scale = _mcQ_frac * _true_michel_Det;
    
    
    _output_tree->Fill();
    
    // don't delete these heap objects ever!!!!
    // delete evt_hits;
    // delete evt_clusters;
    // delete evt_mcshower;
    // delete evt_ass_data;

    delete proj_start;
    
    std::cout << "\n\t == Wrote event: " << _evt << "\n";
    _evt++;

    _large_frac_shower_hits_X.clear();
    _large_frac_shower_hits_Y.clear();
    _ALL_hits_p2_X.clear();
    _ALL_hits_p2_Y.clear();

    return true;
  }
  
  bool Michel2DAna::finalize() {
    
    _output_tree->Write();
    
    return true;
  }
  
  bool Michel2DAna::convert_2d(const event_hit     *evt_hits,
			       const event_cluster *evt_clusters,
			       const event_ass     *evt_ass_data) {
    
    
    AssSet_t cluster_to_hit_ass;
    
    try {
      cluster_to_hit_ass = 
	evt_ass_data->association(evt_clusters->id(),evt_hits->id()); }
    catch(DataFormatException) { return false; }
    
    size_t out_cnt = -1;
    size_t in_cnt  = 0;
    for(auto const& hit_indicies : cluster_to_hit_ass) {
      out_cnt++;
      
      if (evt_clusters->at(out_cnt).View() == 2) {
	
	std::vector<hit> the_hits; the_hits.resize(hit_indicies.size());
	
	if(the_hits.size() < 4) // control the minimum size of clusters, user should set this...
	  continue;
	
	for(unsigned int i = 0; i < the_hits.size(); ++i)
	  the_hits[i] = evt_hits->at(hit_indicies[i]);
	
	_clusters.push_back(new ClusterYPlane(the_hits,
					      evt_clusters->at(out_cnt)));
	in_cnt++;
      }
      
    }
    
    // merge them
    check_cluster_boundaries();
    
    return true;
  }
  
  void Michel2DAna::check_cluster_boundaries() {
    
    //this code will probably leak memory
    //i create new heap ClusterYPlane from addition
    //the old ones get erased from the _clusters array
    //but they surely exist somewhere right?

    //here is vic's bullshit implementation to 
    //check pairwise touching clusters
    //then use old school goto to remove them from _clusters
    //the added cluster get added to the rear with push_back(1+2)
    
    bool j = false;
    size_t a = 0;
    size_t b = 0;

    while(1) {
      for( a = 0; a < _clusters.size(); ++a) {
	for( b = 0; b < _clusters.size(); ++b) {
	  if(a != b && (_clusters[a])->touching(_clusters[b]) ) {
	    auto bb = *_clusters[a] + _clusters[b]; //real object, how do I put this in a reference of pointers??
	    _clusters.push_back(new ClusterYPlane(bb)); //wow
	    goto baka; //ouch
	  }
	}
      }
  
      if(!j) //probably not necessary
	break;
	    
      
    baka: //holy shit this is ghetto, goto past break statement
      
      _clusters.erase(_clusters.begin() + a);
      
      if(a > b) _clusters.erase(_clusters.begin() + b);
      else      _clusters.erase(_clusters.begin() + b - 1);
      
    }
    
  }
  
  void Michel2DAna::clear_all() {
    
    for (std::vector<ClusterYPlane*>::iterator it = _clusters.begin() 
	   ; it != _clusters.end(); ++it) { delete (*it); } 
    _clusters.clear();
    
    _ahits_X_copy.clear();
    _ahits_Y_copy.clear();
    _charges_copy.clear();
    _ordered_pts_copy.clear();
    _mean_charges_copy.clear();
    _dqds_copy.clear();
    _s_copy.clear();
  }

  bool Michel2DAna::find_projected_start(TVector2*& p, 
					 const event_mcshower* evt_mcshower) {
    
    TLorentzVector true_start;
    auto bb = false;
    for(const auto& shower : *evt_mcshower) {
      if (shower.Process() == "muMinusCaptureAtRest" &&
	  shower.Charge(2) > 2.0 )	  {
	true_start = shower.Start().Position();
	_true_michel_E   = shower.Charge(2);
	_true_michel_Det = shower.DetProfile().E();
	bb = true; 
      }
    }
    if(!bb)
      return false;
    
    TVector3 *ttt = new TVector3(true_start.Vect());
    
    ::larutil::PxPoint pxpoint;
    
    try{ 
      pxpoint = ::larutil::GeometryUtilities::GetME()->Get2DPointProjection2(ttt,2);
    } catch(larutil::LArUtilException) { return false; }
    
    p = new TVector2(pxpoint.w * 0.3,
		     pxpoint.t * 0.0802814);
    
    delete ttt;
    return true;
    
  }
  
  std::pair<Double_t,Double_t> Michel2DAna::get_summed_mcshower_other(const ::btutil::MCBTAlg& aho,
								      const std::vector<larlite::hit>& hits,
								      bool plane2hits) {
    Double_t baka1 = 0.0;
    Double_t baka2 = 0.0;

    
    //std::cout << "in get_summed... " << std::endl;
    //for(size_t u = 0; u < c->_ordered_pts.size(); ++u) {
    for(const auto& h : hits) {
      ::btutil::WireRange_t wire_hit(h.Channel(),h.StartTick(),h.EndTick());
      baka1 += aho.MCQ(wire_hit)[0] * _ne2ADC;
      baka2 += aho.MCQ(wire_hit)[1] * _ne2ADC;
    }
    

    if(plane2hits) {
      for(const auto& h : hits) {
	::btutil::WireRange_t wire_hit(h.Channel(),h.StartTick(),h.EndTick());
	if(aho.MCQFrac(wire_hit)[0] > 0.5) {
	  _large_frac_shower_hits_X.push_back(h.WireID().Wire * 0.3);
	  _large_frac_shower_hits_Y.push_back(h.PeakTime()  * 0.0802814);
	}
      }
    }
    
    return std::make_pair(baka1,baka2); //first is MCShower, //second is other!
    
  }

  
  
  
  
}



#endif
