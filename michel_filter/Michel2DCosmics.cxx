#ifndef LARLITE_MICHEL2DCOSMICS_CXX
#define LARLITE_MICHEL2DCOSMICS_CXX

#include "Michel2DCosmics.h"

namespace larlite {

  bool Michel2DCosmics::initialize() {

    _wire2cm   = ::larutil::GeometryUtilities::GetME()->WireToCm();
    _time2cm   = ::larutil::GeometryUtilities::GetME()->TimeToCm();    

    _output_tree = new TTree("out_tree","aho_tree");
    
    _output_tree->Branch("ahits_X"     , "std::vector<std::vector<Double_t> >" 
			 , &_ahits_X_copy);
    _output_tree->Branch("ahits_Y"     , "std::vector<std::vector<Double_t> >" 
			 , &_ahits_Y_copy);
    _output_tree->Branch("ordered_pts" ,  "std::vector<std::vector<size_t> >" 
			 , &_ordered_pts_copy);
    
    _output_tree->Branch("true_X", "std::vector<Double_t>", &_tX);
    _output_tree->Branch("true_Y", "std::vector<Double_t>", &_tY);
    
    return true;
  }
  
  bool Michel2DCosmics::analyze(storage_manager* storage) {
    clear_all();
    
    //we now must do a cluster-wise michel search,
    //cosmic data now have \emph{many} clusters
    
    //True
    auto evt_mcshower = storage->get_data<event_mcshower>("mcreco");
    
    //Reco
    auto evt_hits     = storage->get_data<event_hit>    ("gaushit");
    auto evt_clusters = storage->get_data<event_cluster>(_cluster_producer);
    auto evt_ass_data = storage->get_data<event_ass>    (_cluster_producer);
    
    
    //may have to remove the "cluster merger!"
    if(!convert_2d(evt_hits,
		   evt_clusters,
		   evt_ass_data)) return false;
    
    
    //get the XZ locations of the michels
    std::vector<TVector2*> proj_starts;
    
    if(!find_projected_starts(proj_starts,evt_mcshower))
      return false;
    
    std::cout << "found " << proj_starts.size() << " number of michels " << "\n";

    for(const auto& p : proj_starts) {
      _tX.push_back(p->X());
      _tY.push_back(p->Y());
    }
    
    //Do clusterwise michel find....
    for(auto& c : _clusters) { //not const we will change them...
      std::vector<TVector2*>::iterator itr;
      bool found = false;
      for(itr = proj_starts.begin(); itr != proj_starts.end(); ++itr) {
    	auto r = c->match(*itr);

	if(r == 2) { //i already have michel
	  found = true;
	  break;
	}
	if(r == 1)  {//i matched with accepted parameters
	  found = true;
	  break;
	}
	if(r == 0) { //michel wasn't right for me continue
	  continue;
	}
      }
      
      if(!found)
	continue;
      std::cout << "current size of proj_starts " << proj_starts.size() << "\n";
      proj_starts.erase(itr); //does this leak memory since these are pointers? beyond me
      if(proj_starts.size() == 0)
	break;
    }
    
    std::cout << "I have " << _clusters.size() << " number of clusters, lets dump them";
    for(const auto& c : _clusters)
      if(c->_has_michel) c->dump();
    
    std::vector<Double_t> XX;
    std::vector<Double_t> YY;
    for(const auto* c   : _clusters) {
      for(const auto& hit : c->_ahits) {
	XX.push_back(hit.vec->X());
	YY.push_back(hit.vec->Y());
	//_charges_copy.push_back(hit.hit.Integral());
      }
      _ordered_pts_copy.push_back(c->_ordered_pts);
      _ahits_X_copy.push_back(XX);
      _ahits_Y_copy.push_back(YY);
      XX.clear();
      YY.clear();
    }
    std::cout << "\n\t == Wrote event: " << _evt << "\n";
    _output_tree->Fill();
    for(auto& p : proj_starts)
      delete p;
    
    _ahits_X_copy.clear();
    _ahits_Y_copy.clear();
    _ordered_pts_copy.clear();


    _tX.clear();
    _tY.clear();
    _evt++;
    return true;
  }
  
  bool Michel2DCosmics::finalize() {
    _output_tree->Write();
    return true;
  }

  bool Michel2DCosmics::convert_2d(const event_hit     *evt_hits,
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

    // merge them...
    std::cout << "I have " << _clusters.size() << " before touching\n";
    
    check_cluster_boundaries();
    std::cout << "I have " << _clusters.size() << "AFTER:) touching\n";
    return true;
  }
  void Michel2DCosmics::check_cluster_boundaries() {
    
    std::vector<ClusterYPlane*>::iterator itr1;    
    std::vector<ClusterYPlane*>::iterator itr2;    
    std::vector<ClusterYPlane*>::iterator trueitr;    
    
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
    // std::cout << _clusters.size() << "\n";
    // std::cout << "\tclusters BEFORE merge\n";
    // for(const auto& _c : _clusters)
    //   _c->dump();

    while(1) {
      //std::cout << "While 1\n";
      for( a = 0; a < _clusters.size(); ++a) {
	//for(itr1 = _clusters.begin(); itr1 != _clusters.end(); ++itr1) {
	//std::cout << "a " << a << " ";
	//	for(itr2 = _clusters.begin(); itr2 != _clusters.end(); ++itr2) {
	for( b = 0; b < _clusters.size(); ++b) {
	  // std::cout << "b " << b << " ";
	  //	  if(itr1 != itr2 && (*itr1)->touching(*itr2) ) {
	  if(a != b && (_clusters[a])->touching(_clusters[b]) ) {
	    //std::cout << "\nc\n ";
	    //std::cout << "Found two touchers \n";
	    auto bb = *_clusters[a] + _clusters[b]; //real object, how do I put this in a reference of pointers??
	    _clusters.push_back(new ClusterYPlane(bb)); //wow
	    goto baka; //ouch
	  }
	}
      }
  
      if(!j) //probably not necessary
	break;
	    
      
    baka: //holy shit this is ghetto, goto past break statement

      Int_t p = 2;

      _clusters.erase(_clusters.begin() + a);
      
      if(a > b) _clusters.erase(_clusters.begin() + b);
      else      _clusters.erase(_clusters.begin() + b - 1);
      
    }
    
    // std::cout << "\tclusters after merge\n";
    // for(const auto& _c : _clusters)
    //   _c->dump();
    
    
  }
  
  void Michel2DCosmics::clear_all() {
    
    for (std::vector<ClusterYPlane*>::iterator it = _clusters.begin() 
	   ; it != _clusters.end(); ++it) { delete (*it); } 
    _clusters.clear();
  }
  
  
  bool Michel2DCosmics::find_projected_starts(std::vector<TVector2*>& p, 
					      const event_mcshower* evt_mcshower) {
    
    //std::cout << Form("Afinding the projected start point...\n");
    std::vector<TLorentzVector> true_starts;
    std::vector<Double_t> true_times;
    auto bb = false;
    for(const auto& shower : *evt_mcshower) {
      if (shower.Process() == "muMinusCaptureAtRest" &&
	  shower.Charge(2) > 10.0 )	  {
	true_starts.push_back(shower.Start().Position());
	true_times.push_back(shower.Start().T());
	std::cout << shower.Start().T() << std::endl;
	bb = true; 
      }
    }
    if(!bb)
      return false;
    
        
    for(const auto& starts : true_starts) {
      TVector3* ttt = new TVector3(starts.Vect());
      ::larutil::PxPoint pxpoint;
      try{ 
	pxpoint = ::larutil::GeometryUtilities::GetME()->Get2DPointProjection2(ttt,2);      
      } catch(larutil::LArUtilException) { continue; }
      
      
      p.push_back(new TVector2(pxpoint.w*0.3,
			       pxpoint.t*0.0802814));
      
      delete ttt;
    }
    
    if(p.size() == 0)
      return false;

    return true;
    
  }
  
  
}
#endif
