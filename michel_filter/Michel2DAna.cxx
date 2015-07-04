#ifndef LARLITE_MICHEL2DANA_CXX
#define LARLITE_MICHEL2DANA_CXX

#include "Michel2DAna.h"

namespace larlite {
 
  bool Michel2DAna::initialize() {
    
    _wire2cm   = ::larutil::GeometryUtilities::GetME()->WireToCm();
    _time2cm   = ::larutil::GeometryUtilities::GetME()->TimeToCm();    
    
    _output_tree = new TTree("out_tree","aho_tree");
    _output_tree->Branch("true_X",&_tX,"tX/D");
    _output_tree->Branch("true_Y",&_tY,"tY/D");
    _output_tree->Branch("reco_X",&_rX,"rX/D");
    _output_tree->Branch("reco_Y",&_rY,"rY/D");
    
    return true;
    
  }
  
  bool Michel2DAna::analyze(storage_manager* storage) {
    clear_all(); // clear out class variables
    std::cout << "\n\t\tOn event... " << storage->event_id() << "\n";
    
    //True
    auto evt_mcshower = storage->get_data<event_mcshower>("mcreco");
    
    //Reco
    auto evt_hits     = storage->get_data<event_hit>    ("gaushit");
    auto evt_clusters = storage->get_data<event_cluster>(_cluster_producer);
    auto evt_ass_data = storage->get_data<event_ass>    (_cluster_producer);
	
    
    if(!convert_2d(evt_hits,
		   evt_clusters,
		   evt_ass_data)) return false;
    
    
    //std::cout << "about to do the windowed means.... \n";
    
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
    

    std::vector<Double_t> b;
    std::vector<Double_t> baka;

    //std::cout << "c  size: " << _clusters.size() << "\n";
    
    auto c = _clusters[largest(_clusters)];
    
    if(c->_ordered_pts.size() < 25)
      return false;

    b = r2d->windowed_means(25,0.25,0,
			    c->_ahits, c->_ordered_pts);
      
    //ok lets cut off the bullshit on the edges of mean
    
    size_t soff = 2;
    b.erase(b.begin()   ,b.begin()+soff);
    b.erase(b.end()-soff,b.end()       );
    
    int s = 3;
      
    //std::cout << "smooth_deriving................          \n";

    for(int o = 0; o < s; ++o)
      baka.push_back(0.0);

    for(int i = s; i < b.size() - s + 1; ++i) {
      std::vector<Double_t> f(b.begin() + i - s,b.begin() + i + s);
      std::vector<Double_t> x(c->_s.begin() + i - s + soff,c->_s.begin() + i + s + soff);
      baka.push_back(r2d->smooth_derive(f,x,2*s+1));
    }
      
    for(int o = 0; o < s; ++o)
      baka.push_back(0.0);
    
    // b holds mean charges
    // baka holds dqdx
    // do the reco
    
    
    //std::cout << " the size of _orderdpts is " << c->_ordered_pts.size() << "\n";
    auto mean_michel_vtx = r2d->DetEVtx(b,baka); //should return index in charge with highest cham
    
    
    // std::cout << "mean_michel_vtx : first ~ "
    // 	      << mean_michel_vtx.first  << " second ~"
    // 	      << mean_michel_vtx.second << " \n";
    
    if(mean_michel_vtx.first == 999)
      return false;
    
    auto real_michel_vtx = r2d->REALDetEVtx(c->_ahits,
					    c->_ordered_pts,
					    mean_michel_vtx.first);
    
    auto the_vtx = size_t{0};
    
    if(mean_michel_vtx.first < mean_michel_vtx.second)
      the_vtx = real_michel_vtx + 1;
    else
      the_vtx = real_michel_vtx - 1;
    
    
    // Get the closest reconstructed hit to the start of the mcshower
    TLorentzVector true_start;
    auto bb = false;
    for(const auto& shower : *evt_mcshower) {
      if (shower.Process() == "muMinusCaptureAtRest" &&
	  shower.Charge(2) > 2.0 )	  {
	true_start = shower.Start().Position();
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
    TVector2 *proj_start = new TVector2(pxpoint.w*_wire2cm,
					pxpoint.t*_wire2cm);
    

    
    //get the closest hit in ordered_pts to the projected start...
    auto real_michel = c->find_closest_hit(proj_start);

    // std::cout << "\tthe_vtx : " << the_vtx
    // 	      << "\treal_mic: " << real_michel;
    
    // std::cout << "wire : " << proj_start->X() << " time: " << proj_start->Y() << "\n";
    
    
    // std::cout << "maybe we found the michel point we don't know...";
    // std::cout << real_michel_vtx << "\n";
    auto thit = c->_ahits[c->_ordered_pts[the_vtx]];
    std::cout << "it could be at... " << the_vtx
	       << ": (" << thit.vec->X() << "," << thit.vec->Y() << ")\n";
    
    
    _tX = c->_ahits[c->_ordered_pts[real_michel]].vec->X();
    _tY = c->_ahits[c->_ordered_pts[real_michel]].vec->Y();
    _rX = c->_ahits[c->_ordered_pts[the_vtx]].vec->X();
    _rY = c->_ahits[c->_ordered_pts[the_vtx]].vec->Y();

    
    _output_tree->Fill();
    
    // don't delete these heap objects...
    // delete evt_hits;
    // delete evt_clusters;
    // delete evt_mcshower;
    //delete evt_ass_data;

    delete proj_start;
    delete ttt;

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

    // merge them...
    
    check_cluster_boundaries();

    return true;
  }
  
  void Michel2DAna::check_cluster_boundaries() {
    
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
	    std::cout << "Found two touchers \n";
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
  
  void Michel2DAna::clear_all() {
    
    for (std::vector<ClusterYPlane*>::iterator it = _clusters.begin() 
	   ; it != _clusters.end(); ++it) { delete (*it); } 
    _clusters.clear();
    
  }
}
#endif
