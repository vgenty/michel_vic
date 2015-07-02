#ifndef LARLITE_MICHEL2DANA_CXX
#define LARLITE_MICHEL2DANA_CXX

#include "Michel2DAna.h"

namespace larlite {
 
  bool Michel2DAna::initialize() {
    
    size_t nplanes = 3;
    
    _wire2cm   = ::larutil::GeometryUtilities::GetME()->WireToCm();
    _time2cm   = ::larutil::GeometryUtilities::GetME()->TimeToCm();    
    
    //_output_tree = new TTree("out_tree","aho_tree");
    //_output_tree->Branch("charge_sum",&_sss,"sss/D");
    //std::cout << "initialize with " << _output_tree->GetNbranches() << " branches\n";
    

    tgMeans  = new TGraph();
    tgPoints = new TGraph();
    tgdqdx   = new TGraph();
    
    return true;
    
  }
  
  bool Michel2DAna::analyze(storage_manager* storage) {
    
    auto evt_hits     = storage->get_data<event_hit>    ("gaushit");
    auto evt_clusters = storage->get_data<event_cluster>(_cluster_producer);
    auto evt_ass_data = storage->get_data<event_ass>    (_cluster_producer);
	
    
    if(!convert_2d(evt_hits,
		   evt_clusters,
		   evt_ass_data)) return false;
    

    std::cout << " I created the clusters...   ";
    
    if(_clusters.size() != 1) {
      std::cout << "too many big bab\|n";
      return false;
    }
    
    std::cout << "about to do the windowed means.... \n";
    for(const auto& c : _clusters) {
      auto b = r2d->windowed_means(25,0.25,0,
				   c->_ahits, c->_ordered_pts);
      
      
      std::cout << "windowed means gave me a b............       ";
      std::cout << "do b and _s have same size? b.size " 
		<< b.size() << " and _s.size "  << c->_s.size() << "\n";

      for(int i = 0; i < b.size(); ++i) {
	tgMeans ->SetPoint(i,c->_s[i],b[i]);
      	tgPoints->SetPoint(i,c->_ahits[c->_ordered_pts[i]].vec->X(),
			   c->_ahits[c->_ordered_pts[i]].vec->Y());
      }
      
      tgMeans->SetName("q");
      tgPoints->SetName("points");
      
      int s = 3;
      std::vector<Double_t> baka;
      std::cout << "smooth_deriving................          \n";
      for(int i = s; i < b.size() - s + 1; ++i) {
	std::vector<Double_t> f(b.begin() + i - s,b.begin() + i + s);
	std::vector<Double_t> x(c->_s.begin() + i - s,c->_s.begin() + i + s);
	
	baka.push_back(r2d->smooth_derive(f,x,2*s+1));
      }
      
      
      tgdqdx->SetName("dqds");
      for(int i = 0; i < baka.size(); ++i) {
	tgdqdx->SetPoint(i,i,baka[i]);
      }

    }
    
    
    //_output_tree->Fill();
    clear_all();
    
  }
  
  bool Michel2DAna::finalize() {
    
    tgMeans->Write();
    tgPoints->Write();
    tgdqdx->Write();
    //_output_tree->Write();
      
    //    RecoMethods::getInstance().aho();
    
    
    return true;
  }
  
  bool Michel2DAna::convert_2d(const event_hit     *evt_hits,
			       const event_cluster *evt_clusters,
			       const event_ass     *evt_ass_data) {
    

    AssSet_t cluster_to_hit_ass;

    //std::cout << Form("how big are the clusters %d\n",evt_clusters->size());
    //_clusters.resize(evt_clusters->size());
    
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
	
	//std::cout << "size of the hits: " << the_hits.size() << "\n";
	
	if(the_hits.size() < 4) // control the minimum size of clusters, hard set
	  continue;
	
	for(unsigned int i = 0; i < the_hits.size(); ++i)
	  the_hits[i] = evt_hits->at(hit_indicies[i]);
	
	
	// _clusters[in_cnt] = new ClusterYPlane(the_hits,
	// 				      evt_clusters->at(out_cnt));

	_clusters.push_back(new ClusterYPlane(the_hits,
					      evt_clusters->at(out_cnt)));
	
	in_cnt++;
      }
      
    }
    
    // merge them...
    
    check_cluster_boundaries();
    
    
    //clear_all();
    return true;
  }
  
  void Michel2DAna::check_cluster_boundaries() {
    
    //std::cout << "Checking all cluster boundaries\n";
    bool aho = true;
    
    std::vector<ClusterYPlane*>::iterator itr1;    
    std::vector<ClusterYPlane*>::iterator itr2;    
    std::vector<ClusterYPlane*>::iterator trueitr;    
    
    //this code will probably leak memory
    //i create new heap ClusterYPlane from addition
    //the old ones get erased from the _clusters array
    //but they surely exist somewhere right?

    //here is vic's bullshit implementation to 
    //check pairwise touching clusters
    //then use old school goto to remove them form _clusters
    //the added cluster get added to the rear with push_back(1+2)
    
    bool j = false;
    size_t a = 0;
    size_t b = 0;
    // std::cout << _clusters.size() << "\n";
    std::cout << "\tclusters BEFORE merge\n";
    for(const auto& _c : _clusters)
      _c->dump();

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
      
      //auto p1 = *itr1;
      //auto p2 = *itr2;
      
      auto p1 = a;
      auto p2 = b;

      Int_t p = 2;


      _clusters.erase(_clusters.begin() + a);
      
      if(a > b) _clusters.erase(_clusters.begin() + b);
      else      _clusters.erase(_clusters.begin() + b - 1);
	
      // while(p) {
      // 	for(size_t y = 0; y < _clusters.size(); ++y) {
      // 	  if(a == _clusters[y]) {
      // 	    _clusters.erase(_clusters.begin() + y);
      // 	    p--;
      // 	    break;
      // 	  }
      // 	  if(p2 == _clusters[y]) {
      // 	    _clusters.erase(_clusters.begin() + y);
      // 	    p--;
      // 	    break;
      // 	  }
      // 	}
      //	 }
      
      
      
      
      
    }
    
    std::cout << "\tclusters after merge\n";
    for(const auto& _c : _clusters)
      _c->dump();
    
    
  
  }
  
  void Michel2DAna::clear_all() {
    
    _clusters.clear();
    
  }
}
#endif
