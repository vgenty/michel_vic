#ifndef LARLITE_TRIANGLES_CXX
#define LARLITE_TRIANGLES_CXX

#include "Triangles.h"

namespace larlite {

  bool Triangles::initialize() {
    
    fTg.resize(3);
    //tTf = new TF1("")
    
    return true;
  }
  
  bool Triangles::analyze(storage_manager* storage) {
    std::cout << "Event: " << storage->get_entries_read() << "\n";
    
    auto evt_hits = storage->get_data<event_hit>("gaushit"); 

    // >> Initialize only once
    const auto wire2cm   = ::larutil::GeometryUtilities::GetME()->WireToCm();
    const auto time2cm   = ::larutil::GeometryUtilities::GetME()->TimeToCm();

    int c = 0;
    for(auto const& h : *evt_hits) {
      // Locations of hits
      double x      = h.WireID().Wire   * wire2cm;
      double y      = h.PeakTime()      * time2cm;
      double errx   = 0.3 * scale 
      double erry   = h.SigmaPeakTime() * time2cm * scale;
      UChar_t plane = ::larutil::Geometry::GetME()->ChannelToPlane(h.Channel());

      fHits_xy[plane].push_back(std::pair<double,double>(x,y));
      fHits_xy_err[plane].push_back(std::pair<double,double>(x,y));
      // print out the plane
      // printf("%u",plane);
    }
    
    
    for(int i = 0; i < 3; ++i) {
      fTg[i] = new TGraphErrors();
      for(int j = 0; j < fHits_xy[i].size(); ++j) {
	fTg[i]->SetPoint(j,
			 fHits_xy[i][j].first,
			 fHits_xy[i][j].second);
	fTg[i]->SetPointError(j,
			      fHits_xy_err[i][j].first,
			      fHits_xy_err[i][j].second);
	
      }
    }
    
    return true;
  }
  
  bool Triangles::finalize() {
    std::cout << "In finalize... \n";
    for(auto k : fTg)
      k->Write();
    
    return true;
  }

}
#endif
