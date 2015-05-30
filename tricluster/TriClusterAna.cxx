#ifndef LARLITE_TRICLUSTERANA_CXX
#define LARLITE_TRICLUSTERANA_CXX

#include "TriClusterAna.h"

namespace larlite {

  bool TriClusterAna::initialize() {
    
    size_t nplanes = 3;
    fArtist.set_num_triangles(nplanes);

    fWire2cm   = ::larutil::GeometryUtilities::GetME()->WireToCm();
    fTime2cm   = ::larutil::GeometryUtilities::GetME()->TimeToCm();    

    return true;
  }
  
  bool TriClusterAna::analyze(storage_manager* storage) {
    //Event loop
    fArtist.erase();

    std::cout << "Event: " << storage->get_entries_read() << "\n";
    auto evt_hits = storage->get_data<event_hit>("gaushit");
    
    if(evt_hits->size() < 2) return true;
    
    for(auto const& h: *evt_hits) {
      //Which plane
      int plane = ::larutil::Geometry::GetME()->ChannelToPlane(h.Channel());
      
      
      // xy coordinates
      double x = h.WireID().Wire * fWire2cm;
      double y = h.PeakTime()    * fTime2cm;
      
      
      //Arbitrary scaling of x/y error bars. Hit's with low charge
      //deposition are ``least important" to the fitter
      double scale = 50.0/h.PeakAmplitude();
      double errx  = scale;
      double erry  = scale;
      
      fArtist.set_hit_point    (plane,std::make_pair(x,y));
      fArtist.set_hit_point_err(plane,std::make_pair(errx,erry));
      fArtist.set_hit_charge   (plane,h.Integral());
    }
    
    for(int i = 0; i < 3; ++i) { // For each plane, artist create triangle
      fArtist.fit_hits(i);
      fArtist.create_triangle(i);
    }
    
    fArtist.choose_best_triangle();
    
    
    return true;
  }

  bool TriClusterAna::finalize() {
    std::cout << "In finalize\n";
    
    return true;
  }

}
#endif
