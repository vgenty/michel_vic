#ifndef LARLITE_MICHELFILTER_CXX
#define LARLITE_MICHELFILTER_CXX

#include "MichelFilter.h"

namespace larlite {

  bool MichelFilter::initialize() {

    return true;
  }
  
  bool MichelFilter::analyze(storage_manager* storage) {
    auto evt_pars = storage->get_data<event_>("gaushit");    
    
    return true;
  }

  bool MichelFilter::finalize() {
  
    return true;
  }

}
#endif
