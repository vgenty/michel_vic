#ifndef MICHEL_CXX
#define MICHEL_CXX

#include "Michel.h"

void Michel::dump() {
  std::cout << "\n\t -- Michel --\n";
  std::cout << "\tenergy : " << _charge << " length: " << _length << "\n";
  std::cout << "\tat position : (" << _start.X() << "," << _start.Y() << ")" << std::endl;
  std::cout << "\twith number of hits " << _num_hits << std::endl;
}

#endif
