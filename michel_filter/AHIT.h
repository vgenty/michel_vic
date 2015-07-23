#ifndef AHIT_H
#define AHIT_H

#include "DataFormat/hit.h"
#include "TVector2.h"

struct ahit {
  larlite::hit hit;
  TVector2 vec;
};


#endif
