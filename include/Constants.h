#ifndef GEP_CONST_H
#define GEP_CONST_H

#include "TMath.h"

namespace constant {

  // math
  static const double pi = TMath::Pi();
  static const double deg2rad = TMath::DegToRad();

  // electron
  static const double Me = 0.5109989461E-03;  // GeV

  // proton
  static const double Mp  = 0.938272081;      // GeV
  static const double mup = 2.7928473446;     // nuclear magnetons
  static const double kappap = mup - 1.;      // anomalous magnetic moment

}

#endif
