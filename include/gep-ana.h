/*
   gep-ana.h includes all the gep-ana libraries. Including this one file in an
   analysis macro is enough to get access to everything (same pattern as
   gmn-ana/include/gmn-ana.h by P. Datta).
*/
#ifndef GEP_ANA_H
#define GEP_ANA_H

#include "Constants.h"                // namespace constant
#include "JSONManager.h"              // class JSONManager (header-only)
#include "../src/GEpConfig.cpp"       // class GEpConfig
#include "../src/SetROOTVar.cpp"      // namespace setrootvar
#include "../src/KinematicVar.cpp"    // namespace kine
#include "../src/Cut.cpp"             // namespace cut
#include "../src/FPP.cpp"             // namespace fpp
#include "../src/Polarization.cpp"    // namespace pol   <-- polarization observables (to implement)
#include "../src/Utilities.cpp"       // namespace util

/* --- List of gep-ana libraries ---
   Constants.h    : namespace constant   (masses, magnetic moments)
   JSONManager.h  : class JSONManager    (reads the JSON config)
   GEpConfig.h    : class GEpConfig      (kinematic setting: beam energy, arm angles, distances)
   SetROOTVar.h   : namespace setrootvar (branch address helpers)
   KinematicVar.h : namespace kine       (elastic ep kinematics, projection to ECAL)
   Cut.h          : namespace cut        (elliptical, window and fiducial cuts)
   FPP.h          : namespace fpp        (FPP scattering angles, closest approach)
   Polarization.h : namespace pol        (helicity asymmetries, P_t/P_l, GE/GM) -- stubs
   Utilities.h    : namespace util       (histogram helpers)
*/

#endif
