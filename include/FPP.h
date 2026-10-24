#ifndef GEP_FPP_H
#define GEP_FPP_H

#include "TVector3.h"

/*
   Focal-plane polarimeter geometry: the proton's secondary scattering in the
   analyzer, from the incoming (front tracker) and outgoing (FPP tracker) tracks.

   Coordinates: transport-like frame, x = dispersive (vertical, down),
   y = non-dispersive, z = along the central ray. Slopes xp = dx/dz, yp = dy/dz.
   Both tracks MUST be expressed in the same frame before they are compared;
   that is what ToCommonFrame() is for.
*/
namespace fpp {

  struct Track {
    double x = 0, y = 0;     // position at z = z0 (m)
    double xp = 0, yp = 0;   // slopes
    double z0 = 0;           // z of the reference plane in the common frame (m)
    bool valid = false;
    TVector3 Dir()   const { return TVector3(xp, yp, 1.).Unit(); }
    TVector3 Point() const { return TVector3(x, y, z0); }
  };

  // TODO: apply the real tracker -> common-frame transformation (offsets,
  // rotations from the GEM alignment). For now it only sets the z offset.
  void ToCommonFrame(Track &trk, double zoffset);

  // polar (theta) and azimuthal (phi) angle of the outgoing track relative to
  // the incoming one. phi = 0 along +x (dispersive), phi in (-pi, pi].
  void ScatteringAngles(const Track &in, const Track &out, double &theta, double &phi);

  // distance of closest approach between the two tracks (sclose) and the z of
  // that point (zclose, should sit inside the analyzer).
  void ClosestApproach(const Track &in, const Track &out, double &sclose, double &zclose);

}

#endif
