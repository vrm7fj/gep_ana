#include <cmath>
#include "../include/FPP.h"

namespace fpp {

  //--------------------------------------------
  void ToCommonFrame(Track &trk, double zoffset) {
    // TODO: replace with the full transformation from your alignment
    // (position offsets + rotation of the FPP tracker relative to the front tracker).
    trk.z0 = zoffset;
  }
  //--------------------------------------------
  void ScatteringAngles(const Track &in, const Track &out, double &theta, double &phi) {
    TVector3 zh = in.Dir();
    TVector3 xd(1., 0., 0.);                    // dispersive direction
    TVector3 yh = zh.Cross(xd).Unit();          // perpendicular to incoming track
    TVector3 xh = yh.Cross(zh).Unit();          // ~ +x, perpendicular to incoming track
    TVector3 o  = out.Dir();
    theta = std::acos(std::max(-1., std::min(1., o.Dot(zh))));
    phi   = std::atan2(o.Dot(yh), o.Dot(xh));
    // NOTE: check this phi convention (sign and zero) against the one used for
    // the analyzing power you will use, before fitting asymmetries.
  }
  //--------------------------------------------
  void ClosestApproach(const Track &in, const Track &out, double &sclose, double &zclose) {
    TVector3 p1 = in.Point(), d1 = in.Dir();
    TVector3 p2 = out.Point(), d2 = out.Dir();
    TVector3 w0 = p1 - p2;
    double a = d1.Dot(d1), b = d1.Dot(d2), c = d2.Dot(d2);
    double d = d1.Dot(w0), e = d2.Dot(w0);
    double den = a * c - b * b;
    if (std::fabs(den) < 1e-12) { sclose = (w0 - d * d1).Mag(); zclose = p1.Z(); return; } // parallel
    double s = (b * e - c * d) / den;
    double t = (a * e - b * d) / den;
    TVector3 c1 = p1 + s * d1, c2 = p2 + t * d2;
    sclose = (c1 - c2).Mag();
    zclose = 0.5 * (c1.Z() + c2.Z());
  }

}
