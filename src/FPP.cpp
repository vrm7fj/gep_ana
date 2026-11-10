#include <cmath>
#include "TMath.h"
#include "../include/FPP.h"
#include "../include/Constants.h"

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

  bool ConeTest(const Track &in, double zc, double zback, double theta, const std::vector<double> &area)   {

    double x_close = in.x + in.xp * (zc - in.z0);
    double y_close = in.y + in.yp * (zc - in.z0);

    double xp_scat_plus = TMath::Tan( TMath::ATan(in.xp) + theta);
    double xp_scat_minus = TMath::Tan( TMath::ATan(in.xp) - theta);
    double yp_scat_plus = TMath::Tan( TMath::ATan(in.yp) + theta);
    double yp_scat_minus = TMath::Tan( TMath::ATan(in.yp) - theta);

    double x_max = x_close + xp_scat_plus * (zback - zc);
    double x_min = x_close + xp_scat_minus * (zback - zc);
    double y_max = y_close + yp_scat_plus * (zback - zc);
    double y_min = y_close + yp_scat_minus * (zback - zc);

    return x_min > area[0] && x_max < area[1] && y_min > area[2] && y_max < area[3];

  }

  //--------------------------------------------
  // focal-plane direction (optics-frame slopes) expressed in the SBS frame, as in SBSEArm::InitOpticsAxes()
  static TVector3 FpDirSBS(double xpfp_r, double ypfp_r, double opticsangle) {
    const TVector3 X(std::cos(opticsangle), 0., std::sin(opticsangle));
    const TVector3 Y(0., 1., 0.);
    const TVector3 Z(-std::sin(opticsangle), 0., std::cos(opticsangle));
    TVector3 n = TVector3(xpfp_r, ypfp_r, 1.).Unit();
    return n.X() * X + n.Y() * Y + n.Z() * Z;
  }

  //--------------------------------------------
  double BendAngle(double xptar, double yptar, double xpfp_r, double ypfp_r, double opticsangle) {
    TVector3 ntg = TVector3(xptar, yptar, 1.).Unit();
    double c = FpDirSBS(xpfp_r, ypfp_r, opticsangle).Dot(ntg);
    if (c > 1.) c = 1.; else if (c < -1.) c = -1.;
    return std::acos(c);
  }

  //--------------------------------------------
  void BendComponents(double xptar, double yptar, double xpfp_r, double ypfp_r, double opticsangle,
                      double &dth, double &dph) {
    TVector3 nfp = FpDirSBS(xpfp_r, ypfp_r, opticsangle);
    dth = std::atan2(nfp.X(), nfp.Z()) - std::atan(xptar);   // dispersive (vertical) plane
    dph = std::atan2(nfp.Y(), nfp.Z()) - std::atan(yptar);   // non-dispersive (horizontal) plane
  }

  //--------------------------------------------
  double PrecessionAngle(double p, double thetabend) {
    if (p <= 0.) return -1.;
    double gamma = std::sqrt(p * p + constant::Mp * constant::Mp) / constant::Mp;
    return gamma * constant::kappap * thetabend;
  }

}
