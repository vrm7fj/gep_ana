#include <cmath>
#include "TVector2.h"
#include "../include/Constants.h"
#include "../include/KinematicVar.h"

namespace kine {

  using constant::Mp;

  //--------------------------------------------
  double pp_elastic(double E, double thp) {
    // p_p = 2 M E (E+M) cos(theta_p) / (M^2 + 2 M E + E^2 sin^2(theta_p))
    double s = std::sin(thp), c = std::cos(thp);
    return 2. * Mp * E * (E + Mp) * c / (Mp * Mp + 2. * Mp * E + E * E * s * s);
  }
  //--------------------------------------------
  double Q2_from_pp(double pp) {
    double Tp = std::sqrt(pp * pp + Mp * Mp) - Mp;
    return 2. * Mp * Tp;
  }
  //--------------------------------------------
  double eprime_from_pp(double E, double pp) {
    double Tp = std::sqrt(pp * pp + Mp * Mp) - Mp;
    return E - Tp;
  }
  //--------------------------------------------
  double etheta_from_Q2(double E, double Ep, double Q2) {
    // Q2 = 4 E E' sin^2(theta_e/2)
    double arg = Q2 / (4. * E * Ep);
    if (arg < 0. || arg > 1.) return -1.;
    return 2. * std::asin(std::sqrt(arg));
  }
  //--------------------------------------------
  double tau(double Q2) { return Q2 / (4. * Mp * Mp); }
  //--------------------------------------------
  double epsilon(double Q2, double the) {
    double t = std::tan(the / 2.);
    return 1. / (1. + 2. * (1. + tau(Q2)) * t * t);
  }
  //--------------------------------------------
  TVector3 eprime_dir(double E, const TVector3 &pp_vec) {
    TVector3 k(0., 0., E);
    return (k - pp_vec).Unit();
  }
  //--------------------------------------------
  void SetDetAxes(double theta, int side, std::vector<TVector3> &axes) {
    // same convention as gmn-ana SetHCALaxes (side = -1 reproduces it)
    TVector3 zaxis(side * std::sin(theta), 0., std::cos(theta));
    TVector3 xaxis(0., -1., 0.);                 // -Y of Hall = X (vertical, down) of detector
    TVector3 yaxis = zaxis.Cross(xaxis).Unit();
    axes = {xaxis, yaxis, zaxis};
  }
  //--------------------------------------------
  TVector3 DetOrigin(double dist, double voff, double hoff, const std::vector<TVector3> &axes) {
    return dist * axes[2] + voff * axes[0] + hoff * axes[1];
  }
  //--------------------------------------------
  void GetxyExpect(const TVector3 &vertex, const TVector3 &dir, const TVector3 &origin,
                   const std::vector<TVector3> &axes, double &xexp, double &yexp) {
    // ray-plane intersection
    double s = (origin - vertex).Dot(axes[2]) / dir.Dot(axes[2]);
    TVector3 hit = vertex + s * dir;
    xexp = (hit - origin).Dot(axes[0]);
    yexp = (hit - origin).Dot(axes[1]);
  }
  //--------------------------------------------
  TVector3 DetToHall(double x, double y, const TVector3 &origin, const std::vector<TVector3> &axes) {
    return origin + x * axes[0] + y * axes[1];
  }
  //--------------------------------------------
  double Coplanarity(double phi_e, double phi_p) {
    return TVector2::Phi_mpi_pi(phi_e - phi_p - constant::pi);
  }

}
