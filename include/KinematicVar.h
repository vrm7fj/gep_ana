#ifndef GEP_KINEMATIC_VAR_H
#define GEP_KINEMATIC_VAR_H

#include <vector>
#include "TVector3.h"

/*
   Elastic ep kinematics for GEp. The proton is measured in SBS (beam right)
   and the electron in ECAL (beam left). Everything is driven by the proton
   angle, which SBS measures best.
   Hall CoS: z along beam, y up, x to beam left.
*/
namespace kine {

  // --- elastic ep from the proton angle ---
  double pp_elastic(double ebeam, double theta_p);           // expected proton momentum p(theta_p)
  double Q2_from_pp(double pp);                              // Q2 = 2 M T_p
  double eprime_from_pp(double ebeam, double pp);            // E' = E - T_p
  double etheta_from_Q2(double ebeam, double eprime, double Q2);
  double tau(double Q2);
  double epsilon(double Q2, double theta_e);                 // virtual photon polarization
  TVector3 eprime_dir(double ebeam, const TVector3 &pp_vec); // unit vector of e' = k - p_p

  // --- detector coordinate systems ---
  // side = +1 for beam left (ECAL), -1 for beam right (SBS/HCAL).
  // axes[0] = x (vertical, down), axes[1] = y (horizontal), axes[2] = z (central ray)
  void SetDetAxes(double theta_rad, int side, std::vector<TVector3> &axes);
  TVector3 DetOrigin(double dist, double voffset, double hoffset,
                     const std::vector<TVector3> &axes);

  // expected (x,y) where a ray from vertex along dir hits the detector plane
  void GetxyExpect(const TVector3 &vertex, const TVector3 &dir, const TVector3 &origin,
                   const std::vector<TVector3> &axes, double &xexp, double &yexp);
  // detector (x,y) -> Hall CoS point
  TVector3 DetToHall(double x, double y, const TVector3 &origin, const std::vector<TVector3> &axes);

  // coplanarity: phi_e - phi_p - pi wrapped into [-pi, pi]
  double Coplanarity(double phi_e, double phi_p);

}

#endif
