#ifndef GEP_POLARIZATION_H
#define GEP_POLARIZATION_H

#include <string>
#include <vector>
#include "TH1D.h"
#include "TTree.h"

/*
   Polarization observables. Two routes are laid out:

   (A) Binned:  phi_fpp histograms per helicity and theta_fpp bin
                -> helicity asymmetry -> fit -> focal-plane polarization
                -> spin transport -> P_t, P_l at the target -> mu GE/GM
   (B) Event-by-event maximum likelihood on the output tree Tout, which keeps
       theta_fpp, phi_fpp, helicity, the proton target variables and kinematics.

   Functions marked TODO are stubs: they print a message and return an
   "invalid" result so the macro still runs end to end.
*/
namespace pol {

  //------------------------------------------------------------------
  // Bookkeeping: phi_fpp distributions for h = +1 and h = -1, per theta_fpp bin
  class FPPPhiHistos {
  public:
    void Book(const std::string &name, int nphi, const std::vector<double> &theta_edges);
    void Fill(double theta_fpp, double phi_fpp, int helicity);   // ignores helicity == 0
    void Write() const;
    int  NThetaBins() const { return int(fHp.size()); }
    TH1D *Hplus(int i)  const { return fHp[i]; }
    TH1D *Hminus(int i) const { return fHm[i]; }
    double ThetaLow(int i)  const { return fEdges[i]; }
    double ThetaHigh(int i) const { return fEdges[i + 1]; }
  private:
    std::vector<double> fEdges;
    std::vector<TH1D*> fHp, fHm;
  };

  //------------------------------------------------------------------
  struct FocalPlanePol {   // transverse components at the focal plane (x = dispersive)
    double Px = 0, Px_err = 0, Py = 0, Py_err = 0;
    bool ok = false;
  };
  struct TargetPol {       // transferred polarization at the target
    double Pt = 0, Pt_err = 0, Pl = 0, Pl_err = 0;
    bool ok = false;
  };
  struct Ratio {           // mu_p GE/GM
    double R = 0, R_err = 0;
    bool ok = false;
  };

  //------------------------------------------------------------------
  // (A) binned route
  TH1D *HelicityAsymmetry(const TH1D *hplus, const TH1D *hminus, const std::string &name); // TODO
  double AnalyzingPower(double theta_fpp, double pp);                                      // TODO
  FocalPlanePol FitAsymmetry(const FPPPhiHistos &h, double pp_mean, double Pbeam);         // TODO
  TargetPol SpinTransport(const FocalPlanePol &fp, double chi);                            // TODO

  // spin precession angle in a dipole: chi = gamma * kappa_p * bend_angle (approximation)
  double SpinPrecessionAngle(double pp, double bend_angle);

  // mu_p GE/GM = -mu_p (Pt/Pl) (E+E')/(2M) tan(theta_e/2)
  Ratio FormFactorRatio(const TargetPol &tp, double ebeam, double eprime, double theta_e);

  //------------------------------------------------------------------
  // (B) event-by-event route
  TargetPol FitMaxLikelihood(TTree *T, double Pbeam);                                       // TODO

}

#endif
