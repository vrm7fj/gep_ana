#include <cmath>
#include <iostream>
#include "TMath.h"
#include "../include/Constants.h"
#include "../include/Polarization.h"

namespace pol {

  static void todo(const char *fn) {
    std::cout << "[pol::" << fn << "] TODO: not implemented yet, returning an invalid result." << std::endl;
  }

  /* ##################################################################
     ##                 Bookkeeping (implemented)                    ##
     ################################################################## */
  void FPPPhiHistos::Book(const std::string &name, int nphi, const std::vector<double> &edges) {
    fEdges = edges;
    for (size_t i = 0; i + 1 < edges.size(); i++) {
      std::string tag = Form("%s_th%zu", name.c_str(), i);
      std::string ttl = Form("%.1f < #theta_{fpp} < %.1f deg;#phi_{fpp} (rad);counts",
                             edges[i] * TMath::RadToDeg(), edges[i + 1] * TMath::RadToDeg());
      fHp.push_back(new TH1D((tag + "_hp").c_str(), ("h=+1, " + ttl).c_str(), nphi, -TMath::Pi(), TMath::Pi()));
      fHm.push_back(new TH1D((tag + "_hm").c_str(), ("h=-1, " + ttl).c_str(), nphi, -TMath::Pi(), TMath::Pi()));
      fHp.back()->Sumw2(); fHm.back()->Sumw2();
    }
  }
  //--------------------------------------------
  void FPPPhiHistos::Fill(double theta, double phi, int hel) {
    if (hel == 0) return;
    for (size_t i = 0; i + 1 < fEdges.size(); i++) {
      if (theta >= fEdges[i] && theta < fEdges[i + 1]) {
        (hel > 0 ? fHp[i] : fHm[i])->Fill(phi);
        return;
      }
    }
  }
  //--------------------------------------------
  void FPPPhiHistos::Write() const {
    for (auto *h : fHp) h->Write();
    for (auto *h : fHm) h->Write();
  }

  /* ##################################################################
     ##             (A) Binned route  -- TO IMPLEMENT                ##
     ################################################################## */
  TH1D *HelicityAsymmetry(const TH1D *hp, const TH1D *hm, const std::string &name) {
    // TODO: build A(phi) = (N+ - N-)/(N+ + N-) bin by bin (with charge/livetime
    //       normalization of each helicity state if needed) and return it.
    todo("HelicityAsymmetry");
    return nullptr;
  }
  //--------------------------------------------
  double AnalyzingPower(double theta_fpp, double pp) {
    // TODO: analyzing power parametrization A_y(theta_fpp, p) for the analyzer
    //       (CH2/CH or Fe), e.g. from previous GEp measurements or calibration.
    todo("AnalyzingPower");
    return 0.;
  }
  //--------------------------------------------
  FocalPlanePol FitAsymmetry(const FPPPhiHistos &h, double pp_mean, double Pbeam) {
    // TODO: for each theta_fpp bin, fit A(phi) with
    //         A(phi) = Pbeam * A_y(theta) * (Px sin(phi) - Py cos(phi))   (check signs!)
    //       then combine bins (or do one global fit) to get Px, Py at the focal plane.
    todo("FitAsymmetry");
    return FocalPlanePol{};
  }
  //--------------------------------------------
  TargetPol SpinTransport(const FocalPlanePol &fp, double chi) {
    // TODO: rotate focal-plane components back to the target frame.
    //       Simple dipole picture: P_t ~ P_y^fpp, P_l ~ P_x^fpp / sin(chi)
    //       Full analysis: spin transport matrix (e.g. COSY) per event or per bin.
    todo("SpinTransport");
    return TargetPol{};
  }

  /* ##################################################################
     ##                 Kinematic helpers (implemented)              ##
     ################################################################## */
  double SpinPrecessionAngle(double pp, double bend) {
    double gamma = std::sqrt(pp * pp + constant::Mp * constant::Mp) / constant::Mp;
    return gamma * constant::kappap * bend;
  }
  //--------------------------------------------
  Ratio FormFactorRatio(const TargetPol &tp, double E, double Ep, double the) {
    Ratio r;
    if (!tp.ok || tp.Pl == 0.) return r;
    double k = (E + Ep) / (2. * constant::Mp) * std::tan(the / 2.);
    r.R = -constant::mup * (tp.Pt / tp.Pl) * k;
    // statistical error from Pt and Pl only (ignores their correlation)
    double rel2 = std::pow(tp.Pt_err / tp.Pt, 2) + std::pow(tp.Pl_err / tp.Pl, 2);
    r.R_err = std::fabs(r.R) * std::sqrt(rel2);
    r.ok = true;
    return r;
  }

  /* ##################################################################
     ##          (B) Event-by-event route -- TO IMPLEMENT            ##
     ################################################################## */
  TargetPol FitMaxLikelihood(TTree *T, double Pbeam) {
    // TODO: maximize  L = prod_i [1 + h_i Pbeam A_y(theta_i) (c_i^t P_t + c_i^l P_l)]
    //       where c_i^t, c_i^l come from the per-event spin transport and phi_i.
    //       Use only events with elasticCut && fppCut && helicity != 0
    //       (all stored as branches in Tout).
    todo("FitMaxLikelihood");
    return TargetPol{};
  }

}
