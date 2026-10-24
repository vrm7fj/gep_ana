/*
   GEp elastic ep analysis: event selection + FPP angles, feeding the
   polarization-observable stage. Structure follows gmn-ana's qelas_ana_data.cpp:

     config (JSON) -> TChain of runs -> global cut (TTreeFormula)
     -> per-event variables -> cut flags -> histograms + output tree Tout
     -> polarization stage (src/Polarization.cpp)

   Reads the SBS-replay GEp full-replay tree (checked on GEP5 run 5969, Aug-2026 pass):
     heep.*                 elastic kinematics + ECAL residuals computed by the replay
     sbs.tr.*               SBS (proton) track
     sbs.gemFPP.track.*     FPP tracks; theta/phi/sclose/zclose are relative to the
                            front-tracker track and are filled only when one exists,
                            so the count is Ndata.sbs.gemFPP.track.theta (not .ntrack)
     earm.ecal.*            ECAL best cluster
     scalhel.*              helicity (-1, 0 = undefined, +1)

   Usage (ROOT):
     root -l -b -q 'scripts/gep_elastic_ana.cpp("config/gep5/conf_gep5_5969.json","out/gep5")'
*/
#include <map>
#include <chrono>
#include <vector>
#include <iostream>

#include "TSystem.h"

#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"

#include "../include/gep-ana.h"

int gep_elastic_ana(const char *configfilename, std::string filebase = "out/gep_elastic_ana")
{
  gErrorIgnoreLevel = kError;
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // ================================================================
  // 1. Config and input
  // ================================================================
  JSONManager *jmgr = new JSONManager(configfilename);
  GEpConfig gepconf(*jmgr);
  gepconf.Print();

  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::string file_pattern = jmgr->GetValueFromKey<std::string>("file_pattern", "*%d*.root");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums", runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana", -1);
  if (nruns < 1 || nruns > int(runnums.size())) nruns = runnums.size();
  int max_files_per_run = jmgr->GetValueFromKey<int>("max_files_per_run", -1);
  TChain *C = new TChain("T");
  for (int i = 0; i < nruns; i++) {
    std::string rfname = rootfile_dir + "/" + Form(file_pattern.c_str(), runnums[i]);
    int nf = 0;
    if (max_files_per_run > 0) {
      TChain tmp("T"); tmp.Add(rfname.c_str());
      TObjArray *fl = tmp.GetListOfFiles();
      for (int k = 0; k < fl->GetEntries() && k < max_files_per_run; k++) nf += C->Add(fl->At(k)->GetTitle());
    } else {
      nf = C->Add(rfname.c_str());
    }
    std::cout << "Run " << runnums[i] << ": " << nf << " files" << std::endl;
  }
  std::cout << "Chained " << C->GetNtrees() << " files from " << nruns << " runs." << std::endl;
  if (C->GetNtrees() == 0) { std::cerr << "Nothing to analyze, check rootfile_dir / file_pattern **!**" << std::endl; return 1; }

  // global cut (string in config, evaluated per event)
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", TCut(gcut.c_str()), C);

  // elastic selection parameters
  std::vector<double> dpp_cut;  jmgr->GetVectorFromKey<double>("dpp_cut", dpp_cut);       // {lo, hi}
  std::vector<double> dx_ecal;  jmgr->GetVectorFromKey<double>("dx_ecal", dx_ecal);       // {mean, sigma}
  std::vector<double> dy_ecal;  jmgr->GetVectorFromKey<double>("dy_ecal", dy_ecal);       // {mean, sigma}
  double nsig_dx = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_ecal");
  double nsig_dy = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_ecal");
  std::vector<double> dphi_cut; jmgr->GetVectorFromKey<double>("dphi_cut", dphi_cut);     // {lo, hi} rad
  bool use_dphi = jmgr->GetValueFromKey<bool>("use_dphi_in_elastic_cut", true);
  std::vector<double> ecal_area; jmgr->GetVectorFromKey<double>("ecal_active_area", ecal_area);     // {xmin,xmax,ymin,ymax}
  std::vector<double> ecal_margin; jmgr->GetVectorFromKey<double>("ecal_fiducial_margin", ecal_margin); // {mx, my}
  std::vector<double> ecal_fidu = cut::Shrink(ecal_area, ecal_margin[0], ecal_margin[1]);

  // FPP selection parameters
  std::vector<double> thfpp_cut;  jmgr->GetVectorFromKey<double>("theta_fpp_cut_deg", thfpp_cut);
  for (auto &v : thfpp_cut) v *= constant::deg2rad;
  std::vector<double> zclose_cut; jmgr->GetVectorFromKey<double>("zclose_cut", zclose_cut);
  double sclose_max = jmgr->GetValueFromKey<double>("sclose_max");

  // helicity: global sign, optionally overridden per run (IHWP state)
  int hel_sign = jmgr->GetValueFromKey<int>("helicity_sign", 1);
  std::map<std::string, int> hel_sign_by_run =
    jmgr->GetValueFromKey<std::map<std::string, int>>("helicity_sign_by_run", {});
  double Pbeam = jmgr->GetValueFromKey<double>("beam_polarization", 0.);

  // ================================================================
  // 2. Branches
  // ================================================================
  const int maxNtr = 1000;
  C->SetBranchStatus("*", 0);

  // elastic kinematics from the replay (1e38 when not valid -> require heep.datavalid in global cut)
  double hv_valid, hv_dpp, hv_dx, hv_dy, hv_dphi, hv_Q2, hv_eps, hv_eprime, hv_eth, hv_ephi;
  double hv_pp, hv_pth, hv_pph, hv_ppel, hv_dt;
  setrootvar::setbranch(C, "heep",
    {"datavalid","dpp","dxECAL","dyECAL","dphi","Q2_pth","eps_pth","eprime_pth","etheta","ephi",
     "pp","ptheta","pphi","pp_pth","dt_ADC"},
    {&hv_valid,&hv_dpp,&hv_dx,&hv_dy,&hv_dphi,&hv_Q2,&hv_eps,&hv_eprime,&hv_eth,&hv_ephi,
     &hv_pp,&hv_pth,&hv_pph,&hv_ppel,&hv_dt});

  // SBS (proton) track: target and focal-plane variables, needed for spin transport
  double ntrack = 0, vz[maxNtr], tgy[maxNtr], tgth[maxNtr], tgph[maxNtr], tgdp[maxNtr];
  double fpx[maxNtr], fpy[maxNtr], fpth[maxNtr], fpph[maxNtr];
  setrootvar::setbranch(C, "sbs.tr",
    {"n","vz","tg_y","tg_th","tg_ph","tg_dp","x","y","th","ph"},
    {&ntrack,&vz,&tgy,&tgth,&tgph,&tgdp,&fpx,&fpy,&fpth,&fpph});

  // FPP tracks (angles relative to the front-tracker track, computed in the replay)
  int nfpp = 0;
  double fpp_th[maxNtr], fpp_ph[maxNtr], fpp_sc[maxNtr], fpp_zc[maxNtr], fpp_chi2[maxNtr];
  setrootvar::setbranch(C, "sbs.gemFPP.track",
    {"theta","theta","phi","sclose","zclose","chi2ndf"},
    {&nfpp,&fpp_th,&fpp_ph,&fpp_sc,&fpp_zc,&fpp_chi2}, std::vector<int>{0});

  // ECAL best cluster
  double eECAL = 0, xECAL = 0, yECAL = 0;
  setrootvar::setbranch(C, "earm.ecal", {"e","x","y"}, {&eECAL,&xECAL,&yECAL});

  // helicity and run number
  double helicity_raw = 0, hel_err = 0, runnum_d = 0;
  setrootvar::setbranch(C, "scalhel", {"hel","errcode"}, {&helicity_raw,&hel_err});
  setrootvar::setbranch(C, "g.runnum", "", &runnum_d);

  // leaves used only by the global cut
  for (auto &b : jmgr->GetValueFromKey<std::vector<std::string>>("extra_branches_for_global_cut", {}))
    C->SetBranchStatus(b.c_str(), 1);

  // ================================================================
  // 3. Output: histograms and tree
  // ================================================================
  TString outFile = Form("%s_%s.root", filebase.c_str(), gepconf.GetLabel().c_str());
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  std::vector<double> l_dpp, l_dx, l_dy, l_dphi, l_Q2;
  jmgr->GetVectorFromKey<double>("h_dpp_lims", l_dpp);
  jmgr->GetVectorFromKey<double>("h_dx_ecal_lims", l_dx);
  jmgr->GetVectorFromKey<double>("h_dy_ecal_lims", l_dy);
  jmgr->GetVectorFromKey<double>("h_dphi_lims", l_dphi);
  jmgr->GetVectorFromKey<double>("h_Q2_lims", l_Q2);

  TH1F *h_dpp      = util::TH1Flims("h_dpp",      ";p/p_{el}(#theta_{p}) - 1;", l_dpp);
  TH1F *h_dpp_cut  = util::TH1Flims("h_dpp_cut",  "ECAL spot cut;p/p_{el}(#theta_{p}) - 1;", l_dpp);
  TH1F *h_dpp_acut = util::TH1Flims("h_dpp_acut", "ECAL spot anti-cut;p/p_{el}(#theta_{p}) - 1;", l_dpp);
  TH1F *h_dx       = util::TH1Flims("h_dx_ecal",  ";x_{ECAL} - x_{exp} (m);", l_dx);
  TH1F *h_dy       = util::TH1Flims("h_dy_ecal",  ";y_{ECAL} - y_{exp} (m);", l_dy);
  TH2F *h2_dxdy    = util::TH2Flims("h2_dxdy_ecal", ";y_{ECAL} - y_{exp} (m);x_{ECAL} - x_{exp} (m)", l_dy, l_dx);
  TH1F *h_dphi     = util::TH1Flims("h_dphi",     ";#phi_{e} - #phi_{p} - #pi (rad);", l_dphi);
  TH1F *h_Q2       = util::TH1Flims("h_Q2",       ";Q^{2} (GeV^{2});", l_Q2);
  TH2F *h2_xyECAL  = new TH2F("h2_xyECAL_exp", "expected e' at ECAL;y_{exp} (m);x_{exp} (m)",
                              100, ecal_area[2] - 0.2, ecal_area[3] + 0.2, 100, ecal_area[0] - 0.2, ecal_area[1] + 0.2);
  TH1F *h_thfpp    = new TH1F("h_theta_fpp", ";#theta_{fpp} (deg);", 200, 0, 25);
  TH1F *h_phfpp    = new TH1F("h_phi_fpp",   ";#phi_{fpp} (rad);", 72, -constant::pi, constant::pi);
  TH1F *h_zclose   = new TH1F("h_zclose",    ";z_{close} (m);", 200, -1, 5);
  TH1F *h_sclose   = new TH1F("h_sclose",    ";s_{close} (m);", 200, 0, 0.1);
  TH1F *h_hel      = new TH1F("h_helicity",  ";helicity;", 3, -1.5, 1.5);

  // phi_fpp per helicity and theta_fpp bin -> input to the polarization stage
  std::vector<double> th_bins; jmgr->GetVectorFromKey<double>("theta_fpp_bins_deg", th_bins);
  for (auto &v : th_bins) v *= constant::deg2rad;
  pol::FPPPhiHistos fppPhi;
  fppPhi.Book("h_phi_fpp", jmgr->GetValueFromKey<int>("n_phi_bins", 24), th_bins);

  TTree *Tout = new TTree("Tout", "GEp elastic events");
  bool dppCut, ecalCut, dphiCut, fiduCut, elasticCut, fppCut;
  Tout->Branch("dppCut", &dppCut, "dppCut/O");
  Tout->Branch("ecalCut", &ecalCut, "ecalCut/O");
  Tout->Branch("dphiCut", &dphiCut, "dphiCut/O");
  Tout->Branch("fiduCut", &fiduCut, "fiduCut/O");
  Tout->Branch("elasticCut", &elasticCut, "elasticCut/O");
  Tout->Branch("fppCut", &fppCut, "fppCut/O");
  int T_run, T_hel, T_nfpp;
  double T_Q2, T_eps, T_pp, T_pp_el, T_dpp, T_thp, T_php, T_eprime, T_the, T_dphi, T_dt;
  double T_vz, T_tgy, T_tgth, T_tgph, T_tgdp, T_fpx, T_fpy, T_fpth, T_fpph;
  double T_eECAL, T_xECAL, T_yECAL, T_xECAL_exp, T_yECAL_exp, T_dx, T_dy;
  double T_thfpp, T_phfpp, T_sclose, T_zclose, T_chi2fpp;
  Tout->Branch("run", &T_run, "run/I");
  Tout->Branch("helicity", &T_hel, "helicity/I");
  Tout->Branch("Q2", &T_Q2, "Q2/D");
  Tout->Branch("epsilon", &T_eps, "epsilon/D");
  Tout->Branch("pp", &T_pp, "pp/D");
  Tout->Branch("pp_el", &T_pp_el, "pp_el/D");
  Tout->Branch("dpp", &T_dpp, "dpp/D");
  Tout->Branch("theta_p", &T_thp, "theta_p/D");
  Tout->Branch("phi_p", &T_php, "phi_p/D");
  Tout->Branch("eprime", &T_eprime, "eprime/D");
  Tout->Branch("theta_e", &T_the, "theta_e/D");
  Tout->Branch("dphi", &T_dphi, "dphi/D");
  Tout->Branch("dt_ADC", &T_dt, "dt_ADC/D");
  Tout->Branch("vz", &T_vz, "vz/D");
  Tout->Branch("tg_y", &T_tgy, "tg_y/D");     // target variables: needed for
  Tout->Branch("tg_th", &T_tgth, "tg_th/D");  // event-by-event spin transport
  Tout->Branch("tg_ph", &T_tgph, "tg_ph/D");
  Tout->Branch("tg_dp", &T_tgdp, "tg_dp/D");
  Tout->Branch("fp_x", &T_fpx, "fp_x/D");
  Tout->Branch("fp_y", &T_fpy, "fp_y/D");
  Tout->Branch("fp_th", &T_fpth, "fp_th/D");
  Tout->Branch("fp_ph", &T_fpph, "fp_ph/D");
  Tout->Branch("eECAL", &T_eECAL, "eECAL/D");
  Tout->Branch("xECAL", &T_xECAL, "xECAL/D");
  Tout->Branch("yECAL", &T_yECAL, "yECAL/D");
  Tout->Branch("xECAL_exp", &T_xECAL_exp, "xECAL_exp/D");
  Tout->Branch("yECAL_exp", &T_yECAL_exp, "yECAL_exp/D");
  Tout->Branch("dx", &T_dx, "dx/D");
  Tout->Branch("dy", &T_dy, "dy/D");
  Tout->Branch("nfpp", &T_nfpp, "nfpp/I");
  Tout->Branch("theta_fpp", &T_thfpp, "theta_fpp/D");
  Tout->Branch("phi_fpp", &T_phfpp, "phi_fpp/D");
  Tout->Branch("sclose", &T_sclose, "sclose/D");
  Tout->Branch("zclose", &T_zclose, "zclose/D");
  Tout->Branch("chi2ndf_fpp", &T_chi2fpp, "chi2ndf_fpp/D");

  // ================================================================
  // 4. Event loop
  // ================================================================
  std::cout << "Counting entries in " << C->GetNtrees() << " files..." << std::endl;
  long nevent = 0, nevents = C->GetEntries(), npass = 0, nelas = 0, nfppok = 0;
  std::cout << nevents << " events to process." << std::endl;
  int treenum = 0, currenttreenum = 0;

  // progress: one line every print_every events (plain lines, so it works in batch / log files)
  long print_every = jmgr->GetValueFromKey<long>("print_every", 200000);
  auto t_start = std::chrono::steady_clock::now();
  auto print_progress = [&](long n) {
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    double frac = nevents > 0 ? double(n) / nevents : 0.;
    double eta = frac > 0 ? el * (1. - frac) / frac : 0.;
    TString fname = C->GetCurrentFile() ? gSystem->BaseName(C->GetCurrentFile()->GetName()) : "";
    printf("[%5.1f%%] %ld/%ld events | file %d/%d %s | passed %ld, elastic %ld, FPP ok %ld | %.0f s, ETA %.0f s\n",
           100. * frac, n, nevents, C->GetTreeNumber() + 1, C->GetNtrees(), fname.Data(),
           npass, nelas, nfppok, el, eta);
    fflush(stdout);
  };

  while (C->GetEntry(nevent++)) {
    if (nevent % print_every == 0) print_progress(nevent);

    // ---- global cut (refresh formula leaves when the chain moves to a new file)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }
    if (GlobalCut->EvalInstance(0) == 0) continue;
    if (hv_valid != 1 || ntrack < 1) continue;     // heep.* is 1e38 when not valid
    npass++;

    // ---- elastic variables (from the replay's heep module)
    double dpp = hv_dpp, dx = hv_dx, dy = hv_dy, dphi = hv_dphi;
    double xexp = xECAL - dx, yexp = yECAL - dy;   // expected e' position at ECAL

    dppCut  = cut::InWindow(dpp, dpp_cut);
    ecalCut = cut::InEllipse(dx, dy, dx_ecal, dy_ecal, nsig_dx, nsig_dy);
    dphiCut = cut::InWindow(dphi, dphi_cut);
    fiduCut = cut::InArea(xexp, yexp, ecal_fidu);    // expected e' well inside ECAL
    elasticCut = dppCut && ecalCut && fiduCut && (!use_dphi || dphiCut);

    // ---- FPP: choose the track with the smallest sclose inside the theta/zclose windows
    double thfpp = -1, phfpp = -10, sclose = -1, zclose = -1000, chi2fpp = -1;
    fppCut = false;
    int ibest = -1; double best = 1e9;
    for (int it = 0; it < nfpp && it < maxNtr; it++) {
      if (!cut::InWindow(fpp_th[it], thfpp_cut) || !cut::InWindow(fpp_zc[it], zclose_cut)) continue;
      if (fpp_sc[it] < best) { best = fpp_sc[it]; ibest = it; }
    }
    if (ibest < 0 && nfpp > 0) ibest = 0;            // keep something for the plots
    if (ibest >= 0) {
      thfpp = fpp_th[ibest]; phfpp = fpp_ph[ibest]; sclose = fpp_sc[ibest];
      zclose = fpp_zc[ibest]; chi2fpp = fpp_chi2[ibest];
      fppCut = cut::InWindow(thfpp, thfpp_cut) && cut::InWindow(zclose, zclose_cut) && sclose < sclose_max;
    }

    // ---- helicity (0 if undefined or decoder error; sign corrected per run for IHWP)
    int run = int(runnum_d);
    int sign = hel_sign;
    auto its = hel_sign_by_run.find(std::to_string(run));
    if (its != hel_sign_by_run.end()) sign = its->second;
    int hel = (hel_err != 0) ? 0 : (helicity_raw > 0 ? sign : (helicity_raw < 0 ? -sign : 0));

    // ---- fill histograms
    h_dpp->Fill(dpp);
    if (fiduCut) (ecalCut ? h_dpp_cut : h_dpp_acut)->Fill(dpp);
    if (dppCut && fiduCut) { h_dx->Fill(dx); h_dy->Fill(dy); h2_dxdy->Fill(dy, dx); }
    if (dppCut && ecalCut) h_dphi->Fill(dphi);
    h2_xyECAL->Fill(yexp, xexp);
    if (elasticCut) {
      nelas++;
      h_Q2->Fill(hv_Q2);
      h_hel->Fill(hel);
      if (thfpp >= 0) {
        h_thfpp->Fill(thfpp / constant::deg2rad);
        h_zclose->Fill(zclose);
        h_sclose->Fill(sclose);
      }
      if (fppCut) {
        nfppok++;
        h_phfpp->Fill(phfpp);
        fppPhi.Fill(thfpp, phfpp, hel);
      }
    }

    // ---- output tree
    T_run = run; T_hel = hel; T_nfpp = nfpp;
    T_Q2 = hv_Q2; T_eps = hv_eps;
    T_pp = hv_pp; T_pp_el = hv_ppel; T_dpp = dpp; T_thp = hv_pth; T_php = hv_pph;
    T_eprime = hv_eprime; T_the = hv_eth; T_dphi = dphi; T_dt = hv_dt;
    T_vz = vz[0]; T_tgy = tgy[0]; T_tgth = tgth[0]; T_tgph = tgph[0]; T_tgdp = tgdp[0];
    T_fpx = fpx[0]; T_fpy = fpy[0]; T_fpth = fpth[0]; T_fpph = fpph[0];
    T_eECAL = eECAL; T_xECAL = xECAL; T_yECAL = yECAL;
    T_xECAL_exp = xexp; T_yECAL_exp = yexp; T_dx = dx; T_dy = dy;
    T_thfpp = thfpp; T_phfpp = phfpp; T_sclose = sclose; T_zclose = zclose; T_chi2fpp = chi2fpp;
    Tout->Fill();
  }
  print_progress(nevents);
  std::cout << nevents << " events read, " << npass << " passed the global cut, "
            << nelas << " elastic, " << nfppok << " elastic with a good FPP scattering." << std::endl;

  // ================================================================
  // 5. Polarization stage (implement in src/Polarization.cpp)
  // ================================================================
  double pp_mean = 0, eprime_mean = 0, the_mean = 0, Q2_mean = h_Q2->GetMean();
  if (Tout->GetEntries("elasticCut") > 0) {
    Tout->Draw("pp:eprime:theta_e", "elasticCut", "goff");
    int n = Tout->GetSelectedRows();
    for (int i = 0; i < n; i++) {
      pp_mean += Tout->GetV1()[i]; eprime_mean += Tout->GetV2()[i]; the_mean += Tout->GetV3()[i];
    }
    pp_mean /= n; eprime_mean /= n; the_mean /= n;
  }
  std::cout << "Elastic sample: <Q2> = " << Q2_mean << " GeV2, <p_p> = " << pp_mean
            << " GeV, <theta_e> = " << the_mean / constant::deg2rad << " deg" << std::endl;

  // (A) binned route
  for (int i = 0; i < fppPhi.NThetaBins(); i++) {
    TH1D *hA = pol::HelicityAsymmetry(fppPhi.Hplus(i), fppPhi.Hminus(i), Form("h_asym_th%d", i));
    if (hA) hA->Write();
  }
  pol::FocalPlanePol fpPol = pol::FitAsymmetry(fppPhi, pp_mean, Pbeam);
  double bend = jmgr->GetValueFromKey<double>("sbs_bend_angle_deg", 0.) * constant::deg2rad;
  double chi = pol::SpinPrecessionAngle(pp_mean, bend);
  pol::TargetPol tgtPol = pol::SpinTransport(fpPol, chi);
  pol::Ratio R = pol::FormFactorRatio(tgtPol, gepconf.GetEbeam(), eprime_mean, the_mean);
  // (B) event-by-event route
  pol::TargetPol tgtPolML = pol::FitMaxLikelihood(Tout, Pbeam);

  std::cout << "chi (dipole approx.) = " << chi / constant::deg2rad << " deg" << std::endl;
  if (R.ok) std::cout << "mu_p GE/GM = " << R.R << " +/- " << R.R_err << std::endl;
  else      std::cout << "mu_p GE/GM: not available yet (polarization stage not implemented)." << std::endl;

  // ================================================================
  // 6. Plots and output
  // ================================================================
  TCanvas *c1 = new TCanvas("c1", "elastic selection", 1200, 1000);
  c1->Divide(2, 2);
  c1->cd(1); h2_dxdy->Draw("colz");
  TEllipse el; el.SetFillStyle(0); el.SetLineColor(2); el.SetLineWidth(2);
  el.DrawEllipse(dy_ecal[0], dx_ecal[0], nsig_dy * dy_ecal[1], nsig_dx * dx_ecal[1], 0, 360, 0);
  c1->cd(2); h_dpp->SetLineColor(1); h_dpp->Draw();
  h_dpp_cut->SetLineColor(2); h_dpp_cut->Draw("same");
  h_dpp_acut->SetLineColor(4); h_dpp_acut->Draw("same");
  c1->cd(3); h_dphi->Draw();
  c1->cd(4); h2_xyECAL->Draw("colz"); util::DrawArea(ecal_area); util::DrawArea(ecal_fidu, 4);
  c1->Write();

  TCanvas *c2 = new TCanvas("c2", "FPP", 1200, 1000);
  c2->Divide(2, 2);
  c2->cd(1); h_thfpp->Draw();
  c2->cd(2); h_phfpp->Draw();
  c2->cd(3); h_zclose->Draw();
  c2->cd(4); h_sclose->Draw();
  c2->Write();

  if (jmgr->GetValueFromKey<bool>("save_pdf", true)) {
    TString pdf = Form("%s_%s.pdf", filebase.c_str(), gepconf.GetLabel().c_str());
    c1->Print(pdf + "(");
    c2->Print(pdf + ")");
    std::cout << "Plots written to " << pdf << std::endl;
  }

  fppPhi.Write();
  fout->Write();
  std::cout << "Output written to " << outFile << std::endl;

  sw->Stop();
  std::cout << "CPU time " << sw->CpuTime() << " s, real time " << sw->RealTime() << " s" << std::endl;
  return 0;
}
