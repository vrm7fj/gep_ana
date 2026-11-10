/*
   GEp elastic ep analysis: event selection + FPP angles -> output tree Tout.
   Structure follows gmn-ana's qelas_ana_data.cpp:

     config (JSON) -> TChain of runs -> global cut (TTreeFormula)
     -> per-event variables -> cut flags -> output tree Tout

   Plotting and the polarization stage are done afterwards from Tout.

   Reads the SBS-replay GEp full-replay tree (checked on GEP5 run 5969, Aug-2026 pass):
     heep.*                 elastic kinematics + ECAL residuals computed by the replay
     sbs.tr.*               SBS (proton) track
     sbs.gemFPP.track.*     FPP tracks; theta/phi/sclose/zclose are relative to the
                            front-tracker track and are filled only when one exists,
                            so the count is Ndata.sbs.gemFPP.track.theta (not .ntrack)
     earm.ecal.*            ECAL best cluster
     scalhel.*              helicity (-1, 0 = undefined, +1)
     IGL1I00OD16_16         IHWP state from EPICS (+1 IN, -1 OUT); per-run values can be
                            given instead with "ihwp_file" (output of scripts/ihwp_runs.C)

   Cut flags in Tout:
     dppCut, ecalCut, dphiCut, fiduCut, elasticCut   elastic selection
     fppNoConeCut   theta_fpp, zclose and sclose windows on the chosen FPP track
     coneCut        cone test: the cone of half-angle theta_fpp (this event's angle)
                    around the incoming track, from zclose, fits inside the FPP
                    active area at fpp_cone_zback (true when use_cone_test is false)
     fppCut         fppNoConeCut && coneCut

   Spin precession (ideal dipole, per event):
     fp_rx, fp_ry, fp_rth, fp_rph   focal-plane track in the optics frame (sbs.tr.r_*), the
                                    coordinates SBS-offline's optics expansion uses
     p_tr           sbs.tr.p
     thetabend      bend angle, computed as in SBSEArm::CalcTargetCoords() from tg_th, tg_ph,
                    fp_rth, fp_rph and "sbs_opticstheta_deg" (= sbs.opticstheta in the
                    SBS-offline DB; 0 when the DB does not set it, as for GEp)
     chi            gamma * kappa_p * thetabend (rad)
     dthbend        dispersive (vertical) bend, atan(x'_fp) - atan(x'_tg), signed (rad)
     dphbend        non-dispersive (horizontal) bend, atan(y'_fp) - atan(y'_tg), signed (rad)
     chi_th, chi_ph gamma * kappa_p * dthbend, gamma * kappa_p * dphbend: precession about the
                    horizontal and vertical axes (two-rotation approximation), signed (rad)
     Check: p_tr*thetabend reproduces the optics (p*thetabend) term event by event.

   Usage (ROOT):
     root -l -b -q 'scripts/gep_elastic_ana.cpp("config/gep5/conf_gep5_5969.json","out/gep5")'
*/
#include <map>
#include <chrono>
#include <vector>
#include <iostream>

#include "TSystem.h"

#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"

#include "../include/gep-ana.h"

int gep_elastic_ana(const char *configfilename, std::string filebase = "../gep_out", int onlyrun = -1)
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
  if (onlyrun > 0 ) runnums = {onlyrun};
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

  // cone test (event by event: cone half-angle = that event's theta_fpp)
  bool use_cone = jmgr->GetValueFromKey<bool>("use_cone_test", false);
  double cone_zback = jmgr->GetValueFromKey<double>("fpp_cone_zback", 0.);
  std::vector<double> fpp_area; jmgr->GetVectorFromKey<double>("fpp_active_area", fpp_area);   // {xmin,xmax,ymin,ymax}
  if (use_cone && fpp_area.size() != 4) { std::cerr << "fpp_active_area needs 4 values **!**" << std::endl; return 1; }

  // helicity: overall sign (Wien flip / decoder convention), optionally overridden per run
  int hel_sign = jmgr->GetValueFromKey<int>("helicity_sign", 1);
  std::map<std::string, int> hel_sign_by_run =
    jmgr->GetValueFromKey<std::map<std::string, int>>("helicity_sign_by_run", {});

  // IHWP state per run, from scripts/ihwp_runs.C (optional; EPICS per event is the fallback)
  std::map<std::string, int> ihwp_by_run;
  std::string ihwp_file = jmgr->GetValueFromKey<std::string>("ihwp_file", "");
  if (!ihwp_file.empty()) {
    // a relative path is taken relative to the config file's folder, not the working directory
    if (ihwp_file[0] != '/') ihwp_file = std::string(gSystem->GetDirName(configfilename).Data()) + "/" + ihwp_file;
    JSONManager jih(ihwp_file);
    ihwp_by_run = jih.GetValueFromKey<std::map<std::string, int>>("ihwp_by_run", {});
    std::cout << "IHWP states for " << ihwp_by_run.size() << " runs read from " << ihwp_file << std::endl;
  }

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

  // SBS (proton) track: target and focal-plane variables, needed for spin transport and the cone test
  double ntrack = 0, vz[maxNtr], tgy[maxNtr], tgth[maxNtr], tgph[maxNtr], tgdp[maxNtr];
  double fpx[maxNtr], fpy[maxNtr], fpth[maxNtr], fpph[maxNtr];
  // optics-frame focal-plane coordinates (what SBS-offline's optics uses) and track momentum:
  // needed for the bend angle / spin precession. Zero-initialised in case a leaf is missing.
  double fprx[maxNtr] = {0}, fpry[maxNtr] = {0}, fprth[maxNtr] = {0}, fprph[maxNtr] = {0}, ptr[maxNtr] = {0};
  setrootvar::setbranch(C, "sbs.tr",
    {"n","vz","tg_y","tg_th","tg_ph","tg_dp","x","y","th","ph","r_x","r_y","r_th","r_ph","p"},
    {&ntrack,&vz,&tgy,&tgth,&tgph,&tgdp,&fpx,&fpy,&fpth,&fpph,&fprx,&fpry,&fprth,&fprph,&ptr});

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

  // IHWP from EPICS (+1 IN, -1 OUT, 0/1e38 before the first EPICS read)
  double ihwp_raw = 0;
  setrootvar::setbranch(C, "IGL1I00OD16_16", "", &ihwp_raw);

  // optics angle of the SBS-offline optics frame (sbs.opticstheta; not set for GEp -> 0)
  double opticsangle = jmgr->GetValueFromKey<double>("sbs_opticstheta_deg", 0.) * constant::deg2rad;

  // ================================================================
  // 3. Output tree
  // ================================================================
  TString outFile = Form("%s_%s.root", filebase.c_str(), gepconf.GetLabel().c_str());
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  TTree *Tout = new TTree("Tout", "GEp elastic events");
  bool dppCut, ecalCut, dphiCut, fiduCut, elasticCut, fppNoConeCut, coneCut, fppCut;
  Tout->Branch("dppCut", &dppCut, "dppCut/O");
  Tout->Branch("ecalCut", &ecalCut, "ecalCut/O");
  Tout->Branch("dphiCut", &dphiCut, "dphiCut/O");
  Tout->Branch("fiduCut", &fiduCut, "fiduCut/O");
  Tout->Branch("elasticCut", &elasticCut, "elasticCut/O");
  Tout->Branch("fppNoConeCut", &fppNoConeCut, "fppNoConeCut/O");
  Tout->Branch("coneCut", &coneCut, "coneCut/O");
  Tout->Branch("fppCut", &fppCut, "fppCut/O");
  int T_run, T_hel, T_hel_raw, T_ihwp, T_nfpp;
  double T_Q2, T_eps, T_pp, T_pp_el, T_dpp, T_thp, T_php, T_eprime, T_the, T_dphi, T_dt;
  double T_vz, T_tgy, T_tgth, T_tgph, T_tgdp, T_fpx, T_fpy, T_fpth, T_fpph;
  double T_eECAL, T_xECAL, T_yECAL, T_xECAL_exp, T_yECAL_exp, T_dx, T_dy;
  double T_thfpp, T_phfpp, T_sclose, T_zclose, T_chi2fpp;
  double T_fprx, T_fpry, T_fprth, T_fprph, T_ptr, T_thb, T_chi;
  double T_dthb, T_dphb, T_chith, T_chiph;
  Tout->Branch("run", &T_run, "run/I");
  Tout->Branch("helicity", &T_hel, "helicity/I");             // corrected: raw * sign * IHWP
  Tout->Branch("helicity_raw", &T_hel_raw, "helicity_raw/I");  // as decoded (scalhel.hel)
  Tout->Branch("IHWP", &T_ihwp, "IHWP/I");                     // +1 IN, -1 OUT, 0 unknown
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
  Tout->Branch("fp_rx", &T_fprx, "fp_rx/D");           // optics frame (sbs.tr.r_*):
  Tout->Branch("fp_ry", &T_fpry, "fp_ry/D");           // used by SBS-offline's optics
  Tout->Branch("fp_rth", &T_fprth, "fp_rth/D");
  Tout->Branch("fp_rph", &T_fprph, "fp_rph/D");
  Tout->Branch("p_tr", &T_ptr, "p_tr/D");              // sbs.tr.p
  Tout->Branch("thetabend", &T_thb, "thetabend/D");    // SBS bend angle (rad)
  Tout->Branch("chi", &T_chi, "chi/D");                // gamma*kappa_p*thetabend (rad)
  Tout->Branch("dthbend", &T_dthb, "dthbend/D");       // dispersive bend, signed (rad)
  Tout->Branch("dphbend", &T_dphb, "dphbend/D");       // non-dispersive bend, signed (rad)
  Tout->Branch("chi_th", &T_chith, "chi_th/D");        // gamma*kappa_p*dthbend (rad)
  Tout->Branch("chi_ph", &T_chiph, "chi_ph/D");        // gamma*kappa_p*dphbend (rad)
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
  long nevent = 0, npass = 0, nelas = 0, nconeok = 0, nfppok = 0;
  std::cout << "Processing " << C->GetNtrees() << " files..." << std::endl;
  int treenum = 0, currenttreenum = 0;

  // progress: one line every print_every events (plain lines, so it works in batch / log files)
  long print_every = jmgr->GetValueFromKey<long>("print_every", 200000);
  auto t_start = std::chrono::steady_clock::now();
  auto print_progress = [&](long n) {
    double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
    TString fname = C->GetCurrentFile() ? gSystem->BaseName(C->GetCurrentFile()->GetName()) : "";
    printf("%ld events | file %d/%d %s | passed %ld, elastic %ld, cone ok %ld, FPP ok %ld | %.0f s, %.0f ev/s\n",
           n, C->GetTreeNumber() + 1, C->GetNtrees(), fname.Data(),
           npass, nelas, nconeok, nfppok, el, el > 0 ? n / el : 0.);
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
    fppNoConeCut = false; coneCut = false; fppCut = false;
    int ibest = -1; double best = 1e9;
    for (int it = 0; it < nfpp && it < maxNtr; it++) {
      if (!cut::InWindow(fpp_th[it], thfpp_cut) || !cut::InWindow(fpp_zc[it], zclose_cut)) continue;
      if (fpp_sc[it] < best) { best = fpp_sc[it]; ibest = it; }
    }
    if (ibest < 0 && nfpp > 0) ibest = 0;            // keep something in the tree
    if (ibest >= 0) {
      thfpp = fpp_th[ibest]; phfpp = fpp_ph[ibest]; sclose = fpp_sc[ibest];
      zclose = fpp_zc[ibest]; chi2fpp = fpp_chi2[ibest];

      fppNoConeCut = cut::InWindow(thfpp, thfpp_cut) && cut::InWindow(zclose, zclose_cut) && sclose < sclose_max;

      // ---- cone test: incoming proton track at the focal plane (z0 = 0, front-tracker frame)
      fpp::Track in;
      in.x  = fpx[0];    // sbs.tr.x
      in.y  = fpy[0];    // sbs.tr.y
      in.xp = fpth[0];   // sbs.tr.th (dx/dz)
      in.yp = fpph[0];   // sbs.tr.ph (dy/dz)
      in.z0 = 0.;
      in.valid = true;

      coneCut = !use_cone || (cut::InWindow(zclose, zclose_cut) &&
                              fpp::ConeTest(in, zclose, cone_zback, thfpp, fpp_area));

      fppCut = fppNoConeCut && coneCut;
    }

    // ---- helicity: raw decoded value (0 if undefined or decoder error),
    //      times the overall sign, times -1 when the IHWP is IN; 0 if the IHWP state is unknown
    int run = int(runnum_d);
    int sign = hel_sign;
    auto its = hel_sign_by_run.find(std::to_string(run));
    if (its != hel_sign_by_run.end()) sign = its->second;

    int hel_raw = (hel_err != 0) ? 0 : (helicity_raw > 0 ? 1 : (helicity_raw < 0 ? -1 : 0));

    int ihwp = 0;
    auto iti = ihwp_by_run.find(std::to_string(run));
    if (iti != ihwp_by_run.end() && iti->second != 0) ihwp = iti->second;   // per-run value
    else ihwp = (ihwp_raw > 0.5) ? 1 : (ihwp_raw < -0.5 ? -1 : 0);          // per-event EPICS

    int hel = (ihwp == 0) ? 0 : hel_raw * sign * (ihwp == 1 ? -1 : 1);

    // ---- counters (printout only)
    if (elasticCut) {
      nelas++;
      if (ibest >= 0 && coneCut) nconeok++;
      if (fppCut) nfppok++;
    }

    // ---- output tree
    T_run = run; T_hel = hel; T_hel_raw = hel_raw; T_ihwp = ihwp; T_nfpp = nfpp;
    T_Q2 = hv_Q2; T_eps = hv_eps;
    T_pp = hv_pp; T_pp_el = hv_ppel; T_dpp = dpp; T_thp = hv_pth; T_php = hv_pph;
    T_eprime = hv_eprime; T_the = hv_eth; T_dphi = dphi; T_dt = hv_dt;
    T_vz = vz[0]; T_tgy = tgy[0]; T_tgth = tgth[0]; T_tgph = tgph[0]; T_tgdp = tgdp[0];
    T_fpx = fpx[0]; T_fpy = fpy[0]; T_fpth = fpth[0]; T_fpph = fpph[0];
    T_fprx = fprx[0]; T_fpry = fpry[0]; T_fprth = fprth[0]; T_fprph = fprph[0]; T_ptr = ptr[0];
    T_thb = fpp::BendAngle(tgth[0], tgph[0], fprth[0], fprph[0], opticsangle);
    T_chi = fpp::PrecessionAngle(ptr[0], T_thb);
    fpp::BendComponents(tgth[0], tgph[0], fprth[0], fprph[0], opticsangle, T_dthb, T_dphb);
    double gk = (ptr[0] > 0.) ? std::sqrt(ptr[0]*ptr[0] + constant::Mp*constant::Mp) / constant::Mp * constant::kappap : 0.;
    T_chith = gk * T_dthb;                               // 0 if p is not valid
    T_chiph = gk * T_dphb;
    T_eECAL = eECAL; T_xECAL = xECAL; T_yECAL = yECAL;
    T_xECAL_exp = xexp; T_yECAL_exp = yexp; T_dx = dx; T_dy = dy;
    T_thfpp = thfpp; T_phfpp = phfpp; T_sclose = sclose; T_zclose = zclose; T_chi2fpp = chi2fpp;
    Tout->Fill();
  }

  long nread = nevent - 1;
  print_progress(nread);
  std::cout << nread << " events read, " << npass << " passed the global cut, "
            << nelas << " elastic, " << nconeok << " elastic passing the cone test, "
            << nfppok << " elastic with a good FPP scattering." << std::endl;

  // ================================================================
  // 5. Output
  // ================================================================
  fout->cd();
  Tout->Write();
  fout->Close();
  std::cout << "Output written to " << outFile << std::endl;

  sw->Stop();
  std::cout << "CPU time " << sw->CpuTime() << " s, real time " << sw->RealTime() << " s" << std::endl;
  return 0;
}
