/*
   draw_fpp_phi_tree.C
   Usage 
   root -l -b -q 'draw_fpp_phi_tree.C("../out/runs/gep5_*_GEP5.root","plots/gep5_all_phi.pdf")'
*/
#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

#include "TH1D.h"
#include "TKey.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TTreeFormula.h"
#include "THLimitsFinder.h"
#include "TF1.h"
#include "TPaveText.h"

// histogram copy with the x axis converted from rad to deg
static TH1D *ToDegrees(const TH1D *h) {
  int n = h->GetNbinsX();
  TH1D *hd = new TH1D(Form("%s_deg", h->GetName()), h->GetTitle(), n,
                      h->GetXaxis()->GetXmin() * TMath::RadToDeg(),
                      h->GetXaxis()->GetXmax() * TMath::RadToDeg());
  hd->SetDirectory(nullptr);
  for (int i = 0; i <= n + 1; i++) {
    hd->SetBinContent(i, h->GetBinContent(i));
    hd->SetBinError(i, h->GetBinError(i));
  }
  return hd;
}

// "3.0 < #theta_{fpp} < 5.0 deg" from the histogram title, or "th<N>" if not found
static TString ThetaLabel(const TH1D *h, int ibin) {
  TString t = h->GetTitle();
  int i = t.Index("<");
  if (i > 0) {
    int j = t.Last(',');                       // titles look like "N^{+} + N^{-}, 3.0 < #theta... deg"
    TString s = (j >= 0 && j < i) ? TString(t(j + 1, t.Length())) : t;
    return s.Strip(TString::kBoth);
  }
  return Form("#theta bin %d", ibin);
}

static void DrawPage(TCanvas *c, const std::map<int, TH1D*> &hs, const char *ytitle,
                     const char *pagetitle, bool deg, bool zero_line) {
  c->Clear();
  int n = hs.size();
  if (n == 0) return;
  int ncol = (n > 1) ? 2 : 1;
  int nrow = (n + ncol - 1) / ncol;

  // page layout (fractions of the canvas)
  const double xl = 0.07, xr = 0.98, yb = 0.07, yt = 0.94, gap = 0.04;
  double colw = (xr - xl - (ncol - 1) * gap) / ncol;
  double rowh = (yt - yb) / nrow;

  c->cd();
  TLatex tt; tt.SetNDC(); tt.SetTextFont(43); tt.SetTextSize(20); tt.SetTextAlign(22);
  tt.DrawLatex(0.5, 0.975, pagetitle);
  TLatex ty; ty.SetNDC(); ty.SetTextFont(43); ty.SetTextSize(18); ty.SetTextAngle(90); ty.SetTextAlign(22);
  ty.DrawLatex(0.02, 0.5 * (yb + yt), ytitle);

  int k = 0;
  for (auto &kv : hs) {
    int col = k / nrow, row = k % nrow;            // fill column 1 top->bottom, then column 2
    bool bottom = (row == nrow - 1) || (k == n - 1);
    double x1 = xl + col * (colw + gap), x2 = x1 + colw;
    double y2 = yt - row * rowh, y1 = y2 - rowh;
    // the bottom pad of each column is taller to make room for the x labels
    double ext = bottom ? 0.04 : 0.;
    c->cd();
    TPad *p = new TPad(Form("p_%s_%d", pagetitle, k), "", x1, y1 - ext, x2, y2);
    p->SetTopMargin(0.);
    p->SetBottomMargin(bottom ? ext / (rowh + ext) : 0.);
    p->SetLeftMargin(0.12); p->SetRightMargin(0.02);
    p->SetTicks(1, 1);
    p->Draw(); p->cd();

    TH1D *h = deg ? ToDegrees(kv.second) : (TH1D*)kv.second->Clone(Form("%s_draw", kv.second->GetName()));
    h->SetDirectory(nullptr);
    h->SetStats(0);
    h->SetTitle("");
    h->SetMarkerStyle(20); h->SetMarkerSize(0.6); h->SetLineColor(kBlack);
    TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
    ax->SetLabelFont(43); ax->SetLabelSize(bottom ? 14 : 0); ax->SetTitleFont(43); ax->SetTitleSize(16);
    ax->SetTitle(bottom ? (deg ? "#phi_{fpp} (deg)" : "#phi_{fpp} (rad)") : "");
    ax->SetTitleOffset(3.2); ax->SetTickLength(0.06);
    ay->SetLabelFont(43); ay->SetLabelSize(12); ay->SetNdivisions(504); ay->SetTitle("");
    // y range: data +- errors with some headroom, rounded to "nice" values so the
    // axis labels land on round numbers and do not collide with the panel below
    double lo = h->GetMinimum() - h->GetBinError(h->GetMinimumBin());
    double hi = h->GetMaximum() + h->GetBinError(h->GetMaximumBin());
    double pad_ = 0.12 * (hi - lo > 0 ? hi - lo : 1.);
    double olo, ohi, bw; int nb;
    THLimitsFinder::Optimize(lo - pad_, hi + pad_, 4, olo, ohi, nb, bw, "");   // only use the nice step bw
    double ylo = std::floor((lo - pad_) / bw) * bw;                            // round OUTWARD to that step
    double yhi = std::ceil ((hi + pad_) / bw) * bw;
    h->SetMinimum(ylo); h->SetMaximum(yhi);
    h->Draw("E1 X0 P");

    if (zero_line && ylo < 0 && yhi > 0) {
      TLine *l = new TLine(ax->GetXmin(), 0., ax->GetXmax(), 0.);
      l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
    }
    TLatex lab; lab.SetNDC(); lab.SetTextFont(43); lab.SetTextSize(13); lab.SetTextAlign(23);
    lab.DrawLatex(0.55, bottom ? 0.97 : 0.95, ThetaLabel(kv.second, kv.first));
    k++;
  }
  c->cd();
}

// overall page: sum, diff, asym (keys 0, 1, 2) stacked in one column, each with its own y title
static void DrawOverall(TCanvas *c, const std::map<int, TH1D*> &hs, bool deg) {
  c->Clear();
  const char *ytit[3] = {"N^{+} + N^{-}", "N^{+} - N^{-}", "(N^{+} - N^{-})/(N^{+} + N^{-})"};
  const double xl = 0.12, xr = 0.97, yb = 0.07, yt = 0.94;
  int n = hs.size();
  if (n == 0) return;
  double rowh = (yt - yb) / n;

  c->cd();
  TLatex tt; tt.SetNDC(); tt.SetTextFont(43); tt.SetTextSize(20); tt.SetTextAlign(22);
  TString rng = ThetaLabel(hs.begin()->second, 0);
  tt.DrawLatex(0.5, 0.975, "All #theta_{fpp}: " + rng);

  int k = 0;
  for (auto &kv : hs) {
    bool bottom = (k == n - 1);
    double y2 = yt - k * rowh, y1 = y2 - rowh, ext = bottom ? 0.04 : 0.;
    c->cd();
    TPad *p = new TPad(Form("p_all_%d", k), "", xl - 0.08, y1 - ext, xr, y2);
    p->SetTopMargin(0.); p->SetBottomMargin(bottom ? ext / (rowh + ext) : 0.);
    p->SetLeftMargin(0.12); p->SetRightMargin(0.02); p->SetTicks(1, 1);
    p->Draw(); p->cd();

    TH1D *h = deg ? ToDegrees(kv.second) : (TH1D*)kv.second->Clone(Form("%s_draw", kv.second->GetName()));
    h->SetDirectory(nullptr); h->SetStats(0); h->SetTitle("");
    h->SetMarkerStyle(20); h->SetMarkerSize(0.6); h->SetLineColor(kBlack);
    TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
    ax->SetLabelFont(43); ax->SetLabelSize(bottom ? 14 : 0); ax->SetTitleFont(43); ax->SetTitleSize(16);
    ax->SetTitle(bottom ? (deg ? "#phi_{fpp} (deg)" : "#phi_{fpp} (rad)") : ""); ax->SetTitleOffset(3.2);
    ay->SetLabelFont(43); ay->SetLabelSize(12); ay->SetNdivisions(505);
    ay->SetTitleFont(43); ay->SetTitleSize(14); ay->SetTitleOffset(2.6);
    ay->SetTitle(ytit[kv.first]);

    double lo = h->GetMinimum() - h->GetBinError(h->GetMinimumBin());
    double hi = h->GetMaximum() + h->GetBinError(h->GetMaximumBin());
    double pad_ = 0.12 * (hi - lo > 0 ? hi - lo : 1.);
    double olo, ohi, bw; int nb;
    THLimitsFinder::Optimize(lo - pad_, hi + pad_, 4, olo, ohi, nb, bw, "");
    double ylo = std::floor((lo - pad_) / bw) * bw, yhi = std::ceil((hi + pad_) / bw) * bw;
    h->SetMinimum(ylo); h->SetMaximum(yhi);
    h->Draw("E1 X0 P");

    if (kv.first > 0 && ylo < 0 && yhi > 0) {          // zero line for diff and asym
      TLine *l = new TLine(ax->GetXmin(), 0., ax->GetXmax(), 0.);
      l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
    }

    // asymmetry panel: fit A(phi) = c1 cos(phi) + s1 sin(phi), blue line, results in a box
    if (kv.first == 2) {
      double k2rad = deg ? TMath::DegToRad() : 1.;     // fit in the units of the drawn axis
      TF1 *fa = new TF1("fasym", Form("[0]*cos(x*%.12g) - [1]*sin(x*%.12g)", k2rad, k2rad),
                        ax->GetXmin(), ax->GetXmax());
      fa->SetParNames("c1", "s1");
      fa->SetParameters(0., 0.);
      fa->SetLineColor(kBlue); fa->SetLineWidth(2); fa->SetNpx(400);
      h->Fit(fa, "Q0 R");                               // Q: quiet, 0: do not draw yet, R: full range
      fa->Draw("same");

      double c1 = fa->GetParameter(0), ec1 = fa->GetParError(0);
      double s1 = fa->GetParameter(1), es1 = fa->GetParError(1);
      double chi2 = fa->GetChisquare(); int ndf = fa->GetNDF();
      TPaveText *pt = new TPaveText(0.135, 0.66, 0.36, 0.97, "NDC");
      pt->SetFillColor(kWhite); pt->SetBorderSize(1); pt->SetTextFont(43); pt->SetTextSize(14);
      pt->SetTextAlign(12); pt->SetTextColor(kBlue);
      pt->AddText(Form("c_{1} = %.4f #pm %.4f", c1, ec1));
      pt->AddText(Form("s_{1} = %.4f #pm %.4f", s1, es1));
      pt->AddText(Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, ndf > 0 ? chi2 / ndf : 0.));
      pt->Draw();
      printf("asymmetry fit c1*cos(phi) + s1*sin(phi):  c1 = %.4f +- %.4f   s1 = %.4f +- %.4f   chi2/ndf = %.1f/%d = %.2f\n",
             c1, ec1, s1, es1, chi2, ndf, ndf > 0 ? chi2 / ndf : 0.);
    }
    k++;
  }
  c->cd();
}

// N+ and N- for one theta range, and the sum / difference / asymmetry built from them
struct PhiSet {
  TH1D *hp = nullptr, *hm = nullptr;
  void Book(const TString &base, int nphi, double lo_deg, double hi_deg) {
    TString rng = Form("%.1f < #theta_{fpp} < %.1f deg", lo_deg, hi_deg);
    hp = new TH1D(base + "_hp", "h=+1, " + rng + ";#phi_{fpp} (rad);counts", nphi, -TMath::Pi(), TMath::Pi());
    hm = new TH1D(base + "_hm", "h=-1, " + rng + ";#phi_{fpp} (rad);counts", nphi, -TMath::Pi(), TMath::Pi());
    hp->Sumw2(); hm->Sumw2();
    hp->SetDirectory(nullptr); hm->SetDirectory(nullptr);
    range = rng; name = base;
  }
  TH1D *Sum() {
    TH1D *h = (TH1D*)hp->Clone(name + "_sum");  h->Add(hm);
    h->SetTitle("N^{+} + N^{-}, " + range);  return h;
  }
  TH1D *Diff() {
    TH1D *h = (TH1D*)hp->Clone(name + "_diff"); h->Add(hm, -1.);
    h->SetTitle("N^{+} - N^{-}, " + range);  return h;
  }
  TH1D *Asym() {                                   // correct errors for (N+ - N-)/(N+ + N-)
    int lvl = gErrorIgnoreLevel; gErrorIgnoreLevel = kError;   // hide harmless "Sumw2 already created"
    TH1D *h = (TH1D*)hp->GetAsymmetry(hm);
    gErrorIgnoreLevel = lvl;
    h->SetName(name + "_asym");
    h->SetTitle("(N^{+} - N^{-})/(N^{+} + N^{-}), " + range);  return h;
  }
  TString range, name;
};

int draw_fpp_phi_tree(const char *rootfile, const char *pdf = "plots/fpp_phi_tree.pdf",
                      const char *thbins = "1,2,3,4,5,6,7,8,9", int nphi = 8,
                      double thmin = 1., double thmax = 10.,
                      const char *cut = "fppCut", bool deg = false, const char *outroot = "")
{
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(26., 19.5);                 // landscape PDF page
  TH1::AddDirectory(kFALSE);

  // chain Tout over all matching files (wildcards allowed; several patterns separated by spaces)
  TChain *T = new TChain("Tout");
  TObjArray *pats = TString(rootfile).Tokenize(" ");
  for (int i = 0; i < pats->GetEntries(); i++) {
    TString p = ((TObjString*)pats->At(i))->GetString();
    int nf = T->Add(p);
    if (nf == 0) std::cerr << "WARNING: no files match " << p << std::endl;
  }
  delete pats;
  if (T->GetNtrees() == 0) { std::cerr << "No input files for: " << rootfile << std::endl; return 1; }
  std::cout << "Chained Tout from " << T->GetNtrees() << " files" << std::endl;

  // theta bin edges (deg) from the comma-separated string
  std::vector<double> edges;
  TObjArray *tok = TString(thbins).Tokenize(",");
  for (int i = 0; i < tok->GetEntries(); i++) edges.push_back(((TObjString*)tok->At(i))->GetString().Atof());
  delete tok;
  if (edges.size() < 2) { std::cerr << "need at least two theta bin edges" << std::endl; return 1; }

  // book: one PhiSet per theta bin + one for the overall window
  std::vector<PhiSet> bins(edges.size() - 1);
  for (size_t i = 0; i + 1 < edges.size(); i++) bins[i].Book(Form("h_phi_fpp_th%zu", i), nphi, edges[i], edges[i + 1]);
  PhiSet all; all.Book("h_phi_fpp_all_th0", nphi, thmin, thmax);

  // branches + optional cut
  double theta = 0, phi = 0; int hel = 0;
  T->SetBranchStatus("*", 0);
  for (auto b : {"theta_fpp", "phi_fpp", "helicity"}) T->SetBranchStatus(b, 1);
  T->SetBranchAddress("theta_fpp", &theta);
  T->SetBranchAddress("phi_fpp", &phi);
  T->SetBranchAddress("helicity", &hel);
  TTreeFormula *fcut = nullptr;
  if (cut && strlen(cut) > 0) {
    T->SetBranchStatus("*", 1);                    // the cut may use any branch
    fcut = new TTreeFormula("fcut", cut, T);
    if (!fcut->GetTree() || fcut->GetNdim() == 0) { std::cerr << "bad cut: " << cut << std::endl; return 1; }
  }

  // one pass over the tree
  const double d2r = TMath::DegToRad();
  long n = T->GetEntries(), nused = 0;
  int treenum = -1;
  for (long i = 0; i < n; i++) {
    T->GetEntry(i);
    if (fcut && T->GetTreeNumber() != treenum) {   // chain moved to the next file
      treenum = T->GetTreeNumber();
      fcut->UpdateFormulaLeaves();
    }
    if (hel != 1 && hel != -1) continue;
    if (fcut && fcut->EvalInstance(0) == 0) continue;
    double thdeg = theta / d2r;
    for (size_t b = 0; b + 1 < edges.size(); b++) {
      if (thdeg >= edges[b] && thdeg < edges[b + 1]) { (hel > 0 ? bins[b].hp : bins[b].hm)->Fill(phi); break; }
    }
    if (thdeg >= thmin && thdeg < thmax) { (hel > 0 ? all.hp : all.hm)->Fill(phi); nused++; }
  }
  std::cout << n << " Tout entries, " << nused << " used in " << thmin << "-" << thmax
            << " deg (cut: \"" << (cut ? cut : "") << "\")" << std::endl;

  // counts per theta bin
  printf("%-26s %9s %9s %9s\n", "theta_fpp bin", "N+", "N-", "(N+-N-)/sum");
  for (size_t i = 0; i < bins.size(); i++) {
    double p = bins[i].hp->Integral(), m = bins[i].hm->Integral();
    printf("%-26s %9.0f %9.0f %9.4f\n", Form("%.1f < theta < %.1f deg", edges[i], edges[i + 1]),
           p, m, p + m > 0 ? (p - m) / (p + m) : 0.);
  }
  {
    double p = all.hp->Integral(), m = all.hm->Integral();
    printf("%-26s %9.0f %9.0f %9.4f\n", Form("all %.1f-%.1f deg", thmin, thmax), p, m, p + m > 0 ? (p - m) / (p + m) : 0.);
  }

  // build sum / diff / asym from the totals
  std::map<std::string, std::map<int, TH1D*>> H;
  for (size_t i = 0; i < bins.size(); i++) {
    H["sum"][i] = bins[i].Sum(); H["diff"][i] = bins[i].Diff(); H["asym"][i] = bins[i].Asym();
  }
  std::map<int, TH1D*> Hall = {{0, all.Sum()}, {1, all.Diff()}, {2, all.Asym()}};

  // draw: same layout as draw_fpp_phi.C
  const std::vector<std::string> types = {"sum", "diff", "asym"};
  const std::map<std::string, std::pair<const char*, const char*>> info = {
    {"sum",  {"N^{+} + N^{-}",                        "Helicity sum"}},
    {"diff", {"N^{+} - N^{-}",                        "Helicity difference"}},
    {"asym", {"(N^{+} - N^{-}) / (N^{+} + N^{-})",    "Helicity asymmetry"}}};

  TCanvas *c = new TCanvas("c_fpp_phi_tree", "phi_fpp from Tout", 1200, 900);
  int npages = 4, ip = 0;
  auto print_page = [&]() {
    ip++;
    TString opt = ip == 1 ? "(" : (ip == npages ? ")" : "");
    c->Print(TString(pdf) + opt, "pdf Landscape");
  };
  for (auto &t : types) {
    DrawPage(c, H[t], info.at(t).first, info.at(t).second, deg, t != "sum");
    print_page();
  }
  DrawOverall(c, Hall, deg);
  print_page();
  std::cout << "Written " << pdf << " (" << npages << " pages)" << std::endl;

  // optional: keep the histograms
  if (outroot && strlen(outroot) > 0) {
    TFile fo(outroot, "RECREATE");
    for (auto &b : bins) { b.hp->Write(); b.hm->Write(); }
    all.hp->Write(); all.hm->Write();
    for (auto &t : types) for (auto &kv : H[t]) kv.second->Write();
    for (auto &kv : Hall) kv.second->Write();
    std::cout << "Histograms written to " << outroot << std::endl;
  }
  return 0;
}
