/*
   draw_fpp_phi.C
   Draws the per-theta_fpp phi_fpp histograms written by gep_elastic_ana:
       h_phi_fpp_th<N>_sum    N+ + N-
       h_phi_fpp_th<N>_diff   N+ - N-
       h_phi_fpp_th<N>_asym   (N+ - N-)/(N+ + N-)
   One page per type, theta bins stacked in two columns (panels touching,
   x axis labelled only on the bottom panel of each column).
   Last page: the overall (all theta_fpp) histograms from fppPhiAll,
       h_phi_fpp_all_th0_sum / _diff / _asym, stacked in one column.

   Usage:
     root -l -b -q 'draw_fpp_phi.C("../out/gep5_GEP5.root","../out/gep5_phi_fpp.pdf")'
   Optional 3rd argument: true -> phi axis in degrees (-180..180) instead of radians.
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
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "THLimitsFinder.h"

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
    k++;
  }
  c->cd();
}

int draw_fpp_phi(const char *rootfile, const char *pdf = "plots/fpp_phi.pdf", bool deg = false)
{
  gStyle->SetOptStat(0);
  gStyle->SetPaperSize(26., 19.5);     // landscape PDF page (cm), same 4:3 aspect as the 1200x900 canvas
  TFile *f = TFile::Open(rootfile, "READ");
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << rootfile << std::endl; return 1; }

  // collect h_phi_fpp_th<N>_<type>; Get() by name returns the highest cycle (;2 over ;1)
  const std::vector<std::string> types = {"sum", "diff", "asym"};
  std::map<std::string, std::map<int, TH1D*>> H;
  TPRegexp re("^h_phi_fpp_th([0-9]+)_(sum|diff|asym)$");
  TIter next(f->GetListOfKeys());
  while (TKey *key = (TKey*)next()) {
    TString name = key->GetName();
    TObjArray *m = re.MatchS(name);
    if (m->GetEntries() == 3) {
      int ith = ((TObjString*)m->At(1))->GetString().Atoi();
      std::string type = ((TObjString*)m->At(2))->GetString().Data();
      if (!H[type].count(ith)) {
        TH1D *h = dynamic_cast<TH1D*>(f->Get(name));
        if (h) H[type][ith] = h;
      }
    }
    delete m;
  }

  // overall (single theta bin) histograms from fppPhiAll
  std::map<int, TH1D*> Hall;          // 0 = sum, 1 = diff, 2 = asym
  const char *alltypes[3] = {"sum", "diff", "asym"};
  for (int i = 0; i < 3; i++) {
    TH1D *h = dynamic_cast<TH1D*>(f->Get(Form("h_phi_fpp_all_th0_%s", alltypes[i])));
    if (h) Hall[i] = h;
  }

  const std::map<std::string, std::pair<const char*, const char*>> info = {
    {"sum",  {"N^{+} + N^{-}",                        "Helicity sum"}},
    {"diff", {"N^{+} - N^{-}",                        "Helicity difference"}},
    {"asym", {"(N^{+} - N^{-}) / (N^{+} + N^{-})",    "Helicity asymmetry"}}};

  TCanvas *c = new TCanvas("c_fpp_phi", "phi_fpp per theta bin", 1200, 900);
  int npages = 0;
  for (auto &t : types) { if (!H[t].empty()) npages++; }
  if (!Hall.empty()) npages++;                                  // overall page
  if (npages == 0) { std::cerr << "No h_phi_fpp_th*/all histograms in " << rootfile << std::endl; return 1; }

  // multi-page PDF: "(" on the first page, ")" on the last, whichever pages exist
  int ip = 0;
  auto print_page = [&]() {
    ip++;
    TString opt = (npages == 1) ? "" : (ip == 1 ? "(" : (ip == npages ? ")" : ""));
    c->Print(TString(pdf) + opt, "pdf Landscape");
  };

  for (auto &t : types) {
    if (H[t].empty()) { std::cout << "No '" << t << "' histograms, skipping that page." << std::endl; continue; }
    std::cout << t << ": " << H[t].size() << " theta bins" << std::endl;
    DrawPage(c, H[t], info.at(t).first, info.at(t).second, deg, t != "sum");
    print_page();
  }
  if (!Hall.empty()) {
    std::cout << "overall: " << Hall.size() << " histograms" << std::endl;
    DrawOverall(c, Hall, deg);
    print_page();
  }
  std::cout << "Written " << pdf << " (" << npages << " pages)" << std::endl;
  return 0;
}
