/*
   draw_fpp_phi_pub.C  --  publication-style phi_fpp figures from Tout

   Usage:
     root -l -b -q 'draw_fpp_phi_pub.C("../out/runs/gep5_*_GEP5.root","plots/gep5", "Q^{2} = 11 GeV^{2}")'

   Writes one vector PDF per figure (journals want single figures):
     <prefix>_phi_overall.pdf   N+ and N- (top) and the helicity asymmetry with fit (bottom)
     <prefix>_asym_thbins.pdf   asymmetry in each theta_fpp bin, common y range, with fits
     <prefix>_fit_vs_theta.pdf  fitted c1 and s1 versus theta_fpp

   Fit (same sign convention as before):
     A(phi) = c1 cos(phi) - s1 sin(phi)
*/
#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>

#include "TH1D.h"
#include "TPad.h"
#include "TFile.h"
#include "TChain.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TError.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TTreeFormula.h"
#include "THLimitsFinder.h"
#include "TGraphErrors.h"
#include "TF1.h"

// ---------------------------------------------------------------- style
static const int kFont = 43;           // Helvetica, size in pixels -> identical text size in every pad
static const int kColP = kRed + 1;     // h = +1
static const int kColM = kAzure + 2;   // h = -1
static const int kColFit = kRed + 1;

static void SetPubStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetCanvasColor(0); gStyle->SetPadColor(0); gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0); gStyle->SetPadBorderMode(0); gStyle->SetFrameBorderMode(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetEndErrorSize(0);
  gStyle->SetErrorX(0.);               // no horizontal bars: bin width is implied by the point spacing
  gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
  gStyle->SetLineScalePS(1.5);
  TH1::AddDirectory(kFALSE);
}

static void StyleAxes(TH1 *h, int labsz, int titsz, bool showx, double xoff, double yoff, int ylabsz = -1) {
  TAxis *ax = h->GetXaxis(), *ay = h->GetYaxis();
  for (TAxis *a : {ax, ay}) { a->SetLabelFont(kFont); a->SetTitleFont(kFont); a->SetTitleSize(titsz); a->CenterTitle(false); }
  ax->SetLabelSize(showx ? labsz : 0); ax->SetTitleSize(showx ? titsz : 0);
  ay->SetLabelSize(ylabsz < 0 ? labsz : ylabsz);
  ax->SetTitleOffset(xoff); ay->SetTitleOffset(yoff);
  ay->SetMaxDigits(4);
}

// phi axis: -pi ... pi with labels in units of pi (radians), or plain degrees
static void PhiAxis(TH1 *h, bool deg, bool hideFirst = false, bool hideLast = false) {
  TAxis *ax = h->GetXaxis();
  ax->SetTitle(deg ? "#phi_{fpp} (deg)" : "#phi_{fpp} (rad)");
  ax->SetNdivisions(-404);
  if (hideFirst) ax->ChangeLabel(1, -1, 0);      // ChangeLabel: the first call for an index wins
  if (hideLast)  ax->ChangeLabel(5, -1, 0);
  if (!deg) {
    ax->ChangeLabel(1, -1, -1, -1, -1, -1, "#minus#pi");
    ax->ChangeLabel(2, -1, -1, -1, -1, -1, "#minus#pi/2");
    ax->ChangeLabel(3, -1, -1, -1, -1, -1, "0");
    ax->ChangeLabel(4, -1, -1, -1, -1, -1, "#pi/2");
    ax->ChangeLabel(5, -1, -1, -1, -1, -1, "#pi");
  }
}

// ---------------------------------------------------------------- helpers
static TH1D *ToDegrees(const TH1D *h) {
  int n = h->GetNbinsX();
  TH1D *hd = new TH1D(Form("%s_deg", h->GetName()), "", n,
                      h->GetXaxis()->GetXmin() * TMath::RadToDeg(), h->GetXaxis()->GetXmax() * TMath::RadToDeg());
  hd->SetDirectory(nullptr);
  for (int i = 0; i <= n + 1; i++) { hd->SetBinContent(i, h->GetBinContent(i)); hd->SetBinError(i, h->GetBinError(i)); }
  return hd;
}
static TH1D *ForDrawing(const TH1D *h, bool deg) {
  TH1D *d = deg ? ToDegrees(h) : (TH1D*)h->Clone(Form("%s_draw", h->GetName()));
  d->SetDirectory(nullptr); d->SetTitle(""); d->SetStats(0);
  return d;
}

// min/max of content -+ error over filled bins
static void DataRange(const TH1 *h, double &lo, double &hi) {
  lo = 1e300; hi = -1e300;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double v = h->GetBinContent(i), e = h->GetBinError(i);
    if (v == 0 && e == 0) continue;
    lo = std::min(lo, v - e); hi = std::max(hi, v + e);
  }
  if (lo > hi) { lo = -1; hi = 1; }
}
// round a range outward to a "nice" step so the end labels are round numbers
static void NiceRange(double lo, double hi, double &ylo, double &yhi, int ndiv = 4) {
  double olo, ohi, bw; int nb;
  THLimitsFinder::Optimize(lo, hi, ndiv, olo, ohi, nb, bw, "");
  ylo = std::floor(lo / bw + 1e-9) * bw; yhi = std::ceil(hi / bw - 1e-9) * bw;
}

struct FitRes { double c1 = 0, ec1 = 0, s1 = 0, es1 = 0, chi2 = 0; int ndf = 0; };

// fit is always done on the radian histogram
static FitRes FitAsym(TH1D *hrad) {
  TF1 f("fasym_tmp", "[0]*cos(x) - [1]*sin(x)", -TMath::Pi(), TMath::Pi());
  f.SetParameters(0., 0.);
  hrad->Fit(&f, "Q N 0 R");
  FitRes r;
  r.c1 = f.GetParameter(0); r.ec1 = f.GetParError(0);
  r.s1 = f.GetParameter(1); r.es1 = f.GetParError(1);
  r.chi2 = f.GetChisquare(); r.ndf = f.GetNDF();
  return r;
}
static TF1 *FitCurve(const FitRes &r, bool deg, const char *name) {
  double k = deg ? TMath::DegToRad() : 1.;
  double xm = deg ? 180. : TMath::Pi();
  TF1 *f = new TF1(name, Form("%.10g*cos(x*%.12g) - %.10g*sin(x*%.12g)", r.c1, k, r.s1, k), -xm, xm);
  f->SetLineColor(kColFit); f->SetLineWidth(2); f->SetNpx(400);
  return f;
}
static void ZeroLine(const TH1 *h) {
  TLine *l = new TLine(h->GetXaxis()->GetXmin(), 0., h->GetXaxis()->GetXmax(), 0.);
  l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->SetLineWidth(1); l->Draw();
}
static void StyleMarker(TH1 *h, int col, int style, double size) {
  h->SetMarkerStyle(style); h->SetMarkerSize(size); h->SetMarkerColor(col); h->SetLineColor(col); h->SetLineWidth(1);
}
static TString ThRange(double a, double b) { return Form("%g^{#circ} < #theta_{fpp} < %g^{#circ}", a, b); }

// PDF page = canvas aspect ratio, so the figure fills the page with no white border
static void PrintFig(TCanvas *c, const TString &out) {
  double w = c->GetWw(), h = c->GetWh();
  gStyle->SetPaperSize(20., 20. * h / w);
  c->Print(out);
}

// N+ and N- for one theta range
struct PhiSet {
  TH1D *hp = nullptr, *hm = nullptr; TString name;
  double lo = 0, hi = 0;
  void Book(const TString &base, int nphi, double lo_deg, double hi_deg) {
    name = base; lo = lo_deg; hi = hi_deg;
    hp = new TH1D(base + "_hp", ";#phi_{fpp} (rad);counts", nphi, -TMath::Pi(), TMath::Pi());
    hm = new TH1D(base + "_hm", ";#phi_{fpp} (rad);counts", nphi, -TMath::Pi(), TMath::Pi());
    hp->Sumw2(); hm->Sumw2();
  }
  TH1D *Asym() const {
    int lvl = gErrorIgnoreLevel; gErrorIgnoreLevel = kError;
    TH1D *h = (TH1D*)hp->GetAsymmetry(hm);
    gErrorIgnoreLevel = lvl;
    h->SetName(name + "_asym"); h->SetDirectory(nullptr);
    return h;
  }
};

// ---------------------------------------------------------------- figure 1
static void FigOverall(const PhiSet &all, const char *label, bool deg, const TString &out, FitRes &fr) {
  const int W = 800, H = 900, labsz = 24, titsz = 28;
  TCanvas *c = new TCanvas("c_overall", "", W, H);
  const double split = 0.47;                                   // bottom pad height fraction

  // top: counts
  c->cd();
  TPad *pt = new TPad("pt", "", 0, split, 1, 1);
  pt->SetLeftMargin(0.17); pt->SetRightMargin(0.04); pt->SetTopMargin(0.08); pt->SetBottomMargin(0.0);
  pt->Draw(); pt->cd();
  TH1D *hp = ForDrawing(all.hp, deg), *hm = ForDrawing(all.hm, deg);
  StyleMarker(hp, kColP, 20, 1.3); StyleMarker(hm, kColM, 25, 1.3);
  double lo1, hi1, lo2, hi2; DataRange(hp, lo1, hi1); DataRange(hm, lo2, hi2);
  double ylo, yhi; NiceRange(0., 1.45 * std::max(hi1, hi2), ylo, yhi, 5);
  hp->SetMinimum(ylo); hp->SetMaximum(yhi);
  hp->GetYaxis()->SetTitle("Counts");
  PhiAxis(hp, deg);
  StyleAxes(hp, labsz, titsz, false, 0, 1.9);
  hp->GetYaxis()->SetNdivisions(505);
  hp->GetYaxis()->ChangeLabel(1, -1, 0);                      // the "0" label would be cut by the panel edge
  hp->Draw("E1 P"); hm->Draw("E1 P SAME");

  TLegend *lg = new TLegend(0.62, 0.70, 0.93, 0.90);
  lg->SetTextFont(kFont); lg->SetTextSize(labsz);
  lg->AddEntry(hp, "h = +1", "p"); lg->AddEntry(hm, "h = #minus1", "p");
  lg->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(kFont); t.SetTextSize(labsz); t.SetTextAlign(13);
  if (label && strlen(label)) t.DrawLatex(0.21, 0.90, label);
  t.DrawLatex(0.21, 0.90 - (label && strlen(label) ? 0.09 : 0.), ThRange(all.lo, all.hi));

  // bottom: asymmetry
  c->cd();
  TPad *pb = new TPad("pb", "", 0, 0, 1, split);
  pb->SetLeftMargin(0.17); pb->SetRightMargin(0.04); pb->SetTopMargin(0.0); pb->SetBottomMargin(0.24);
  pb->Draw(); pb->cd();
  TH1D *arad = all.Asym();
  fr = FitAsym(arad);
  TH1D *a = ForDrawing(arad, deg);
  StyleMarker(a, kBlack, 20, 1.3);
  double alo, ahi; DataRange(a, alo, ahi);
  double amax = 1.6 * std::max(std::fabs(alo), std::fabs(ahi));  // headroom for the fit results
  NiceRange(-amax, amax, ylo, yhi, 4);
  a->SetMinimum(ylo); a->SetMaximum(yhi);
  a->GetYaxis()->SetTitle("A_{#phi}");
  PhiAxis(a, deg);
  StyleAxes(a, labsz, titsz, true, 1.15, 1.9);
  a->GetYaxis()->SetNdivisions(504);
  a->GetYaxis()->ChangeLabel(-1, -1, 0);                       // top label would touch the panel above
  a->Draw("E1 P");
  ZeroLine(a);
  FitCurve(fr, deg, "fasym_overall")->Draw("SAME");
  a->Draw("E1 P SAME");                                         // points on top of the curve

  TLatex f; f.SetNDC(); f.SetTextFont(kFont); f.SetTextSize(labsz - 2); f.SetTextAlign(13);
  f.DrawLatex(0.21, 0.94, Form("s_{1} = %.4f #pm %.4f", fr.s1, fr.es1));
  f.DrawLatex(0.21, 0.84, Form("c_{1} = %.4f #pm %.4f", fr.c1, fr.ec1));
  f.SetTextAlign(33);
  f.DrawLatex(0.93, 0.94, Form("#chi^{2}/ndf = %.1f/%d", fr.chi2, fr.ndf));

  PrintFig(c, out);
}

// ---------------------------------------------------------------- figure 2
static void FigThetaBins(const std::vector<PhiSet> &bins, const char *label, bool deg, const TString &out,
                         std::vector<FitRes> &fits) {
  int n = bins.size();
  int ncol = n <= 3 ? n : (n <= 4 ? 2 : (n <= 9 ? 3 : 4));
  int nrow = (n + ncol - 1) / ncol;
  const int W = 1200, H = 260 * nrow + 140, labsz = 22, titsz = 26;
  TCanvas *c = new TCanvas("c_thbins", "", W, H);

  std::vector<TH1D*> hr(n), hd(n);
  double gmax = 0;
  fits.assign(n, FitRes());
  for (int i = 0; i < n; i++) {
    hr[i] = bins[i].Asym(); fits[i] = FitAsym(hr[i]); hd[i] = ForDrawing(hr[i], deg);
    double lo, hi; DataRange(hd[i], lo, hi); gmax = std::max(gmax, std::max(std::fabs(lo), std::fabs(hi)));
  }
  double ylo, yhi; NiceRange(-1.1 * gmax, 1.1 * gmax, ylo, yhi, 4);   // one y range for all pads

  // pad grid: outer margins in canvas fractions, panels touch each other
  const double L = 0.09, R = 0.015, B = 110. / H, T = 30. / H;
  double pw = (1 - L - R) / ncol, ph = (1 - B - T) / nrow;
  for (int i = 0; i < n; i++) {
    int col = i % ncol, row = i / ncol;
    bool left = (col == 0), bottom = (row == nrow - 1) || (i + ncol >= n);
    double x1 = L + col * pw, x2 = x1 + pw, y2 = 1 - T - row * ph, y1 = y2 - ph;
    double ex1 = left ? L : 0., ey1 = bottom ? B : 0.;          // extend outer pads to hold the labels
    c->cd();
    TPad *p = new TPad(Form("pth%d", i), "", x1 - ex1, y1 - ey1, x2, y2);
    p->SetLeftMargin(ex1 / (pw + ex1)); p->SetRightMargin(0.);
    p->SetBottomMargin(ey1 / (ph + ey1)); p->SetTopMargin(0.);
    p->Draw(); p->cd();

    TH1D *h = hd[i];
    StyleMarker(h, kBlack, 20, 1.0);
    h->SetMinimum(ylo); h->SetMaximum(yhi);
    PhiAxis(h, deg, col > 0, col < ncol - 1);   // no clipped labels where panels touch
    h->GetYaxis()->SetTitle("A_{#phi}");
    StyleAxes(h, labsz, titsz, bottom, 0, 2.2, left ? labsz : 0);
    h->GetXaxis()->SetTitleSize(0);                            // common x title drawn on the canvas
    if (!left) h->GetYaxis()->SetTitleSize(0);
    h->GetYaxis()->SetNdivisions(504);
    if (row > 0)  h->GetYaxis()->ChangeLabel(-1, -1, 0);
    if (!bottom)  h->GetYaxis()->ChangeLabel(1, -1, 0);
    h->GetXaxis()->SetTickLength(0.04);
    h->Draw("E1 P");
    ZeroLine(h);
    FitCurve(fits[i], deg, Form("fth%d", i))->Draw("SAME");
    h->Draw("E1 P SAME");

    TLatex t; t.SetNDC(); t.SetTextFont(kFont); t.SetTextSize(labsz - 2); t.SetTextAlign(13);
    double xl = p->GetLeftMargin() + 0.04, yt = 1 - p->GetTopMargin() - 0.05;
    t.DrawLatex(xl, yt, ThRange(bins[i].lo, bins[i].hi));
  }
  c->cd();
  TLatex t; t.SetNDC(); t.SetTextFont(kFont);
  t.SetTextSize(titsz); t.SetTextAlign(22);
  t.DrawLatex(L + 0.5 * (1 - L - R), 0.35 * B, deg ? "#phi_{fpp} (deg)" : "#phi_{fpp} (rad)");
  if (label && strlen(label)) {
    t.SetTextSize(labsz); t.SetTextAlign(31);
    t.DrawLatex(1 - R, 1 - T + 0.006, label);
  }
  PrintFig(c, out);
}

// ---------------------------------------------------------------- figure 3
static void FigFitVsTheta(const std::vector<PhiSet> &bins, const std::vector<FitRes> &fits,
                          const char *label, const TString &out) {
  int n = bins.size();
  TGraphErrors *gs = new TGraphErrors(n), *gc = new TGraphErrors(n);
  double lo = 1e9, hi = -1e9;
  for (int i = 0; i < n; i++) {
    double x = 0.5 * (bins[i].lo + bins[i].hi), dx = 0.06 * (bins[i].hi - bins[i].lo);   // small offset, no overlap
    gs->SetPoint(i, x - dx, fits[i].s1); gs->SetPointError(i, 0, fits[i].es1);
    gc->SetPoint(i, x + dx, fits[i].c1); gc->SetPointError(i, 0, fits[i].ec1);
    lo = std::min({lo, fits[i].s1 - fits[i].es1, fits[i].c1 - fits[i].ec1});
    hi = std::max({hi, fits[i].s1 + fits[i].es1, fits[i].c1 + fits[i].ec1});
  }
  const int labsz = 24, titsz = 28;
  TCanvas *c = new TCanvas("c_fitth", "", 800, 600);
  c->SetLeftMargin(0.16); c->SetRightMargin(0.04); c->SetTopMargin(0.05); c->SetBottomMargin(0.15);
  double xlo = bins.front().lo, xhi = bins.back().hi;
  double span = hi - lo, ylo, yhi;
  NiceRange(std::min(lo - 0.15 * span, 0.), std::max(hi + 0.45 * span, 0.), ylo, yhi, 5);
  TH1D *fr = new TH1D("frame_fitth", "", 1, xlo, xhi);
  fr->SetMinimum(ylo); fr->SetMaximum(yhi);
  fr->GetXaxis()->SetTitle("#theta_{fpp} (deg)");
  fr->GetYaxis()->SetTitle("Fit amplitude");
  StyleAxes(fr, labsz, titsz, true, 1.3, 1.9);
  fr->GetYaxis()->SetNdivisions(505);
  fr->Draw("AXIS");
  ZeroLine(fr);
  gs->SetMarkerStyle(20); gs->SetMarkerSize(1.4); gs->SetMarkerColor(kBlack); gs->SetLineColor(kBlack);
  gc->SetMarkerStyle(24); gc->SetMarkerSize(1.4); gc->SetMarkerColor(kGray + 2); gc->SetLineColor(kGray + 2);
  gc->Draw("P Z SAME"); gs->Draw("P Z SAME");
  TLegend *lg = new TLegend(0.66, 0.76, 0.93, 0.92);
  lg->SetTextFont(kFont); lg->SetTextSize(labsz);
  lg->AddEntry(gs, "s_{1}", "p"); lg->AddEntry(gc, "c_{1}", "p");
  lg->Draw();
  if (label && strlen(label)) {
    TLatex t; t.SetNDC(); t.SetTextFont(kFont); t.SetTextSize(labsz); t.SetTextAlign(13);
    t.DrawLatex(0.20, 0.91, label);
  }
  PrintFig(c, out);
}

// ---------------------------------------------------------------- main
int draw_fpp_phi_pub(const char *rootfile, const char *prefix = "plots/fpp",
                     const char *label = "", const char *thbins = "1,2,3,4,5,6,7,8,9", int nphi = 18,
                     double thmin = 1., double thmax = 10., const char *cut = "fppCut",
                     bool deg = false, const char *outroot = "")
{
  SetPubStyle();

  TChain *T = new TChain("Tout");
  TObjArray *pats = TString(rootfile).Tokenize(" ");
  for (int i = 0; i < pats->GetEntries(); i++) {
    TString p = ((TObjString*)pats->At(i))->GetString();
    if (T->Add(p) == 0) std::cerr << "WARNING: no files match " << p << std::endl;
  }
  delete pats;
  if (T->GetNtrees() == 0) { std::cerr << "No input files for: " << rootfile << std::endl; return 1; }
  std::cout << "Chained Tout from " << T->GetNtrees() << " files" << std::endl;

  std::vector<double> edges;
  TObjArray *tok = TString(thbins).Tokenize(",");
  for (int i = 0; i < tok->GetEntries(); i++) edges.push_back(((TObjString*)tok->At(i))->GetString().Atof());
  delete tok;
  if (edges.size() < 2) { std::cerr << "need at least two theta bin edges" << std::endl; return 1; }

  std::vector<PhiSet> bins(edges.size() - 1);
  for (size_t i = 0; i + 1 < edges.size(); i++) bins[i].Book(Form("h_phi_fpp_th%zu", i), nphi, edges[i], edges[i + 1]);
  PhiSet all; all.Book("h_phi_fpp_all", nphi, thmin, thmax);

  double theta = 0, phi = 0; int hel = 0;
  T->SetBranchStatus("*", 0);
  for (auto b : {"theta_fpp", "phi_fpp", "helicity"}) T->SetBranchStatus(b, 1);
  T->SetBranchAddress("theta_fpp", &theta);
  T->SetBranchAddress("phi_fpp", &phi);
  T->SetBranchAddress("helicity", &hel);
  TTreeFormula *fcut = nullptr;
  if (cut && strlen(cut) > 0) {
    T->SetBranchStatus("*", 1);
    fcut = new TTreeFormula("fcut", cut, T);
    if (!fcut->GetTree() || fcut->GetNdim() == 0) { std::cerr << "bad cut: " << cut << std::endl; return 1; }
  }

  const double r2d = TMath::RadToDeg();
  long n = T->GetEntries(), nused = 0;
  int treenum = -1;
  for (long i = 0; i < n; i++) {
    T->GetEntry(i);
    if (fcut && T->GetTreeNumber() != treenum) { treenum = T->GetTreeNumber(); fcut->UpdateFormulaLeaves(); }
    if (hel != 1 && hel != -1) continue;
    if (fcut && fcut->EvalInstance(0) == 0) continue;
    double thdeg = theta * r2d;
    for (size_t b = 0; b + 1 < edges.size(); b++)
      if (thdeg >= edges[b] && thdeg < edges[b + 1]) { (hel > 0 ? bins[b].hp : bins[b].hm)->Fill(phi); break; }
    if (thdeg >= thmin && thdeg < thmax) { (hel > 0 ? all.hp : all.hm)->Fill(phi); nused++; }
  }
  std::cout << n << " Tout entries, " << nused << " used in " << thmin << "-" << thmax
            << " deg (cut: \"" << (cut ? cut : "") << "\")" << std::endl;

  TString pre(prefix);
  FitRes fr; std::vector<FitRes> fits;
  FigOverall(all, label, deg, pre + "_phi_overall.pdf", fr);
  FigThetaBins(bins, label, deg, pre + "_asym_thbins.pdf", fits);
  FigFitVsTheta(bins, fits, label, pre + "_fit_vs_theta.pdf");

  printf("\nfit A(phi) = c1 cos(phi) - s1 sin(phi)\n");
  printf("%-22s %9s %9s %10s %8s %10s %8s %12s\n", "theta_fpp (deg)", "N+", "N-", "c1", "dc1", "s1", "ds1", "chi2/ndf");
  for (size_t i = 0; i < bins.size(); i++)
    printf("%6.1f - %-13.1f %9.0f %9.0f %10.4f %8.4f %10.4f %8.4f %8.1f/%d\n", bins[i].lo, bins[i].hi,
           bins[i].hp->Integral(), bins[i].hm->Integral(), fits[i].c1, fits[i].ec1, fits[i].s1, fits[i].es1,
           fits[i].chi2, fits[i].ndf);
  printf("%6.1f - %-13.1f %9.0f %9.0f %10.4f %8.4f %10.4f %8.4f %8.1f/%d\n", all.lo, all.hi,
         all.hp->Integral(), all.hm->Integral(), fr.c1, fr.ec1, fr.s1, fr.es1, fr.chi2, fr.ndf);
  std::cout << "Written " << pre << "_{phi_overall,asym_thbins,fit_vs_theta}.pdf" << std::endl;

  if (outroot && strlen(outroot) > 0) {
    TFile fo(outroot, "RECREATE");
    for (auto &b : bins) { b.hp->Write(); b.hm->Write(); b.Asym()->Write(); }
    all.hp->Write(); all.hm->Write(); all.Asym()->Write();
    std::cout << "Histograms written to " << outroot << std::endl;
  }
  return 0;
}
