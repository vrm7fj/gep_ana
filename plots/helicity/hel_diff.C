// hel_diff.C : per-helicity normalized phi_fpp, sum, difference and asymmetry
void hel_diff(double thmin_deg = 3., double thmax_deg = 4., int nbins = 24)
{
  TChain *T = new TChain("Tout");
  T->Add("../../out/runs/gep5_*.root");
  if (T->GetNtrees() == 0) { std::cout << "no files matched" << std::endl; return; }

  TString base = Form("elasticCut && fppCut && theta_fpp>%f && theta_fpp<%f",
                      thmin_deg * TMath::DegToRad(), thmax_deg * TMath::DegToRad());
  TString thlab = Form("%.1f < #theta_{fpp} < %.1f deg", thmin_deg, thmax_deg);

  TH1D *hp = new TH1D("hp", thlab + ";#phi_{fpp} (rad);f^{#pm}", nbins, -TMath::Pi(), TMath::Pi());
  TH1D *hm = new TH1D("hm", "h = -1;#phi_{fpp} (rad);f^{-}", nbins, -TMath::Pi(), TMath::Pi());
  hp->Sumw2(); hm->Sumw2();
  T->Draw("phi_fpp>>hp", base + " && helicity==1",  "goff");
  T->Draw("phi_fpp>>hm", base + " && helicity==-1", "goff");
  double Np = hp->Integral(), Nm = hm->Integral();
  std::cout << "N+ = " << Np << ", N- = " << Nm << std::endl;
  if (Np <= 0 || Nm <= 0) return;

  // normalize each helicity to its own total: average bin = 1
  hp->Scale(nbins / Np);
  hm->Scale(nbins / Nm);

  // S = (f+ + f-)/2  (acceptance)
  TH1D *hsum = (TH1D*)hp->Clone("hsum");
  hsum->Add(hm); hsum->Scale(0.5);
  hsum->SetTitle(thlab + ";#phi_{fpp} (rad);(f^{+} + f^{-})/2");

  // D = (f+ - f-)/2  (signal x acceptance)
  TH1D *hd = (TH1D*)hp->Clone("hd");
  hd->Add(hm, -1.); hd->Scale(0.5);
  hd->SetTitle(thlab + ";#phi_{fpp} (rad);(f^{+} - f^{-})/2");

  // A = (f+ - f-)/(f+ + f-) = D/S  (acceptance cancels to first order)
  TH1D *ha = (TH1D*)hd->Clone("ha");
  ha->Divide(hsum);
  ha->SetTitle(thlab + ";#phi_{fpp} (rad);(f^{+} - f^{-})/(f^{+} + f^{-})");

  // fits
  TF1 *fs = new TF1("fs", "[0]*(1+[1]*cos(x)+[2]*sin(x)+[3]*cos(2*x)+[4]*sin(2*x))", -TMath::Pi(), TMath::Pi());
  fs->SetParameters(1, 0, 0, 0, 0);
  fs->SetParNames("c0", "a1 (cos)", "b1 (sin)", "a2 (cos2)", "b2 (sin2)");
  TF1 *fd = new TF1("fd", "[0] + [1]*cos(x) + [2]*sin(x)", -TMath::Pi(), TMath::Pi());
  fd->SetParNames("const", "c (cos)", "s (sin)");
  TF1 *fa = (TF1*)fd->Clone("fa");

  gStyle->SetOptFit(1111);
  TCanvas *c = new TCanvas("c_heldiff", "helicity sum / difference", 1400, 1000);
  c->Divide(2, 2);

  c->cd(1);
  hp->SetLineColor(2); hm->SetLineColor(4);
  hp->SetMinimum(0); hp->SetMaximum(1.2 * std::max(hp->GetMaximum(), hm->GetMaximum()));
  hp->SetStats(0); hm->SetStats(0);
  hp->Draw("E"); hm->Draw("E same");
  TLegend *leg = new TLegend(0.70, 0.75, 0.88, 0.88);
  leg->AddEntry(hp, "h = +1", "l"); leg->AddEntry(hm, "h = -1", "l"); leg->Draw();

  c->cd(2); hsum->Draw("E"); hsum->Fit(fs, "R");
  c->cd(3); hd->Draw("E");   hd->Fit(fd, "R");
  c->cd(4); ha->Draw("E");   ha->Fit(fa, "R");
}
