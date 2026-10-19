// ============================================================================
//  NoteFigures.C -- figures for the Sigma0 analysis note and slides.
//
//    bdt_score.pdf       PeLEE chain: cross-validated BDTG score at NTracksCut,
//                        background stacked by sample, NuWro signal x100
//    island_display.pdf  Lambda-style chain: the island test on the collection
//                        plane for one signal and one background event (hits
//                        coloured by cluster, seeds at the mu / p / pi starts)
//    lambda_mass.pdf     Lambda-style chain: p-pi mass after the pointing cut and
//                        photon gap, before the mass window and the island test
//
//  The expected-limit figure is written by Sensitivity.C (its plotFile argument).
//
//  Run:  root -l -b -q 'NoteFigures.C+("<outdir>")'
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TColor.h"
#include <vector>
#include <cmath>
#include <cstdio>
#include "Island.h"

namespace {
void style(){
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  gStyle->SetPadLeftMargin(0.13); gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.04); gStyle->SetPadTopMargin(0.06);
  gStyle->SetLabelSize(0.045, "XY"); gStyle->SetTitleSize(0.05, "XY");
  gStyle->SetLegendBorderSize(0); gStyle->SetLegendFillColor(0);
}

// ── BDT score ──────────────────────────────────────────────────────────────
void bdtScore(const char* in, const char* out, double cut){
  TFile f(in);
  TH1D* h[6]; const char* nm[6] = {"sig","hyp","ovl","dirt","data","ext"};
  for (int k = 0; k < 6; ++k) { h[k] = new TH1D(nm[k], ";BDTG score;events / 0.05 (1.15#times10^{21} POT)", 40, -1, 1); h[k]->SetDirectory(0); }
  for (int isS = 0; isS < 2; ++isS) {
    TTree* t = (TTree*)f.Get(isS ? "TreeS" : "TreeB");
    float w, sc, smp; t->SetBranchAddress("Weight", &w); t->SetBranchAddress("BdtScore", &sc); t->SetBranchAddress("Sample", &smp);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); h[isS ? 0 : (int)smp]->Fill(sc, w); }
  }
  const int col[6] = {kBlack, kViolet - 4, kAzure - 9, kOrange - 3, 0, kGreen - 6};
  THStack st("st", ";BDTG score;events / 0.05");
  for (int k : {1, 3, 5, 2}) { h[k]->SetFillColor(col[k]); h[k]->SetLineColor(col[k]); st.Add(h[k]); }
  TH1D* hs = (TH1D*)h[0]->Clone("sigx"); hs->Scale(100); hs->SetLineColor(kRed + 1); hs->SetLineWidth(3);
  TCanvas c("c", "", 800, 600); c.SetLogy();
  st.SetMinimum(0.05); st.SetMaximum(2e6); st.Draw("hist");
  st.GetXaxis()->SetTitle("BDTG score (cross-validated)"); st.GetYaxis()->SetTitle("events / 0.05 at 1.15#times10^{21} POT");
  hs->Draw("hist same");
  TLine l(cut, 0.05, cut, 2e3); l.SetLineStyle(2); l.SetLineWidth(2); l.Draw();
  TLegend lg(0.16, 0.64, 0.56, 0.93); lg.SetTextSize(0.037);
  lg.AddEntry(hs, "#Sigma^{0} signal (NuWro) #times 100", "l");
  lg.AddEntry(h[2], "#nu overlay (GENIE)", "f"); lg.AddEntry(h[5], "beam-off (EXT)", "f");
  lg.AddEntry(h[3], "dirt", "f"); lg.AddEntry(h[1], "other hyperons (NuWro)", "f");
  lg.Draw();
  TLatex tx; tx.SetTextSize(0.038); tx.DrawLatex(cut - 0.13, 3e3, Form("cut %.3f", cut));
  tx.SetNDC(); tx.SetTextSize(0.035); tx.DrawLatex(0.66, 0.89, "RHC, 1.15#times10^{21} POT");
  tx.DrawLatex(0.66, 0.84, "at NTracksCut");
  c.SaveAs(out);
}

// ── island display ─────────────────────────────────────────────────────────
struct Ev { std::vector<int> w; std::vector<float> t; int seedW[3]; float seedT[3]; int evt; double flight, mass; };

// cluster labels as Island::planeTestP builds them (same connectivity)
std::vector<int> labels(const Ev& e, int& nl){
  const int n = e.w.size(); std::vector<int> lab(n, -1); nl = 0;
  for (int i = 0; i < n; ++i) { if (lab[i] >= 0) continue;
    std::vector<int> st{i}; lab[i] = nl;
    while (!st.empty()) { int q = st.back(); st.pop_back();
      for (int j = 0; j < n; ++j)
        if (lab[j] < 0 && std::fabs(e.t[j] - e.t[q]) <= Island::CDT && std::abs(e.w[j] - e.w[q]) <= Island::CDW) { lab[j] = nl; st.push_back(j); } }
    ++nl; }
  return lab;
}
int seedCluster(const Ev& e, const std::vector<int>& lab, int s){
  double best = 1e18; int bi = -1;
  for (size_t i = 0; i < e.w.size(); ++i) {
    const double dw = e.w[i] - e.seedW[s], dt = (e.t[i] - e.seedT[s]) * 0.055 / 0.3, d = dw * dw + dt * dt;
    if (d < best) { best = d; bi = i; } }
  return (bi >= 0 && best < 400.0) ? lab[bi] : -1;
}
// pass = p and pi share a cluster that the muon is not in; mode 1 wants a pass with
// flight in [3,6] cm, mode 0 a fail with the muon in the p/pi cluster
bool pick(const char* fn, int mode, Ev& out){
  TFile f(fn); TTree* t = (TTree*)f.Get("cache");
  std::vector<int>* hw = nullptr; std::vector<float>* ht = nullptr;
  int seedW[3][3], evt; float seedT[3]; double flight, mass, len, cosA;
  t->SetBranchAddress("hw2", &hw); t->SetBranchAddress("ht2", &ht);
  t->SetBranchAddress("seedW", seedW); t->SetBranchAddress("seedT", seedT); t->SetBranchAddress("evt", &evt);
  t->SetBranchAddress("trueFlight", &flight); t->SetBranchAddress("mass", &mass);
  t->SetBranchAddress("decayLen", &len); t->SetBranchAddress("cosAlpha", &cosA);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (!(std::fabs(mass - 1.115683) < 0.020 && len > 0.5 && cosA > 0.95)) continue;
    if (mode == 1 && !(flight > 3 && flight < 6)) continue;
    if (hw->size() < 30 || hw->size() > 400) continue;
    Ev e; e.w = *hw; e.t = *ht; e.evt = evt; e.flight = flight; e.mass = mass;
    for (int k = 0; k < 3; ++k) { e.seedW[k] = seedW[2][k]; e.seedT[k] = seedT[k]; }
    int nl; const auto lab = labels(e, nl);
    const int m = seedCluster(e, lab, 0), p = seedCluster(e, lab, 1), q = seedCluster(e, lab, 2);
    if (m < 0 || p < 0 || q < 0) continue;
    const bool pass = (p == q) && (m != p);
    if ((mode == 1 && pass) || (mode == 0 && !pass && m == p && m == q)) { out = e; return true; }
  }
  return false;
}
void drawEvent(const Ev& e, const char* title){
  int nl; const auto lab = labels(e, nl);
  const int cl[3] = { seedCluster(e, lab, 0), seedCluster(e, lab, 1), seedCluster(e, lab, 2) };
  double w0 = 1e9, w1 = -1e9, t0 = 1e9, t1 = -1e9;
  for (int k = 0; k < 3; ++k) { w0 = std::min(w0, (double)e.seedW[k]); w1 = std::max(w1, (double)e.seedW[k]);
                                t0 = std::min(t0, (double)e.seedT[k]); t1 = std::max(t1, (double)e.seedT[k]); }
  w0 -= 25; w1 += 25; t0 -= 250; t1 += 250;
  TH1D* fr = new TH1D(Form("fr%d", e.evt), Form(";collection-plane wire;time tick"), 10, w0, w1);
  fr->GetYaxis()->SetTitleOffset(1.45);
  fr->SetMinimum(t0); fr->SetMaximum(t1); fr->Draw("axis"); fr->SetDirectory(0);
  // colour: the cluster of each seed, everything else grey
  const int seedCol[3] = { kAzure + 1, kRed + 1, kGreen + 2 };
  for (size_t i = 0; i < e.w.size(); ++i) {
    if (e.w[i] < w0 || e.w[i] > w1 || e.t[i] < t0 || e.t[i] > t1) continue;
    int c = kGray + 1;
    if (lab[i] == cl[0]) c = seedCol[0];
    if (lab[i] == cl[1] || lab[i] == cl[2]) c = (cl[1] == cl[2]) ? kOrange + 7 : (lab[i] == cl[1] ? seedCol[1] : seedCol[2]);
    if (lab[i] == cl[0] && (cl[0] == cl[1] || cl[0] == cl[2])) c = kMagenta + 1;
    TGraph* g = new TGraph(1); g->SetPoint(0, e.w[i], e.t[i]); g->SetMarkerStyle(20); g->SetMarkerSize(0.55); g->SetMarkerColor(c); g->Draw("P");
  }
  const char* sn[3] = {"#mu", "p", "#pi"};
  for (int k = 0; k < 3; ++k) {
    TGraph* g = new TGraph(1); g->SetPoint(0, e.seedW[k], e.seedT[k]); g->SetMarkerStyle(29); g->SetMarkerSize(2.6);
    g->SetMarkerColor(kBlack); g->Draw("P");
    TLatex* tx = new TLatex(e.seedW[k] + 1.2, e.seedT[k] + 25, sn[k]); tx->SetTextSize(0.05); tx->Draw();
  }
  TLatex* tt = new TLatex(); tt->SetNDC(); tt->SetTextSize(0.045); tt->DrawLatex(0.15, 0.95, title);
  // colour key
  struct K { int c; const char* t; };
  const K key[4] = { {kAzure + 1, "#mu cluster"}, {kOrange + 7, "p+#pi cluster"},
                     {kMagenta + 1, "#mu+p+#pi one cluster"}, {kGray + 1, "other activity"} };
  for (int k = 0; k < 4; ++k) {
    TLatex* kt = new TLatex(); kt->SetNDC(); kt->SetTextSize(0.036); kt->SetTextColor(key[k].c);
    kt->DrawLatex(0.62, 0.86 - 0.05 * k, key[k].t); }
}
void islandDisplay(const char* sigCache, const char* bkgCache, const char* out){
  Ev s, b;
  if (!pick(sigCache, 1, s) || !pick(bkgCache, 0, b)) { printf("  no display events found\n"); return; }
  printf("  island display: signal evt %d (true flight %.2f cm, m %.4f), background evt %d (m %.4f)\n",
         s.evt, s.flight, s.mass, b.evt, b.mass);
  TCanvas c("ci", "", 1300, 560); c.Divide(2, 1);
  c.cd(1); gPad->SetTopMargin(0.08); gPad->SetLeftMargin(0.15); drawEvent(s, Form("signal: p+#pi merge, #mu apart (flight %.1f cm)", s.flight));
  c.cd(2); gPad->SetTopMargin(0.08); gPad->SetLeftMargin(0.15); drawEvent(b, "background: p, #pi and #mu one cluster");
  c.SaveAs(out);
}

// ── Lambda mass in the Lambda-style chain ─────────────────────────────────
void lambdaMass(const char* sigFile, const char* bkgFile, const char* out){
  TH1D* h[2];
  for (int k = 0; k < 2; ++k) {
    h[k] = new TH1D(Form("m%d", k), ";m(p#pi^{-}) [GeV];events / 10 MeV", 22, 1.07, 1.29); h[k]->SetDirectory(0);
    TFile f(k ? bkgFile : sigFile); TTree* t = (TTree*)f.Get("isl");
    int hasPair = 0; double w = 0, mass = -1, len = -1, cosA = -2, phDist = -1;
    t->SetBranchAddress("hasPair", &hasPair); t->SetBranchAddress("weight", &w); t->SetBranchAddress("mass", &mass);
    t->SetBranchAddress("decayLen", &len); t->SetBranchAddress("cosAlpha", &cosA); t->SetBranchAddress("phDist", &phDist);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i);
      if (hasPair && len > 0.5 && cosA > 0.95 && phDist > 2.0) h[k]->Fill(mass, w); }
  }
  h[1]->SetFillColor(kAzure - 9); h[1]->SetLineColor(kAzure - 9);
  h[0]->Scale(10); h[0]->SetLineColor(kRed + 1); h[0]->SetLineWidth(3);
  TCanvas c("cm", "", 800, 600);
  h[1]->SetMaximum(1.25 * std::max(h[1]->GetMaximum(), h[0]->GetMaximum())); h[1]->SetMinimum(0);
  h[1]->Draw("hist"); h[0]->Draw("hist same");
  const double mL = 1.115683;
  TLine a(mL - 0.020, 0, mL - 0.020, h[1]->GetMaximum()), b(mL + 0.020, 0, mL + 0.020, h[1]->GetMaximum());
  for (TLine* l : {&a, &b}) { l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
  TLegend lg(0.52, 0.72, 0.93, 0.92); lg.SetTextSize(0.04);
  lg.AddEntry(h[0], "#Sigma^{0} signal (GENIE) #times 10", "l"); lg.AddEntry(h[1], "#nu overlay background", "f"); lg.Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
  tx.DrawLatex(0.52, 0.66, "decay > 0.5 cm, cos#alpha > 0.95,");
  tx.DrawLatex(0.52, 0.61, "photon gap > 2 cm; before the island");
  tx.DrawLatex(0.52, 0.56, "RHC, 1.15#times10^{21} POT");
  c.SaveAs(out);
}
}

void NoteFigures(const char* outDir = "figs/note", double bdtCut = 0.759)
{
  style();
  gSystem->mkdir(outDir, true);
  bdtScore("BDTInput.root", Form("%s/bdt_score.pdf", outDir), bdtCut);
  islandDisplay("islcache_sig.root", "islcache_bkg.root", Form("%s/island_display.pdf", outDir));
  lambdaMass("isl_sig_genie.root", "isl_bkg.root", Form("%s/lambda_mass.pdf", outDir));
}
