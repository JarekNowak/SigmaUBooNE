// ============================================================================
//  PlotSigma.C  —  standalone plotting for the Sigma analysis.
//
//  Reads SigmaHists.root (produced by SigmaSelection.C) and, for every variable
//  at every cut stage, writes two PNGs:
//
//    figs/stack/<Cut>/<Var>.png   POT-scaled stacked MC (EXT+Dirt+Background+Hyperons)
//                                 with the signal (Sigma) scaled x N and overlaid.
//    figs/shape/<Cut>/<Var>.png   unit-normalized shape comparison of all samples.
//
//  Run:  root -b -l -q PlotSigma.C+      (or  bash runPlots.sh)
// ============================================================================
#include "SigmaSelConfig.h"

#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace SigmaCfg;

// ---- global style -----------------------------------------------------------
static void SetStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleFontSize(0.045);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillStyle(0);
}

// ---- fetch one histogram from the per-cut directory -------------------------
static TH1F *GetHist(TFile *f, int v, int s, int c)
{
  TString name = VarName[v] + Names[s] + CutName[c];
  return (TH1F *)f->Get(CutName[c] + "/" + name);
}

// ---- round a raw scale factor up to a "nice" 1/2/5 x 10^n value -------------
static double NiceFactor(double raw)
{
  if (raw <= 1) return 1.0;
  double p = pow(10.0, floor(log10(raw)));
  double m = raw / p;
  double nm = (m <= 1) ? 1 : (m <= 2) ? 2 : (m <= 5) ? 5 : 10;
  return nm * p;
}

// ---- stacked, POT-scaled plot with signal overlaid -------------------------
static void MakeStackPlot(TFile *f, int v, int c, TCanvas *cv)
{
  TH1F *hsig  = GetHist(f, v, 0, c);
  TH1F *hhyp  = GetHist(f, v, 1, c);
  TH1F *hbkg  = GetHist(f, v, 2, c);
  TH1F *hdirt = GetHist(f, v, 3, c);
  TH1F *hdata = GetHist(f, v, 4, c);
  TH1F *hext  = GetHist(f, v, 5, c);   // absent from files made before EXT was added
  if (!hsig || !hhyp || !hbkg || !hdirt) return;

  // background sum (for max + stack ordering)
  double stackMax = 0;
  TH1F *sum = (TH1F *)hdirt->Clone("sum_tmp");
  sum->Add(hbkg); sum->Add(hhyp);
  if (hext) sum->Add(hext);
  stackMax = sum->GetMaximum();
  delete sum;

  int nb; double lo, hi; bool logy; Binning(v, nb, lo, hi, logy);

  // nothing to show at all -> skip (avoids empty-frame draws)
  if (stackMax <= 0 && hsig->GetMaximum() <= 0) return;

  // background fills
  hhyp->SetFillColor(SampleColor[1]);  hhyp->SetLineColor(kBlack);  hhyp->SetLineWidth(1);
  hbkg->SetFillColor(SampleColor[2]);  hbkg->SetLineColor(kBlack);  hbkg->SetLineWidth(1);
  hdirt->SetFillColor(SampleColor[3]); hdirt->SetLineColor(kBlack); hdirt->SetLineWidth(1);
  if (hext) { hext->SetFillColor(SampleColor[5]); hext->SetLineColor(kBlack); hext->SetLineWidth(1); }

  THStack *st = new THStack(Form("st_%d_%d", v, c),
                            VarName[v] + "  -  " + CutName[c] + ";" + VarName[v] +
                            ";Events (1.1526e21 POT)");   // full RHC, matches SigmaSelection.C dataPOT
  if (hext) st->Add(hext);
  st->Add(hdirt);   // smallest typically at bottom
  st->Add(hbkg);
  st->Add(hhyp);

  // signal scaled x N so its peak ~ 0.3 x stack max
  double sigMax = hsig->GetMaximum();
  double factor = 1.0;
  if (sigMax > 0 && stackMax > 0) factor = NiceFactor(0.3 * stackMax / sigMax);
  TH1F *hsigS = (TH1F *)hsig->Clone(Form("sigS_%d_%d", v, c));
  hsigS->Scale(factor);
  hsigS->SetLineColor(SampleColor[0]);
  hsigS->SetLineWidth(3);
  hsigS->SetFillStyle(0);

  cv->Clear();
  cv->SetLogy(logy && stackMax > 0);

  st->Draw("HIST");
  double top = std::max(stackMax, hsigS->GetMaximum());
  st->SetMaximum(top * (logy ? 50.0 : 1.5));
  if (logy) st->SetMinimum(0.5);

  hsigS->Draw("HIST SAME");
  if (hdata && hdata->GetEntries() > 0) {
    hdata->SetMarkerStyle(20);
    hdata->SetMarkerColor(kBlack);
    hdata->SetLineColor(kBlack);
    hdata->Draw("E SAME");
  }

  TLegend *leg = new TLegend(0.60, 0.62, 0.93, 0.90);
  leg->AddEntry(hsigS, Form("Sigma #times%g  (%.2f)", factor, hsig->Integral()), "l");
  leg->AddEntry(hhyp,  Form("Hyperons  (%.1f)",  hhyp->Integral()),  "f");
  leg->AddEntry(hbkg,  Form("Background  (%.0f)", hbkg->Integral()),  "f");
  leg->AddEntry(hdirt, Form("Dirt  (%.1f)",       hdirt->Integral()), "f");
  if (hext) leg->AddEntry(hext, Form("EXT  (%.1f)", hext->Integral()), "f");
  if (hdata && hdata->GetEntries() > 0)
    leg->AddEntry(hdata, Form("Data  (%.0f)", hdata->Integral()), "lep");
  leg->Draw();

  TString dir = "figs/stack/" + CutName[c];
  gSystem->mkdir(dir, kTRUE);
  cv->Print(dir + "/" + VarName[v] + ".png");

  delete hsigS;
  delete st;
  delete leg;
}

// ---- unit-normalized shape comparison --------------------------------------
static void MakeShapePlot(TFile *f, int v, int c, TCanvas *cv)
{
  TH1F *hn[4] = {0, 0, 0, 0};
  double ymax = 0;
  for (int s = 0; s < 4; s++) {
    TH1F *h = GetHist(f, v, s, c);
    if (!h || h->Integral() <= 0) continue;
    hn[s] = (TH1F *)h->Clone(Form("shape_%d_%d_%d", v, s, c));
    hn[s]->Scale(1.0 / hn[s]->Integral());
    hn[s]->SetLineColor(SampleColor[s]);
    hn[s]->SetLineWidth(2);
    hn[s]->SetFillStyle(0);
    hn[s]->SetMarkerColor(SampleColor[s]);
    ymax = std::max(ymax, hn[s]->GetMaximum());
  }

  cv->Clear();
  cv->SetLogy(false);

  TLegend *leg = new TLegend(0.62, 0.70, 0.93, 0.90);
  bool first = true;
  for (int s = 0; s < 4; s++) {
    if (!hn[s]) continue;
    if (first) {
      hn[s]->SetTitle(VarName[v] + " (shape)  -  " + CutName[c] + ";" + VarName[v] + ";a.u.");
      hn[s]->SetMaximum(ymax * 1.3);
      hn[s]->SetMinimum(0.0);
      hn[s]->Draw("HIST");
      first = false;
    } else {
      hn[s]->Draw("HIST SAME");
    }
    leg->AddEntry(hn[s], Names[s], "l");
  }
  if (first) { delete leg; return; }  // nothing had entries
  leg->Draw();

  TString dir = "figs/shape/" + CutName[c];
  gSystem->mkdir(dir, kTRUE);
  cv->Print(dir + "/" + VarName[v] + ".png");

  for (int s = 0; s < 4; s++) if (hn[s]) delete hn[s];
  delete leg;
}

// ---- entry point ------------------------------------------------------------
void PlotSigma(const char *infile = "SigmaHists.root")
{
  SetStyle();

  TFile *f = TFile::Open(infile);
  if (!f || f->IsZombie()) {
    cout << "ERROR: cannot open " << infile << endl;
    return;
  }

  gSystem->mkdir("figs/stack", kTRUE);
  gSystem->mkdir("figs/shape", kTRUE);

  TCanvas *cv = new TCanvas("cv", "", 800, 600);

  for (int v = 0; v < Nv; v++) {
    for (int c = 0; c < Nc; c++) {
      MakeStackPlot(f, v, c, cv);
      MakeShapePlot(f, v, c, cv);
    }
    cout << "\rPlotted variable " << (v + 1) << " / " << Nv << flush;
  }
  cout << endl;

  f->Close();
  cout << "Done. Plots under figs/stack/<cut>/ and figs/shape/<cut>/" << endl;
}
