// ============================================================================
//  validate_closure.C -- closure test of the weights on generated samples: GENIE as distributed
//  (Pais) reweighted with GenieHypXsecWeight.h against GENIE with the corrected cross section,
//  both NuMI RHC on argon (MicroBooNE). Reads the chunk directories of the two productions
//  (chunk.summary: POT; chunk.hepmcev: the Sigma0 events as HepMC-like text: status pdg m1 m2 d1 d2
//  px py pz E m x y z t), and compares numubar CC Sigma0 per POT, in total and in E_nu / Q2 bins.
//  Expected (2026-11-02): old 191.9, old x w 127.9, fixed 127.4 per 1e21 POT.
//  Run:  root -l -b -q 'validate_closure.C("<old chunks dir>", "<fixed chunks dir>", 500)'
// ============================================================================
#include "GenieHypXsecWeight.h"
#include "TSystem.h"
#include <fstream>
#include <sstream>
#include <vector>

static double readChunks(const char* base, int maxChunks, const GenieHypXsecWeight::Weights* W,
                         double& n, double& nw, std::vector<double>& hE, std::vector<double>& hEw)
{
  double pot = 0; int used = 0;
  for (int c = 0; c < 100000 && used < maxChunks; ++c) {
    TString d = Form("%s/c%04d", base, c);
    if (gSystem->AccessPathName(d + "/DONE")) { if (gSystem->AccessPathName(d)) break; else continue; }
    std::ifstream s(d + "/chunk.summary"); std::string k; double v;
    while (s >> k) { if (k == "pot") { s >> v; pot += v; } else s.ignore(1000, '\n'); }
    std::ifstream in(d + "/chunk.hepmcev"); std::string line;
    while (std::getline(in, line)) {
      std::istringstream hs(line); int ie, np; if (!(hs >> ie >> np) || !hs.eof()) continue;
      TLorentzVector nu, nuc, lep; int nupdg = 0; bool hasMu = false, hasY = false;
      for (int i = 0; i < np && std::getline(in, line); ++i) {
        std::istringstream ps(line); int st, pdg, a, b, cc, dd; double px, py, pz, E;
        ps >> st >> pdg >> a >> b >> cc >> dd >> px >> py >> pz >> E;
        if (st == 0) { nu.SetPxPyPzE(px, py, pz, E); nupdg = pdg; }
        else if (st == 11) nuc.SetPxPyPzE(px, py, pz, E);
        else if (st == 1 && pdg == -13) { lep.SetPxPyPzE(px, py, pz, E); hasMu = true; }
        else if (st == 1 && pdg == 3212) hasY = true;
      }
      if (nupdg != -14 || !hasMu || !hasY) continue;
      const double w = W ? W->weight(nupdg, 3212, nu, nuc, lep) : 1.;
      const int b = std::min(5, (int)(nu.E() / 1.0));            // E_nu bins of 1 GeV, last open
      n += 1; nw += w; hE[b] += 1; hEw[b] += w;
    }
    ++used;
  }
  printf("  %s: %d chunks, %.4g POT\n", base, used, pot);
  return pot;
}

void validate_closure(const char* oldDir, const char* newDir, int maxOld = 500)
{
  GenieHypXsecWeight::Weights W; if (!W.load()) return;
  double no = 0, nwo = 0, nn = 0, nwn = 0;
  std::vector<double> eo(6), ewo(6), en(6), ewn(6);
  const double po = readChunks(oldDir, maxOld, &W, no, nwo, eo, ewo);
  const double pn = readChunks(newDir, 100000, nullptr, nn, nwn, en, ewn);
  const double s = 1e21;
  printf("numubar CC Sigma0 per 1e21 POT:  old %.1f   old x w %.1f +- %.1f   fixed %.1f +- %.1f\n",
         no / po * s, nwo / po * s, nwo / po * s / std::sqrt(no), nn / pn * s, nn / pn * s / std::sqrt(nn));
  printf("E_nu bin [GeV]   old x w / fixed\n");
  for (int b = 0; b < 6; ++b)
    if (en[b] > 0 && eo[b] > 0)
      printf("  %d-%s      %.3f +- %.3f\n", b, b == 5 ? "inf" : Form("%d", b + 1), (ewo[b] / po) / (en[b] / pn),
             (ewo[b] / po) / (en[b] / pn) * std::sqrt(1. / eo[b] + 1. / en[b]));
}
