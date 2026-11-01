// ============================================================================
//  ShowerOrigin.C -- the truth behind every shower-like PFP (track score < 0.5) in Sigma0
//  signal events, to see what the 2nd and 3rd showers allowed by 0 < nShowers < 4 are.
//  Per event: nSh and the number of showers of each class, from pfp_true_pdg and
//  pfp_true_origin (HyperonProduction: 1 nu primary, 2 Lambda daughter, 3 other, 5 Sigma0
//  daughter, -1 no MC match = overlaid cosmic data):
//    nGam  the Sigma0 photon (22, origin 5)      nLp / nLpi  Lambda proton / pi- (origin 2)
//    nMu   the muon (|13|, origin 1)             nPrim  other primaries (origin 1)
//    nOth  secondaries (origin 3: Michel e+, delta rays, hadron interactions)   nCos  origin -1
//  The best-matched shower of each class is also stored by hits. Joined to the island cache
//  (run, sub, evt) to take the classes at any stage of the Lambda-style chain.
//  Run: root -l -b -q 'ShowerOrigin.C+("lists/sig_genie_qe.list","sho_genie.root")'
// ============================================================================
#include "HypInput.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <cstdio>

void ShowerOrigin(const char* list, const char* out)
{
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while (std::getline(in, l)) if (l.size()) { ch.Add(l.c_str()); nf++; }
  HypInput H; if (!H.attach(&ch)) { printf("  not a hypana tree\n"); return; }

  TFile fo(out, "RECREATE"); TTree T("sho", "shower truth classes per signal event");
  int run, sub, evt, nSh, nTr, nGam, nLp, nLpi, nMu, nPrim, nOth, nCos, oth3pdg;
  T.Branch("run", &run, "run/I"); T.Branch("sub", &sub, "sub/I"); T.Branch("evt", &evt, "evt/I");
  T.Branch("nSh", &nSh, "nSh/I"); T.Branch("nTr", &nTr, "nTr/I");
  T.Branch("nGam", &nGam, "nGam/I"); T.Branch("nLp", &nLp, "nLp/I"); T.Branch("nLpi", &nLpi, "nLpi/I");
  T.Branch("nMu", &nMu, "nMu/I"); T.Branch("nPrim", &nPrim, "nPrim/I"); T.Branch("nOth", &nOth, "nOth/I");
  T.Branch("nCos", &nCos, "nCos/I"); T.Branch("othPdg", &oth3pdg, "othPdg/I");

  const Long64_t N = ch.GetEntries();
  printf("  %s: %d files, %lld events\n", list, nf, N);
  for (Long64_t ie = 0; ie < N; ++ie) {
    ch.GetEntry(ie); H.convert();
    if (!(H.h_hyperon_pdg == 3212 && H.nu_pdg == -14 && H.ccnc == 0)) continue;
    run = H.run; sub = H.sub; evt = H.evt;
    nSh = nTr = nGam = nLp = nLpi = nMu = nPrim = nOth = nCos = 0; oth3pdg = 0;
    const int n = (int)H.trk_score.size();
    for (int i = 0; i < n; ++i) {
      if (H.trk_score[i] > 0.5) { nTr++; continue; }
      nSh++;
      const int tp = (H.h_pfp_true_pdg && i < (int)H.h_pfp_true_pdg->size()) ? H.h_pfp_true_pdg->at(i) : 0;
      const int to = (H.h_pfp_true_origin && i < (int)H.h_pfp_true_origin->size()) ? H.h_pfp_true_origin->at(i) : 0;
      if (to == 5 && tp == 22) nGam++;
      else if (to == 2 && tp == 2212) nLp++;
      else if (to == 2 && tp == -211) nLpi++;
      else if (to == 1 && std::abs(tp) == 13) nMu++;
      else if (to == 1) nPrim++;
      else if (to == -1) nCos++;
      else { nOth++; if (!oth3pdg) oth3pdg = tp; }
    }
    T.Fill();
  }
  fo.Write();
  printf("  -> %s\n", out);
}
