// ============================================================================
//  example_gst.C -- weight the hyperon QE events of a GENIE gst file (gntpc -f gst).
//  Selection: nubar probe (neu < 0), CC, QEL, a Lambda / Sigma0 / Sigma- among the final-state
//  particles (GENIE does not re-interact hyperons, so the primary hyperon is in pdgf).
//  Prints the weighted and unweighted counts per hyperon and writes a friend tree "hypw"
//  (one entry per gst entry, w = 1 for every other event).
//  Run:  root -l -b -q 'example_gst.C("gntp.0.gst.root", "hypw_friend.root")'
// ============================================================================
#include "GenieHypXsecWeight.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>

void example_gst(const char* gst, const char* out = "hypw_friend.root",
                 const char* tmu = "hyp_xsec_weights_mu.root", const char* te = "hyp_xsec_weights_e.root")
{
  GenieHypXsecWeight::Weights W;
  if (!W.load(tmu, te)) return;
  TFile f(gst); TTree* t = (TTree*)f.Get("gst");
  if (!t) { printf("no gst tree in %s\n", gst); return; }
  int neu, nf, pdgf[250]; bool cc, qel;
  double Ev, pxv, pyv, pzv, En, pxn, pyn, pzn, El, pxl, pyl, pzl;
  t->SetBranchAddress("neu", &neu); t->SetBranchAddress("cc", &cc); t->SetBranchAddress("qel", &qel);
  t->SetBranchAddress("nf", &nf); t->SetBranchAddress("pdgf", pdgf);
  t->SetBranchAddress("Ev", &Ev); t->SetBranchAddress("pxv", &pxv); t->SetBranchAddress("pyv", &pyv); t->SetBranchAddress("pzv", &pzv);
  t->SetBranchAddress("En", &En); t->SetBranchAddress("pxn", &pxn); t->SetBranchAddress("pyn", &pyn); t->SetBranchAddress("pzn", &pzn);
  t->SetBranchAddress("El", &El); t->SetBranchAddress("pxl", &pxl); t->SetBranchAddress("pyl", &pyl); t->SetBranchAddress("pzl", &pzl);

  TFile fo(out, "RECREATE"); TTree tw("hypw", "GENIE hyperon QE cross-section weight (friend of gst)");
  double w; int hyp; tw.Branch("w", &w, "w/D"); tw.Branch("hyperon", &hyp, "hyperon/I");
  std::map<int, double> sw, n;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i); w = 1; hyp = 0;
    if (neu < 0 && cc && qel) {
      for (int k = 0; k < nf; ++k) if (GenieHypXsecWeight::isHyperonQE(pdgf[k])) { hyp = pdgf[k]; break; }
      if (hyp) {
        TLorentzVector nu(pxv, pyv, pzv, Ev), nuc(pxn, pyn, pzn, En), lep(pxl, pyl, pzl, El);
        w = W.weight(neu, hyp, nu, nuc, lep);
        sw[hyp] += w; n[hyp] += 1;
      }
    }
    tw.Fill();
  }
  fo.cd(); tw.Write();
  printf("hyperon   events   weighted   weighted/events\n");
  for (auto& [p, c] : n) printf("%7d %8.0f %10.1f %10.3f\n", p, c, sw[p], sw[p] / c);
  printf("friend tree -> %s (gst->AddFriend(\"hypw\", \"%s\"))\n", out, out);
}
