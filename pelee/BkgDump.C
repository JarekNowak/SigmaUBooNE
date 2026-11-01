// ============================================================================
//  BkgDump.C -- the full truth of the background events an island-cache selection keeps,
//  read back from the hypana files (the cache holds only the candidates' true pdg).
//  One text line per event: run sub evt, neutrino pdg, CC/NC, GENIE mode, hyperon pdg,
//  true vertex in the FV, lepton momentum, primary pdg list with momenta, and every PFP as
//  pdg/origin/track score/hits (HyperonProduction origin: 1 nu primary, 2 Lambda daughter,
//  3 other, 5 Sigma0 daughter, -1 no MC match = overlaid cosmic data).
//  Run: root -l -b -q 'BkgDump.C+("lists/bkg_overlay.list","isl_bkg.root","<sel>","bkgdump.txt")'
// ============================================================================
#include "HypInput.h"
#include "SigmaFV.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <set>
#include <map>
#include <cstdio>

static unsigned long long rseK(unsigned r, unsigned s, unsigned e)
{ return ((unsigned long long)r << 40) | ((unsigned long long)s << 20) | e; }

void BkgDump(const char* list, const char* cache, const char* sel, const char* out)
{
  std::set<unsigned long long> want;
  { TFile c(cache); TTree* t = (TTree*)c.Get("isl");
    t->SetEstimate(t->GetEntries() + 1);
    t->Draw("run:sub:evt", sel, "goff para");
    for (Long64_t i = 0; i < t->GetSelectedRows(); ++i)
      want.insert(rseK((unsigned)t->GetVal(0)[i], (unsigned)t->GetVal(1)[i], (unsigned)t->GetVal(2)[i])); }
  printf("  %zu events selected from %s\n", want.size(), cache);
  FILE* fo = fopen(out, "w");
  std::ifstream lf(list); std::string fn; std::set<unsigned long long> done; long nfile = 0;
  while (std::getline(lf, fn) && done.size() < want.size()) {
    if (fn.empty()) continue;
    TFile* f = TFile::Open(fn.c_str()); if (!f || f->IsZombie()) { delete f; continue; }
    TTree* t = (TTree*)f->Get("ana/OutputTree"); if (!t) { delete f; continue; }
    nfile++;
    UInt_t r, s, e; std::vector<Long64_t> hits;
    t->SetBranchStatus("*", 0); t->SetBranchStatus("run", 1); t->SetBranchStatus("subrun", 1); t->SetBranchStatus("event", 1);
    t->SetBranchAddress("run", &r); t->SetBranchAddress("subrun", &s); t->SetBranchAddress("event", &e);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); const auto k = rseK(r, s, e);
      if (want.count(k) && !done.count(k)) hits.push_back(i);
    }
    t->ResetBranchAddresses(); t->SetBranchStatus("*", 1);
    if (hits.empty()) { delete f; continue; }
    HypInput H; if (!H.attach(t)) { delete f; continue; }
    for (Long64_t i : hits) {
      t->GetEntry(i); H.convert();
      done.insert(rseK(H.h_run, H.h_subrun, H.h_event));
      fprintf(fo, "%d %d %d nu %d %s mode %s hyp %d trueFV %d pl %.3f | prim:", H.run, H.sub, H.evt, H.nu_pdg,
              H.ccnc ? "NC" : "CC", H.h_mode ? H.h_mode->c_str() : "?", H.h_hyperon_pdg,
              (int)SigmaFV::isFiducial(H.true_vtx_x, H.true_vtx_y, H.true_vtx_z), H.h_lepton_mom);
      if (H.h_mc_primary_pdg) for (size_t k = 0; k < H.h_mc_primary_pdg->size(); ++k)
        fprintf(fo, " %d(%.2f)", H.h_mc_primary_pdg->at(k),
                (H.h_mc_primary_mom && k < H.h_mc_primary_mom->size()) ? H.h_mc_primary_mom->at(k) : -1.0);
      fprintf(fo, " | pfp:");
      for (size_t k = 0; k < H.trk_score.size(); ++k)
        fprintf(fo, " %d/%d/%.2f/%u",
                (H.h_pfp_true_pdg && k < H.h_pfp_true_pdg->size()) ? H.h_pfp_true_pdg->at(k) : 0,
                (H.h_pfp_true_origin && k < H.h_pfp_true_origin->size()) ? H.h_pfp_true_origin->at(k) : 0,
                H.trk_score[k], H.pfnhits[k]);
      fprintf(fo, "\n");
    }
    delete f;
  }
  fclose(fo);
  printf("  %zu of %zu found in %ld files -> %s\n", done.size(), want.size(), nfile, out);
}
