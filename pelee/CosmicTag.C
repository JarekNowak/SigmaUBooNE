// ============================================================================
//  CosmicTag.C -- event-level cosmic origin for an IslandProcess.C cache, as a
//  friend tree, read from the hypana truth without touching the wire branches.
//
//  WHY: the hypana EXT is too small to be followed to the Lambda-style operating
//  point (150 data events per EXT event), so its contribution there has to be
//  factorised: EXT rate at an early stage x the pass rate of the later cuts for
//  COSMIC-INDUCED candidates. The overlay MC has such candidates (the selected
//  slice is a cosmic in an event that also holds a neutrino), but the caches do
//  not record it. The slice-truth branches (true_nu_slice_ID, purity) are not
//  usable for this: in the overlay the flash-matched slice never equals
//  true_nu_slice_ID and the purities are ~1%.
//
//  cosFrac = fraction of the selected slice's PFP hits (U+V+Y) whose PFP has no MC
//  match (pfp_true_origin == -1): in the overlay those are the cosmic DATA hits.
//  pfp_true_origin is HyperonProduction's own enum (1 neutrino primary, 2 Lambda
//  daughter, 3 other, 5 Sigma0 daughter, -1 unmatched), not simb::Origin.
//  -1 if the event is not found or has no hits. Entry i of the friend matches
//  entry i of the cache.
//
//  Run:  root -l -b -q 'CosmicTag.C("lists/bkg_overlay.list","isl_bkg.root","isl_bkg_cos.root")'
//  Use:  TFile f("isl_bkg.root"); TTree* t=(TTree*)f.Get("isl"); t->AddFriend("cos","isl_bkg_cos.root");
// ============================================================================
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include <unordered_map>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>

void CosmicTag(const char* list = "lists/bkg_overlay.list", const char* cache = "isl_bkg.root",
               const char* outName = "isl_bkg_cos.root")
{
  auto key = [](long r, long s, long e){ return (unsigned long long)((r << 40) | (s << 20) | e); };

  // cache keys -> entry
  TFile fc(cache); TTree* T = (TTree*)fc.Get("isl");
  if (!T) { printf("  no isl tree in %s\n", cache); return; }
  int run=0, sub=0, evt=0;
  T->SetBranchAddress("run",&run); T->SetBranchAddress("sub",&sub); T->SetBranchAddress("evt",&evt);
  const Long64_t N = T->GetEntries();
  std::unordered_map<unsigned long long, Long64_t> idx; idx.reserve(N);
  long dupCache = 0;
  for (Long64_t i = 0; i < N; ++i) { T->GetEntry(i); if (!idx.emplace(key(run,sub,evt), i).second) dupCache++; }
  std::vector<float> cosFrac(N, -1.f);

  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while (std::getline(in, l)) if (l.size()) { ch.Add(l.c_str()); nf++; }
  ch.SetBranchStatus("*", 0);
  for (const char* b : {"run","subrun","event","pfp_true_origin","pfpnplanehits_U","pfpnplanehits_V","pfpnplanehits_Y"})
    ch.SetBranchStatus(b, 1);
  UInt_t r=0, s=0, e=0;
  std::vector<int> *org=nullptr, *hu=nullptr, *hv=nullptr, *hy=nullptr;
  ch.SetBranchAddress("run",&r); ch.SetBranchAddress("subrun",&s); ch.SetBranchAddress("event",&e);
  ch.SetBranchAddress("pfp_true_origin",&org);
  ch.SetBranchAddress("pfpnplanehits_U",&hu); ch.SetBranchAddress("pfpnplanehits_V",&hv);
  ch.SetBranchAddress("pfpnplanehits_Y",&hy);
  const Long64_t M = ch.GetEntries();
  printf("  %s: %lld cached events (%ld repeated keys); %d files, %lld events\n", cache, N, dupCache, nf, M);
  long found = 0;
  for (Long64_t j = 0; j < M; ++j) {
    ch.GetEntry(j);
    if (j % 1000000 == 0) { printf("  ... %lld / %lld\n", j, M); fflush(stdout); }
    auto it = idx.find(key(r,s,e));
    if (it == idx.end()) continue;
    double all = 0, cos = 0;
    for (size_t k = 0; org && k < org->size(); ++k) {
      const double h = (hu && k<hu->size() ? hu->at(k) : 0) + (hv && k<hv->size() ? hv->at(k) : 0)
                     + (hy && k<hy->size() ? hy->at(k) : 0);
      all += h; if (org->at(k) == -1) cos += h;
    }
    if (cosFrac[it->second] < 0) found++;
    cosFrac[it->second] = all > 0 ? (float)(cos / all) : -1.f;
  }
  TFile fo(outName, "RECREATE"); TTree F("cos", "cosmic fraction of the selected slice, friend of isl");
  float cf = 0; F.Branch("cosFrac", &cf, "cosFrac/F");
  for (Long64_t i = 0; i < N; ++i) { cf = cosFrac[i]; F.Fill(); }
  F.Write(); fo.Close();
  printf("  matched %ld of %lld cached events -> %s\n", found, N, outName);
}
