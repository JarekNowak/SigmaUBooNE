// ============================================================================
//  PreselFlow.C -- the Lambda-style chain's preselection, event by event, to see
//  where signal is lost before the island cache (IslandProcess.C caches only
//  events with a muon candidate) and what relaxing each requirement would recover.
//
//  Same stages and constants as IslandProcess.C: SliceID (flash-matched nu slice),
//  reco vertex in the FV, topological score > 0.20, 0 < nShowers < 4, nTracks < 5,
//  muon candidate (track score > 0.8, within 20 cm, length > 10 cm, LLR > 0.7).
//  Recorded per event, with no cut applied: the stage reached, nTr, nSh, topo, whether
//  a muon candidate exists regardless of the earlier cuts, the number of non-muon
//  tracks (a pair needs two), the true vertex in the FV, and the truth (signal flag).
//  Signal-definition flags as SigmaSelection.C: the true muon momentum (pmuTrue) and any other
//  meson at the primary vertex (pi0, eta, eta', K -- the list of SigmaSelection.C; mesonVeto). The analysis signal is
//  isSig && trueFV && pmuTrue > 0.15 && !mesonVeto.
//
//  Output tree "pf"; run: root -l -b -q 'PreselFlow.C+("lists/sig_genie_qe.list","pf_sig.root",true)'
// ============================================================================
#include "HypInput.h"
#include "SigmaFV.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <fstream>
#include <cstdio>

void PreselFlow(const char* list, const char* out, bool signalOnly = true)
{
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while (std::getline(in, l)) if (l.size()) { ch.Add(l.c_str()); nf++; }
  double topo = 0; ch.SetBranchAddress("topological_score", &topo);
  HypInput H; if (!H.attach(&ch)) { printf("  not a hypana tree\n"); return; }
  int fmSlice = -1; const bool haveSlice = ch.GetBranch("flash_match_nu_slice_ID") != nullptr;
  if (haveSlice) ch.SetBranchAddress("flash_match_nu_slice_ID", &fmSlice);

  TFile fo(out, "RECREATE"); TTree T("pf", "Lambda-style preselection flow");
  int run, sub, evt, isSig, trueFV, slice, recoFV, nTr, nSh, hasMu, nOther, stage, mesonVeto;
  double topoS, pmuTrue;
  T.Branch("run", &run, "run/I"); T.Branch("sub", &sub, "sub/I"); T.Branch("evt", &evt, "evt/I");
  T.Branch("isSig", &isSig, "isSig/I"); T.Branch("trueFV", &trueFV, "trueFV/I");
  T.Branch("slice", &slice, "slice/I"); T.Branch("recoFV", &recoFV, "recoFV/I"); T.Branch("topo", &topoS, "topo/D");
  T.Branch("nTr", &nTr, "nTr/I"); T.Branch("nSh", &nSh, "nSh/I"); T.Branch("hasMu", &hasMu, "hasMu/I");
  T.Branch("nOther", &nOther, "nOther/I"); T.Branch("stage", &stage, "stage/I");
  T.Branch("pmuTrue", &pmuTrue, "pmuTrue/D"); T.Branch("mesonVeto", &mesonVeto, "mesonVeto/I");

  const Long64_t N = ch.GetEntries();
  printf("  %s: %d files, %lld events\n", list, nf, N);
  for (Long64_t ie = 0; ie < N; ++ie) {
    ch.GetEntry(ie); H.convert();
    isSig = (H.h_hyperon_pdg == 3212 && H.nu_pdg == -14 && H.ccnc == 0);
    if (signalOnly != (bool)isSig) continue;
    run = H.run; sub = H.sub; evt = H.evt;
    trueFV = SigmaFV::isFiducial(H.true_vtx_x, H.true_vtx_y, H.true_vtx_z);
    pmuTrue = H.h_lepton_mom;
    mesonVeto = 0;
    if (H.h_mc_primary_pdg) for (int pd : *H.h_mc_primary_pdg) {
      const int ap = std::abs(pd);
      if (ap == 111 || ap == 221 || ap == 331 || ap == 130 || ap == 310 || ap == 311 || ap == 321) mesonVeto = 1;
    }
    slice = !(haveSlice && fmSlice < 0);
    recoFV = SigmaFV::isFiducial(H.vtx_x, H.vtx_y, H.vtx_z);
    topoS = topo;
    const int n = (int)H.trk_score.size();
    const TVector3 vtx(H.vtx_x, H.vtx_y, H.vtx_z);
    nTr = 0; nSh = 0; hasMu = 0; nOther = 0; int mu = -1;
    for (int i = 0; i < n; ++i) {
      const bool trk = H.trk_score[i] > 0.5; (trk ? nTr : nSh)++;
      if (!trk) continue;
      const double d = (TVector3(H.trk_sx[i], H.trk_sy[i], H.trk_sz[i]) - vtx).Mag();
      if (H.trk_score[i] > 0.8 && d < 20.0 && H.trk_len[i] > 10.0 && H.trk_llr[i] > 0.7)
        if (mu == -1 || H.trk_len[i] > H.trk_len[mu]) mu = i;
    }
    hasMu = mu >= 0; nOther = nTr - (hasMu ? 1 : 0);
    // stage reached in IslandProcess.C order: 0 All, 1 Slice, 2 FV, 3 Cosmic, 4 NShowers/NTracks, 5 Muon
    stage = 0;
    if (slice) { stage = 1; if (recoFV) { stage = 2; if (topo >= 0.20) { stage = 3;
      if (nSh > 0 && nSh < 4 && nTr < 5) { stage = 4; if (hasMu) stage = 5; } } } }
    T.Fill();
  }
  fo.Write();
  printf("  -> %s\n", out);
}
