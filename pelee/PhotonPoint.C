// ============================================================================
//  PhotonPoint.C -- photon-pointing variables for the Lambda-style chain, added
//  after the fact for the events of an island cache (IslandProcess.C output) that
//  pass a selection, without re-running the island test.
//
//  The Sigma0 decays at the primary vertex, so its photon points back to it. The
//  background with a real Lambda (GENIE numubar CC QE Lambda) has no such photon:
//  the "photon" the chain picks there is an unrelated shower (66-452 cm away).
//
//  For each selected event: the chain's photon (nearest shower-like PFP to the
//  reco vertex, trk_score <= 0.5, as IslandProcess.C; checked against the cache's
//  phDist) and
//    phCos    = shower direction . unit(shower start - vertex)    (1 = points away from the vertex)
//    phImpact = distance of the vertex from the shower axis [cm]
//  plus the alternative photon choice "best pointing shower" (largest phCos among
//  shower-like PFPs within 100 cm): its distance, energy and phCos.
//
//  Output tree "pp": run sub evt, distChk (|dist - cache phDist|), phCos, phImpact,
//  altIdx, altDist, altE, altCos, altSame (alt == chain photon), and the chain photon's
//  shower properties: hits per plane (nhU/V/Y), length, opening angle, Moliere avg/rms,
//  dE/dx and energy per plane, truth pdg and origin (HyperonProduction: 1 nu primary,
//  2 Lambda daughter, 3 other, 5 Sigma0 daughter, -1 unmatched = cosmic), the number of
//  shower-like PFPs (nShw) and the second-nearest shower (dist2, e2).
//  Muon attachment (2026-10-31): the chain's muon (IslandProcess.C: track score > 0.8, within
//  20 cm of the vertex, length > 10 cm, LLR > 0.7, longest); for the chain photon its distance
//  from the muon axis (phMuPerp) and position along it (phMuAlong, 0 = muon start), muLen; and the
//  photon re-picked as the nearest shower-like PFP NOT attached to the muon (attached: within
//  MU_PERP cm of the axis and between MU_ALONG_LO and length + MU_ALONG_HI along it -- delta rays
//  and the Michel electron): repIdx (-1 = none left), repDist, repE, repPdg, repOrigin.
//  Run:  root -l -b -q 'PhotonPoint.C+("lists/bkg_overlay.list","isl_bkg.root","<sel>","pp_bkg.root")'
// ============================================================================
#include "HypInput.h"
#include "FluxWeights.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <cstdio>

void PhotonPoint(const char* list, const char* cache, const char* sel, const char* out)
{
  // events to process, with the cache's phDist
  std::map<unsigned long long, double> want;
  { TFile c(cache); TTree* t = (TTree*)c.Get("isl");
    t->SetEstimate(t->GetEntries() + 1);
    t->Draw("run:sub:evt:phDist", sel, "goff para");
    for (Long64_t i = 0; i < t->GetSelectedRows(); ++i)
      want[rseKey((unsigned)t->GetVal(0)[i], (unsigned)t->GetVal(1)[i], (unsigned)t->GetVal(2)[i])] = t->GetVal(3)[i];
  }
  printf("  %zu events selected from %s\n", want.size(), cache);

  TFile fo(out, "RECREATE");
  TTree T("pp", "photon pointing for the Lambda-style chain");
  int run, sub, evt, altIdx, altSame;
  double distChk, phCos, phImpact, altDist, altE, altCos;
  int nhU, nhV, nhY, phPdg, phOrigin, nShw;
  double shLen, shOpen, molAvg, molRms, dedx0, dedx1, dedx2, e0, e1, e2p, dist2, e2nd;
  int repIdx, repPdg, repOrigin, phMuAtt; double phMuPerp, phMuAlong, muLen, repDist, repE;
  const double MU_PERP = 5.0, MU_ALONG_LO = -5.0, MU_ALONG_HI = 10.0;
  T.Branch("run", &run, "run/I"); T.Branch("sub", &sub, "sub/I"); T.Branch("evt", &evt, "evt/I");
  T.Branch("distChk", &distChk, "distChk/D"); T.Branch("phCos", &phCos, "phCos/D"); T.Branch("phImpact", &phImpact, "phImpact/D");
  T.Branch("altIdx", &altIdx, "altIdx/I"); T.Branch("altSame", &altSame, "altSame/I");
  T.Branch("altDist", &altDist, "altDist/D"); T.Branch("altE", &altE, "altE/D"); T.Branch("altCos", &altCos, "altCos/D");
  T.Branch("nhU", &nhU, "nhU/I"); T.Branch("nhV", &nhV, "nhV/I"); T.Branch("nhY", &nhY, "nhY/I");
  T.Branch("shLen", &shLen, "shLen/D"); T.Branch("shOpen", &shOpen, "shOpen/D");
  T.Branch("molAvg", &molAvg, "molAvg/D"); T.Branch("molRms", &molRms, "molRms/D");
  T.Branch("dedx0", &dedx0, "dedx0/D"); T.Branch("dedx1", &dedx1, "dedx1/D"); T.Branch("dedx2", &dedx2, "dedx2/D");
  T.Branch("e0", &e0, "e0/D"); T.Branch("e1", &e1, "e1/D"); T.Branch("e2", &e2p, "e2/D");
  T.Branch("phPdg", &phPdg, "phPdg/I"); T.Branch("phOrigin", &phOrigin, "phOrigin/I");
  T.Branch("nShw", &nShw, "nShw/I"); T.Branch("dist2", &dist2, "dist2/D"); T.Branch("e2nd", &e2nd, "e2nd/D");
  T.Branch("phMuPerp", &phMuPerp, "phMuPerp/D"); T.Branch("phMuAlong", &phMuAlong, "phMuAlong/D"); T.Branch("muLen", &muLen, "muLen/D");
  T.Branch("phMuAtt", &phMuAtt, "phMuAtt/I");
  T.Branch("repIdx", &repIdx, "repIdx/I"); T.Branch("repDist", &repDist, "repDist/D"); T.Branch("repE", &repE, "repE/D");
  T.Branch("repPdg", &repPdg, "repPdg/I"); T.Branch("repOrigin", &repOrigin, "repOrigin/I");

  std::ifstream lf(list); std::string fn; long nfile = 0, nfound = 0;
  std::set<unsigned long long> done;
  while (std::getline(lf, fn) && done.size() < want.size()) {
    if (fn.empty()) continue;
    TFile* f = TFile::Open(fn.c_str()); if (!f || f->IsZombie()) { delete f; continue; }
    TTree* t = (TTree*)f->Get("ana/OutputTree"); if (!t) { delete f; continue; }
    nfile++;
    // cheap pass: which entries are wanted
    UInt_t r, s, e; std::vector<Long64_t> hits;
    t->SetBranchStatus("*", 0); t->SetBranchStatus("run", 1); t->SetBranchStatus("subrun", 1); t->SetBranchStatus("event", 1);
    t->SetBranchAddress("run", &r); t->SetBranchAddress("subrun", &s); t->SetBranchAddress("event", &e);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); const auto k = rseKey(r, s, e);
      if (want.count(k) && !done.count(k)) hits.push_back(i);
    }
    t->ResetBranchAddresses();
    if (hits.empty()) { delete f; continue; }
    HypInput H; if (!H.attach(t)) { delete f; continue; }
    std::vector<double> *d0 = nullptr, *d1 = nullptr, *en0 = nullptr, *en1 = nullptr;
    if (t->GetBranch("shr_dedx_plane0"))   t->SetBranchAddress("shr_dedx_plane0", &d0);
    if (t->GetBranch("shr_dedx_plane1"))   t->SetBranchAddress("shr_dedx_plane1", &d1);
    if (t->GetBranch("shr_energy_plane0")) t->SetBranchAddress("shr_energy_plane0", &en0);
    if (t->GetBranch("shr_energy_plane1")) t->SetBranchAddress("shr_energy_plane1", &en1);
    for (Long64_t i : hits) {
      t->GetEntry(i); H.convert();
      const auto k = rseKey(H.h_run, H.h_subrun, H.h_event);
      done.insert(k); nfound++;
      run = H.run; sub = H.sub; evt = H.evt;
      const int n = (int)H.trk_score.size();
      const TVector3 vtx(H.vtx_x, H.vtx_y, H.vtx_z);
      int s1 = -1; double bs = 1e9;
      std::vector<double> dv(n, 1e9);
      for (int j = 0; j < n; ++j) {
        if (H.trk_score[j] > 0.5) continue;
        dv[j] = (TVector3(H.shr_sx[j], H.shr_sy[j], H.shr_sz[j]) - vtx).Mag();
        if (dv[j] < bs) { bs = dv[j]; s1 = j; }
      }
      auto point = [&](int j, double& c, double& b) {
        const TVector3 st(H.shr_sx[j], H.shr_sy[j], H.shr_sz[j]);
        TVector3 d(H.shr_px[j], H.shr_py[j], H.shr_pz[j]); c = -2; b = -1;
        if (d.Mag() <= 0) return; d = d.Unit();
        const TVector3 rv = st - vtx;
        if (rv.Mag() > 0) c = d.Dot(rv.Unit());
        b = (vtx - st).Cross(d).Mag();
      };
      distChk = -1; phCos = -2; phImpact = -1;
      if (s1 >= 0) { distChk = std::fabs(dv[s1] - want[k]); point(s1, phCos, phImpact); }
      altIdx = -1; altCos = -2; altDist = -1; altE = -1;
      for (int j = 0; j < n; ++j) {
        if (H.trk_score[j] > 0.5 || dv[j] > 100) continue;
        double c, b; point(j, c, b);
        if (c > altCos) { altCos = c; altIdx = j; }
      }
      if (altIdx >= 0) { altDist = dv[altIdx]; altE = H.shr_energy_y[altIdx] / 1000.0; }
      altSame = (altIdx == s1);
      auto at = [](const std::vector<double>* v, int j) { return (v && j >= 0 && j < (int)v->size()) ? v->at(j) : -999.; };
      auto ati = [](const std::vector<int>* v, int j) { return (v && j >= 0 && j < (int)v->size()) ? v->at(j) : -999; };
      nhU = ati(H.h_hitsU, s1); nhV = ati(H.h_hitsV, s1); nhY = ati(H.h_hitsY, s1);
      shLen = at(H.h_shr_length, s1); shOpen = at(H.h_shr_open_angle, s1);
      molAvg = at(H.h_shr_moliere_avg, s1); molRms = at(H.h_shr_moliere_rms, s1);
      dedx0 = at(d0, s1); dedx1 = at(d1, s1); dedx2 = at(H.h_shr_dedx_p2, s1);
      e0 = at(en0, s1); e1 = at(en1, s1); e2p = at(H.h_shr_energy_p2, s1);
      phPdg = ati(H.h_pfp_true_pdg, s1); phOrigin = ati(H.h_pfp_true_origin, s1);
      nShw = 0; dist2 = -1; e2nd = -1; double b2 = 1e9;
      for (int j = 0; j < n; ++j) {
        if (H.trk_score[j] > 0.5) continue; nShw++;
        if (j != s1 && dv[j] < b2) { b2 = dv[j]; dist2 = dv[j]; e2nd = H.shr_energy_y[j] / 1000.0; }
      }
      // muon as IslandProcess.C, and attachment of each shower-like PFP to it
      int mu = -1;
      for (int j = 0; j < n; ++j) {
        if (H.trk_score[j] <= 0.8) continue;
        const double dm = (TVector3(H.trk_sx[j], H.trk_sy[j], H.trk_sz[j]) - vtx).Mag();
        if (dm < 20.0 && H.trk_len[j] > 10.0 && H.trk_llr[j] > 0.7 && (mu < 0 || H.trk_len[j] > H.trk_len[mu])) mu = j;
      }
      auto rel = [&](int j, double& perp, double& along) {
        perp = -1; along = -1e9; if (mu < 0 || !H.h_trk_end_x) return;
        const TVector3 ms(H.trk_sx[mu], H.trk_sy[mu], H.trk_sz[mu]);
        const TVector3 me(H.h_trk_end_x->at(mu), H.h_trk_end_y->at(mu), H.h_trk_end_z->at(mu));
        const TVector3 ps(H.shr_sx[j], H.shr_sy[j], H.shr_sz[j]);
        const TVector3 ax = me - ms; if (ax.Mag() <= 0) return; const TVector3 u = ax.Unit();
        along = (ps - ms).Dot(u); perp = (ps - ms).Cross(u).Mag();
      };
      auto attached = [&](int j) { double pe, al; rel(j, pe, al); if (pe < 0) return false;
        const double L = (mu >= 0 && H.h_trk_end_x) ? (TVector3(H.h_trk_end_x->at(mu), H.h_trk_end_y->at(mu), H.h_trk_end_z->at(mu))
                                                        - TVector3(H.trk_sx[mu], H.trk_sy[mu], H.trk_sz[mu])).Mag() : 0;
        return pe < MU_PERP && al > MU_ALONG_LO && al < L + MU_ALONG_HI; };
      muLen = (mu >= 0 && H.h_trk_end_x) ? (TVector3(H.h_trk_end_x->at(mu), H.h_trk_end_y->at(mu), H.h_trk_end_z->at(mu))
                                            - TVector3(H.trk_sx[mu], H.trk_sy[mu], H.trk_sz[mu])).Mag() : -1;
      phMuPerp = -1; phMuAlong = -1e9; phMuAtt = 0;
      if (s1 >= 0) { rel(s1, phMuPerp, phMuAlong); phMuAtt = attached(s1); }
      repIdx = -1; double br = 1e9;
      for (int j = 0; j < n; ++j) {
        if (H.trk_score[j] > 0.5 || attached(j)) continue;
        if (dv[j] < br) { br = dv[j]; repIdx = j; }
      }
      repDist = repIdx >= 0 ? dv[repIdx] : -1; repE = repIdx >= 0 ? H.shr_energy_y[repIdx] / 1000.0 : -1;
      repPdg = ati(H.h_pfp_true_pdg, repIdx); repOrigin = ati(H.h_pfp_true_origin, repIdx);
      T.Fill();
    }
    t->ResetBranchAddresses();
    delete f;
  }
  printf("  %ld files read, %ld of %zu events found -> %s\n", nfile, nfound, want.size(), out);
  fo.Write();
}
