#ifndef HypInput_h
#define HypInput_h
// ============================================================================
//  HypInput -- read HyperonProduction ("hypana") ntuples through the same
//  variable names the PeLEE-based Sigma0 selection already uses.
//
//  WHY: the PeLEE truth record stops at the primary vertex. mc_pdg lists the
//  Sigma0 and the muon but never the Sigma0 -> Lambda gamma -> p pi- gamma
//  chain, so zero of 17721 signal events have a backtracked pi- PFP and the
//  per-candidate truth columns are unusable (measured: pion candidate matches
//  a true pi- 0.0% of the time; 60% carry no match at all and 36.5% leak the
//  muon's truth). HyperonProduction writes the full decay chain plus per-PFP
//  truth matching, so candidate purity, the pion candidate and the photon can
//  finally be studied against truth.
//
//  Layout is identical in both files: every trk_*/shr_*/pfp_* vector is
//  PER-PFP and all have the same length (verified on 396 events), so this is a
//  branch-name remap with a double->float copy, not a structural rewrite.
//
//  Usage:
//      HypInput hyp;
//      if (hyp.attach(tree)) { ... hyp.convert(); ... }
//  then point the selection's existing pointers at hyp's buffers.
// ============================================================================
#include <vector>
#include <string>
#include <cmath>
#include "TTree.h"

struct HypInput {
  // ── hypana-side branch pointers ──────────────────────────────────────────
  std::vector<double> *h_trk_length=nullptr, *h_trk_llrpid=nullptr,
    *h_trk_start_x=nullptr, *h_trk_start_y=nullptr, *h_trk_start_z=nullptr,
    *h_trk_end_x=nullptr, *h_trk_end_y=nullptr, *h_trk_end_z=nullptr,
    *h_trk_dir_x=nullptr, *h_trk_dir_y=nullptr, *h_trk_dir_z=nullptr,
    *h_trk_pid_chipr=nullptr, *h_trk_pid_chipi=nullptr,
    *h_trk_energy_proton=nullptr, *h_trk_range_muon_mom=nullptr,
    *h_trk_mcs_muon_mom=nullptr, *h_trk_rms_deflection=nullptr,
    *h_trk_bragg_p_fwd=nullptr, *h_trk_bragg_p_bwd=nullptr,
    *h_trk_bragg_mu_fwd=nullptr, *h_trk_bragg_mu_bwd=nullptr,
    *h_trk_three_plane_dedx=nullptr,
    *h_shr_start_x=nullptr, *h_shr_start_y=nullptr, *h_shr_start_z=nullptr,
    *h_shr_dir_x=nullptr, *h_shr_dir_y=nullptr, *h_shr_dir_z=nullptr,
    *h_shr_energy_p2=nullptr, *h_shr_dedx_p2=nullptr,
    *h_shr_moliere_avg=nullptr, *h_shr_moliere_rms=nullptr,
    *h_shr_length=nullptr, *h_shr_open_angle=nullptr,
    *h_pfp_trk_shr_score=nullptr, *h_pfp_purity=nullptr, *h_pfp_completeness=nullptr,
    *h_mc_primary_mom=nullptr;
  std::vector<int> *h_pfp_pdg=nullptr, *h_pfp_true_pdg=nullptr,
    *h_pfp_true_origin=nullptr, *h_mc_primary_pdg=nullptr,
    *h_hitsU=nullptr, *h_hitsV=nullptr, *h_hitsY=nullptr;
  std::string *h_ccnc=nullptr, *h_mode=nullptr;
  Double_t h_vx=0,h_vy=0,h_vz=0, h_nx=0,h_ny=0,h_nz=0, h_lepton_mom=0;
  Int_t    h_nu_pdg=0, h_hyperon_pdg=0;
  UInt_t   h_run=0,h_subrun=0,h_event=0;
  Bool_t   h_vtx_inFV=false;

  // ── PeLEE-named buffers the selection body reads ─────────────────────────
  std::vector<float> trk_score, trk_len, trk_llr,
    trk_sx, trk_sy, trk_sz, trk_dx, trk_dy, trk_dz,
    trk_calo_y, trk_chipr, trk_chipi, trk_eproton, trk_rmumom,
    shr_sx, shr_sy, shr_sz, shr_dist, shr_energy_y,
    shr_px, shr_py, shr_pz, shr_dedx_y, shr_moliere,
    bkt_overlay_purity, mc_px, mc_py, mc_pz;
  std::vector<int> pfnplanehits_Y, bkt_pdg, mc_pdg_buf;
  std::vector<unsigned int> pfnhits;

  int   run=0, sub=0, evt=0, nu_pdg=0, ccnc=0, interaction=0;
  float nu_e=-1.f;                   // NOT in hypana; flagged as -1
  float vtx_x=0, vtx_y=0, vtx_z=0, true_vtx_x=0, true_vtx_y=0, true_vtx_z=0;

  bool attach(TTree* t){
    if(!t || !t->GetBranch("pfp_true_pdg") || !t->GetBranch("trk_length")) return false;
    t->SetBranchStatus("*",1);
    auto BV=[&](const char* n, void* p){ if(t->GetBranch(n)) t->SetBranchAddress(n,p); };
    BV("trk_length",&h_trk_length);           BV("trk_llrpid",&h_trk_llrpid);
    BV("trk_start_x",&h_trk_start_x); BV("trk_start_y",&h_trk_start_y); BV("trk_start_z",&h_trk_start_z);
    BV("trk_end_x",&h_trk_end_x);     BV("trk_end_y",&h_trk_end_y);     BV("trk_end_z",&h_trk_end_z);
    BV("trk_dir_x",&h_trk_dir_x);     BV("trk_dir_y",&h_trk_dir_y);     BV("trk_dir_z",&h_trk_dir_z);
    BV("trk_pid_chipr",&h_trk_pid_chipr); BV("trk_pid_chipi",&h_trk_pid_chipi);
    BV("trk_energy_proton",&h_trk_energy_proton); BV("trk_range_muon_mom",&h_trk_range_muon_mom);
    BV("trk_mcs_muon_mom",&h_trk_mcs_muon_mom);   BV("trk_rms_deflection",&h_trk_rms_deflection);
    BV("trk_bragg_p_fwd",&h_trk_bragg_p_fwd);     BV("trk_bragg_p_bwd",&h_trk_bragg_p_bwd);
    BV("trk_bragg_mu_fwd",&h_trk_bragg_mu_fwd);   BV("trk_bragg_mu_bwd",&h_trk_bragg_mu_bwd);
    BV("trk_three_plane_dedx",&h_trk_three_plane_dedx);
    BV("shr_start_x",&h_shr_start_x); BV("shr_start_y",&h_shr_start_y); BV("shr_start_z",&h_shr_start_z);
    BV("shr_dir_x",&h_shr_dir_x);     BV("shr_dir_y",&h_shr_dir_y);     BV("shr_dir_z",&h_shr_dir_z);
    BV("shr_energy_plane2",&h_shr_energy_p2); BV("shr_dedx_plane2",&h_shr_dedx_p2);
    BV("shr_moliere_avg",&h_shr_moliere_avg); BV("shr_moliere_rms",&h_shr_moliere_rms);
    BV("shr_length",&h_shr_length);           BV("shr_open_angle",&h_shr_open_angle);
    BV("pfp_trk_shr_score",&h_pfp_trk_shr_score);
    BV("pfp_purity",&h_pfp_purity); BV("pfp_completeness",&h_pfp_completeness);
    BV("pfp_pdg",&h_pfp_pdg); BV("pfp_true_pdg",&h_pfp_true_pdg);
    BV("pfp_true_origin",&h_pfp_true_origin);
    BV("mc_primary_pdg",&h_mc_primary_pdg); BV("mc_primary_mom",&h_mc_primary_mom);
    BV("pfpnplanehits_U",&h_hitsU); BV("pfpnplanehits_V",&h_hitsV); BV("pfpnplanehits_Y",&h_hitsY);
    BV("mc_ccnc",&h_ccnc); BV("mc_mode",&h_mode);
    t->SetBranchAddress("reco_primary_vtx_x",&h_vx);
    t->SetBranchAddress("reco_primary_vtx_y",&h_vy);
    t->SetBranchAddress("reco_primary_vtx_z",&h_vz);
    t->SetBranchAddress("mc_nu_pos_x",&h_nx);
    t->SetBranchAddress("mc_nu_pos_y",&h_ny);
    t->SetBranchAddress("mc_nu_pos_z",&h_nz);
    t->SetBranchAddress("mc_nu_pdg",&h_nu_pdg);
    t->SetBranchAddress("mc_hyperon_pdg",&h_hyperon_pdg);
    t->SetBranchAddress("mc_lepton_mom",&h_lepton_mom);
    t->SetBranchAddress("run",&h_run); t->SetBranchAddress("subrun",&h_subrun);
    t->SetBranchAddress("event",&h_event);
    if(t->GetBranch("reco_primary_vtx_inFV")) t->SetBranchAddress("reco_primary_vtx_inFV",&h_vtx_inFV);
    return true;
  }

  // Per-event conversion. Call after GetEntry().
  void convert(){
    const size_t n = h_pfp_pdg ? h_pfp_pdg->size() : 0;
    auto F=[&](std::vector<float>& d, std::vector<double>* s){
      d.assign(n, -999.f);
      if(s) for(size_t k=0;k<n&&k<s->size();++k) d[k]=(float)s->at(k); };
    F(trk_score,h_pfp_trk_shr_score); F(trk_len,h_trk_length); F(trk_llr,h_trk_llrpid);
    F(trk_sx,h_trk_start_x); F(trk_sy,h_trk_start_y); F(trk_sz,h_trk_start_z);
    F(trk_dx,h_trk_dir_x);   F(trk_dy,h_trk_dir_y);   F(trk_dz,h_trk_dir_z);
    F(trk_chipr,h_trk_pid_chipr); F(trk_chipi,h_trk_pid_chipi);
    // HyperonProduction's trk_energy_proton is the proton MOMENTUM (GeV/c), not
    // the kinetic energy the PeLEE branch of the same name carries. Verified
    // against pfp_true_p on true protons: 0.594 vs 0.609, 0.516 vs 0.519, 0.653
    // vs 0.674. Treating it as KE inflates the proton momentum ~2.3x (median 930
    // vs 410 MeV/c at 4.8 cm range). Convert to KE here so trk_eproton means the
    // same thing in both formats and every consumer stays unchanged.
    F(trk_eproton,h_trk_energy_proton);
    for(size_t k=0;k<n;++k){ const double pm=trk_eproton[k];
      trk_eproton[k] = (pm>0) ? (float)(std::sqrt(pm*pm+0.938272*0.938272)-0.938272) : trk_eproton[k]; }
    F(trk_rmumom,h_trk_range_muon_mom);
    F(shr_sx,h_shr_start_x); F(shr_sy,h_shr_start_y); F(shr_sz,h_shr_start_z);
    F(shr_px,h_shr_dir_x);   F(shr_py,h_shr_dir_y);   F(shr_pz,h_shr_dir_z);
    F(shr_energy_y,h_shr_energy_p2); F(shr_dedx_y,h_shr_dedx_p2);
    F(shr_moliere,h_shr_moliere_avg);
    // trk_calo_energy_y has no hypana equivalent; flag it rather than fake it.
    trk_calo_y.assign(n,-1.f);
    // shower distance to the primary vertex is derived, not stored.
    shr_dist.assign(n,-1.f);
    for(size_t k=0;k<n;++k){
      if(h_shr_start_x && k<h_shr_start_x->size()){
        const double ddx=h_shr_start_x->at(k)-h_vx, ddy=h_shr_start_y->at(k)-h_vy,
                     ddz=h_shr_start_z->at(k)-h_vz;
        shr_dist[k]=(float)std::sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
      }
    }
    // per-PFP truth: now genuinely usable (pi-, gamma and hyperons all appear).
    bkt_pdg.assign(n,0);
    if(h_pfp_true_pdg) for(size_t k=0;k<n&&k<h_pfp_true_pdg->size();++k) bkt_pdg[k]=h_pfp_true_pdg->at(k);
    // Cosmic tag, PeLEE-style (1 = overlay/cosmic). pfp_true_origin is NOT simb::Origin:
    // HyperonProduction's util::OriginType is 1 neutrino primary, 2 Lambda daughter,
    // 3 other, 5 Sigma0 daughter, and -1 when the PFP has no MC match -- which is what the
    // overlaid cosmic DATA gives. (Until 2026-10-17 this tested ==2 and so flagged the Lambda
    // daughters; no hypana macro read it.)
    bkt_overlay_purity.assign(n,0.f);
    if(h_pfp_true_origin) for(size_t k=0;k<n&&k<h_pfp_true_origin->size();++k)
      bkt_overlay_purity[k] = (h_pfp_true_origin->at(k)==-1) ? 1.f : 0.f;
    pfnplanehits_Y.assign(n,0); pfnhits.assign(n,0u);
    for(size_t k=0;k<n;++k){
      const int u=(h_hitsU&&k<h_hitsU->size())?h_hitsU->at(k):0;
      const int v=(h_hitsV&&k<h_hitsV->size())?h_hitsV->at(k):0;
      const int y=(h_hitsY&&k<h_hitsY->size())?h_hitsY->at(k):0;
      pfnplanehits_Y[k]=y; pfnhits[k]=(unsigned)(u+v+y);
    }
    // truth primaries. The selection only ever takes |p| from mc_px/py/pz, so
    // packing the stored magnitude into x reproduces it exactly.
    const size_t m = h_mc_primary_pdg ? h_mc_primary_pdg->size() : 0;
    mc_pdg_buf.assign(m,0); mc_px.assign(m,0.f); mc_py.assign(m,0.f); mc_pz.assign(m,0.f);
    for(size_t k=0;k<m;++k){
      mc_pdg_buf[k]=h_mc_primary_pdg->at(k);
      if(h_mc_primary_mom && k<h_mc_primary_mom->size()) mc_px[k]=(float)h_mc_primary_mom->at(k);
    }
    run=(int)h_run; sub=(int)h_subrun; evt=(int)h_event;
    nu_pdg=h_nu_pdg;
    ccnc = (h_ccnc && *h_ccnc=="NC") ? 1 : 0;
    interaction = !h_mode ? -1 : (*h_mode=="QEL"?0 : *h_mode=="RES"?1 : *h_mode=="DIS"?2
                               : *h_mode=="COH"?3 : *h_mode=="MEC"?10 : -1);
    nu_e = -1.f;                                  // absent from hypana
    vtx_x=(float)h_vx; vtx_y=(float)h_vy; vtx_z=(float)h_vz;
    true_vtx_x=(float)h_nx; true_vtx_y=(float)h_ny; true_vtx_z=(float)h_nz;
  }
};
#endif
