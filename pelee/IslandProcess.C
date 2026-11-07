// ============================================================================
//  IslandProcess.C -- ONE expensive pass of the island-finding (connectedness)
//  algorithm over a hypana sample, caching everything the later studies need.
//
//  WHY a caching pass: the island test reads the raw wire windows, which is what
//  makes every hypana loop slow. Running it once per sample and writing the
//  per-event result lets the cut placement, the flight dependence, the POT-scaled
//  S/B and the BDT integration all be studied from a small tree without touching
//  the wire branches again. Analysis of the output: IslandAna.C.
//
//  Blind analysis: signal, overlay, EXT and (since 2026-11-06) dirt; never beam-on.
//
//  The chain is the Lambda-style one (SigmaLambdaSel.C) but NOTHING after the
//  muon is a hard cut here: every event reaching MuonCandidate is written, with
//  the stage it would have reached recorded in `stage`, so any later cut order
//  can be evaluated offline. The island test is evaluated for every event that
//  has a proton/pion pair at all -- the best precut-passing pair, or, when none
//  passes, the best pair regardless (precutOK=0) -- so its value on the events
//  the precuts would remove is also known.
//
//  Run:  root -l -b -q 'IslandProcess.C("lists/sig_genie_qe.list",0,true,"isl_sig_genie_qe.root",9.37035e-4)'
//        (the previous hyperon sample: "lists/sig_genie.list" -> "isl_sig_genie.root", 1.9284e-2)
//        root -l -b -q 'IslandProcess.C("lists/bkg_overlay.list",0,false,"isl_bkg.root",1.6352e-1)'
//  The last argument is the POT scale to data (1.1526e21 POT), see HypCutFlow.C.
//  For lists/bkg_overlay.list it became 1.6536e-1 on 2026-09-24, when the overlay
//  extension took that list to 6.9701e21 POT (was 4.1468e-1 at 2.7795e21), and
//  1.6352e-1 on 2026-10-16, when the residue pass took it to 7.0488e21.
//  NuWro signal (lists/sig_nuwro.list, 1.3394e24 POT -> 8.6055e-4) is on the OLD NuMI
//  flux: pass the per-event table as fluxWeights (see FluxWeights.h) and `weight` becomes
//  potScale x flux weight; the table is applied to signal passes only.
//        root -l -b -q 'IslandProcess.C("lists/sig_nuwro.list",0,true,"isl_sig_nuwro.root",8.6055e-4,
//                       "/home/t2k/nowak/MicroBooNE/production/nuwro_prod/flux_rw/nuwro_rhc_run3b_reco2_fluxweights.tsv")'
//  Old-flux background samples (dirt) take the table with bkgFluxWeights=true:
//        root -l -b -q 'IslandProcess.C("lists/bkg_dirt_run3b.list",0,false,"isl_dirt.root",<scale>,"fw_dirt_run3b.tsv",0,true)'
// ============================================================================
#include "HypInput.h"
#include "SigmaFV.h"
#include "Island.h"
#include "LambdaStyle.h"
#include "FluxWeights.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <fstream>
#include <vector>
#include <cstdio>

using namespace LambdaStyle;

static const double MUON_TRK_SCORE = 0.8;
static const double MUON_TRK_DIST  = 20.0;
static const double MUON_TRK_LEN   = 10.0;
static const double MUON_LLR       = 0.7;
static const int    NSHOWERS_MAX   = 4;
static const int    NTRACKS_MAX    = 5;
static const double TRK_SHR_SPLIT  = 0.5;
static const double TOPO_MIN       = 0.20;
static const double LAM_MASS_WIN   = 0.040;
static const double COS_ALPHA_MIN  = 0.90;
static const double DECAY_LEN_MIN  = 0.5;

// stage indices written to the tree
enum { ST_ALL=0, ST_SLICE, ST_FV, ST_COSMIC, ST_NSHOW, ST_MUON, ST_PAIR, ST_MASS, ST_ALPHA, NST };
static const char* STNAME[NST] = { "AllEvents","SliceID","FiducialV","CosmicVeto","NShowers",
                                   "MuonCandidate","PairCandidate","LambdaMass","PointingAlpha" };

void IslandProcess(const char* list = "lists/sig_genie_qe.list", int maxFiles = 0,
                   bool signalOnly = true, const char* outName = "isl_sig_genie_qe.root",
                   double potScale = 1.0, const char* fluxWeights = "",
                   int showerMode = 0,   // 0: 0 < nShowers < 4 (standard); 1: ONLY events with no
                                         // shower (2026-11-01: the Sigma0 photon is unreconstructed in
                                         // 20% of signal at this stage, and the chain uses the photon only
                                         // for the gap cut). Combine a mode-1 cache with the standard one.
                   bool bkgFluxWeights = false)  // apply fluxWeights to a background pass too: for a
                                         // sample whose keys are unique on their own (the old-flux dirt,
                                         // fw_dirt_run3b.tsv; the hyperon detvar CV)
{
  const std::map<unsigned long long,double> FW = (signalOnly || bkgFluxWeights) ? loadFluxWeights(fluxWeights)
                                                            : std::map<unsigned long long,double>();
  double sumW=0; long nW=0, nWmiss=0;
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf = 0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0 || nf<maxFiles)){ ch.Add(l.c_str()); nf++; }

  double topo = 0; ch.SetBranchAddress("topological_score", &topo);
  HypInput H; if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }
  // extras AFTER attach() (TChain bind-order rule)
  int fmSlice = -1;
  const bool haveSlice = ch.GetBranch("flash_match_nu_slice_ID") != nullptr;
  if(haveSlice) ch.SetBranchAddress("flash_match_nu_slice_ID", &fmSlice);
  std::vector<double>* d3 = nullptr;
  const bool haveDedx = ch.GetBranch("trk_three_plane_dedx") != nullptr;
  if(haveDedx) ch.SetBranchAddress("trk_three_plane_dedx", &d3);
  std::vector<std::vector<float>>* win[3] = {0,0,0};
  std::vector<int>* vw = nullptr;
  for(int p=0;p<3;++p) ch.SetBranchAddress(Form("ct_test_window_plane%d",p), &win[p]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires", &vw);
  double lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0;
  auto B=[&](const char* n, void* q){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,q); };
  // v2 additions: muon end point (Michel tag), true lepton charge, true photon momentum
  std::vector<double> *tex=nullptr,*tey=nullptr,*tez=nullptr;
  B("trk_end_x",&tex); B("trk_end_y",&tey); B("trk_end_z",&tez);
  int lepPdg=0; double truePhMom=-1;
  B("mc_lepton_pdg",&lepPdg); B("mc_sigmazero_photon_mom",&truePhMom);
  double trueQ2=-1; B("mc_nu_q2",&trueQ2);   // v3: truth for background-model checks (hyperon samples)
  B("mc_sigmazero_lambda_start_x",&lsx); B("mc_sigmazero_lambda_start_y",&lsy);
  B("mc_sigmazero_lambda_start_z",&lsz); B("mc_sigmazero_lambda_end_x",&lex);
  B("mc_sigmazero_lambda_end_y",&ley);   B("mc_sigmazero_lambda_end_z",&lez);

  // ── output ────────────────────────────────────────────────────────────────
  TFile fout(outName, "RECREATE");
  TTree T("isl", "island-finding pass, Lambda-style chain");
  int    isSig=0, mode=-1, run=0, sub=0, evt=0, stage=0;
  double weight=potScale;
  int    nTr=0, nSh=0, hasPair=0, precutOK=0, hasPhoton=0;
  double mass=-1, decayLen=-1, cosAlpha=-2, vfitChi2=-1, sep=-1, llrP=-9, llrPi=-9,
         pMom=-1, iMom=-1, shrE=-1, trueFlight=-1, pairScore=-1e30;
  int    islPass[3]={0,0,0}, islOK[3]={0,0,0}, islNPlanes=0, islAny=0, islSeeded=0;
  int    muTrue=0, prTrue=0, piTrue=0, s1True=0;
  // v2: photon usage and Michel tag
  double sigmaMass=-1, cosPhLam=-2, phDist=-1, phDedx=-1, phMoliere=-1, truePhE=-1;
  double michelDist=-1, michelE=-1; int michelN5=0, muContained=0, lepTrue=0;
  T.Branch("isSig",&isSig,"isSig/I");       T.Branch("weight",&weight,"weight/D");
  T.Branch("mode",&mode,"mode/I");          T.Branch("run",&run,"run/I");
  T.Branch("sub",&sub,"sub/I");             T.Branch("evt",&evt,"evt/I");
  T.Branch("stage",&stage,"stage/I");
  T.Branch("nTr",&nTr,"nTr/I");             T.Branch("nSh",&nSh,"nSh/I");
  T.Branch("hasPair",&hasPair,"hasPair/I"); T.Branch("precutOK",&precutOK,"precutOK/I");
  T.Branch("hasPhoton",&hasPhoton,"hasPhoton/I");
  T.Branch("mass",&mass,"mass/D");          T.Branch("decayLen",&decayLen,"decayLen/D");
  T.Branch("cosAlpha",&cosAlpha,"cosAlpha/D"); T.Branch("vfitChi2",&vfitChi2,"vfitChi2/D");
  T.Branch("sep",&sep,"sep/D");             T.Branch("llrP",&llrP,"llrP/D");
  T.Branch("llrPi",&llrPi,"llrPi/D");       T.Branch("pMom",&pMom,"pMom/D");
  T.Branch("iMom",&iMom,"iMom/D");          T.Branch("shrE",&shrE,"shrE/D");
  T.Branch("pairScore",&pairScore,"pairScore/D");
  T.Branch("trueFlight",&trueFlight,"trueFlight/D");
  T.Branch("islPass",islPass,"islPass[3]/I"); T.Branch("islOK",islOK,"islOK[3]/I");
  T.Branch("islNPlanes",&islNPlanes,"islNPlanes/I"); T.Branch("islAny",&islAny,"islAny/I");
  T.Branch("islSeeded",&islSeeded,"islSeeded/I");
  T.Branch("muTrue",&muTrue,"muTrue/I");    T.Branch("prTrue",&prTrue,"prTrue/I");
  T.Branch("piTrue",&piTrue,"piTrue/I");    T.Branch("s1True",&s1True,"s1True/I");
  T.Branch("sigmaMass",&sigmaMass,"sigmaMass/D"); T.Branch("cosPhLam",&cosPhLam,"cosPhLam/D");
  T.Branch("phDist",&phDist,"phDist/D");    T.Branch("phDedx",&phDedx,"phDedx/D");
  T.Branch("phMoliere",&phMoliere,"phMoliere/D"); T.Branch("truePhE",&truePhE,"truePhE/D");
  T.Branch("michelDist",&michelDist,"michelDist/D"); T.Branch("michelE",&michelE,"michelE/D");
  T.Branch("michelN5",&michelN5,"michelN5/I"); T.Branch("muContained",&muContained,"muContained/I");
  T.Branch("lepTrue",&lepTrue,"lepTrue/I");
  int hypPdg=0, nuPdg=0, ccnc=0; double q2=-1;
  T.Branch("hypPdg",&hypPdg,"hypPdg/I");    T.Branch("nuPdg",&nuPdg,"nuPdg/I");
  T.Branch("ccnc",&ccnc,"ccnc/I");          T.Branch("q2",&q2,"q2/D");

  long nCut[NST] = {0}; long nIsland = 0;
  const PairCuts pc;
  const Long64_t N = ch.GetEntries();
  printf("\n  %s : %d files, %lld events, POT scale %.4e\n",
         signalOnly?"SIGNAL":"BACKGROUND", nf, N, potScale);
  fflush(stdout);

  for(Long64_t ie = 0; ie < N; ++ie){
    ch.GetEntry(ie); H.convert();
    isSig = (H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
    if(signalOnly != (bool)isSig) continue;
    if(ie % 20000 == 0){ printf("  ... %lld / %lld\n", ie, N); fflush(stdout); }
    weight = potScale;
    if(!FW.empty()){
      auto it=FW.find(rseKey(H.h_run,H.h_subrun,H.h_event));
      if(it!=FW.end()){ weight*=it->second; sumW+=it->second; ++nW; } else ++nWmiss;
    }

    stage = ST_ALL; nCut[ST_ALL]++;
    if(haveSlice && fmSlice < 0) continue;                        stage=ST_SLICE;  nCut[ST_SLICE]++;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;   stage=ST_FV;     nCut[ST_FV]++;
    if(topo < TOPO_MIN) continue;                                 stage=ST_COSMIC; nCut[ST_COSMIC]++;

    const int n = (int)H.trk_score.size(); if(n < 1) continue;
    const TVector3 vtx(H.vtx_x, H.vtx_y, H.vtx_z);
    auto isTrk = [&](int i){ return H.trk_score[i] > TRK_SHR_SPLIT; };
    std::vector<TVector3> st(n), dir(n); std::vector<double> dv(n, 1e9);
    for(int i=0;i<n;++i){
      if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else         st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
      dir[i].SetXYZ(H.trk_dx[i],H.trk_dy[i],H.trk_dz[i]);
      dv[i] = (st[i]-vtx).Mag();
    }
    nTr=0; nSh=0; for(int i=0;i<n;++i)(isTrk(i)?nTr:nSh)++;
    if(showerMode == 0 && !(nSh>0 && nSh<NSHOWERS_MAX)) continue;
    if(showerMode == 1 && nSh != 0) continue;
    if(!(nTr < NTRACKS_MAX)) continue;                            stage=ST_NSHOW;  nCut[ST_NSHOW]++;

    int mu = -1;
    for(int i=0;i<n;++i)
      if(isTrk(i) && H.trk_score[i]>MUON_TRK_SCORE && dv[i]<MUON_TRK_DIST &&
         H.trk_len[i]>MUON_TRK_LEN && H.trk_llr[i]>MUON_LLR)
        if(mu==-1 || H.trk_len[i]>H.trk_len[mu]) mu = i;
    if(mu < 0) continue;                                          stage=ST_MUON;   nCut[ST_MUON]++;

    // ── everything below is RECORDED, not cut ───────────────────────────────
    mode = H.interaction; run=H.run; sub=H.sub; evt=H.evt;
    hypPdg = H.h_hyperon_pdg; nuPdg = H.nu_pdg; ccnc = H.ccnc; q2 = trueQ2;
    hasPair=0; precutOK=0; hasPhoton=0;
    mass=-1; decayLen=-1; cosAlpha=-2; vfitChi2=-1; sep=-1; llrP=-9; llrPi=-9;
    pMom=-1; iMom=-1; shrE=-1; pairScore=-1e30;
    for(int p=0;p<3;++p){ islPass[p]=0; islOK[p]=0; }
    islNPlanes=0; islAny=0; islSeeded=0;
    sigmaMass=-1; cosPhLam=-2; phDist=-1; phDedx=-1; phMoliere=-1; truePhE = isSig ? truePhMom : -1;
    michelDist=-1; michelE=-1; michelN5=0; muContained=0; lepTrue=lepPdg;
    muTrue = (mu<(int)H.bkt_pdg.size()) ? H.bkt_pdg[mu] : 0; prTrue=piTrue=s1True=0;
    trueFlight = isSig ? std::sqrt((lex-lsx)*(lex-lsx)+(ley-lsy)*(ley-lsy)+(lez-lsz)*(lez-lsz)) : -1;

    // best pair passing precuts; else best pair of all
    int pr=-1, pi=-1, prAny=-1, piAny=-1; double best=-1e30, bestAny=-1e30;
    double sepB=-1, llrPB=-9, llrPiB=-9, sepA=-1, llrPA=-9, llrPiA=-9;
    for(int a=0;a<n;++a){
      if(a==mu || !isTrk(a)) continue;
      for(int b=0;b<n;++b){
        if(b==mu || b==a || !isTrk(b)) continue;
        PairFeat f;
        f.sep=(st[a]-st[b]).Mag(); f.scoreP=H.trk_score[a]; f.scorePi=H.trk_score[b];
        f.llrP=H.trk_llr[a]; f.llrPi=H.trk_llr[b];
        if(haveDedx && d3 && a<(int)d3->size() && std::isfinite(d3->at(a))) f.dedxP =d3->at(a);
        if(haveDedx && d3 && b<(int)d3->size() && std::isfinite(d3->at(b))) f.dedxPi=d3->at(b);
        const double s = provisionalPairScore(f);
        if(s > bestAny){ bestAny=s; prAny=a; piAny=b; sepA=f.sep; llrPA=f.llrP; llrPiA=f.llrPi; }
        if(!passPrecuts(f, pc)) continue;
        if(s > best){ best=s; pr=a; pi=b; sepB=f.sep; llrPB=f.llrP; llrPiB=f.llrPi; }
      }
    }
    if(pr >= 0){ precutOK=1; hasPair=1; pairScore=best; sep=sepB; llrP=llrPB; llrPi=llrPiB; }
    else if(prAny >= 0){ hasPair=1; pr=prAny; pi=piAny; pairScore=bestAny; sep=sepA; llrP=llrPA; llrPi=llrPiA; }

    // photon: nearest shower to the vertex, excluding the candidates
    int s1=-1; double bs=1e9;
    for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi||isTrk(i)) continue;
      if(dv[i] < bs){ bs=dv[i]; s1=i; } }
    if(s1 >= 0){ hasPhoton=1; shrE=H.shr_energy_y[s1]/1000.0;
                 s1True = (s1<(int)H.bkt_pdg.size()) ? H.bkt_pdg[s1] : 0;
                 phDist=dv[s1]; phDedx=H.shr_dedx_y[s1]; phMoliere=H.shr_moliere[s1]; }

    // Michel tag: the muon candidate's end point against the start of every other
    // PFP. A mu+ always decays (Michel e+), a mu- captures 75% of the time in argon;
    // 68% of the surviving background muons are mu-.
    if(tex && mu<(int)tex->size()){
      const TVector3 mEnd(tex->at(mu), tey->at(mu), tez->at(mu));
      muContained = SigmaFV::isFiducial(mEnd.X(), mEnd.Y(), mEnd.Z());
      double bestD=1e9; int bi=-1;
      for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi) continue;
        const double d=(st[i]-mEnd).Mag(); if(d<5.0) michelN5++;
        if(d<bestD){ bestD=d; bi=i; } }
      if(bi>=0){ michelDist=bestD; michelE = isTrk(bi) ? -1 : H.shr_energy_y[bi]/1000.0; }
    }

    if(hasPair){
      prTrue = (pr<(int)H.bkt_pdg.size()) ? H.bkt_pdg[pr] : 0;
      piTrue = (pi<(int)H.bkt_pdg.size()) ? H.bkt_pdg[pi] : 0;
      pMom = protonMomFromRange(H.trk_len[pr]);
      iMom = H.trk_rmumom[pi];
      const DecayObs obs = decayObservables(vtx, st[pr],dir[pr],H.trk_len[pr],
                                                 st[pi],dir[pi],H.trk_len[pi], pMom, iMom);
      if(obs.ok){ mass=obs.mass; decayLen=obs.decayLength; cosAlpha=obs.cosAlpha; vfitChi2=obs.vfitChi2;
        // Sigma0 -> Lambda gamma: the photon's only use in this chain so far
        if(s1>=0 && mass>0 && shrE>0 && shrE<10){
          TVector3 ds(H.shr_px[s1],H.shr_py[s1],H.shr_pz[s1]); if(ds.Mag()>0) ds=ds.Unit();
          const TVector3 pL=obs.lambdaMom; const double EL=std::sqrt(mass*mass+pL.Mag2());
          if(pL.Mag()>0) cosPhLam = ds.Dot(pL.Unit());
          const TVector3 pS=pL+ds*shrE; const double ES=EL+shrE; const double m2=ES*ES-pS.Mag2();
          if(m2>0) sigmaMass=std::sqrt(m2); } }

      // island test on every plane, recorded per plane
      if(vw && vw->size() > 2)
        for(int p=0;p<3;++p){
          if(!win[p] || win[p]->empty() || vw->at(p) < 0) continue;
          bool ok=false;
          const bool pass = Island::planeTest(*win[p], vw->at(p)-(int)win[p]->size()/2, p,
                               st[mu].X(),st[mu].Y(),st[mu].Z(),
                               st[pr].X(),st[pr].Y(),st[pr].Z(),
                               st[pi].X(),st[pi].Y(),st[pi].Z(), ok);
          islPass[p]=pass; islOK[p]=ok;
          if(pass) islNPlanes++; if(ok) islSeeded++;
        }
      islAny = islNPlanes>0;

      // the hard-cut chain, for the printed table only
      if(precutOK){ stage=ST_PAIR; nCut[ST_PAIR]++;
        if(obs.ok && mass>0 && std::fabs(mass-M_LAMBDA)<LAM_MASS_WIN){ stage=ST_MASS; nCut[ST_MASS]++;
          if(decayLen>DECAY_LEN_MIN && cosAlpha>COS_ALPHA_MIN){ stage=ST_ALPHA; nCut[ST_ALPHA]++;
            if(hasPhoton && islAny) nIsland++; } } }
    }
    T.Fill();
  }
  fout.cd(); T.Write();
  printf("\n  %-16s %10s %10s\n", "cut", "events", "rel [%]");
  for(int c=0;c<NST;++c)
    printf("  %-16s %10ld %9.1f\n", STNAME[c], nCut[c], (c>0&&nCut[c-1])?100.0*nCut[c]/nCut[c-1]:100.0);
  printf("  %-16s %10ld %9.1f   (photon + island, hard-cut chain)\n", "Island", nIsland,
         nCut[ST_ALPHA]?100.0*nIsland/nCut[ST_ALPHA]:0.0);
  if(!FW.empty())
    printf("  flux weight: %ld matched (<w>=%.4f), %ld missing (left at 1)\n",
           nW, nW?sumW/nW:0.0, nWmiss);
  printf("  wrote %lld events to %s\n", T.GetEntries(), outName);
  fflush(stdout);
}
