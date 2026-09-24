#ifndef DisplayUtil_h
#define DisplayUtil_h
// ============================================================================
//  DisplayUtil -- shared drawing helpers for the candidate-level event displays.
//  Neither ntuple format carries wire hits, so a "display" here is the set of
//  reconstructed candidates: the vertex, the muon/proton/pion track sticks and
//  the photon shower, drawn in two orthogonal projections.
// ============================================================================
#include "TVector3.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TVirtualPad.h"
#include "TString.h"
#include <vector>

struct Seg { TVector3 a,b; int col,wid,sty; };

// PDG code -> ROOT latex label.
inline const char* pdgName(int p){
  switch(p){ case 13:return "#mu^{-}"; case -13:return "#mu^{+}"; case 2212:return "p";
    case 2112:return "n"; case 211:return "#pi^{+}"; case -211:return "#pi^{-}";
    case 111:return "#pi^{0}"; case 11:return "e^{-}"; case -11:return "e^{+}";
    case 22:return "#gamma"; case 321:return "K^{+}"; case -321:return "K^{-}";
    case 3122:return "#Lambda"; case 3212:return "#Sigma^{0}";
    case 0: case -1: return "unmatched"; default:return Form("%d",p); }
}
// GENIE interaction code -> short name.
inline const char* modeName(int c){
  switch(c){ case 0:return "QE"; case 1:return "RES"; case 2:return "DIS";
             case 3:return "COH"; case 10:return "MEC"; default:return "Other"; }
}

// Draw one projection. ia/ib pick the axes (0=X,1=Y,2=Z).
inline void drawView(TVirtualPad* pad, const std::vector<Seg>& segs, const TVector3& vtx,
                     int ia, int ib, const char* xt, const char* yt)
{
  pad->cd();
  double lo[3]={1e9,1e9,1e9}, hi[3]={-1e9,-1e9,-1e9};
  auto upd=[&](const TVector3& v){ double c[3]={v.X(),v.Y(),v.Z()};
    for(int k=0;k<3;k++){ if(c[k]<lo[k])lo[k]=c[k]; if(c[k]>hi[k])hi[k]=c[k]; } };
  upd(vtx); for(const auto& s:segs){ upd(s.a); upd(s.b); }
  const double xa[3]={vtx.X(),vtx.Y(),vtx.Z()};
  double xmin=lo[ia],xmax=hi[ia],ymin=lo[ib],ymax=hi[ib];
  const double mx=0.15*(xmax-xmin)+5, my=0.15*(ymax-ymin)+5;
  TH2F* fr=new TH2F(Form("fr%p",(void*)pad),"",10,xmin-mx,xmax+mx,10,ymin-my,ymax+my);
  fr->SetStats(0); fr->GetXaxis()->SetTitle(xt); fr->GetYaxis()->SetTitle(yt); fr->Draw();
  for(const auto& s:segs){
    const double A[3]={s.a.X(),s.a.Y(),s.a.Z()}, B[3]={s.b.X(),s.b.Y(),s.b.Z()};
    TLine* l=new TLine(A[ia],A[ib],B[ia],B[ib]);
    l->SetLineColor(s.col); l->SetLineWidth(s.wid); l->SetLineStyle(s.sty); l->Draw();
  }
  TMarker* m=new TMarker(xa[ia],xa[ib],29);
  m->SetMarkerSize(2.2); m->SetMarkerColor(kBlack); m->Draw();
}
#endif
