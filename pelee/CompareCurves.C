// CompareCurves.C -- overlay the BDT curves before and after the 2026-09-09 retrain
// on the completed hypana file set. Each curve is the one produced from its own tree,
// so this is a like-for-like comparison of the two published working points.
#include <fstream>
#include <sstream>
struct Row{ double eff,cut,S,B,effb,pur,fom,xb; };
static std::vector<Row> parse(const char* fn){
  std::vector<Row> v; std::ifstream in(fn); std::string l;
  while(std::getline(in,l)){
    std::istringstream ss(l); Row r;
    if(!(ss>>r.eff>>r.cut>>r.S>>r.B>>r.effb>>r.pur>>r.fom>>r.xb)) continue;
    if(r.eff<1||r.eff>100) continue;          // skip the -1e9 baseline row
    v.push_back(r);
  }
  return v;
}
void CompareCurves(){
  auto A=parse("curve_prelists.log");   // before
  auto B=parse("curve_full.log");       // after
  printf("  parsed before=%zu points  after=%zu points\n",A.size(),B.size());
  auto mk=[](std::vector<Row>&v,int which,std::vector<double>&x,std::vector<double>&y){
    x.clear(); y.clear();
    for(auto&r:v){ x.push_back(r.eff);
      y.push_back(which==0? 100-r.effb : which==1? r.xb : r.pur); } };
  gStyle->SetOptStat(0);
  TCanvas c("c","",1600,520); c.Divide(3,1);
  const char* tit[3]={"ROC;signal efficiency [%];background rejection [%]",
                      "figure of merit;signal efficiency [%];FoM / baseline",
                      "purity;signal efficiency [%];purity [%]"};
  std::vector<TGraph*> keep;
  for(int p=0;p<3;++p){
    c.cd(p+1); gPad->SetGrid();
    std::vector<double> ax,ay,bx,by;
    mk(A,p,ax,ay); mk(B,p,bx,by);
    TGraph* ga=new TGraph(ax.size(),&ax[0],&ay[0]);
    TGraph* gb=new TGraph(bx.size(),&bx[0],&by[0]);
    ga->SetTitle(tit[p]);
    ga->SetLineColor(kGray+2); ga->SetLineWidth(3); ga->SetLineStyle(2);
    ga->SetMarkerStyle(24); ga->SetMarkerColor(kGray+2);
    gb->SetLineColor(kAzure+2); gb->SetLineWidth(3);
    gb->SetMarkerStyle(20); gb->SetMarkerColor(kAzure+2);
    double lo=1e9,hi=-1e9;
    for(double y:ay){lo=std::min(lo,y);hi=std::max(hi,y);}
    for(double y:by){lo=std::min(lo,y);hi=std::max(hi,y);}
    double pad=0.08*(hi-lo); ga->GetYaxis()->SetRangeUser(lo-pad,hi+pad);
    ga->Draw("ALP"); gb->Draw("LP same");
    keep.push_back(ga); keep.push_back(gb);
    if(p==1){ TLine* l=new TLine(0,1,100,1); l->SetLineStyle(3); l->Draw(); }
    // park the legend in each panel's empty corner: ROC and FoM fall to the
    // right, so bottom-left is clear; purity decays, so top-right is clear.
    TLegend* lg = (p==2) ? new TLegend(0.40,0.72,0.89,0.87)
                         : new TLegend(0.14,0.14,0.63,0.29);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
    lg->AddEntry(ga,"before  (3806 sig / 1689 bkg files)","lp");
    lg->AddEntry(gb,"after   (4851 sig / 3491 bkg files)","lp");
    lg->Draw();
  }
  c.SaveAs("bdt_curve_before_after.png");
  printf("  wrote bdt_curve_before_after.png\n");
}
