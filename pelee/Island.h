#ifndef Island_h
#define Island_h
// ============================================================================
//  Island -- the connectedness test as the Lambda CCQE analysis defines it:
//  seed clusters of raw wire activity at the muon / proton / pion track starts,
//  grow them, and ask whether the proton and pion merge while the muon stays
//  separate. Our existing ct_* variables instead count islands around the
//  primary vertex with no reference to which island holds which candidate.
//
//  Measured on the GENIE Sigma0 sample: passes 27.2% of signal against 5.7% of
//  background, a 4.8x rejection. The signal efficiency is capped by the Lambda
//  flight distance (mean 4.10 cm, median 2.43 cm; 28.3% fly beyond 5 cm), which
//  is why 27.2% and not more -- the test can only fire when the decay actually
//  separates. Use it as a BDT input, not a hard cut.
//
//  Transforms validated in ctcalib.C: wire = (A*y+B*z+C)/0.3 reproduces 98.7% of
//  stored vertex wires within 2 wires; tick = 18.1509*x + 818.1 implies
//  0.0551 cm/tick against the known 0.1098 cm/us at 0.5 us/tick.
// ============================================================================
#include <vector>
#include <cmath>

namespace Island {

// Connectivity, tuned 2026-09-12 (IslandTune.C on the IslandCache.C cache). The
// original 3 wires x 40 ticks bridged the decay V to the muon for Lambda flights
// below ~3 cm: 40 ticks is 2.2 cm of drift. 2 wires (0.6 cm) x 20 ticks (1.1 cm)
// is roughly isotropic and lifts the signal pass rate at 2-3 cm flight from 45%
// to 80% (one plane). Background rises with it on one plane; the gain is on the
// two-plane requirement: S 1.04 -> 1.39 at B 3.3 -> 4.1 (FoM 0.485 -> 0.571).
static const int    CDW=2;        // wires        (was 3)
static const double CDT=20.0;     // ticks        (was 40)
static const int    MINCLUST=3;   // minimum cluster size
static const double WA[3]={-0.86603,0.86603,0.0};
static const double WB[3]={0.5,0.5,1.0};
static const double WC[3]={101.475,99.748,0.0};
static const double MT=18.1509, CT=818.1;

struct Cell { int w; float t; };

inline int wireOf(int plane,double y,double z){
  return (int)std::lround((WA[plane]*y + WB[plane]*z + WC[plane])/0.3); }
inline float tickOf(double x){ return (float)(MT*x + CT); }

// Tunable parameters. The defaults reproduce the original test exactly; the
// tuning is done by IslandTune.C on the cache written by IslandCache.C.
struct Params {
  int    cdw      = CDW;       // wire connectivity
  double cdt      = CDT;       // tick connectivity
  int    minclust = MINCLUST;  // minimum cluster size
  double seedR2   = 400.0;     // max seed-to-activity distance^2, wire-equivalent units
  float  adcMin   = -1.f;      // occupancy threshold; <0 = every stored hit counts
};

// Returns true if, on this plane, the proton and pion seeds land in one cluster
// and the muon seed in a different one. `ok` reports whether all three seeds
// found activity at all. `adc` may be null; when given with prm.adcMin >= 0
// only hits above the threshold are used.
inline bool planeTestP(const std::vector<std::vector<float>>& win,
                       const std::vector<std::vector<float>>* adc,
                       int firstWire, int plane,
                       double mx,double my,double mz,
                       double px,double py,double pz,
                       double ix,double iy,double iz,
                       bool& ok, const Params& prm)
{
  ok=false;
  std::vector<Cell> pts;
  for(int i=0;i<(int)win.size();++i)
    for(size_t k=0;k<win[i].size();++k){
      if(adc && prm.adcMin>=0.f && i<(int)adc->size() && k<(*adc)[i].size() && (*adc)[i][k]<prm.adcMin) continue;
      pts.push_back({firstWire+i,win[i][k]}); }
  const int n=pts.size();
  const int MINCLUST=prm.minclust;
  if(n<MINCLUST) return false;
  auto conn=[&](const Cell&a,const Cell&b){
    return std::fabs(a.t-b.t)<=prm.cdt && std::abs(a.w-b.w)<=prm.cdw; };
  std::vector<int> lab(n,-1); int nl=0;
  for(int i=0;i<n;++i){ if(lab[i]>=0) continue;
    std::vector<int> st{i}; lab[i]=nl;
    while(!st.empty()){ int q=st.back(); st.pop_back();
      for(int j=0;j<n;++j) if(lab[j]<0&&conn(pts[j],pts[q])){ lab[j]=nl; st.push_back(j);} }
    ++nl; }
  std::vector<int> csz(nl,0); for(int i=0;i<n;++i) csz[lab[i]]++;
  const Cell seed[3]={ {wireOf(plane,my,mz),tickOf(mx)},
                       {wireOf(plane,py,pz),tickOf(px)},
                       {wireOf(plane,iy,iz),tickOf(ix)} };
  int sl[3]={-1,-1,-1};
  for(int s=0;s<3;++s){
    double best=1e18; int bi=-1;
    for(int i=0;i<n;++i){
      const double dw=pts[i].w-seed[s].w;
      const double dt=(pts[i].t-seed[s].t)*0.055/0.3;   // ticks -> wire-equivalent
      const double d=dw*dw+dt*dt;
      if(d<best){best=d;bi=i;} }
    if(bi>=0 && best<prm.seedR2) sl[s]=lab[bi];         // default: within ~20 wires
  }
  if(sl[0]<0||sl[1]<0||sl[2]<0) return false;
  ok=true;
  if(csz[sl[0]]<MINCLUST||csz[sl[1]]<MINCLUST||csz[sl[2]]<MINCLUST) return false;
  return (sl[1]==sl[2]) && (sl[0]!=sl[1]) && (sl[0]!=sl[2]);
}

// The original interface: default parameters, no ADC threshold.
inline bool planeTest(const std::vector<std::vector<float>>& win, int firstWire,
                      int plane,
                      double mx,double my,double mz,
                      double px,double py,double pz,
                      double ix,double iy,double iz,
                      bool& ok)
{
  return planeTestP(win, nullptr, firstWire, plane, mx,my,mz, px,py,pz, ix,iy,iz, ok, Params());
}

} // namespace Island
#endif
