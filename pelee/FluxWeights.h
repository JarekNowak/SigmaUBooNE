// ============================================================================
//  FluxWeights.h -- old -> current NuMI flux weights for the NuWro signal,
//  keyed on (run,subrun,event). Shared by HypCutFlow.C and IslandProcess.C.
//
//  NuWro was generated on the old (bugged) NuMI flux while the overlay background and
//  the data are on the current g4 10.4 flux, so NuWro yields must be weighted down
//  (<w> = 0.677, and the same for Sigma0 events despite their higher energy). Table:
//  ~/MicroBooNE/production/nuwro_prod/flux_rw/nuwro_rhc_run3b_reco2_fluxweights.tsv.
//  Apply to the SIGNAL pass only -- the NuWro and overlay samples are built from the
//  same run3b EXT stream, so (run,subrun,event) is not unique across the two.
// ============================================================================
#ifndef FLUXWEIGHTS_H
#define FLUXWEIGHTS_H
#include <fstream>
#include <map>
#include <string>
#include <cstdio>

static inline unsigned long long rseKey(unsigned r, unsigned s, unsigned e){
  return ((unsigned long long)r*100000ULL + s)*10000000ULL + e;
}

static std::map<unsigned long long,double> loadFluxWeights(const char* p){
  std::map<unsigned long long,double> m;
  if(!p || !*p) return m;
  std::ifstream in(p);
  if(!in){ printf("  WARNING: flux weight file %s not readable -- NuWro left on the OLD flux\n",p); return m; }
  std::string l;
  while(std::getline(in,l)){
    if(l.empty() || l[0]=='#') continue;
    unsigned r,su,e; double w;
    if(sscanf(l.c_str(),"%u %u %u %lf",&r,&su,&e,&w)==4) m[rseKey(r,su,e)]=w;
  }
  printf("  flux weights: %zu events from %s\n",m.size(),p);
  return m;
}
#endif
