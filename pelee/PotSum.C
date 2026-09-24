// Sum generated POT over a file list (ana/MetaTree has exactly one entry/file).
#include <fstream>
void PotSum(const char* list, int maxFiles=-1){
  std::ifstream in(list); std::string l; int n=0; double pot=0; long bad=0;
  while(std::getline(in,l)){ if(l.empty()) continue; if(maxFiles>0&&n>=maxFiles) break;
    TFile* f=TFile::Open(l.c_str(),"READ");
    if(!f||f->IsZombie()){ bad++; if(f) delete f; continue; }
    TTree* mt=(TTree*)f->Get("ana/MetaTree");
    if(mt&&mt->GetBranch("POT")){ double p=0; mt->SetBranchAddress("POT",&p); mt->GetEntry(0); pot+=p; n++; }
    else bad++;
    f->Close(); delete f; }
  printf("  %s : %d files, POT %.4e  (bad %ld)\n", list, n, pot, bad);
  if(n) printf("  mean POT/file %.4e ; full-list estimate %.4e\n", pot/n, pot/n*std::count(std::istreambuf_iterator<char>(*(new std::ifstream(list))), std::istreambuf_iterator<char>(), '\n'));
  fflush(stdout);
}
