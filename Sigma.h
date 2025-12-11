//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  3 13:33:43 2025 by ROOT version 6.36.04
// from TChain OutputTree/
//////////////////////////////////////////////////////////

#ifndef Sigma_h
#define Sigma_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"
#include "vector"
#include "vector"

class Sigma {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          subrun;
   UInt_t          event;
   UInt_t          n_mctruths;
   Int_t           mc_nu_pdg;
   Double_t        mc_nu_q2;
   string          *mc_ccnc;
   string          *mc_mode;
   Double_t        mc_nu_pos_x;
   Double_t        mc_nu_pos_y;
   Double_t        mc_nu_pos_z;
   Double_t        mc_lepton_start_x;
   Double_t        mc_lepton_start_y;
   Double_t        mc_lepton_start_z;
   Double_t        mc_lepton_end_x;
   Double_t        mc_lepton_end_y;
   Double_t        mc_lepton_end_z;
   Int_t           mc_lepton_pdg;
   Double_t        mc_lepton_mom;
   Double_t        mc_hyperon_start_x;
   Double_t        mc_hyperon_start_y;
   Double_t        mc_hyperon_start_z;
   Double_t        mc_hyperon_end_x;
   Double_t        mc_hyperon_end_y;
   Double_t        mc_hyperon_end_z;
   Int_t           mc_hyperon_pdg;
   Double_t        mc_hyperon_mom;
   Double_t        mc_sigmazero_photon_start_x;
   Double_t        mc_sigmazero_photon_start_y;
   Double_t        mc_sigmazero_photon_start_z;
   Double_t        mc_sigmazero_photon_end_x;
   Double_t        mc_sigmazero_photon_end_y;
   Double_t        mc_sigmazero_photon_end_z;
   Int_t           mc_sigmazero_photon_pdg;
   Double_t        mc_sigmazero_photon_mom;
   Double_t        mc_sigmazero_lambda_start_x;
   Double_t        mc_sigmazero_lambda_start_y;
   Double_t        mc_sigmazero_lambda_start_z;
   Double_t        mc_sigmazero_lambda_end_x;
   Double_t        mc_sigmazero_lambda_end_y;
   Double_t        mc_sigmazero_lambda_end_z;
   Int_t           mc_sigmazero_lambda_pdg;
   Double_t        mc_sigmazero_lambda_mom;
   vector<int>     *mc_decay_pdg;
   vector<double>  *mc_decay_start_x;
   vector<double>  *mc_decay_start_y;
   vector<double>  *mc_decay_start_z;
   vector<double>  *mc_decay_end_x;
   vector<double>  *mc_decay_end_y;
   vector<double>  *mc_decay_end_z;
   vector<double>  *mc_decay_mom;
   Int_t           true_nu_slice_ID;
   Double_t        true_nu_slice_completeness;
   Double_t        true_nu_slice_purity;
   UInt_t          n_slices;
   Int_t           flash_match_nu_slice_ID;
   Int_t           ct_nu_slice_ID;
   Int_t           pandora_nu_slice_ID;
   vector<double>  *pfp_purity;
   vector<double>  *pfp_completeness;
   vector<bool>    *pfp_has_truth;
   vector<int>     *pfp_trackID;
   vector<int>     *pfp_true_pdg;
   vector<double>  *pfp_true_energy;
   vector<double>  *pfp_true_ke;
   vector<double>  *pfp_true_px;
   vector<double>  *pfp_true_py;
   vector<double>  *pfp_true_pz;
   vector<double>  *pfp_true_length;
   vector<int>     *pfp_true_origin;
   vector<int>     *pfp_pdg;
   vector<double>  *pfp_trk_shr_score;
   vector<double>  *pfp_x;
   vector<double>  *pfp_y;
   vector<double>  *pfp_z;
   vector<int>     *pfpnplanehits_U;
   vector<int>     *pfpnplanehits_V;
   vector<int>     *pfpnplanehits_Y;
   Bool_t          reco_primary_vtx_inFV;
   Double_t        reco_primary_vtx_x;
   Double_t        reco_primary_vtx_y;
   Double_t        reco_primary_vtx_z;
   vector<double>  *trk_length;
   vector<double>  *trk_start_x;
   vector<double>  *trk_start_y;
   vector<double>  *trk_start_z;
   vector<double>  *trk_end_x;
   vector<double>  *trk_end_y;
   vector<double>  *trk_end_z;
   vector<double>  *trk_dir_x;
   vector<double>  *trk_dir_y;
   vector<double>  *trk_dir_z;
   vector<double>  *trk_mean_dedx_plane0;
   vector<double>  *trk_mean_dedx_plane1;
   vector<double>  *trk_mean_dedx_plane2;
   vector<double>  *trk_three_plane_dedx;
   vector<double>  *trk_llrpid;
   vector<double>  *shr_length;
   vector<double>  *shr_open_angle;
   vector<double>  *shr_start_x;
   vector<double>  *shr_start_y;
   vector<double>  *shr_start_z;
   vector<double>  *shr_dir_x;
   vector<double>  *shr_dir_y;
   vector<double>  *shr_dir_z;
   vector<double>  *shr_energy_plane0;
   vector<double>  *shr_energy_plane1;
   vector<double>  *shr_energy_plane2;
   vector<double>  *shr_dedx_plane0;
   vector<double>  *shr_dedx_plane1;
   vector<double>  *shr_dedx_plane2;
   vector<int>     *ct_test_primary_vtx_wires;
   vector<vector<float> > *ct_test_window_plane0;
   vector<vector<float> > *ct_test_window_plane1;
   vector<vector<float> > *ct_test_window_plane2;
   vector<double>  *pfp_photon_bdt_score;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_event;   //!
   TBranch        *b_n_mctruths;   //!
   TBranch        *b_mc_nu_pdg;   //!
   TBranch        *b_mc_nu_q2;   //!
   TBranch        *b_mc_ccnc;   //!
   TBranch        *b_mc_mode;   //!
   TBranch        *b_mc_nu_pos_x;   //!
   TBranch        *b_mc_nu_pos_y;   //!
   TBranch        *b_mc_nu_pos_z;   //!
   TBranch        *b_mc_lepton_start_x;   //!
   TBranch        *b_mc_lepton_start_y;   //!
   TBranch        *b_mc_lepton_start_z;   //!
   TBranch        *b_mc_lepton_end_x;   //!
   TBranch        *b_mc_lepton_end_y;   //!
   TBranch        *b_mc_lepton_end_z;   //!
   TBranch        *b_mc_lepton_pdg;   //!
   TBranch        *b_mc_lepton_mom;   //!
   TBranch        *b_mc_hyperon_start_x;   //!
   TBranch        *b_mc_hyperon_start_y;   //!
   TBranch        *b_mc_hyperon_start_z;   //!
   TBranch        *b_mc_hyperon_end_x;   //!
   TBranch        *b_mc_hyperon_end_y;   //!
   TBranch        *b_mc_hyperon_end_z;   //!
   TBranch        *b_mc_hyperon_pdg;   //!
   TBranch        *b_mc_hyperon_mom;   //!
   TBranch        *b_mc_sigmazero_photon_start_x;   //!
   TBranch        *b_mc_sigmazero_photon_start_y;   //!
   TBranch        *b_mc_sigmazero_photon_start_z;   //!
   TBranch        *b_mc_sigmazero_photon_end_x;   //!
   TBranch        *b_mc_sigmazero_photon_end_y;   //!
   TBranch        *b_mc_sigmazero_photon_end_z;   //!
   TBranch        *b_mc_sigmazero_photon_pdg;   //!
   TBranch        *b_mc_sigmazero_photon_mom;   //!
   TBranch        *b_mc_sigmazero_lambda_start_x;   //!
   TBranch        *b_mc_sigmazero_lambda_start_y;   //!
   TBranch        *b_mc_sigmazero_lambda_start_z;   //!
   TBranch        *b_mc_sigmazero_lambda_end_x;   //!
   TBranch        *b_mc_sigmazero_lambda_end_y;   //!
   TBranch        *b_mc_sigmazero_lambda_end_z;   //!
   TBranch        *b_mc_sigmazero_lambda_pdg;   //!
   TBranch        *b_mc_sigmazero_lambda_mom;   //!
   TBranch        *b_mc_decay_pdg;   //!
   TBranch        *b_mc_decay_start_x;   //!
   TBranch        *b_mc_decay_start_y;   //!
   TBranch        *b_mc_decay_start_z;   //!
   TBranch        *b_mc_decay_end_x;   //!
   TBranch        *b_mc_decay_end_y;   //!
   TBranch        *b_mc_decay_end_z;   //!
   TBranch        *b_mc_decay_mom;   //!
   TBranch        *b_true_nu_slice_ID;   //!
   TBranch        *b_true_nu_slice_completeness;   //!
   TBranch        *b_true_nu_slice_purity;   //!
   TBranch        *b_n_slices;   //!
   TBranch        *b_flash_match_nu_slice_ID;   //!
   TBranch        *b_ct_nu_slice_ID;   //!
   TBranch        *b_pandora_nu_slice_ID;   //!
   TBranch        *b_pfp_purity;   //!
   TBranch        *b_pfp_completeness;   //!
   TBranch        *b_pfp_has_truth;   //!
   TBranch        *b_pfp_trackID;   //!
   TBranch        *b_pfp_true_pdg;   //!
   TBranch        *b_pfp_true_energy;   //!
   TBranch        *b_pfp_true_ke;   //!
   TBranch        *b_pfp_true_px;   //!
   TBranch        *b_pfp_true_py;   //!
   TBranch        *b_pfp_true_pz;   //!
   TBranch        *b_pfp_true_length;   //!
   TBranch        *b_pfp_true_origin;   //!
   TBranch        *b_pfp_pdg;   //!
   TBranch        *b_pfp_trk_shr_score;   //!
   TBranch        *b_pfp_x;   //!
   TBranch        *b_pfp_y;   //!
   TBranch        *b_pfp_z;   //!
   TBranch        *b_pfpnplanehits_U;   //!
   TBranch        *b_pfpnplanehits_V;   //!
   TBranch        *b_pfpnplanehits_Y;   //!
   TBranch        *b_reco_primary_vtx_inFV;   //!
   TBranch        *b_reco_primary_vtx_x;   //!
   TBranch        *b_reco_primary_vtx_y;   //!
   TBranch        *b_reco_primary_vtx_z;   //!
   TBranch        *b_trk_length;   //!
   TBranch        *b_trk_start_x;   //!
   TBranch        *b_trk_start_y;   //!
   TBranch        *b_trk_start_z;   //!
   TBranch        *b_trk_end_x;   //!
   TBranch        *b_trk_end_y;   //!
   TBranch        *b_trk_end_z;   //!
   TBranch        *b_trk_dir_x;   //!
   TBranch        *b_trk_dir_y;   //!
   TBranch        *b_trk_dir_z;   //!
   TBranch        *b_trk_mean_dedx_plane0;   //!
   TBranch        *b_trk_mean_dedx_plane1;   //!
   TBranch        *b_trk_mean_dedx_plane2;   //!
   TBranch        *b_trk_three_plane_dedx;   //!
   TBranch        *b_trk_llrpid;   //!
   TBranch        *b_shr_length;   //!
   TBranch        *b_shr_open_angle;   //!
   TBranch        *b_shr_start_x;   //!
   TBranch        *b_shr_start_y;   //!
   TBranch        *b_shr_start_z;   //!
   TBranch        *b_shr_dir_x;   //!
   TBranch        *b_shr_dir_y;   //!
   TBranch        *b_shr_dir_z;   //!
   TBranch        *b_shr_energy_plane0;   //!
   TBranch        *b_shr_energy_plane1;   //!
   TBranch        *b_shr_energy_plane2;   //!
   TBranch        *b_shr_dedx_plane0;   //!
   TBranch        *b_shr_dedx_plane1;   //!
   TBranch        *b_shr_dedx_plane2;   //!
   TBranch        *b_ct_test_primary_vtx_wires;   //!
   TBranch        *b_ct_test_window_plane0;   //!
   TBranch        *b_ct_test_window_plane1;   //!
   TBranch        *b_ct_test_window_plane2;   //!
   TBranch        *b_pfp_photon_bdt_score;   //!

   Sigma(TTree *tree=0);
   virtual ~Sigma();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Sigma_cxx
Sigma::Sigma(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("OutputTree",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("OutputTree","");
     chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Hyperon_ana_hyperon_sigmazero_mSlice_prodgenie_run3b.root/OutputTree");
     chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Background_prodgenie_numi_overlay_rhc_run3b_reco2.root/OutputTree");
      chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Dirt_prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b.root/OutputTree");

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

Sigma::~Sigma()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Sigma::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Sigma::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Sigma::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   mc_ccnc = 0;
   mc_mode = 0;
   mc_decay_pdg = 0;
   mc_decay_start_x = 0;
   mc_decay_start_y = 0;
   mc_decay_start_z = 0;
   mc_decay_end_x = 0;
   mc_decay_end_y = 0;
   mc_decay_end_z = 0;
   mc_decay_mom = 0;
   pfp_purity = 0;
   pfp_completeness = 0;
   pfp_has_truth = 0;
   pfp_trackID = 0;
   pfp_true_pdg = 0;
   pfp_true_energy = 0;
   pfp_true_ke = 0;
   pfp_true_px = 0;
   pfp_true_py = 0;
   pfp_true_pz = 0;
   pfp_true_length = 0;
   pfp_true_origin = 0;
   pfp_pdg = 0;
   pfp_trk_shr_score = 0;
   pfp_x = 0;
   pfp_y = 0;
   pfp_z = 0;
   pfpnplanehits_U = 0;
   pfpnplanehits_V = 0;
   pfpnplanehits_Y = 0;
   trk_length = 0;
   trk_start_x = 0;
   trk_start_y = 0;
   trk_start_z = 0;
   trk_end_x = 0;
   trk_end_y = 0;
   trk_end_z = 0;
   trk_dir_x = 0;
   trk_dir_y = 0;
   trk_dir_z = 0;
   trk_mean_dedx_plane0 = 0;
   trk_mean_dedx_plane1 = 0;
   trk_mean_dedx_plane2 = 0;
   trk_three_plane_dedx = 0;
   trk_llrpid = 0;
   shr_length = 0;
   shr_open_angle = 0;
   shr_start_x = 0;
   shr_start_y = 0;
   shr_start_z = 0;
   shr_dir_x = 0;
   shr_dir_y = 0;
   shr_dir_z = 0;
   shr_energy_plane0 = 0;
   shr_energy_plane1 = 0;
   shr_energy_plane2 = 0;
   shr_dedx_plane0 = 0;
   shr_dedx_plane1 = 0;
   shr_dedx_plane2 = 0;
   ct_test_primary_vtx_wires = 0;
   ct_test_window_plane0 = 0;
   ct_test_window_plane1 = 0;
   ct_test_window_plane2 = 0;
   pfp_photon_bdt_score = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("n_mctruths", &n_mctruths, &b_n_mctruths);
   fChain->SetBranchAddress("mc_nu_pdg", &mc_nu_pdg, &b_mc_nu_pdg);
   fChain->SetBranchAddress("mc_nu_q2", &mc_nu_q2, &b_mc_nu_q2);
   fChain->SetBranchAddress("mc_ccnc", &mc_ccnc, &b_mc_ccnc);
   fChain->SetBranchAddress("mc_mode", &mc_mode, &b_mc_mode);
   fChain->SetBranchAddress("mc_nu_pos_x", &mc_nu_pos_x, &b_mc_nu_pos_x);
   fChain->SetBranchAddress("mc_nu_pos_y", &mc_nu_pos_y, &b_mc_nu_pos_y);
   fChain->SetBranchAddress("mc_nu_pos_z", &mc_nu_pos_z, &b_mc_nu_pos_z);
   fChain->SetBranchAddress("mc_lepton_start_x", &mc_lepton_start_x, &b_mc_lepton_start_x);
   fChain->SetBranchAddress("mc_lepton_start_y", &mc_lepton_start_y, &b_mc_lepton_start_y);
   fChain->SetBranchAddress("mc_lepton_start_z", &mc_lepton_start_z, &b_mc_lepton_start_z);
   fChain->SetBranchAddress("mc_lepton_end_x", &mc_lepton_end_x, &b_mc_lepton_end_x);
   fChain->SetBranchAddress("mc_lepton_end_y", &mc_lepton_end_y, &b_mc_lepton_end_y);
   fChain->SetBranchAddress("mc_lepton_end_z", &mc_lepton_end_z, &b_mc_lepton_end_z);
   fChain->SetBranchAddress("mc_lepton_pdg", &mc_lepton_pdg, &b_mc_lepton_pdg);
   fChain->SetBranchAddress("mc_lepton_mom", &mc_lepton_mom, &b_mc_lepton_mom);
   fChain->SetBranchAddress("mc_hyperon_start_x", &mc_hyperon_start_x, &b_mc_hyperon_start_x);
   fChain->SetBranchAddress("mc_hyperon_start_y", &mc_hyperon_start_y, &b_mc_hyperon_start_y);
   fChain->SetBranchAddress("mc_hyperon_start_z", &mc_hyperon_start_z, &b_mc_hyperon_start_z);
   fChain->SetBranchAddress("mc_hyperon_end_x", &mc_hyperon_end_x, &b_mc_hyperon_end_x);
   fChain->SetBranchAddress("mc_hyperon_end_y", &mc_hyperon_end_y, &b_mc_hyperon_end_y);
   fChain->SetBranchAddress("mc_hyperon_end_z", &mc_hyperon_end_z, &b_mc_hyperon_end_z);
   fChain->SetBranchAddress("mc_hyperon_pdg", &mc_hyperon_pdg, &b_mc_hyperon_pdg);
   fChain->SetBranchAddress("mc_hyperon_mom", &mc_hyperon_mom, &b_mc_hyperon_mom);
   fChain->SetBranchAddress("mc_sigmazero_photon_start_x", &mc_sigmazero_photon_start_x, &b_mc_sigmazero_photon_start_x);
   fChain->SetBranchAddress("mc_sigmazero_photon_start_y", &mc_sigmazero_photon_start_y, &b_mc_sigmazero_photon_start_y);
   fChain->SetBranchAddress("mc_sigmazero_photon_start_z", &mc_sigmazero_photon_start_z, &b_mc_sigmazero_photon_start_z);
   fChain->SetBranchAddress("mc_sigmazero_photon_end_x", &mc_sigmazero_photon_end_x, &b_mc_sigmazero_photon_end_x);
   fChain->SetBranchAddress("mc_sigmazero_photon_end_y", &mc_sigmazero_photon_end_y, &b_mc_sigmazero_photon_end_y);
   fChain->SetBranchAddress("mc_sigmazero_photon_end_z", &mc_sigmazero_photon_end_z, &b_mc_sigmazero_photon_end_z);
   fChain->SetBranchAddress("mc_sigmazero_photon_pdg", &mc_sigmazero_photon_pdg, &b_mc_sigmazero_photon_pdg);
   fChain->SetBranchAddress("mc_sigmazero_photon_mom", &mc_sigmazero_photon_mom, &b_mc_sigmazero_photon_mom);
   fChain->SetBranchAddress("mc_sigmazero_lambda_start_x", &mc_sigmazero_lambda_start_x, &b_mc_sigmazero_lambda_start_x);
   fChain->SetBranchAddress("mc_sigmazero_lambda_start_y", &mc_sigmazero_lambda_start_y, &b_mc_sigmazero_lambda_start_y);
   fChain->SetBranchAddress("mc_sigmazero_lambda_start_z", &mc_sigmazero_lambda_start_z, &b_mc_sigmazero_lambda_start_z);
   fChain->SetBranchAddress("mc_sigmazero_lambda_end_x", &mc_sigmazero_lambda_end_x, &b_mc_sigmazero_lambda_end_x);
   fChain->SetBranchAddress("mc_sigmazero_lambda_end_y", &mc_sigmazero_lambda_end_y, &b_mc_sigmazero_lambda_end_y);
   fChain->SetBranchAddress("mc_sigmazero_lambda_end_z", &mc_sigmazero_lambda_end_z, &b_mc_sigmazero_lambda_end_z);
   fChain->SetBranchAddress("mc_sigmazero_lambda_pdg", &mc_sigmazero_lambda_pdg, &b_mc_sigmazero_lambda_pdg);
   fChain->SetBranchAddress("mc_sigmazero_lambda_mom", &mc_sigmazero_lambda_mom, &b_mc_sigmazero_lambda_mom);
   fChain->SetBranchAddress("mc_decay_pdg", &mc_decay_pdg, &b_mc_decay_pdg);
   fChain->SetBranchAddress("mc_decay_start_x", &mc_decay_start_x, &b_mc_decay_start_x);
   fChain->SetBranchAddress("mc_decay_start_y", &mc_decay_start_y, &b_mc_decay_start_y);
   fChain->SetBranchAddress("mc_decay_start_z", &mc_decay_start_z, &b_mc_decay_start_z);
   fChain->SetBranchAddress("mc_decay_end_x", &mc_decay_end_x, &b_mc_decay_end_x);
   fChain->SetBranchAddress("mc_decay_end_y", &mc_decay_end_y, &b_mc_decay_end_y);
   fChain->SetBranchAddress("mc_decay_end_z", &mc_decay_end_z, &b_mc_decay_end_z);
   fChain->SetBranchAddress("mc_decay_mom", &mc_decay_mom, &b_mc_decay_mom);
   fChain->SetBranchAddress("true_nu_slice_ID", &true_nu_slice_ID, &b_true_nu_slice_ID);
   fChain->SetBranchAddress("true_nu_slice_completeness", &true_nu_slice_completeness, &b_true_nu_slice_completeness);
   fChain->SetBranchAddress("true_nu_slice_purity", &true_nu_slice_purity, &b_true_nu_slice_purity);
   fChain->SetBranchAddress("n_slices", &n_slices, &b_n_slices);
   fChain->SetBranchAddress("flash_match_nu_slice_ID", &flash_match_nu_slice_ID, &b_flash_match_nu_slice_ID);
   fChain->SetBranchAddress("ct_nu_slice_ID", &ct_nu_slice_ID, &b_ct_nu_slice_ID);
   fChain->SetBranchAddress("pandora_nu_slice_ID", &pandora_nu_slice_ID, &b_pandora_nu_slice_ID);
   fChain->SetBranchAddress("pfp_purity", &pfp_purity, &b_pfp_purity);
   fChain->SetBranchAddress("pfp_completeness", &pfp_completeness, &b_pfp_completeness);
   fChain->SetBranchAddress("pfp_has_truth", &pfp_has_truth, &b_pfp_has_truth);
   fChain->SetBranchAddress("pfp_trackID", &pfp_trackID, &b_pfp_trackID);
   fChain->SetBranchAddress("pfp_true_pdg", &pfp_true_pdg, &b_pfp_true_pdg);
   fChain->SetBranchAddress("pfp_true_energy", &pfp_true_energy, &b_pfp_true_energy);
   fChain->SetBranchAddress("pfp_true_ke", &pfp_true_ke, &b_pfp_true_ke);
   fChain->SetBranchAddress("pfp_true_px", &pfp_true_px, &b_pfp_true_px);
   fChain->SetBranchAddress("pfp_true_py", &pfp_true_py, &b_pfp_true_py);
   fChain->SetBranchAddress("pfp_true_pz", &pfp_true_pz, &b_pfp_true_pz);
   fChain->SetBranchAddress("pfp_true_length", &pfp_true_length, &b_pfp_true_length);
   fChain->SetBranchAddress("pfp_true_origin", &pfp_true_origin, &b_pfp_true_origin);
   fChain->SetBranchAddress("pfp_pdg", &pfp_pdg, &b_pfp_pdg);
   fChain->SetBranchAddress("pfp_trk_shr_score", &pfp_trk_shr_score, &b_pfp_trk_shr_score);
   fChain->SetBranchAddress("pfp_x", &pfp_x, &b_pfp_x);
   fChain->SetBranchAddress("pfp_y", &pfp_y, &b_pfp_y);
   fChain->SetBranchAddress("pfp_z", &pfp_z, &b_pfp_z);
   fChain->SetBranchAddress("pfpnplanehits_U", &pfpnplanehits_U, &b_pfpnplanehits_U);
   fChain->SetBranchAddress("pfpnplanehits_V", &pfpnplanehits_V, &b_pfpnplanehits_V);
   fChain->SetBranchAddress("pfpnplanehits_Y", &pfpnplanehits_Y, &b_pfpnplanehits_Y);
   fChain->SetBranchAddress("reco_primary_vtx_inFV", &reco_primary_vtx_inFV, &b_reco_primary_vtx_inFV);
   fChain->SetBranchAddress("reco_primary_vtx_x", &reco_primary_vtx_x, &b_reco_primary_vtx_x);
   fChain->SetBranchAddress("reco_primary_vtx_y", &reco_primary_vtx_y, &b_reco_primary_vtx_y);
   fChain->SetBranchAddress("reco_primary_vtx_z", &reco_primary_vtx_z, &b_reco_primary_vtx_z);
   fChain->SetBranchAddress("trk_length", &trk_length, &b_trk_length);
   fChain->SetBranchAddress("trk_start_x", &trk_start_x, &b_trk_start_x);
   fChain->SetBranchAddress("trk_start_y", &trk_start_y, &b_trk_start_y);
   fChain->SetBranchAddress("trk_start_z", &trk_start_z, &b_trk_start_z);
   fChain->SetBranchAddress("trk_end_x", &trk_end_x, &b_trk_end_x);
   fChain->SetBranchAddress("trk_end_y", &trk_end_y, &b_trk_end_y);
   fChain->SetBranchAddress("trk_end_z", &trk_end_z, &b_trk_end_z);
   fChain->SetBranchAddress("trk_dir_x", &trk_dir_x, &b_trk_dir_x);
   fChain->SetBranchAddress("trk_dir_y", &trk_dir_y, &b_trk_dir_y);
   fChain->SetBranchAddress("trk_dir_z", &trk_dir_z, &b_trk_dir_z);
   fChain->SetBranchAddress("trk_mean_dedx_plane0", &trk_mean_dedx_plane0, &b_trk_mean_dedx_plane0);
   fChain->SetBranchAddress("trk_mean_dedx_plane1", &trk_mean_dedx_plane1, &b_trk_mean_dedx_plane1);
   fChain->SetBranchAddress("trk_mean_dedx_plane2", &trk_mean_dedx_plane2, &b_trk_mean_dedx_plane2);
   fChain->SetBranchAddress("trk_three_plane_dedx", &trk_three_plane_dedx, &b_trk_three_plane_dedx);
   fChain->SetBranchAddress("trk_llrpid", &trk_llrpid, &b_trk_llrpid);
   fChain->SetBranchAddress("shr_length", &shr_length, &b_shr_length);
   fChain->SetBranchAddress("shr_open_angle", &shr_open_angle, &b_shr_open_angle);
   fChain->SetBranchAddress("shr_start_x", &shr_start_x, &b_shr_start_x);
   fChain->SetBranchAddress("shr_start_y", &shr_start_y, &b_shr_start_y);
   fChain->SetBranchAddress("shr_start_z", &shr_start_z, &b_shr_start_z);
   fChain->SetBranchAddress("shr_dir_x", &shr_dir_x, &b_shr_dir_x);
   fChain->SetBranchAddress("shr_dir_y", &shr_dir_y, &b_shr_dir_y);
   fChain->SetBranchAddress("shr_dir_z", &shr_dir_z, &b_shr_dir_z);
   fChain->SetBranchAddress("shr_energy_plane0", &shr_energy_plane0, &b_shr_energy_plane0);
   fChain->SetBranchAddress("shr_energy_plane1", &shr_energy_plane1, &b_shr_energy_plane1);
   fChain->SetBranchAddress("shr_energy_plane2", &shr_energy_plane2, &b_shr_energy_plane2);
   fChain->SetBranchAddress("shr_dedx_plane0", &shr_dedx_plane0, &b_shr_dedx_plane0);
   fChain->SetBranchAddress("shr_dedx_plane1", &shr_dedx_plane1, &b_shr_dedx_plane1);
   fChain->SetBranchAddress("shr_dedx_plane2", &shr_dedx_plane2, &b_shr_dedx_plane2);
   fChain->SetBranchAddress("ct_test_primary_vtx_wires", &ct_test_primary_vtx_wires, &b_ct_test_primary_vtx_wires);
   fChain->SetBranchAddress("ct_test_window_plane0", &ct_test_window_plane0, &b_ct_test_window_plane0);
   fChain->SetBranchAddress("ct_test_window_plane1", &ct_test_window_plane1, &b_ct_test_window_plane1);
   fChain->SetBranchAddress("ct_test_window_plane2", &ct_test_window_plane2, &b_ct_test_window_plane2);
   fChain->SetBranchAddress("pfp_photon_bdt_score", &pfp_photon_bdt_score, &b_pfp_photon_bdt_score);
   Notify();
}

bool Sigma::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void Sigma::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Sigma::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Sigma_cxx
