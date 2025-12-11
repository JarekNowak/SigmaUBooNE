#define Sigma_cxx
#include "Sigma.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Sigma::Loop()
{
  gStyle->SetOptStat(1001111);
  TString Names[5];
  Names[0]="Sigma";
  Names[1]="Hyperons";
  Names[2]="Background";
  Names[3]="Dirt";
  Names[4]="Data";
  
  
  int NEvents[5]={0};
  NEvents[0]=262576;//Sigma
  NEvents[1]=NEvents[0];
  NEvents[2]=651830+NEvents[0];//Background
  NEvents[3]=349879+NEvents[2];//Dirt
  NEvents[4]=0 +NEvents[3]; //data
  
  for(int i=0;i<5;i++) cout<<NEvents[i]<<endl;
  double PoT[5]={0};
  PoT[0] = 552.0; //*e20
  PoT[1] = PoT[0];
  PoT[2] = 9.280;  // *e20
  PoT[3] = 17.30;//*e20
  PoT[4] = 4.90; //*e20 - pot for run 3b	
  
  double scale[5]={0};
  scale[0]=PoT[4]/PoT[0];
  scale[1]=PoT[4]/PoT[1];
  scale[2]=PoT[4]/PoT[2];
  scale[3]=PoT[4]/PoT[3];
  scale[4]=1.0;
  
  for(int i=0;i<5;i++)cout<<"Scale="<<scale[i]<<endl;
  const int Nc=11;
  TString CutName[Nc];
  CutName[0]="AllEvents";
  CutName[1]="FducialV";
  CutName[2]="NShowers"; 
  CutName[3]="NoBigShowers";
  CutName[4]="MuonCandidate";
  CutName[5]="PionCandidate";
  CutName[6]="ProtonCandidate";
  CutName[7]="NTracksCut";
  CutName[8]="2ndShower";
  CutName[9]="1stShower";
  CutName[10]="Distances";
 
	  
  const int Nv=80;
  TString VarName[Nv];
  VarName[0]="NTracks";
  VarName[1]="NShowers";
  VarName[2]="NHitsY";
  VarName[3]="MuonPID";
  VarName[4]="MuonTrkLenght";
  VarName[5]="MuonTrkDistVtx";
  VarName[6]="MuonTrkScore";
  
  VarName[7]="PionPID";
  VarName[8]="PionTrkLenght";
  VarName[9]="PionTrkDistVtx";
  VarName[10]="PionTrkScore";
  
  VarName[11]="ProtonPID";
  VarName[12]="ProtonTrkLenght";
  VarName[13]="ProtonTrkDistVtx";
  VarName[14]="ProtonTrkScore";
  
  
  VarName[15]="GammaPID";
  VarName[16]="GammaTrkLenght";
  VarName[17]="GammaTrkDistVtx";
  VarName[18]="GammaTrkScore";
  
  VarName[19]="MuonNHitsY";
  VarName[20]="PionNHitsY";
  VarName[21]="ProtonNHitsY";
  VarName[22]="GammaNHitsY";
  
  
  VarName[23]="MuonNHitsU";
  VarName[24]="PionNHitsU";
  VarName[25]="ProtonNHitsU";
  VarName[26]="GammaNHitsU";
  
  
  VarName[27]="MuonNHitsV";
  VarName[28]="PionNHitsV";
  VarName[29]="ProtonNHitsV";
  VarName[30]="GammaNHitsV";
  
  
  VarName[31]="MuonMeanDeDxPlane0";//trk_mean_dedx_plane0
  VarName[32]="MuonMeanDeDxPlane1";
  VarName[33]="MuonMeanDeDxPlane2";
  
  VarName[34]="PionMeanDeDxPlane0";//trk_mean_dedx_plane0
  VarName[35]="PionMeanDeDxPlane1";
  VarName[36]="PionMeanDeDxPlane2";
  
  VarName[37]="ProtonMeanDeDxPlane0";//trk_mean_dedx_plane0
  VarName[38]="ProtonMeanDeDxPlane1";
  VarName[39]="ProtonMeanDeDxPlane2";

  VarName[40]="GammaMeanDeDxPlane0";//trk_mean_dedx_plane0
  VarName[41]="GammaMeanDeDxPlane1";
  VarName[42]="GammaMeanDeDxPlane2";
  
  VarName[43]="MuonPhotonPID";
  VarName[44]="PionPhotonPID";
  VarName[45]="ProtonPhotonPID";
  VarName[46]="GammaPhotonPID";
  
 
  VarName[47]="MuonCandShrOpenAngle";
  VarName[48]="PionCandShrOpenAngle";
  VarName[49]="ProtonCandShrOpenAngle";
  VarName[50]="Shower1CandShrOpenAngle";
  VarName[51]="Shower2CandShrOpenAngle";

  VarName[52]="Shower1CandDeDxPlan0";//shr_dedx_plane1
  VarName[53]="Shower1CandDeDxPlan1";
  VarName[54]="Shower1CandDeDxPlan2";

  VarName[55]="Shower2CandDeDxPlan0";//shr_dedx_plane1
  VarName[56]="Shower3CandDeDxPlan1";
  VarName[57]="Shower4CandDeDxPlan2";

 VarName[58]="OpeningAngleMuonPion";
 VarName[59]="OpeningAngleProtonPion";
 VarName[60]="OpeningAngleProtonMuon";

 VarName[61]="OpeningAngleShower1Muon";
 VarName[62]="OpeningAngleShower1Pion";
 VarName[63]="OpeningAngleShower1Proton";

 VarName[64]="OpeningAngleShower2Muon";
 VarName[65]="OpeningAngleShower2Pion";
 VarName[66]="OpeningAngleShower2Proton";


 VarName[67]="Shower1NHitsY";
 VarName[68]="Shower2NHitsY";

 VarName[69]="Shower1VtxDist";
 VarName[70]="Shower2VtxDist";

 VarName[71]="DistMuonPion";
 VarName[72]="DistProtonPion";
 VarName[73]="DistProtonMuon";

 VarName[74]="DistShower1Muon";
 VarName[75]="DistShower1Pion";
 VarName[76]="DistShower1Proton";

 VarName[77]="DistShower2Muon";
 VarName[78]="DistShower2Pion";
 VarName[79]="DistShower2Proton";

  TH1F *Hist[Nv][5][Nc];


  for(int s =0;s<5;s++){//Samples
    for(int v=0;v<Nv;v++){//variables
      for(int c=0;c<Nc;c++){
	Hist[v][s][c] = new TH1F(VarName[v]+Names[s]+CutName[c],"",100,-1,-1);
	//cout<<VarName[v]+Names[s]<<endl;
      }
    }
  }

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    bool IsSignal=false;
    int mc_decay = 0;
    int s=-1;//sample
    int muon_index=-1;
    int pion_index=-1;
    int proton_index=-1;
    int shower1_index=-1;
    int shower2_index=-1;
    bool Cuts[Nc]={false};

    bool sel_FV=false;
    bool sel_ntracks=false;
    bool sel_nshowers=false;
    bool sel_muon_candidate=false;
    bool sel_pion_candidate=false;
    bool sel_proton_candidate=false;

    for(int i=0; i< mc_decay_pdg->size();i++){
      //     cout<<i<<"\t"<<mc_decay_pdg->at(i)<<endl;
      
      if(mc_decay_pdg->at(i) == 2212) mc_decay++;
      if(mc_decay_pdg->at(i) == -211) mc_decay++;
    }
    
    if(  mc_nu_pos_x<3.00    && mc_nu_pos_x > 253.35
	 && mc_nu_pos_y<-112.53 && mc_nu_pos_y > 114.47
	 && mc_nu_pos_z<-3.1    && mc_nu_pos_y > 1033.0) continue; //FV
    if(mc_nu_pos_z > 675.1   && mc_nu_pos_y <775.1 ) continue; // dead reagion
    
    //  if(!reco_primary_vtx_inFV) continue; 
    if( (mc_nu_pdg == -14) && (mc_hyperon_pdg ==3212) && (mc_decay==2)) IsSignal=true; 
    
    
    if(IsSignal && jentry<NEvents[0]+1)  s=0;
    if(!IsSignal && jentry<NEvents[0]+1) s=1;
    if(jentry>NEvents[0] && jentry<NEvents[2]+1) s=2;
    if(jentry>NEvents[2] && jentry<NEvents[3]+4) s=3;
    
    
    
    TVector3 vtx(reco_primary_vtx_x, reco_primary_vtx_y, reco_primary_vtx_z);
   
    if(reco_primary_vtx_inFV) sel_FV=true;;

    
    int ntracks=0, nshowers=0;
    bool NoBigShowers=true; 

    for (int i =0; i< pfp_trk_shr_score->size();i++){
	 TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
          double dist = (start-vtx).Mag();

        if(pfp_trk_shr_score->at(i)<.5){
          nshowers++;
          if(pfpnplanehits_Y->at(i)>600) NoBigShowers=false;
               
     
	}
        else{	
          ntracks++;
        }


    }

//Muon Candidate
     for (int i =0; i< pfp_trk_shr_score->size();i++){
         TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
          double dist = (start-vtx).Mag();

         if(pfp_trk_shr_score->at(i)>0.8 && dist<20 && trk_length->at(i)>10 && trk_llrpid->at(i)>0.7)
         {
          sel_muon_candidate=true;
          if(muon_index==-1){
                  muon_index=i;
          }
          else if(trk_length->at(i)>trk_length->at(muon_index)) {
                  muon_index=i;
          }
         }


     }

//Pion Candidate


        for (int i =0; i< pfp_trk_shr_score->size();i++){	
	  if(i==muon_index) continue;
         TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
          double dist = (start-vtx).Mag();

         if(pfp_trk_shr_score->at(i)>0.5/* && dist<20 && trk_length->at(i)>1 && trk_length->at(i)<50 && trk_llrpid->at(i)>-0.6*/)
         {
          sel_pion_candidate=true;
          if(pion_index==-1){
                  pion_index=i;
          }
          else if(trk_llrpid->at(i)>trk_llrpid->at(pion_index)) {
                  pion_index=i;
          }
         }


     }


//Proton Candidate
//
        for (int i =0; i< pfp_trk_shr_score->size();i++){
          if(i==muon_index || i==pion_index) continue;
//	  cout<<"test"<<pfp_trk_shr_score->at(i)<<endl;
         TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
          double dist = (start-vtx).Mag();

         if(pfp_trk_shr_score->at(i)>0.5 /*&& dist<20 && trk_length->at(i)>1 && trk_length->at(i)<50 && trk_llrpid->at(i)<0.1*/)
         {
		           //cout<<"test"<<pfp_trk_shr_score->at(i)<<endl;

          sel_proton_candidate=true;
          if(proton_index==-1){
                  proton_index=i;
          }
          else if(trk_llrpid->at(i)<trk_llrpid->at(proton_index)) {
                  proton_index=i;
          }
         }


     }
	

//Showers
         double shr_dist1=1000;
	 for (int i =0; i< pfp_trk_shr_score->size();i++){
          if(i==muon_index || i==pion_index || i==proton_index) continue;
         
	  TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
          double dist = (start-vtx).Mag();

         if(pfp_trk_shr_score->at(i)<0.5 /*&& dist<20 && trk_length->at(i)>1 && trk_length->at(i)<30 && trk_llrpid->at(i)<0.1*/ )
         {

          if(shower1_index==-1){
                  shower1_index=i;
          }
          else if(dist<shr_dist1 ) {
		  shower2_index=shower1_index;
                  shower1_index=i;

          }
         }


     }
	
   bool sel_more_showers=true;

   if(shower2_index>-1){
	  if (pfpnplanehits_Y->at(shower2_index)>100) sel_more_showers=false;  
   }
    if(nshowers>0 && nshowers<3)sel_nshowers=true;
    if(ntracks<5)sel_ntracks=true;

    bool sel_1st_shower= false;
  if(shower1_index>-1) {
  if (pfpnplanehits_Y->at(shower1_index)>10 && pfpnplanehits_Y->at(shower1_index)<300) sel_1st_shower=true;
  }


    Cuts[0]=true;
    Cuts[1]=Cuts[0] && sel_FV;
    Cuts[2]=Cuts[1] && sel_nshowers;
    Cuts[3]=Cuts[2] && NoBigShowers;
    Cuts[4]=Cuts[3] && sel_muon_candidate;
    Cuts[5]=Cuts[4] && sel_pion_candidate;
    Cuts[6]=Cuts[5] && sel_proton_candidate;
    Cuts[7]=Cuts[6] && sel_ntracks;
    Cuts[8]=Cuts[7] && sel_more_showers;
    Cuts[9]=Cuts[8] && sel_1st_shower;


         TVector3 muonStart;
         TVector3 pionStart;
         TVector3 protonStart;
         TVector3 shower1Start;
         TVector3 shower2Start;

         if(muon_index>-1)    muonStart.SetXYZ(trk_start_x->at(muon_index), trk_start_y->at(muon_index),trk_start_z->at(muon_index));
  	 if(pion_index>-1)    pionStart.SetXYZ(trk_start_x->at(pion_index), trk_start_y->at(pion_index),trk_start_z->at(pion_index));
	 if(proton_index>-1)  protonStart.SetXYZ(trk_start_x->at(proton_index), trk_start_y->at(proton_index),trk_start_z->at(proton_index));
    	 if(shower1_index>-1) shower1Start.SetXYZ(trk_start_x->at(shower1_index), trk_start_y->at(shower1_index),trk_start_z->at(shower1_index));
         if(shower2_index>-1) shower2Start.SetXYZ(trk_start_x->at(shower2_index), trk_start_y->at(shower2_index),trk_start_z->at(shower2_index));

 	 double distMuonPion=99999;
	 double distProtonPion=99999;
	 double distProtonMuon=99999;
	 double distShower1Muon=99999;
	 double distShower1Pion=99999;
	 double distShower1Proton=99999;
   	 double distShower2Muon=99999;
         double distShower2Pion=99999;
         double distShower2Proton=99999;

         if(muon_index>-1    && pion_index>-1)  distMuonPion      = (muonStart-pionStart).Mag();
         if(proton_index>-1  && pion_index>-1)  distProtonPion    = (protonStart  - pionStart).Mag();
         if(proton_index>-1  && muon_index>-1)  distProtonMuon    = (protonStart  - muonStart).Mag();
         if(shower1_index>-1 && muon_index>-1)  distShower1Muon   = (shower1Start - muonStart).Mag();
         if(shower1_index>-1 && pion_index>-1)  distShower1Pion   = (shower1Start - pionStart).Mag();
         if(shower1_index>-1 && proton_index>-1)distShower1Proton = (shower1Start - protonStart).Mag();
         if(shower2_index>-1 && muon_index>-1)  distShower2Muon   = (shower2Start - muonStart).Mag();
         if(shower2_index>-1 && pion_index>-1)  distShower2Pion   = (shower2Start - pionStart).Mag();
         if(shower2_index>-1 && proton_index>-1)distShower2Proton = (shower2Start - protonStart).Mag();

	 //distance cut
	bool sel_distances=false;
	 if(distMuonPion<30 && distProtonMuon<20 && distProtonPion<20
           && distShower1Muon<100 && distShower1Pion<100 & distShower1Proton<100) sel_distances=true;

	 Cuts[10]=Cuts[9]&& sel_distances;

//if(Cuts[10]){}

    for (int c=0;c<Nc;c++){//CHANGE INTO A LOOP over cuts
      
    if(Cuts[c]){
      
      for (int i =0; i< pfp_trk_shr_score->size();i++){
	
         TVector3 muonDir; 
	 TVector3 pionDir;
	 TVector3 protonDir;
	 TVector3 shower1Dir;
	 TVector3 shower2Dir;	
       
	 TVector3 muonStart;
         TVector3 pionStart;
         TVector3 protonStart;
         TVector3 shower1Start;
         TVector3 shower2Start;

	if(muon_index>-1){
		Hist[47][s][c]->Fill(shr_open_angle->at(muon_index), scale[s]);
		muonDir.SetXYZ(trk_dir_x->at(muon_index), trk_dir_y->at(muon_index), trk_dir_z->at(muon_index));
	        muonStart.SetXYZ(trk_start_x->at(muon_index), trk_start_y->at(muon_index),trk_start_z->at(muon_index));

	}
        if(pion_index>-1){
                Hist[48][s][c]->Fill(shr_open_angle->at(pion_index), scale[s]);
                pionDir.SetXYZ(trk_dir_x->at(pion_index), trk_dir_y->at(pion_index), trk_dir_z->at(pion_index));
                pionStart.SetXYZ(trk_start_x->at(pion_index), trk_start_y->at(pion_index),trk_start_z->at(pion_index));

        }
  	if(proton_index>-1){
                Hist[49][s][c]->Fill(shr_open_angle->at(proton_index), scale[s]);
                protonDir.SetXYZ(trk_dir_x->at(proton_index), trk_dir_y->at(proton_index), trk_dir_z->at(proton_index));
                protonStart.SetXYZ(trk_start_x->at(proton_index), trk_start_y->at(proton_index),trk_start_z->at(proton_index));

	}
  	if(shower1_index>-1){
                Hist[50][s][c]->Fill(shr_open_angle ->at(shower1_index), scale[s]);
    		Hist[52][s][c]->Fill(shr_dedx_plane0->at(shower1_index), scale[s]);
                Hist[53][s][c]->Fill(shr_dedx_plane1->at(shower1_index), scale[s]);
                Hist[54][s][c]->Fill(shr_dedx_plane2->at(shower1_index), scale[s]);
                shower1Dir.SetXYZ(trk_dir_x->at(shower1_index), trk_dir_y->at(shower1_index), trk_dir_z->at(shower1_index));
                Hist[67][s][c]->Fill(pfpnplanehits_Y->at(shower1_index), scale[s]);

          TVector3 start(trk_start_x->at(shower1_index), trk_start_y->at(shower1_index),trk_start_z->at(shower1_index));
          double dist = (start-vtx).Mag();
          Hist[69][s][c]->Fill(dist, scale[s]);  //VarName[5]="MuonTrkDistVtx";

	  shower1Start.SetXYZ(trk_start_x->at(shower1_index), trk_start_y->at(shower1_index),trk_start_z->at(shower1_index));

    	}
  	if(shower2_index>-1){
                Hist[51][s][c]->Fill(shr_open_angle ->at(shower2_index), scale[s]);
                Hist[55][s][c]->Fill(shr_dedx_plane0->at(shower2_index), scale[s]);
                Hist[56][s][c]->Fill(shr_dedx_plane1->at(shower2_index), scale[s]);
                Hist[57][s][c]->Fill(shr_dedx_plane2->at(shower2_index), scale[s]);

                shower2Dir.SetXYZ(trk_dir_x->at(shower2_index), trk_dir_y->at(shower2_index), trk_dir_z->at(shower2_index));
                Hist[68][s][c]->Fill(pfpnplanehits_Y->at(shower2_index), scale[s]);
          TVector3 start(trk_start_x->at(shower2_index), trk_start_y->at(shower2_index),trk_start_z->at(shower2_index));
          double dist = (start-vtx).Mag();
          Hist[70][s][c]->Fill(dist, scale[s]);  //VarName[5]="MuonTrkDistVtx";
                shower2Start.SetXYZ(trk_start_x->at(shower2_index), trk_start_y->at(shower2_index),trk_start_z->at(shower2_index));

	}


	 if(muon_index>-1 && pion_index>-1)     Hist[58][s][c]->Fill(muonDir.Angle(pionDir), scale[s]);
         if(proton_index>-1 && pion_index>-1)   Hist[59][s][c]->Fill(protonDir.Angle(pionDir), scale[s]);
         if(proton_index>-1 && muon_index>-1)   Hist[60][s][c]->Fill(protonDir.Angle(muonDir), scale[s]);
         if(shower1_index>-1 && muon_index>-1)  Hist[61][s][c]->Fill(shower1Dir.Angle(muonDir), scale[s]);
         if(shower1_index>-1 && pion_index>-1)  Hist[62][s][c]->Fill(shower1Dir.Angle(pionDir), scale[s]);
         if(shower1_index>-1 && proton_index>-1)Hist[63][s][c]->Fill(shower1Dir.Angle(protonDir), scale[s]);
         if(shower2_index>-1 && muon_index>-1)  Hist[64][s][c]->Fill(shower2Dir.Angle(muonDir), scale[s]);
         if(shower2_index>-1 && pion_index>-1)  Hist[65][s][c]->Fill(shower2Dir.Angle(pionDir), scale[s]);
         if(shower2_index>-1 && proton_index>-1)Hist[66][s][c]->Fill(shower2Dir.Angle(protonDir), scale[s]);

         if(muon_index>-1 && pion_index>-1)     Hist[71][s][c]->Fill((muonStart-pionStart).Mag(), scale[s]);
         if(proton_index>-1 && pion_index>-1)   Hist[72][s][c]->Fill((protonStart  - pionStart).Mag(), scale[s]);
         if(proton_index>-1 && muon_index>-1)   Hist[73][s][c]->Fill((protonStart  - muonStart).Mag(), scale[s]);
         if(shower1_index>-1 && muon_index>-1)  Hist[74][s][c]->Fill((shower1Start - muonStart).Mag(), scale[s]);
         if(shower1_index>-1 && pion_index>-1)  Hist[75][s][c]->Fill((shower1Start - pionStart).Mag(), scale[s]);
         if(shower1_index>-1 && proton_index>-1)Hist[76][s][c]->Fill((shower1Start - protonStart).Mag(), scale[s]);
         if(shower2_index>-1 && muon_index>-1)  Hist[77][s][c]->Fill((shower2Start - muonStart).Mag(), scale[s]);
         if(shower2_index>-1 && pion_index>-1)  Hist[78][s][c]->Fill((shower2Start - pionStart).Mag(), scale[s]);
         if(shower2_index>-1 && proton_index>-1)Hist[79][s][c]->Fill((shower2Start - protonStart).Mag(), scale[s]);


	
	
	
	if(pfp_true_pdg->at(i) ==-13 && trk_length->at(i)>0){ 
	  
	  TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
	  double dist = (start-vtx).Mag();
	  Hist[3][s][c]->Fill(trk_llrpid->at(i), scale[s]); //VarName[3]="MuonPID";
	  Hist[4][s][c]->Fill(trk_length->at(i), scale[s]);	//VarName[4]="MuonTrkLenght";
	  Hist[5][s][c]->Fill(dist, scale[s]);	//VarName[5]="MuonTrkDistVtx";
	  Hist[6][s][c]->Fill(pfp_trk_shr_score->at(i), scale[s]);	//VarName[6]="MuonTrkScore";
	  Hist[19][s][c]->Fill(pfpnplanehits_Y->at(i), scale[s]);
	  Hist[23][s][c]->Fill(pfpnplanehits_U->at(i), scale[s]);
	  Hist[27][s][c]->Fill(pfpnplanehits_V->at(i), scale[s]);      
	  Hist[31][s][c]->Fill(trk_mean_dedx_plane0->at(i),scale[s]);
	  Hist[32][s][c]->Fill(trk_mean_dedx_plane1->at(i),scale[s]);
	  Hist[33][s][c]->Fill(trk_mean_dedx_plane2->at(i),scale[s]);
	  Hist[43][s][c]->Fill(pfp_photon_bdt_score->at(i), scale[s]);
	  
	}
	
	if(pfp_true_pdg->at(i) ==-211 && trk_length->at(i)>0){
	  
	  TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
	  double dist = (start-vtx).Mag();
	  Hist[7][s][c]->Fill(trk_llrpid->at(i), scale[s]); //VarName[3]="MuonPID";
	  Hist[8][s][c]->Fill(trk_length->at(i), scale[s]);      //VarName[4]="MuonTrkLenght";
	  Hist[9][s][c]->Fill(dist, scale[s]);   //VarName[5]="MuonTrkDistVtx";
	  Hist[10][s][c]->Fill(pfp_trk_shr_score->at(i), scale[s]);       //VarName[6]="MuonTrkScore";
	  Hist[20][s][c]->Fill(pfpnplanehits_Y->at(i), scale[s]);       
	  Hist[24][s][c]->Fill(pfpnplanehits_U->at(i), scale[s]);
	  Hist[28][s][c]->Fill(pfpnplanehits_V->at(i), scale[s]);
	  Hist[34][s][c]->Fill(trk_mean_dedx_plane0->at(i),scale[s]);
	  Hist[35][s][c]->Fill(trk_mean_dedx_plane1->at(i),scale[s]);
	  Hist[36][s][c]->Fill(trk_mean_dedx_plane2->at(i),scale[s]);
	  Hist[44][s][c]->Fill(pfp_photon_bdt_score->at(i), scale[s]);
	  
	  
	}
	
	if(pfp_true_pdg->at(i) ==2212 && trk_length->at(i)>0){
	  
	  TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
	  double dist = (start-vtx).Mag();
	  Hist[11][s][c]->Fill(trk_llrpid->at(i), scale[s]); //VarName[3]="MuonPID";
	  Hist[12][s][c]->Fill(trk_length->at(i), scale[s]);      //VarName[4]="MuonTrkLenght";
	  Hist[13][s][c]->Fill(dist, scale[s]);   //VarName[5]="MuonTrkDistVtx";
	  Hist[14][s][c]->Fill(pfp_trk_shr_score->at(i), scale[s]);       //VarName[6]="MuonTrkScore";
	  Hist[21][s][c]->Fill(pfpnplanehits_Y->at(i), scale[s]);       
	  Hist[25][s][c]->Fill(pfpnplanehits_U->at(i), scale[s]);
	  Hist[29][s][c]->Fill(pfpnplanehits_V->at(i), scale[s]);
	  Hist[37][s][c]->Fill(trk_mean_dedx_plane0->at(i),scale[s]);
	  Hist[38][s][c]->Fill(trk_mean_dedx_plane1->at(i),scale[s]);
	  Hist[39][s][c]->Fill(trk_mean_dedx_plane2->at(i),scale[s]);
	  Hist[45][s][c]->Fill(pfp_photon_bdt_score->at(i), scale[s]);
	  
	  
	}
	
	if(pfp_true_pdg->at(i) ==22 && trk_length->at(i)>0){
	  
	  TVector3 start(trk_start_x->at(i), trk_start_y->at(i),trk_start_z->at(i));
	  double dist = (start-vtx).Mag();
	  Hist[15][s][c]->Fill(trk_llrpid->at(i), scale[s]); //VarName[3]="MuonPID";
	  Hist[16][s][c]->Fill(trk_length->at(i), scale[s]);      //VarName[4]="MuonTrkLenght";
	  Hist[17][s][c]->Fill(dist, scale[s]);   //VarName[5]="MuonTrkDistVtx";
	  Hist[18][s][c]->Fill(pfp_trk_shr_score->at(i), scale[s]);       //VarName[6]="MuonTrkScore";
	  Hist[22][s][c]->Fill(pfpnplanehits_Y->at(i), scale[s]);       
	  Hist[26][s][c]->Fill(pfpnplanehits_U->at(i), scale[s]);
	  Hist[30][s][c]->Fill(pfpnplanehits_V->at(i), scale[s]);
	  Hist[40][s][c]->Fill(trk_mean_dedx_plane0->at(i),scale[s]);
	  Hist[41][s][c]->Fill(trk_mean_dedx_plane1->at(i),scale[s]);
	  Hist[42][s][c]->Fill(trk_mean_dedx_plane2->at(i),scale[s]);
	  Hist[46][s][c]->Fill(pfp_photon_bdt_score->at(i), scale[s]);
	  
	}
	
	
	
	
      }
      
    // No shower with number of hits in Y plan more that 300
      if(nshowers>0){
	for (int i =0; i< pfp_trk_shr_score->size();i++){
	  if(pfp_trk_shr_score->at(i)<.5) Hist[2][s][c]->Fill(pfpnplanehits_Y->at(i),scale[s]);
	  
	}
	
      }
      
      Hist[0][s][c]->Fill(ntracks,scale[s]); //NTracks
      Hist[1][s][c]->Fill(nshowers,scale[s]); //NShowers
      
    }
    }//END OF Loop over CUTS   
    
  }
  
cout<<"cut \t\t\t Sig \t hyp \t backg \t dirt \t eff[%] \t purity[%] "<<endl;
 
for (int c=0;c<Nc;c++){
double sig=Hist[0][0][c]->Integral();
double back=Hist[0][1][c]->Integral() + Hist[0][2][c]->Integral()+ Hist[0][3][c]->Integral();

cout<<CutName[c]<<"\t \t"<<Hist[0][0][c]->Integral()<<"\t"<<Hist[0][1][c]->Integral()<< "\t"<<Hist[0][2][c]->Integral()<< "\t"<<Hist[0][3][c]->Integral()
	<<"\t"<<sig/Hist[0][0][0]->Integral()*100 <<"\t"<< sig/back*100<<endl;
}


  TCanvas *c2 = new TCanvas("c2","",5000,4000);
  
  for( int v=0;v<Nv;v++){
    c2->Divide(Nc,4);// 4 rows and 5 columns
    
    for (int c=0; c<Nc; c++){
      for(int s=0;s<4;s++){
	c2->cd(c+1+s*Nc);
	Hist[v][s][c]->Draw();
      }
    }
    c2->Print("figs/"+VarName[v]+".png");
    c2->Clear();
  }
  
  

  
}
