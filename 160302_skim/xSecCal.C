#include "TROOT.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TSystem.h"
#include "TPad.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "../phoRaaCuts/phoRaaCuts_v2.h"
#include <stdlib.h>
using namespace std;

void xSecCal(sampleType colli=kPPMC){
//void xSecCal(string coll = "pp", string channel = "AllQCDPhoton"){
    int nPthat = 5; 
    const int nFile = 2;
    const char *fileName[nFile];

    if(colli==kPPMCEmEnr) nPthat=4;
    const char* lowestPthatFileName;
    float pthatCut[nPthat];

    if(colli==kPPMC) {
        float temp[] = {15,30,50,80,120,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = temp[j];     
        }
//        fileName[0] = "/home/goyeonju/CMS/Files/photon2016/officialMC/HINppWinter16DR/Pythia8_Photon15_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1.root";
//        fileName[1]= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINppWinter16DR/Pythia8_Photon15_30_50_80_120_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1.root";
        fileName[0] = "/home/goyeonju/CMS/Files/photon2016/gsfs-Pythia8_Photon_pp_RECO_forest_v28/gsfs-Pythia8_Photon15_pp_RECO_forest_v28.root";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/gsfs-Pythia8_Photon_pp_RECO_forest_v28/gsfs-Pythia8_Photon15_30_50_80_120_pp_RECO_forest_v28.root";
   } else if(colli==kHIMC) {
        float temp[] = {15,30,50,80,120,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = temp[j];     
        }
        //fileName[0] = "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_Photon15_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
       // fileName[1]= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_Photon15_30_50_80_120_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        fileName[0] ="/home/goyeonju/CMS/Files/photon2016/Pythia8_Photon_Hydjet_RECO_20160306_forest_v28_2/Pythia8_Photon15_Hydjet_RECO_20160306_forest_v28_2.root";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/Pythia8_Photon_Hydjet_RECO_20160306_forest_v28_2/Pythia8_Photon15_30_50_80_120_Hydjet_RECO_20160306_forest_v28_2.root";
 
    } else if(colli==kPPMCEmEnr) {
        float temp[] = {50,80,120,170,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = temp[j];         
        }
        //fileName[0] = "Pythia8_EmEnrDijet50_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_forest_v1.root";
        //fileName[1]= "Pythia8_EmEnrDijet50_80_120_170_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_forest_v1.root";
        fileName[0] ="/home/goyeonju/CMS/Files/photon2016/";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/";
    } else if(colli==kHIMCEmEnr){
        float tmpPthat[] = {30,50,80,120,170,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = tmpPthat[j];
        }
        //fileName[0] = "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_EmEnrDijet30_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        //fileName[1]= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_EmEnrDijet30_50_80_120_170_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        fileName[0] ="/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC_Pythia8_EmEnrichedDijet/2015-PbPb-MC_Pythia8_EmEnrichedDijet30.root";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC_Pythia8_EmEnrichedDijet/2015-PbPb-MC_Pythia8_EmEnrichedDijet30_50_80_120_170.root";
    }
/*
    if(coll =="pp" && channel =="AllQCDPhoton"){
        pthatCut[nPthat] = {15,30,50,80,120,9999};
        fileName[0] = "/home/goyeonju/CMS/Files/photon2016/gsfs-Pythia8_Photon_pp_RECO_forest_v28/gsfs-Pythia8_Photon15_pp_RECO_forest_v28.root";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/gsfs-Pythia8_Photon_pp_RECO_forest_v28/gsfs-Pythia8_Photon15_30_50_80_120_pp_RECO_forest_v28.root";
    } else if(coll =="pbpb" && channel =="AllQCDPhoton"){
        fileName[0] ="/home/goyeonju/CMS/Files/photon2016/Pythia8_Photon_Hydjet_RECO_20160306_forest_v28_2/Pythia8_Photon15_Hydjet_RECO_20160306_forest_v28_2.root";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/Pythia8_Photon_Hydjet_RECO_20160306_forest_v28_2/Pythia8_Photon15_30_50_80_120_Hydjet_RECO_20160306_forest_v28_2.root";
    } else if(coll =="pp" && channel =="EmEnrichedDijet"){
        fileName[0] ="/home/goyeonju/CMS/Files/photon2016/";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/";
    } else if(coll =="pbpb" && channel =="EmEnrichedDijet"){
        fileName[0] ="/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC_Pythia8_EmEnrichedDijet/2015-PbPb-MC_Pythia8_EmEnrichedDijet30.root";
        fileName[1] ="/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC_Pythia8_EmEnrichedDijet/2015-PbPb-MC_Pythia8_EmEnrichedDijet30_50_80_120_170.root";
    }
*/
//    cout << "pthatWeight calculation of " << getSampleTypeName(coll) << " " << channel << endl;

    TFile *fin[nFile];
    TTree *t[nFile];
    float entries[nFile][nPthat];
    for(int i=0; i<nFile ; i++){
        for(int j=0; j<nPthat ; j++){
            entries[i][j]=0.0;
        }
    }
    
    for(int ifile=0; ifile<nFile; ifile++){
        fin[ifile] = new TFile(fileName[ifile]);
        t[ifile] = (TTree*) fin[ifile] -> Get("hiEvtAnalyzer/HiTree");
        Float_t pthat;
        TBranch *b_pthat;
        t[ifile]->SetBranchAddress("pthat",&pthat, &b_pthat);
        for(int j=0;j<nPthat;j++){
            entries[ifile][j] = t[ifile]->GetEntries(Form("pthat>= %f && pthat< %f", pthatCut[j], pthatCut[j+1]));
        }
    }

    //
    // Calculate weighting factor
    //
    double weight[nPthat];
    for(int j=0;j<nPthat;j++){
        weight[j] = entries[0][j]/entries[1][j];
        cout << "pthat "<< pthatCut[j] << " weight = " << weight[j] << endl;
    }
   /* 
    TCanvas *can = new TCanvas("can", "can", 300,300);
    TH1D *hpthat[nFile];
    TH1D *hpthat_total;
    for(int ifile=0; ifile<nFile; ifile++){
        hpthat[ifile] = new TH1D(Form("hpthat%d",ifile),"", 210, 10, 1000);
        //t[ifile]-> SetWeight(weight[ifile]);
        t[ifile]-> Draw("pthat>>+hpthat[ifile]");
        hpthat[ifile] = (TH1D*)gDirectory->Get("hpthat[ifile]");
    }
    hpthat_total = (TH1D*) hpthat[0]->Clone("hpthat_total");
    hpthat_total->Add(hpthat[1]);
    hpthat_total->Add(hpthat[2]);
    hpthat_total->Add(hpthat[3]);
    hpthat_total->Add(hpthat[4]);
    
    hpthat_total->Draw();
*/
}
