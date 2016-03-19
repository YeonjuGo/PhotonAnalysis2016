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
#include <stdlib.h>
using namespace std;

void xSecCal(){
    int pthatCut[] = {15,30,50,80,120,9999};
    const int nFile = sizeof(pthatCut)/sizeof(int)-1;
    
    //pbpb
    const char *fileName[nFile] = {
        "/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC_AllQCDPhoton15_v1/0.root",
        "/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC_AllQCDPhoton30_v1/0.root",
        "/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC_AllQCDPhoton50_v1/0.root",
        "/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC_AllQCDPhoton80_v1/0.root",
        "/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC_AllQCDPhoton120_v1/0.root",
   };

    TFile *fin[nFile];
    TTree *t[nFile];
    int evtN[nFile];
    float pEnt[nFile][nFile];
    float totalEnt[nFile];
    for(int i=0; i<nFile ; i++){
        totalEnt[i]=0.0;
        for(int j=0; j<nFile ; j++){
            pEnt[i][j]=0.0;
        }
    }

    for(int ifile=0; ifile<nFile; ifile++){
        fin[ifile] = new TFile(fileName[ifile]);
        t[ifile] = (TTree*) fin[ifile] -> Get("hiEvtAnalyzer/HiTree");
        Float_t pthat;
        TBranch *b_pthat;
        t[ifile]->SetBranchAddress("pthat",&pthat, &b_pthat);
        for(int jj=0;jj<nFile;jj++){
            pEnt[ifile][jj] = t[ifile]->GetEntries(Form("pthat>= %d && pthat< %d", pthatCut[jj], pthatCut[jj+1]));
        }
    }

    //
    // Calculate weighting factor
    //
    double weight[nFile];

    weight[0] = 1.0;
    weight[1] = pEnt[0][1]/(pEnt[0][1]*weight[0]+pEnt[1][1]);
    weight[2] = pEnt[0][2]/(pEnt[0][2]*weight[0]+pEnt[1][2]*weight[1]+pEnt[2][2]);
    weight[3] = pEnt[0][3]/(pEnt[0][3]*weight[0]+pEnt[1][3]*weight[1]+pEnt[2][3]*weight[2]+pEnt[3][3]);
    weight[4] = pEnt[0][4]/(pEnt[0][4]*weight[0]+pEnt[1][4]*weight[1]+pEnt[2][4]*weight[2]+pEnt[3][3]*weight[3]+pEnt[4][4]);
     
    
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
}
