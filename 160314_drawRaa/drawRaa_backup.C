// drawRaa.C
// Author: Yeonju Go 
// 

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCut.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "stdio.h"
#include <iostream>
#include "../yjUtility.h"
#include "../phoRaaCuts_v1.h"
#include "../ElectroWeak-Jet-Track-Analyses/Histogramming/PhotonPurity.h"
#include "../ElectroWeak-Jet-Track-Analyses/Utilities/interface/CutConfigurationParser.h"
#include "../ElectroWeak-Jet-Track-Analyses/TreeHeaders/CutConfigurationTree.h"
//#include "../ElectroWeak-Jet-Track-Analyses/Plotting/commonUtility.h"

//////////////////////// Cuts & Bins ///////////////////////////
const TString LABEL = "PbPb #sqrt{s}_{_{NN}}=5.02 TeV";
const TCut sampleIsolation = "(pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20) < 1.0 && phoHoverE<0.1";

// last entry is upper bound on last bin
const Int_t CENTBINS[] = {0,60, 200};
//const Int_t CENTBINS[] = {0, 200};
const Int_t nCENTBINS = sizeof(CENTBINS)/sizeof(Int_t) -1;

const Double_t PTBINS[] = {40,50,60,70,80};
const Int_t nPTBINS = sizeof(PTBINS)/sizeof(Double_t) -1;

const Double_t ETABINS[] = {-1.44, 1.44};
//const Double_t ETABINS[] = {-1.44, -1, -0.5, 0, 0.5, 1, 1.44};
const Int_t nETABINS = sizeof(ETABINS)/sizeof(Double_t) -1;



void drawRaa()
    //void drawRaa(const TString configFile, const TString inputData, const TString inputMC, const TString outputName, const TString coll="pbpb")
{
    TH1::SetDefaultSumw2();
    //CutConfiguration config = CutConfigurationParser::Parse(configFile.Data());
    //TTree *configTree = setupConfigurationTreeForWriting(config);
    TString pbpbEff_fname[nCentBin], ppEff_fname;
    pbpbEff_fname[0] = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pbpb_cent0-60_v1_nohotspot.root";
    pbpbEff_fname[1] = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pbpb_cent60-200_v1_nohotspot.root";
    ppEff_fname = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pp_v1_nohotspot.root";
    ////// values
    double lumi_pp = 27.87 * 1000; //nb-1
    double lumi_pbpb = 0.404; // 0.56 nb-1
    double nColl[nCentBin] = {(1626.0+1005.0+606.4)/3.,  (348.3+186.2+90.69+40.14+15.87+5.502+1.642)/7.}; 
    TFile *f1 = TFile::Open(pbpbDatafname);
    TTree *tpbpb = (TTree*)f1->Get("EventTree");
    TTree *tpbpb_hi = (TTree*)f1->Get("HiTree");
    TTree *tpbpb_skim = (TTree*)f1->Get("HltTree");
    tpbpb->AddFriend(tpbpb_hi);
    tpbpb->AddFriend(tpbpb_skim);

    TFile *f2 = TFile::Open(ppDatafname);
    TTree *tpp = (TTree*)f2->Get("ggHiNtuplizerGED/EventTree");
    TTree *tpp_hi = (TTree*)f2->Get("hiEvtAnalyzer/HiTree");
    TTree *tpp_skim = (TTree*)f2->Get("skimanalysis/HltTree");
    tpp->AddFriend(tpp_hi);
    tpp->AddFriend(tpp_skim);


    TFile* eff[nCentBin];
    TFile* effpp;

    TH1D* h1D_eff[nCentBin];
    TH1D* h1D_effpp;
    TH1D* h1D_pur[nCentBin];
    TH1D* h1D_purpp;
    TH1D* h1D_raw[nCentBin];
    TH1D* h1D_rawpp;

    TH1D* h1D_test = new TH1D("hhh","",4,40,80);
    TH1D* h1D_corr[3][nCentBin]; // corrected yield 1. purity, 2. efficiency 
    TH1D* h1D_corrpp[3]; // corrected yield 1. purity, 2. efficiency 
    TH1D* h1D_Raa[nCentBin]; // total Raa!!

    ////// Set purity values. 
    cout << "nPtBin = " << nPtBin << endl;
    double val_purpp[nPtBin] = {0.91,0.89,0.88,0.84};
    double val_purpbpb[nCentBin][nPtBin] = {{0.53,0.51,0.47,0.38},{0.53,0.41,0.24,0.22}};
    /*
    for(int j=0;j<nCentBin;j++){
        for(int i=0;i<nPtBin;i++){
            cout << "pbpb purity i =" << i << ", val = " << val_purpbpb[j][i]<<endl;
            if(j==0) cout << "pp purity i =" << i << ", val = " << val_purpp[i]<<endl;
        }
    }
    */
    ////// Get efficiency & purity hist. 
    ////// Define RAW, Raa histogram.
    for(int j=0;j<nCentBin;j++){
        eff[j] = new TFile(pbpbEff_fname[j]);
        h1D_eff[j] = (TH1D*) eff[j]->Get("h1D_Eff_pt_total");
        h1D_eff[j] ->SetName(Form("h1D_eff_cent%d",j)); 
        //h1D_raw[j] = (TH1D*) h1D_eff[j]->Clone(Form("h1D_raw_cent%d",j));
        h1D_raw[j] = new TH1D(Form("h1D_raw_cent%d",j),"hraw",4,40,80);
        //h1D_raw[j] = new TH1D(Form("h1D_raw_cent%d",j),"hraw",nPTBINS,PTBINS);
        //h1D_raw[j] = new TH1D(Form("h1D_raw_cent%d",j),"hraw",nPtBin,ptBins);
        h1D_pur[j] = (TH1D*) h1D_eff[j]->Clone(Form("h1D_pur_cent%d",j));
        h1D_Raa[j] = (TH1D*) h1D_eff[j]->Clone(Form("h1D_Raa_cent%d",j));
        if(j==0){
            effpp = new TFile(ppEff_fname); 
            h1D_effpp = (TH1D*) effpp->Get("h1D_Eff_pt_total");
            h1D_effpp -> SetName("h1D_eff_pp"); 
            h1D_rawpp = (TH1D*) h1D_effpp ->Clone(Form("h1D_raw_pp"));
            h1D_rawpp->Reset();
            h1D_purpp = (TH1D*) h1D_effpp ->Clone(Form("h1D_pur_pp"));
        }
        for(int i=0;i<nPtBin;i++){
            cout << "pbpb purity i =" << i << ", val = " << val_purpbpb[j][i]<<endl;
            h1D_pur[j]->SetBinContent(i+1,val_purpbpb[j][i]);
            cout << "pp purity i =" << i << ", val = " << val_purpp[i]<<endl;
            if(j==0) h1D_purpp->SetBinContent(i+1,val_purpp[i]);
        }
    }
    /*
       TCanvas *cpur = new TCanvas();
       h1D_pur[0]->Draw();
       */  


    ////// DRAW RAW DATA
    for(int j=0;j<1;j++){
        //for(int j=0;j<nCentBin;j++){
        TCut etaCut = "abs(phoEta) <=1.44";
        TCut ptCut = "phoEt>=40"; 
        cout << "################ pbpb filling in cent "<< j <<  endl;
        //tpbpb->Draw(Form("phoEt>>htest(4,40,80)"),dataCut && phoSignalCut);
        //tpbpb->Draw(Form("phoEt>>htest(4,40,80)"));
        tpbpb->Draw(Form("phoEt>>%s",h1D_test->GetName()));
        //tpbpb->Draw(Form("phoEt>>%s",h1D_raw[j]->GetName()));
        //tpbpb->Draw(Form("phoEt>>%s",h1D_raw[j]->GetName()));
        //tpbpb->Draw(Form("phoEt>>%s",h1D_raw[j]->GetName()),dataCut && phoSignalCut);
        //h1D_raw[j]=(TH1D*)gDirectory->Get(Form("htest")); 
        //h1D_raw[j]=(TH1D*)gDirectory->Get(h1D_raw[j]->GetName()); 
        //h1D_raw[j]=(TH1D*)gDirectory->Get(h1D_raw[j]->GetName()); 
        cout << "################ pbpb done in cent "<< j<< endl;
        if(j==10){
            cout << "################ pp filling"<< endl;
            tpp->Draw(Form("phoEt>>h1D_raw_pp"),dataCut_pp && phoSignalCut_ppGED); 
            cout << "################ pp done"<< endl;
        }
    }
    TCanvas *craw = new TCanvas("craw","",400,400);
    h1D_raw[0]->Draw();
#if 0 
    TCanvas *ccorr0 = new TCanvas("ccorr0","",400,400);


    ////// Raa calculation. 
    for(int j=0;j<nCentBin;j++){
        h1D_raw[j]->Scale(1.,"width");
        cout << "zzkkk" << endl;
        h1D_corr[0][j]-> Multiply(h1D_raw[j], h1D_pur[j]);
        ccorr0->cd();
        h1D_corr[0][j]->Draw();
        cout << "zzkkk" << endl;
        h1D_corr[1][j]-> Divide(h1D_corr[0][j], h1D_eff[j]);
        cout << "zzkk" << endl;
        if(j==0) {
            h1D_rawpp->Scale(1.,"width");
            h1D_corrpp[0]-> Multiply(h1D_rawpp, h1D_purpp);
            h1D_corrpp[1]-> Divide(h1D_corrpp[0], h1D_effpp);
        }
    }

    ///// Draw Raa
    cout << "cc" << endl;
    TCanvas* c1 = new TCanvas("craa","",800,800);

    for(int j=0;j<nCentBin;j++){
        cout << "dd" << endl;
        h1D_Raa[j]->Divide(h1D_corr[1][j],h1D_corrpp[1]);
        h1D_Raa[j]->Scale(lumi_pp/lumi_pbpb/(double)(208*208));
        h1D_Raa[j]->Scale(1./nColl[j]);
        if(j==0) h1D_Raa[j]->SetMarkerColor(2);
        //h1D_Raa[j]->SetMarkerColor(4);
        h1D_Raa[j]->SetMarkerStyle(20);

        h1D_Raa[j]->SetTitle(Form("%d-%d %s;isolated p_{T}^{#gamma} (GeV);R_{AA}",centBins[j]/2,centBins[j+1]/2,"%"));
        if(j==0) h1D_Raa[j]->Draw();
        else h1D_Raa[j]->Draw("same");
    }
    c1->SaveAs("pdf/phoRaa_v1.pdf");
#endif
}
