// drawRaa.C
// Author: Yeonju Go 

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
#include "../phoRaaCuts/phoRaaCuts_v1.h"
#include "../ElectroWeak-Jet-Track-Analyses/Histogramming/PhotonPurity.h"
#include "../ElectroWeak-Jet-Track-Analyses/Utilities/interface/CutConfigurationParser.h"
#include "../ElectroWeak-Jet-Track-Analyses/TreeHeaders/CutConfigurationTree.h"

const TString LABEL = "PbPb #sqrt{s}_{_{NN}}=5.02 TeV";

void drawRaa()
{
    TH1::SetDefaultSumw2();
/*    TString pbpbEff_fname[nCentBin], ppEff_fname;
    pbpbEff_fname[0] = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pbpb_cent0-60_v1_nohotspot.root";
    pbpbEff_fname[1] = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pbpb_cent60-200_v1_nohotspot.root";
    ppEff_fname = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pp_v1_nohotspot.root";
  */  ////// values
    double lumi_pp = 27.87*1e12; //27.87 pb-1
    double lumi_pbpb = 404*1e6; // 0.56 nb-1
    double nColl[nCentBin] = {(1626.0+1005.0+606.4)/3.,  (348.3+186.2+90.69+40.14+15.87+5.502+1.642)/7.};//version1 
    double nMB[nCentBin]; // ?*pbpb lumi?
    double TA[nCentBin]; //TAA=<ncoll>/sigma_{inel,pp} =<ncoll> /(70 mb)
    for(int j=0;j<nCentBin;j++){
        TA[j] = nColl[j]/(70e-3);  
    }
    nMB[0] = 7.75*lumi_pbpb*0.3;//*0.3; // ?*pbpb lumi?
    nMB[1] = 7.75*lumi_pbpb*0.7;//*0.7; // ?*pbpb lumi?
    //http://dde.web.cern.ch/dde/glauber_lhc.htm
    //double nColl[nCentBin] = {(1626.0+1005.0+606.4)*(0.1/392.5),  (348.3+186.2+90.69+40.14+15.87+5.502+1.642)*(0.1/392.5)};//version2
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

    TH1D* h1D_corr[3][nCentBin]; // corrected yield 1. purity, 2. efficiency 
    TH1D* h1D_corrpp[3]; // corrected yield 1. purity, 2. efficiency 
    TH1D* h1D_Raa[nCentBin]; // total Raa!!

    ////// Set purity values. 
    cout << "nPtBin = " << nPtBin << endl;
    double val_purpp[nPtBin] = {0.91,0.89,0.88,0.84};
    double val_purpbpb[nCentBin][nPtBin] = {{0.53,0.51,0.47,0.38},{0.53,0.41,0.24,0.22}};

    ////// Get efficiency & purity hist. 
    ////// Define Raa histogram.
    for(int j=0;j<nCentBin;j++){
        if(j==0) eff[j] = new TFile(pbpbEff_fname_cent1);
        else if(j==1) eff[j] = new TFile(pbpbEff_fname_cent2);
        h1D_eff[j] = (TH1D*) eff[j]->Get("h1D_Eff_pt_total");
        h1D_eff[j] ->SetName(Form("h1D_eff_cent%d",j)); 
        h1D_pur[j] = (TH1D*) h1D_eff[j]->Clone(Form("h1D_pur_cent%d",j));
        h1D_Raa[j] = (TH1D*) h1D_eff[j]->Clone(Form("h1D_Raa_cent%d",j));
        if(j==0){
            effpp = new TFile(ppEff_fname); 
            h1D_effpp = (TH1D*) effpp->Get("h1D_Eff_pt_total");
            h1D_effpp -> SetName("h1D_eff_pp"); 
            h1D_purpp = (TH1D*) h1D_effpp ->Clone(Form("h1D_pur_pp"));
        }
        for(int i=0;i<nPtBin;i++){
            cout << "pbpb purity i =" << i << ", val = " << val_purpbpb[j][i]<<endl;
            h1D_pur[j]->SetBinContent(i+1,val_purpbpb[j][i]);
            cout << "pp purity i =" << i << ", val = " << val_purpp[i]<<endl;
            if(j==0) h1D_purpp->SetBinContent(i+1,val_purpp[i]);
        }
    }

    ////// Define RAW histogram.
    for(int j=0;j<nCentBin;j++){
        h1D_raw[j] = new TH1D(Form("h1D_raw_cent%d",j),"",nPtBin, ptBins);
        if(j==0) h1D_rawpp = (TH1D*) h1D_effpp ->Clone(Form("h1D_raw_pp"));
    } 

    ////// DRAW RAW DATA
    for(int j=0;j<nCentBin;j++){
        TCut etaCut = "abs(phoEta) <=1.44";
        TCut ptCut = "phoEt>=40"; 
        TCut centCut = Form("hiBin>=%d && hiBin < %d",centBins[j],centBins[j+1]); 
        TCut kineCut = etaCut && ptCut;
        cout << "################ pbpb filling in cent "<< j <<  endl;
        tpbpb->Draw(Form("phoEt>>%s",h1D_raw[j]->GetName()),dataCut && phoSignalCut && kineCut && centCut);
        //tpbpb->Draw(Form("phoEt>>%s",h1D_raw[j]->GetName()));
        //h1D_raw[j]=(TH1D*)gDirectory->Get(h1D_raw[j]->GetName()); 
        cout << "################ pbpb done in cent "<< j<< endl;
        if(j==0){
            cout << "################ pp filling"<< endl;
            tpp->Draw(Form("phoEt>>h1D_raw_pp"),dataCut_pp && phoSignalCut_ppGED && kineCut); 
            cout << "################ pp done"<< endl;
        }
    }
    ////// Raa calculation. 
    for(int j=0;j<nCentBin;j++){
        h1D_raw[j]->Scale(1.,"width");
        for(int ii=0;ii<3;ii++){
            h1D_corr[ii][j] = new TH1D(Form("h1D_dNdpt_corr%d_cent%d",ii,j),";p_{T}^{#gamma} (GeV); dN/dp_{T} (GeV^{-1})",nPtBin,ptBins);
            if(j==0) h1D_corrpp[ii] = new TH1D(Form("h1D_dNdpt_corr%d_pp",ii),";p_{T}^{#gamma} (GeV); dN/dp_{T} (GeV^{-1})",nPtBin,ptBins);
        }
        h1D_corr[0][j]-> Multiply(h1D_raw[j], h1D_pur[j]);
        h1D_corr[1][j]-> Divide(h1D_corr[0][j], h1D_eff[j]);
        if(j==0) {
            h1D_rawpp->Scale(1.,"width");
            h1D_corrpp[0]-> Multiply(h1D_rawpp, h1D_purpp);
            h1D_corrpp[1]-> Divide(h1D_corrpp[0], h1D_effpp);
        }
    }

    ///// Draw Raa
    cout << "!!DRAW RAA!!" << endl;
    TCanvas* c1 = new TCanvas("craa","",800,400);
    c1->Divide(2,1);
    for(int j=0;j<nCentBin;j++){
        c1->cd(j+1);
        h1D_Raa[j]->Divide(h1D_corr[1][j],h1D_corrpp[1]);
        h1D_Raa[j]->Scale(lumi_pp/(TA[j]*nMB[j]));
        //h1D_Raa[j]->Scale(lumi_pp/lumi_pbpb/(double)(208*208));
        h1D_Raa[j]->Scale(1./nColl[j]);
        if(j==0) h1D_Raa[j]->SetMarkerColor(2);
        h1D_Raa[j]->SetMarkerStyle(20);

        h1D_Raa[j]->SetTitle(Form("%d-%d %s;isolated p_{T}^{#gamma} (GeV);R_{AA}",centBins[j]/2,centBins[j+1]/2,"%"));
        h1D_Raa[j]->SetAxisRange(0.0,1.2,"Y");
        h1D_Raa[j]->Draw();
        //else h1D_Raa[j]->Draw("same");
    }
    c1->SaveAs("pdf/phoRaa_v1.pdf");

    TFile* outf = new TFile("phoRaa_v1.root","recreate");
    outf->cd();
    for(int j=0;j<nCentBin;j++){
        h1D_raw[j]->Write();
        h1D_eff[j]->Write();
        h1D_pur[j]->Write();
        h1D_Raa[j]->Write();
        for(int ii=0;ii<2;ii++){
            h1D_corr[ii][j]->Write();
        }
    }
    h1D_rawpp->Write();
    h1D_effpp->Write();
    h1D_purpp->Write();
    for(int ii=0;ii<2;ii++){
        h1D_corrpp[ii]->Write();
    }

}
