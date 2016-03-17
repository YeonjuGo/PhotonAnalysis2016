// Author : Yeonju Go

#include "../yjUtility.h"

const int colhere[] = {1,2,4};
TH1D* cutoffHist(TH1* h=0, float xLow=40.0, float xHigh=200.0, int nbinnew=16);

void compare_phoPt_withJetphox(TString coll="pbpb"){
    TH1::SetDefaultSumw2();

    const char* fname_data;
    const char* fname_hydjet;
    const char* fname_jetphox;
    fname_data = "/home/goyeonju/CMS/Files/photon2016/forestSkimed_photonSkim_pbpb_2015data.root";
    fname_hydjet = "/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC-AllQCDPhoton_v1_ptHat_15_30_50_80_120_with_pthatWeight.root";
    fname_jetphox = "/home/goyeonju/CMS/Files/photon2016/ggdpbpb5tev_150314.root";

    const TCut spikeRejection = "(phoSigmaIEtaIEta_2012>=0.002) && (pho_swissCrx<=0.9) && (abs(pho_seedTime)<=3)";
    const TCut isoCut = "(pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20) < 1.0 && phoHoverE<0.1";
    const TCut kineCut = "phoEt>40 && abs(phoEta)< 1.44";
    const TCut hotspotCut = "!(((phoSCEta > -1.194 && phoSCEta < -1.169) && (phoSCPhi > -1.211 && phoSCPhi < -1.162)) || ((phoSCEta > -0.935 && phoSCEta < -0.910) && (phoSCPhi > 2.036 && phoSCPhi < 2.065)) || ((phoSCEta > -0.850 && phoSCEta < -0.819) && (phoSCPhi > -1.723 && phoSCPhi < -1.676)) || ((phoSCEta > -0.764 && phoSCEta < -0.723) && (phoSCPhi > -0.786 && phoSCPhi < -0.720)) || ((phoSCEta > -0.237 && phoSCEta < -0.210) && (phoSCPhi > 0.106 && phoSCPhi < 0.147)) ||((phoSCEta > -0.062 && phoSCEta < -0.036) && (phoSCPhi > 0.860 && phoSCPhi< 0.931)) || ((phoSCEta > 0.973 && phoSCEta < 1.035) && (phoSCPhi > 0.623 && phoSCPhi < 0.697)) || ((phoSCEta > 1.075 && phoSCEta < 1.108) && (phoSCPhi > -3.114 && phoSCPhi < -3.054)) || ((phoSCEta > 1.426 && phoSCEta < 1.447) && (phoSCPhi > 2.744 && phoSCPhi < 2.774)) || ((phoSCEta > 0.915 && phoSCEta < 0.930) && (phoSCPhi > 1.69 && phoSCPhi < 1.72)))";
    const TCut dataTotCut = spikeRejection && isoCut && kineCut && hotspotCut &&"pcollisionEventSelection";
    /*
    TFile* f1 = new TFile(fname_data);
    TTree* t1 = (TTree*) f1 -> Get("EventTree");
    TTree* t1_hlt = (TTree*) f1 -> Get("HltTree");
    t1->AddFriend(t1_hlt);
    t1->Draw("phoEt>>+hdata(17,30,200)",dataTotCut);
    TH1D* hdata=(TH1D*)gDirectory->Get("hdata");
    hdata->Scale(1./(hdata->Integral()));
    hdata->SetMarkerStyle(20);
*/

    const TCut mcIsolation = "mcPt>40 && abs(mcEta)<1.44 && mcCalIsoDR04<5 && abs(mcPID)==22 && abs(mcMomPID)<=22";
    TFile* f2 = new TFile(fname_hydjet);
    TTree* t2 = (TTree*) f2 -> Get("EventTree");
    TTree* t2_hi = (TTree*) f2 -> Get("HiTree");
    t2->AddFriend(t2_hi);
    t2->Draw("mcPt>>hmc(16,40,200)","pthatWeight"*mcIsolation);
    TH1D* hmc =(TH1D*)gDirectory->Get("hmc");
    hmc->Scale(1./hmc->Integral());
    hmc->SetLineColor(2);

    TFile* f3 = new TFile(fname_jetphox);
    TH1D* hpt_jetphox = (TH1D*) f3->Get("hp20"); 
    hpt_jetphox->Rebin(10);
    TH1D* hjetphox = cutoffHist(hpt_jetphox,40.0,200.0,16);
    hjetphox->Scale(1./hjetphox->Integral());
    hjetphox->SetLineColor(4);

    ////////////////////////////////////////
    // DRAW 
    ////////////////////////////////////////
    TH1D* htemp = new TH1D(Form("htemp"), ";p_{T}^{#gamma};Arbitrary Normalisation",17,30,200);
    TCanvas* c = new TCanvas("c","", 400,400);
    htemp->Draw();
    gPad->SetLogy();
    //hdata->Draw("pe same");
    hmc->Draw("hist e same");
    hjetphox->Draw("hist e same");

} // main function

TH1D* cutoffHist(TH1* h, float xLow, float xHigh, int nbinNew){
    int nbins = h->GetNbinsX();
    cout << "nbinsx = " << nbins << endl;
    float val[nbinNew];
    float valerr[nbinNew];
    int tmpnewbin=1;
    TH1D* hnew = new TH1D("hnew","",nbinNew,xLow,xHigh);
    for(int i=0;i<nbins;i++){
        if((h->GetBinLowEdge(i))<xLow) continue;
        if((h->GetBinLowEdge(i))>=xHigh) continue;
        val[tmpnewbin] = h->GetBinContent(i);
        valerr[tmpnewbin] = h->GetBinError(i);
        cout <<"i = " << i <<  ", tempNewBin = " << tmpnewbin << ", h BinLowEdge = " << h->GetBinLowEdge(i) << ", val = " << val[tmpnewbin]<< endl;
        hnew->SetBinContent(tmpnewbin,val[tmpnewbin]);
        hnew->SetBinError(tmpnewbin,valerr[tmpnewbin]);
        tmpnewbin++;
    }
    return hnew;
}
