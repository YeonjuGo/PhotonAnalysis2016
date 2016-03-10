// 2016 Mar 3. comparing two datasets (ex. allqcdphoton and emenricheddijet) 
// Author : Yeonju Go

#include "../yjUtility.h"


void getRatioPlot(const char* fname_1="", const char* fname_2="", const char* cap="");

void compare_ratio(){
   getRatioPlot("output/hist_efficiency_total_pbpb_cent0-60.root","output/hist_efficiency_total_pbpb_cent0-60_nohotspot.root","pbpb_0-60");
   getRatioPlot("output/hist_efficiency_total_pbpb_cent60-200.root","output/hist_efficiency_total_pbpb_cent60-200_nohotspot.root","pbpb_60-200");
   
}

void getRatioPlot(const char* fname_1, const char* fname_2, const char* cap){
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    static int i = 1;
    TFile* f1 = new TFile(fname_1);
    TH1D* h1 = (TH1D*) f1 -> Get("h1D_Eff_pt_total");
    h1->SetName(Form("h1D_Eff_pt_total_withHotSpot_%d",i));
    TFile* f2 = new TFile(fname_2);
    TH1D* h2 = (TH1D*) f2 -> Get("h1D_Eff_pt_total");
    h2->SetName(Form("h1D_Eff_pt_total_noHotSpot_%d",i));
 
    TCanvas* c =new TCanvas(Form("c_%d",i),"", 400,800);
    c->Divide(1,2);
    c->cd(1);
    //gPad->SetLogy();

    h1->SetAxisRange(0.4,1.0,"Y");
    //h1->Sumw2();
    //h2->Sumw2();

    h1->SetMarkerStyle(20);
    h1->SetMarkerSize(0.8);
    h1->SetMarkerColor(1);
 
    h2->SetMarkerStyle(20);
    h2->SetMarkerSize(0.8);
    h2->SetMarkerColor(2);
 
    h1->DrawCopy("pe");
    h2->DrawCopy("pe same");
    TLegend* l1 = new TLegend(0.7,0.75,0.95,0.9);
    legStyle(l1);
    l1->AddEntry(h1,"with hotspot","pl");
    l1->AddEntry(h2,"without hotspot","pl");
    l1->Draw();
    double xMin=h1->GetBinLowEdge(1);
    double xMax=h1->GetBinLowEdge(h1->GetNbinsX())+h1->GetBinWidth(h1->GetNbinsX());
    jumSun(xMin,1,xMax,1);
    drawText(cap,0.7,0.2);

    c->cd(2);
    h2->Divide(h1);
    h2->SetYTitle("without hotspot / with hotspot Ratio");
    h2->SetAxisRange(0.9,1.1,"Y");
    h2->DrawCopy("le1");
    jumSun(xMin,1,xMax,1);
    c-> SaveAs(Form("figures/compare_%s.pdf",cap));
    i++;
}
