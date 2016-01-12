// To get (( Reconstruction Efficiency = recoMatched generated photons / generated isolated photons )) : reco_efficiency2()
// (( Isolation Efficiency = recoMatched isolated generated photons / recoMatched generated photons )) : iso_efficiency()
// Using Event, mc, photon For loops(3 kinds of), Fill the gen photons one by one. 
// Data : 2016 Jan 07
// Creator : Yeonju Go 

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "TCut.h"
#include "TEfficiency.h"
#include "stdio.h"
#include "../yjUtility.h"
#include "/cms/home/ygo/photon_run2/ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"
const float delta2 = 0.15*0.15;
const float delta = 0.15;

const double ptBins[] = {30,40,50,60,70,80,90,100,120,140,180,220};
const int nPtBin = sizeof(ptBins)/sizeof(double) - 1;
const double etaBins[] = {0.0,1.44,2.0,2.5};
//const double etaBins[] = {0.0,0.5,1.0,1.44,2.0,2.5};
const int nEtaBin = sizeof(etaBins)/sizeof(double) - 1;

void reco_efficiency()
{
    //TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleXSize(0.05);
    TFile *f;
    TTree* inTree;
    f = new TFile(Form("/cms/scratch/ygo/photons/Pyquen_AllQCDPhoton30_PhotonFilter20GeV_eta24-HiForest.root"));
    inTree = (TTree*) f->Get("ggHiNtuplizerGED/EventTree");
    ggHiNtuplizer pho;
    pho.setupTreeForReading(inTree);
    TTree *eventTree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    ULong64_t event;
    unsigned run;
    unsigned lumi;
    int hiBin;
    eventTree->SetBranchAddress("evt",&event);
    eventTree->SetBranchAddress("run",&run);
    eventTree->SetBranchAddress("lumi",&lumi);
    eventTree->SetBranchAddress("hiBin",&hiBin);

    //reconstructed photons in generated photons (reco pt & reco eta)
    TH2D* h2D_Num = new TH2D("h2D_Num", ";|#eta|;p_{T} (GeV)", nEtaBin, etaBins, nPtBin, ptBins);
    //generated photons (gen pt & gen eta)
    TH2D* h2D_Den = (TH2D*) h2D_Num->Clone("h2D_Den");
    //reconstructed photons in generated photons / generated photons
    TH2D* h2D_Eff = (TH2D*) h2D_Num->Clone("h2D_Eff");
    //inTree->Draw("phoEt:abs(phoEta)>>h2D_Num","pho_genMatchedIndex!=-1","colz");
    //inTree->Draw("mcEt[pho_genMatchedIndex]:abs(mcEta[pho_genMatchedIndex])>>h2D_Den","pho_genMatchedIndex!=-1","colz");
    //inTree->Draw("mcE:abs(mcEta)>>h2D_Den","(mcPID==22||mcPID==-22)&&(mcStatus==1)","colz");

    Long64_t nentries = inTree->GetEntries();
    for(int ientries = 0; ientries < nentries; ++ientries){
        inTree->GetEntry(ientries);
        eventTree->GetEntry(ientries);

        for (int igen = 0; igen < pho.nMC; ++igen){
            if (pho.mcStatus->at(igen) != 1 || (pho.mcPID->at(igen)) != 22 ) continue;
            if (pho.mcPt->at(igen)<30) continue;
            if (pho.mcCalIsoDR04->at(igen)>5) continue;

            h2D_Den->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));

            Float_t currentDiffPt(1000);
            Float_t currentRecoPt(-1);
            int matchedIndex = -1;
            for(int ipho = 0; ipho < pho.nPho; ipho++){
                float tempDR = getDR(pho.phoEta->at(ipho), pho.phoPhi->at(ipho), pho.mcEta->at(igen), pho.mcPhi->at(igen));
                float tempDiffPt = abs( pho.mcPt->at(igen) - pho.phoEt->at(ipho) );
                if( (tempDR*tempDR < delta2) && (tempDiffPt < currentDiffPt) ) {
                    currentDiffPt = tempDiffPt;
                    currentRecoPt = pho.phoEt->at(ipho);
                    matchedIndex = ipho;
                }
            }
/*
            for(int ipho = 0; ipho < pho.nPho; ipho++){
                if(pho.pho_genMatchedIndex->at(ipho)!=-1) 
                    cout << "genMatchedIndexPt="<<pho.mcPt->at(pho.pho_genMatchedIndex->at(ipho))<<", genPt="<<pho.mcPt->at(igen)<<", recoPt="<<pho.phoEt->at(ipho)<<endl;
            }
            */
            if(matchedIndex!=-1) {
               // cout << "genMatchedIndexPt="<<pho.mcPt->at(pho.pho_genMatchedIndex->at(matchedIndex))<<", genPt="<<pho.mcPt->at(igen)<<", recoPt="<<pho.phoEt->at(matchedIndex)<<endl;
                h2D_Num->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));
            }

        }
    }
    h2D_Eff->Divide(h2D_Num,h2D_Den,1,1,"B");
    
    TH1D* h1D_Num_pt = (TH1D*) h2D_Num->ProjectionY("h1D_Num_pt",0,1); 
    TH1D* h1D_Den_pt = (TH1D*) h2D_Den->ProjectionY("h1D_Den_pt",0,1);
    TH1D* h1D_Eff_pt = (TH1D*) h1D_Den_pt->Clone("h1D_Eff_pt");
    h1D_Eff_pt->Reset();
    h1D_Eff_pt->Divide(h1D_Num_pt,h1D_Den_pt,1,1,"B"); 
    //TEfficiency* hEff = new TEfficiency(h1D_Num_pt,h1D_Den_pt);

    TFile* outf = new TFile("hist_efficiency_reco.root", "RECREATE");
    outf->cd();
    h2D_Num->Write();
    h2D_Den->Write();
    h2D_Eff->Write();
    h1D_Num_pt->Write();
    h1D_Den_pt->Write();
    h1D_Eff_pt->Write();
    outf->Close();
}


void iso_efficiency()
{
    //TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleXSize(0.05);
    TFile *f;
    TTree* inTree;
    f = new TFile(Form("/cms/scratch/ygo/photons/Pyquen_AllQCDPhoton30_PhotonFilter20GeV_eta24-HiForest.root"));
    inTree = (TTree*) f->Get("ggHiNtuplizerGED/EventTree");
    ggHiNtuplizer pho;
    pho.setupTreeForReading(inTree);

    TTree *eventTree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    ULong64_t event;
    unsigned run;
    unsigned lumi;
    int hiBin;
    eventTree->SetBranchAddress("evt",&event);
    eventTree->SetBranchAddress("run",&run);
    eventTree->SetBranchAddress("lumi",&lumi);
    eventTree->SetBranchAddress("hiBin",&hiBin);

    //reconstructed photons in generated photons (reco pt & reco eta)
    TH2D* h2D_Num = new TH2D("h2D_Num", ";|#eta|;p_{T} (GeV)", nEtaBin, etaBins, nPtBin, ptBins);
    //generated photons (gen pt & gen eta)
    TH2D* h2D_Den = (TH2D*) h2D_Num->Clone("h2D_Den");
    //reconstructed photons in generated photons / generated photons
    TH2D* h2D_Eff = (TH2D*) h2D_Num->Clone("h2D_Eff");

    Long64_t nentries = inTree->GetEntries();
    for(int ientries = 0; ientries < nentries; ++ientries){
        inTree->GetEntry(ientries);
        eventTree->GetEntry(ientries);

        for (int igen = 0; igen < pho.nMC; ++igen){
            if (pho.mcStatus->at(igen) != 1 || (pho.mcPID->at(igen)) != 22 ) continue;
            if (pho.mcPt->at(igen)<30) continue;
            if (pho.mcCalIsoDR04->at(igen)>5) continue;

            Float_t currentDiffPt(1000);
            Float_t currentRecoPt(-1);
            int matchedIndex = -1;
            for(int ipho = 0; ipho < pho.nPho; ipho++){
                float tempDR = getDR(pho.phoEta->at(ipho), pho.phoPhi->at(ipho), pho.mcEta->at(igen), pho.mcPhi->at(igen));
                float tempDiffPt = abs( pho.mcPt->at(igen) - pho.phoEt->at(ipho) );
                if( (tempDR*tempDR < delta2) && (tempDiffPt < currentDiffPt) ) {
                    currentDiffPt = tempDiffPt;
                    currentRecoPt = pho.phoEt->at(ipho);
                    matchedIndex = ipho;
                }
            }

            if(matchedIndex!=-1) {
                h2D_Den->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));
                float sumIso = pho.pho_ecalClusterIsoR4->at(matchedIndex) + pho.pho_hcalRechitIsoR4->at(matchedIndex) + pho.pho_trackIsoR4PtCut20->at(matchedIndex);
                if(sumIso<1 && pho.phoHoverE->at(matchedIndex)<0.1) h2D_Num->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));
            }
        }
    }
    h2D_Eff->Divide(h2D_Num,h2D_Den,1,1,"B");

    TH1D* h1D_Num_pt = (TH1D*) h2D_Num->ProjectionY("h1D_Num_pt",0,1); 
    TH1D* h1D_Den_pt = (TH1D*) h2D_Den->ProjectionY("h1D_Den_pt",0,1);
    TH1D* h1D_Eff_pt = (TH1D*) h1D_Den_pt->Clone("h1D_Eff_pt");
    h1D_Eff_pt->Reset();
    h1D_Eff_pt->Divide(h1D_Num_pt,h1D_Den_pt,1,1,"B"); 
//    TEfficiency* hEff = new TEfficiency(h1D_Num_pt,h1D_Den_pt);

    TFile* outf = new TFile("hist_efficiency_iso.root", "RECREATE");
    outf->cd();
    h2D_Num->Write();
    h2D_Den->Write();
    h2D_Eff->Write();
    h1D_Num_pt->Write();
    h1D_Den_pt->Write();
    h1D_Eff_pt->Write();
    outf->Close();
}

int main(){
    reco_efficiency();
    iso_efficiency();
    return 0;
}





