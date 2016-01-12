// To get (( Isolation Efficiency = reconstructed photons passing through isolation cuts / reconstructed photons ))
// Using Draw() function. 
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

    double ptBins[] = {30,40,50,60,70,80,90,100,120,140,180,220};
    int nPtBin = sizeof(ptBins)/sizeof(double) - 1;
    double etaBins[] = {0.0,1.44,2.0,2.5};
    //double etaBins[] = {0.0,0.5,1.0,1.44,2.0,2.5};
    int nEtaBin = sizeof(etaBins)/sizeof(double) - 1;

    //pass the isolation condition
    TH2D* h2D_Num = new TH2D("h2D_Num", ";|#eta|;p_{T} (GeV)", nEtaBin, etaBins, nPtBin, ptBins);
    //all reconstructed photons 
    TH2D* h2D_Den = (TH2D*) h2D_Num->Clone("h2D_Den");
    //reconstructed photons passing through isolation condition / all reconstructed photons
    TH2D* h2D_Eff = (TH2D*) h2D_Num->Clone("h2D_Eff");


    inTree->Draw("phoEt:abs(phoEta)>>h2D_Num","","colz");
    inTree->Draw("phoEt:abs(phoEta)>>h2D_Den","pho_genMatchedIndex!=-1","colz");

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
/*    
    TH1D* h1D_Num_pt = (TH1D*) h2D_Num->ProjectionY("h1D_Num_pt",0,2), ";;p_{T} (GeV)", nPtBin, ptBins);
    TH1D* h1D_Den_pt = (TH1D*) h1D_Num_pt->Clone("h1D_Den_pt");
    
    TEfficiency* hEff = new TEfficiency(h1D_Num,h1D_Den);
*/
    TFile* outf = new TFile("hist_efficiency2.root", "RECREATE");
    outf->cd();
    h2D_Num->Write();
    h2D_Den->Write();
    h2D_Eff->Write();
    hEff->Write();
    outf->Close();
}

