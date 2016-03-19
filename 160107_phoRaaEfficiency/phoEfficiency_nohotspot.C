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
#include <vector>
#include "../yjUtility.h"
#include "../phoRaaCuts_v1.h"
#include "../ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"
//using namespace std;
const float delta2 = 0.15*0.15;
const float delta = 0.15;
/*
const double ptBins[] = {20,30,40,50,60,70,80,90,100,120,140,180,220};
//const double ptBins[] = {30,40,50,60,70,80,90,100,120,140,180,220};
const int nPtBin = sizeof(ptBins)/sizeof(double) - 1;
const double etaBins[] = {0.0,1.44,2.0,2.5};
const int nEtaBin = sizeof(etaBins)/sizeof(double) - 1;
const int centBins[] = {0,60,200};
const int nCentBin = sizeof(centBins)/sizeof(int) -1;
///// KNU
const TString pbpbfname = "/d3/scratch/goyeonju/files/photons2016/mergedpthat_2015-PbPb-MC-AllQCDPhoton_v1.root";
const TString ppfname = "/d3/scratch/goyeonju/files/photons2016/2015-PP-MC_Pythia8_Photon30_pp502_TuneCUETP8M1.root";
///// MIT
//const TString pbpbfname = "/cms/scratch/ygo/photons/Pyquen_AllQCDPhoton30_PhotonFilter20GeV_eta24-HiForest.root";
//const TString ppfname = "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/partial_PbPb_gammaJet_MC/HiForest_QCDPhoton30.root";
*/
void reco_efficiency(const TString coll="pbpb", int cent_i=0, int cent_f=200, const char* ver="v0")
{
    TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleXSize(0.05);
    TFile *f;
    TTree* inTree, *eventTree;
    TString fname;
    TString outfname; 
    const char* cap="";
    if(coll=="pbpb"){ 
        fname = pbpbMCfname;
        cap = Form("pbpb_cent%d-%d_%s",cent_i,cent_f,ver);
        f = TFile::Open(fname);
        inTree = (TTree*) f->Get("EventTree");
        eventTree = (TTree*)f->Get("HiTree");
    }
    else if(coll=="pp"){
        fname = ppMCfname;
        cap=Form("pp_%s",ver);
        f = TFile::Open(fname);
        inTree = (TTree*) f->Get("ggHiNtuplizerGED/EventTree");
        eventTree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    }
    outfname = Form("hist_efficiency_reco_%s_nohotspot.root",cap);
    cout << "fname= " << fname << endl;
    cout << "outfname= " << outfname << endl;

    ggHiNtuplizer pho;
    pho.setupTreeForReading(inTree);


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
    int nhotspot = 0;
    Long64_t nentries = inTree->GetEntries();
    for(int ientries = 0; ientries < nentries; ++ientries){
        inTree->GetEntry(ientries);
        eventTree->GetEntry(ientries);
        if(coll=="pbpb" && !(hiBin>=cent_i && hiBin<cent_f)) continue;
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
                if( (tempDR < delta) && (tempDiffPt < currentDiffPt) ) {
                    currentDiffPt = tempDiffPt;
                    currentRecoPt = pho.phoEt->at(ipho);
                    matchedIndex = ipho;
                }
            }
            // removing hotspot
            if((matchedIndex!=-1) && (((pho.phoSCEta->at(matchedIndex) > -1.194 && pho.phoSCEta->at(matchedIndex) < -1.169) && (pho.phoSCPhi->at(matchedIndex) > -1.211 && pho.phoSCPhi->at(matchedIndex) < -1.162)) || ((pho.phoSCEta->at(matchedIndex) > -0.935 && pho.phoSCEta->at(matchedIndex) < -0.910) && (pho.phoSCPhi->at(matchedIndex) > 2.036 && pho.phoSCPhi->at(matchedIndex) < 2.065)) || ((pho.phoSCEta->at(matchedIndex) > -0.850 && pho.phoSCEta->at(matchedIndex) < -0.819) && (pho.phoSCPhi->at(matchedIndex) > -1.723 && pho.phoSCPhi->at(matchedIndex) < -1.676)) || ((pho.phoSCEta->at(matchedIndex) > -0.764 && pho.phoSCEta->at(matchedIndex) < -0.723) && (pho.phoSCPhi->at(matchedIndex) > -0.786 && pho.phoSCPhi->at(matchedIndex) < -0.720)) || ((pho.phoSCEta->at(matchedIndex) > -0.237 && pho.phoSCEta->at(matchedIndex) < -0.210) && (pho.phoSCPhi->at(matchedIndex) > 0.106 && pho.phoSCPhi->at(matchedIndex) < 0.147)) ||((pho.phoSCEta->at(matchedIndex) > -0.062 && pho.phoSCEta->at(matchedIndex) < -0.036) && (pho.phoSCPhi->at(matchedIndex) > 0.860 && pho.phoSCPhi->at(matchedIndex)< 0.931)) || ((pho.phoSCEta->at(matchedIndex) > 0.973 && pho.phoSCEta->at(matchedIndex) < 1.035) && (pho.phoSCPhi->at(matchedIndex) > 0.623 && pho.phoSCPhi->at(matchedIndex) < 0.697)) || ((pho.phoSCEta->at(matchedIndex) > 1.075 && pho.phoSCEta->at(matchedIndex) < 1.108) && (pho.phoSCPhi->at(matchedIndex) > -3.114 && pho.phoSCPhi->at(matchedIndex) < -3.054)) || ((pho.phoSCEta->at(matchedIndex) > 1.426 && pho.phoSCEta->at(matchedIndex) < 1.447) && (pho.phoSCPhi->at(matchedIndex) > 2.744 && pho.phoSCPhi->at(matchedIndex) < 2.774)) || ((pho.phoSCEta->at(matchedIndex) > 0.915 && pho.phoSCEta->at(matchedIndex) < 0.930) && (pho.phoSCPhi->at(matchedIndex) > 1.69 && pho.phoSCPhi->at(matchedIndex) < 1.72)))) { 
                matchedIndex=-1; 
                nhotspot++;
                cout << "this is " << nhotspot << "th hotspot" << endl;
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

    saveHistogramsToPicture(h1D_Eff_pt, "pdf",Form("reco_%s_nohotspot",cap));
    saveHistogramsToPicture(h2D_Eff, "pdf",Form("reco_%s_nohotspot",cap));
    TFile* outf = new TFile(Form("output/%s",outfname.Data()), "RECREATE");
    outf->cd();
    h2D_Num->Write();
    h2D_Den->Write();
    h2D_Eff->Write();
    h1D_Num_pt->Write();
    h1D_Den_pt->Write();
    h1D_Eff_pt->Write();
    outf->Close();
}


void iso_efficiency(const TString coll="pbpb",
        int cent_i=0, int cent_f=200, const char* ver="v0")
{
    //TH1::SetDefaultSumw2();
    gStyle->SetLabelSize(0.03,"Y");
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleXSize(0.05);

    TFile *f;
    TTree* inTree, *eventTree;
    TString fname;
    TString outfname; 
    const char* cap="";
    if(coll=="pbpb"){ 
        fname = pbpbMCfname;
        cap = Form("pbpb_cent%d-%d_%s",cent_i,cent_f,ver);
        f = TFile::Open(fname);
        inTree = (TTree*) f->Get("EventTree");
        eventTree = (TTree*)f->Get("HiTree");
    }
    else if(coll=="pp"){
        fname = ppMCfname;
        cap=Form("pp_%s",ver);
        f = TFile::Open(fname);
        inTree = (TTree*) f->Get("ggHiNtuplizerGED/EventTree");
        eventTree = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
    }
    outfname = Form("hist_efficiency_iso_%s_nohotspot.root",cap);
    cout << "fname= " << fname << endl;
    cout << "outfname= " << outfname << endl;

    ggHiNtuplizer pho;
    pho.setupTreeForReading(inTree);

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

        if(coll=="pbpb" && !(hiBin>=cent_i && hiBin<cent_f)) continue;
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
                if( (tempDR < delta) && (tempDiffPt < currentDiffPt) ) {
                    currentDiffPt = tempDiffPt;
                    currentRecoPt = pho.phoEt->at(ipho);
                    matchedIndex = ipho;
                }
            }
            // removing hotspot
            if((matchedIndex!=-1) && (((pho.phoSCEta->at(matchedIndex) > -1.194 && pho.phoSCEta->at(matchedIndex) < -1.169) && (pho.phoSCPhi->at(matchedIndex) > -1.211 && pho.phoSCPhi->at(matchedIndex) < -1.162)) || ((pho.phoSCEta->at(matchedIndex) > -0.935 && pho.phoSCEta->at(matchedIndex) < -0.910) && (pho.phoSCPhi->at(matchedIndex) > 2.036 && pho.phoSCPhi->at(matchedIndex) < 2.065)) || ((pho.phoSCEta->at(matchedIndex) > -0.850 && pho.phoSCEta->at(matchedIndex) < -0.819) && (pho.phoSCPhi->at(matchedIndex) > -1.723 && pho.phoSCPhi->at(matchedIndex) < -1.676)) || ((pho.phoSCEta->at(matchedIndex) > -0.764 && pho.phoSCEta->at(matchedIndex) < -0.723) && (pho.phoSCPhi->at(matchedIndex) > -0.786 && pho.phoSCPhi->at(matchedIndex) < -0.720)) || ((pho.phoSCEta->at(matchedIndex) > -0.237 && pho.phoSCEta->at(matchedIndex) < -0.210) && (pho.phoSCPhi->at(matchedIndex) > 0.106 && pho.phoSCPhi->at(matchedIndex) < 0.147)) ||((pho.phoSCEta->at(matchedIndex) > -0.062 && pho.phoSCEta->at(matchedIndex) < -0.036) && (pho.phoSCPhi->at(matchedIndex) > 0.860 && pho.phoSCPhi->at(matchedIndex)< 0.931)) || ((pho.phoSCEta->at(matchedIndex) > 0.973 && pho.phoSCEta->at(matchedIndex) < 1.035) && (pho.phoSCPhi->at(matchedIndex) > 0.623 && pho.phoSCPhi->at(matchedIndex) < 0.697)) || ((pho.phoSCEta->at(matchedIndex) > 1.075 && pho.phoSCEta->at(matchedIndex) < 1.108) && (pho.phoSCPhi->at(matchedIndex) > -3.114 && pho.phoSCPhi->at(matchedIndex) < -3.054)) || ((pho.phoSCEta->at(matchedIndex) > 1.426 && pho.phoSCEta->at(matchedIndex) < 1.447) && (pho.phoSCPhi->at(matchedIndex) > 2.744 && pho.phoSCPhi->at(matchedIndex) < 2.774)) || ((pho.phoSCEta->at(matchedIndex) > 0.915 && pho.phoSCEta->at(matchedIndex) < 0.930) && (pho.phoSCPhi->at(matchedIndex) > 1.69 && pho.phoSCPhi->at(matchedIndex) < 1.72)))) { 
                matchedIndex=-1; 
            }

            if(matchedIndex!=-1) {
                h2D_Den->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));
                float sumIso = pho.pho_ecalClusterIsoR4->at(matchedIndex) + pho.pho_hcalRechitIsoR4->at(matchedIndex) + pho.pho_trackIsoR4PtCut20->at(matchedIndex);
                if(sumIso<1 && pho.phoHoverE->at(matchedIndex)<0.1 && pho.phoSigmaIEtaIEta->at(matchedIndex)<0.0102) h2D_Num->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));
                //if(sumIso<1 && pho.phoHoverE->at(matchedIndex)<0.1) h2D_Num->Fill(abs(pho.mcEta->at(igen)), pho.mcEt->at(igen));
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

    saveHistogramsToPicture(h1D_Eff_pt, "pdf",Form("iso_%s_nohotspot",cap));
    saveHistogramsToPicture(h2D_Eff, "pdf",Form("iso_%s_nohotspot",cap));

    TFile* outf = new TFile(Form("output/%s",outfname.Data()), "RECREATE");
    outf->cd();
    h2D_Num->Write();
    h2D_Den->Write();
    h2D_Eff->Write();
    h1D_Num_pt->Write();
    h1D_Den_pt->Write();
    h1D_Eff_pt->Write();
    outf->Close();
}

int main(int argc, char **argv){
    gROOT->ProcessLine("#include <vector>");
    //    gROOT->ProcessLine("#include <vector>");
    gStyle -> SetOptStat(0);
    if(argc==3){
        for(int jj=0;jj<nCentBin;++jj){
            TString collst; 
            if(std::string(argv[1])=="pbpb"){ 
                collst = Form("pbpb_cent%d-%d_%s",centBins[jj],centBins[jj+1],argv[2]);
            }
            else if(std::string(argv[1])=="pp"){
                collst = Form("pp_%s",argv[2]);
                if(jj<nCentBin-1) continue;
            }
            cout << "collst = " << collst << endl;
            reco_efficiency(argv[1],centBins[jj],centBins[jj+1], argv[2]);
            iso_efficiency(argv[1],centBins[jj],centBins[jj+1], argv[2]);

            TH1::SetDefaultSumw2();
            TFile* freco = new TFile(Form("output/hist_efficiency_reco_%s_nohotspot.root",collst.Data()));
            TH1D* hreco = (TH1D*) freco->Get("h1D_Eff_pt");
            hreco->SetName("h1D_Eff_pt_reco");
            TFile* fiso = new TFile(Form("output/hist_efficiency_iso_%s_nohotspot.root",collst.Data()));
            TH1D* hiso = (TH1D*) fiso->Get("h1D_Eff_pt");
            hiso->SetName("h1D_Eff_pt_iso");
            TH1D htotal = (*hreco)*(*hiso); 
            htotal.SetName("h1D_Eff_pt_total");

            //saveHistogramsToPicture(htotal, "pdf");
            TFile* fout = new TFile(Form("output/hist_efficiency_total_%s_nohotspot.root",collst.Data()),"RECREATE");
            fout->cd();
            htotal.Write();
            fout->Close();
        }
        return 0;
    } else{
        return 1;
    }
}





