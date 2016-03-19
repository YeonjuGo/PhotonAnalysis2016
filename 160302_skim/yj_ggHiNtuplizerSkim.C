#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <iostream>

#include "../../CorrelationTuple/EventMatcher.h"
#include "../../TreeHeaders/ggHiNtuplizerTree.h"
#include "../../TreeHeaders/CutConfigurationTree.h"
#include "../../TreeHeaders/ggHiNtuplizerPhotonSkim.h"
#include "../../Utilities/interface/CutConfigurationParser.h"
#include "../../Utilities/interface/InputConfigurationParser.h"

const long MAXTREESIZE = 500000000000; // set maximum tree size from 10 GB to 100 GB, so that the code does not switch to a new file after 10 GB7

void ggHiNtuplizerSkim(const TString configFile, const TString inputFile, const TString outputSkimFile="photonSkim.root")
{
    bool isMC = true;
    TStopwatch timer;
    timer.Start();

    std::cout<<"running ggHiNtuplizerSkim()"   <<std::endl;
    std::cout<<"configFile  = "<< configFile.Data() <<std::endl;
    std::cout<<"inputFile   = "<< inputFile.Data() <<std::endl;
    std::cout<<"outputFile  = "<< outputSkimFile.Data() <<std::endl;

    // input configuration
    CutConfiguration configCuts = CutConfigurationParser::Parse(configFile.Data());
    TTree *configTree = setupConfigurationTreeForWriting(configCuts);

    int montecarlo = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].i[CUTS::PHO::k_MonteCarlo];

    std::string treePath = "ggHiNtuplizer/EventTree";

    // cut configuration
    float cut_phoSigmaIEtaIEta_spike;
    float cut_pho_swissCrx;
    float cut_pho_seedTime;
    float cut_pho_et; 
    if (configCuts.isValid) {
        cut_phoSigmaIEtaIEta_spike = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].f[CUTS::PHO::k_phoSigmaIEtaIEta_spike];
        cut_pho_swissCrx = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].f[CUTS::PHO::k_pho_swissCrx];
        cut_pho_seedTime = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].f[CUTS::PHO::k_pho_seedTime];
        cut_pho_et = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].f[CUTS::PHO::k_et];
    }
    else {
        cut_phoSigmaIEtaIEta_spike = 0.002;
        cut_pho_swissCrx = 0.9;
        cut_pho_seedTime = 3;
        cut_pho_et = 20.0; 
    }

    // verbose about cut configuration
    std::cout<<"Cut Configuration :"<<std::endl;
    std::cout<<"cut_phoSigmaIEtaIEta_spike = "<<cut_phoSigmaIEtaIEta_spike<<std::endl;   // not used
    std::cout<<"cut_pho_swissCrx           = "<<cut_pho_swissCrx<<std::endl;             // not used
    std::cout<<"cut_pho_seedTime           = "<<cut_pho_seedTime<<std::endl;             // not used


    std::vector<std::string> inputFiles = InputConfigurationParser::ParseFiles(inputFile.Data());

    std::cout<<"input ROOT files : num = "<<inputFiles.size()<< std::endl;
    std::cout<<"#####"<< std::endl;
    for (std::vector<std::string>::iterator it = inputFiles.begin() ; it != inputFiles.end(); ++it) {
        std::cout<<(*it).c_str()<< std::endl;
    }
    std::cout<<"##### END #####"<< std::endl;

    TChain* treeHLT   = new TChain("hltanalysis/HltTree");
    TChain* treeggHiNtuplizer  = new TChain("ggHiNtuplizer/EventTree");
    TChain* treeHiEvt = new TChain("hiEvtAnalyzer/HiTree");
    TChain* skimTree = new TChain("skimanalysis/HltTree");
    TChain* genTree = new TChain("HiGenParticleAna/hi");
    int fnum = 0;
    for (std::vector<std::string>::iterator it = inputFiles.begin() ; it != inputFiles.end(); ++it) {
        fnum++;
        if(fnum%200==0) 
            std::cout <<fnum<<" / "<<inputFiles.size()<<"       "<<(double)fnum/(inputFiles.size())*100<<std::endl;
            //std::cout <<fnum<<" / "<<inputFiles.size()<<" "<<std::setprecision(2)<<(double)fnum/(inputFiles.size())*100<<std::endl;
        treeHLT->Add((*it).c_str());
        treeggHiNtuplizer->Add((*it).c_str());
        treeHiEvt->Add((*it).c_str());
        skimTree->Add((*it).c_str());
        if(isMC) genTree->Add((*it).c_str()); 
    }
/*
    treeggHiNtuplizer->SetBranchStatus("*",0);     // disable all branches
    treeggHiNtuplizer->SetBranchStatus("pho*",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("tower*",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("pf*",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("nPho",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("run",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("event",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("lumis",1);     // enable photon branches
    treeggHiNtuplizer->SetBranchStatus("isData",1);     // enable photon branches
*/
    treeHLT->SetBranchStatus("*",0);     // disable all branches
    treeHLT->SetBranchStatus("HLT_HI*SinglePhoton*Eta*v1*",1);     // enable photon branches
    treeHLT->SetBranchStatus("HLT_HI*DoublePhoton*Eta*v1*",1);     // enable photon branches
    if(isMC) 
    ggHiNtuplizer pho;
    pho.setupTreeForReading(treeggHiNtuplizer);

    Int_t HBHENoiseFilterResult;
    Int_t pcollisionEventSelection;
    skimTree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
    skimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);

    ULong64_t event;
    unsigned run;
    unsigned lumi;
    float pthat;
    int hiBin;
    treeHiEvt->SetBranchAddress("evt",&event);
    treeHiEvt->SetBranchAddress("run",&run);
    treeHiEvt->SetBranchAddress("lumi",&lumi);
    treeHiEvt->SetBranchAddress("hiBin",&hiBin);
    treeHiEvt->SetBranchAddress("pthat",&pthat);

    TFile *outFile = TFile::Open(outputSkimFile,"RECREATE");
    outFile->cd();
    PhotonSkim phoSkim(montecarlo);
    TTree *outTree = new TTree("photonSkimTree","photonSkimTree");
    phoSkim.setupTreeForWriting(outTree);

    // output tree variables
    TTree *outputTreeHLT           = treeHLT->CloneTree(0);
    outputTreeHLT->SetName("hltTree");
    outputTreeHLT->SetTitle("subbranches of hltanalysis/HltTree");
    TTree *outputTreeggHiNtuplizer = treeggHiNtuplizer->CloneTree(0);
    TTree *outputTreeHiEvt = treeHiEvt->CloneTree(0);
    TTree *outputTreeSkim = skimTree->CloneTree(0);

    outputTreeHLT->SetMaxTreeSize(MAXTREESIZE);
    outputTreeggHiNtuplizer->SetMaxTreeSize(MAXTREESIZE);
    outputTreeHiEvt->SetMaxTreeSize(MAXTREESIZE);
    outputTreeSkim->SetMaxTreeSize(MAXTREESIZE);


    EventMatcher* em = new EventMatcher();
    Long64_t duplicateEntries = 0;

    Long64_t nentries = treeggHiNtuplizer->GetEntries();
    Long64_t entriesAnalyzed = 0;
    std::cout << "entries = " << nentries << std::endl;
    std::cout<< "Loop : " << treePath.c_str() <<std::endl;


    for(int ientries = 0; ientries < nentries; ++ientries)
    {
        if (ientries% 5000 == 0)  {
            std::cout <<ientries<<" / "<<nentries<<"       "<<(double)ientries/nentries*100<<std::endl;
        }
        treeggHiNtuplizer->GetEntry(ientries);
        skimTree->GetEntry(ientries);

        bool eventAdded = em->addEvent(pho.run,pho.lumis,pho.event,ientries);
        if(!eventAdded) // this event is duplicate, skip this one.
        {
            duplicateEntries++;
            continue;
        }


        if(pho.nPho < 1) continue;
        Float_t maxPho =-1;
        for(int i = 0; i < pho.nPho; ++i)
        {
            //bool failedEtaCut = (TMath::Abs(pho.phoEta->at(i)) > cutPhoEta) ;
            bool failedSpikeRejection;
           // if (isHI) {
                failedSpikeRejection = (pho.phoSigmaIEtaIEta->at(i) < 0.002 ||
                        pho.pho_swissCrx->at(i)     > 0.9   ||
                        TMath::Abs(pho.pho_seedTime->at(i)) > 3);
           // }
           // else {
           //     failedSpikeRejection = (pho.phoSigmaIEtaIEta->at(i) < 0.002);
           // }

            bool failedHoverE = (pho.phoHoverE->at(i) > 0.2);      // <0.1 cut is applied after corrections
            //               bool failedEnergyRatio = ((float)pho.phoSCRawE->at(i)/pho.phoE->at(i) < 0.5);

            //if (failedEtaCut)         continue;
            if (failedSpikeRejection) continue;
            if (failedHoverE)         continue; 

            if(maxPho < pho.phoEt->at(i)) maxPho = pho.phoEt->at(i);
        }
        if(maxPho < cut_pho_et) continue;

        //std::cout << "pass" << std::endl;
        phoSkim.reset();

        treeHLT->GetEntry(ientries);
        treeHiEvt->GetEntry(ientries);
        phoSkim.run = run;
        phoSkim.lumis = lumi;
        phoSkim.event = event;
        phoSkim.hiBin = hiBin;

        phoSkim.HBHENoiseFilterResult = HBHENoiseFilterResult;
        phoSkim.pcollisionEventSelection = pcollisionEventSelection;

        phoSkim.nPho = pho.nPho;

        phoSkim.phoE = * pho.phoE;
        phoSkim.phoEt = * pho.phoEt;
        phoSkim.phoEta = * pho.phoEta;
        phoSkim.phoPhi = * pho.phoPhi;
        phoSkim.phoSCE = * pho.phoSCE;
        phoSkim.phoSCRawE = * pho.phoSCRawE;
        phoSkim.phoESEn = * pho.phoESEn;
        phoSkim.phoSCEta = * pho.phoSCEta;
        phoSkim.phoSCPhi = * pho.phoSCPhi;
        phoSkim.phoSCEtaWidth = * pho.phoSCEtaWidth;
        phoSkim.phoSCPhiWidth = * pho.phoSCPhiWidth;
        phoSkim.phoSCBrem = * pho.phoSCBrem;
        phoSkim.phohasPixelSeed = * pho.phohasPixelSeed;
        phoSkim.phoR9 = * pho.phoR9;
        phoSkim.phoHoverE = * pho.phoHoverE;
        phoSkim.phoSigmaIEtaIEta = * pho.phoSigmaIEtaIEta;
        phoSkim.phoE1x3 = * pho.phoE1x3;
        phoSkim.phoE2x2 = * pho.phoE2x2;
        //phoSkim.phoE3x3 = * pho.phoE3x3;
        phoSkim.phoE2x5Max = * pho.phoE2x5Max;
        //phoSkim.phoE1x5 = * pho.phoE1x5;
        //phoSkim.phoE2x5 = * pho.phoE2x5;
        phoSkim.phoE5x5 = * pho.phoE5x5;
        //phoSkim.phoMaxEnergyXtal = * pho.phoMaxEnergyXtal;
        //phoSkim.phoSigmaEtaEta = * pho.phoSigmaEtaEta;
        //phoSkim.phoR1x5 = * pho.phoR1x5;
        //phoSkim.phoR2x5 = * pho.phoR2x5;
        phoSkim.phoESEffSigmaRR = * pho.phoESEffSigmaRR;
        phoSkim.phoSigmaIEtaIEta_2012 = * pho.phoSigmaIEtaIEta_2012;
        phoSkim.phoSigmaIEtaIPhi_2012 = * pho.phoSigmaIEtaIPhi_2012;
        phoSkim.phoSigmaIPhiIPhi_2012 = * pho.phoSigmaIPhiIPhi_2012;
        phoSkim.phoE1x3_2012 = * pho.phoE1x3_2012;
        phoSkim.phoE2x2_2012 = * pho.phoE2x2_2012;
        //phoSkim.phoE3x3_2012 = * pho.phoE3x3_2012;
        phoSkim.phoE2x5Max_2012 = * pho.phoE2x5Max_2012;
        phoSkim.phoE5x5_2012 = * pho.phoE5x5_2012;
        phoSkim.phoBC1E = * pho.phoBC1E;
        phoSkim.phoBC1Eta = * pho.phoBC1Eta;
        phoSkim.phoBC2E = * pho.phoBC2E;
        phoSkim.phoBC2Eta = * pho.phoBC2Eta;
        phoSkim.pho_ecalClusterIsoR2 = * pho.pho_ecalClusterIsoR2;
        phoSkim.pho_ecalClusterIsoR3 = * pho.pho_ecalClusterIsoR3;
        phoSkim.pho_ecalClusterIsoR4 = * pho.pho_ecalClusterIsoR4;
        phoSkim.pho_ecalClusterIsoR5 = * pho.pho_ecalClusterIsoR5;
        phoSkim.pho_hcalRechitIsoR1 = * pho.pho_hcalRechitIsoR1;
        phoSkim.pho_hcalRechitIsoR2 = * pho.pho_hcalRechitIsoR2;
        phoSkim.pho_hcalRechitIsoR3 = * pho.pho_hcalRechitIsoR3;
        phoSkim.pho_hcalRechitIsoR4 = * pho.pho_hcalRechitIsoR4;
        phoSkim.pho_hcalRechitIsoR5 = * pho.pho_hcalRechitIsoR5;
        phoSkim.pho_trackIsoR1PtCut20 = * pho.pho_trackIsoR1PtCut20;
        phoSkim.pho_trackIsoR2PtCut20 = * pho.pho_trackIsoR2PtCut20;
        phoSkim.pho_trackIsoR3PtCut20 = * pho.pho_trackIsoR3PtCut20;
        phoSkim.pho_trackIsoR4PtCut20 = * pho.pho_trackIsoR4PtCut20;
        phoSkim.pho_trackIsoR5PtCut20 = * pho.pho_trackIsoR5PtCut20;
        phoSkim.pho_swissCrx = * pho.pho_swissCrx;
        phoSkim.pho_seedTime = * pho.pho_seedTime;

        phoSkim.pfcIso1 = * pho.pfcIso1;
        phoSkim.pfcIso2 = * pho.pfcIso2;
        phoSkim.pfcIso3 = * pho.pfcIso3;
        phoSkim.pfcIso4 = * pho.pfcIso4;
        phoSkim.pfcIso5 = * pho.pfcIso5;

        phoSkim.pfpIso1 = * pho.pfpIso1;
        phoSkim.pfpIso2 = * pho.pfpIso2;
        phoSkim.pfpIso3 = * pho.pfpIso3;
        phoSkim.pfpIso4 = * pho.pfpIso4;
        phoSkim.pfpIso5 = * pho.pfpIso5;

        phoSkim.pfnIso1 = * pho.pfnIso1;
        phoSkim.pfnIso2 = * pho.pfnIso2;
        phoSkim.pfnIso3 = * pho.pfnIso3;
        phoSkim.pfnIso4 = * pho.pfnIso4;
        phoSkim.pfnIso5 = * pho.pfnIso5;

        phoSkim.pfcVsIso1 = * pho.pfcVsIso1;
        phoSkim.pfcVsIso2 = * pho.pfcVsIso2;
        phoSkim.pfcVsIso3 = * pho.pfcVsIso3;
        phoSkim.pfcVsIso4 = * pho.pfcVsIso4;
        phoSkim.pfcVsIso5 = * pho.pfcVsIso5;
        phoSkim.pfcVsIso1th1 = * pho.pfcVsIso1th1;
        phoSkim.pfcVsIso2th1 = * pho.pfcVsIso2th1;
        phoSkim.pfcVsIso3th1 = * pho.pfcVsIso3th1;
        phoSkim.pfcVsIso4th1 = * pho.pfcVsIso4th1;
        phoSkim.pfcVsIso5th1 = * pho.pfcVsIso5th1;
        phoSkim.pfcVsIso1th2 = * pho.pfcVsIso1th2;
        phoSkim.pfcVsIso2th2 = * pho.pfcVsIso2th2;
        phoSkim.pfcVsIso3th2 = * pho.pfcVsIso3th2;
        phoSkim.pfcVsIso4th2 = * pho.pfcVsIso4th2;
        phoSkim.pfcVsIso5th2 = * pho.pfcVsIso5th2;
        phoSkim.pfnVsIso1 = * pho.pfnVsIso1;
        phoSkim.pfnVsIso2 = * pho.pfnVsIso2;
        phoSkim.pfnVsIso3 = * pho.pfnVsIso3;
        phoSkim.pfnVsIso4 = * pho.pfnVsIso4;
        phoSkim.pfnVsIso5 = * pho.pfnVsIso5;
        phoSkim.pfnVsIso1th1 = * pho.pfnVsIso1th1;
        phoSkim.pfnVsIso2th1 = * pho.pfnVsIso2th1;
        phoSkim.pfnVsIso3th1 = * pho.pfnVsIso3th1;
        phoSkim.pfnVsIso4th1 = * pho.pfnVsIso4th1;
        phoSkim.pfnVsIso5th1 = * pho.pfnVsIso5th1;
        phoSkim.pfnVsIso1th2 = * pho.pfnVsIso1th2;
        phoSkim.pfnVsIso2th2 = * pho.pfnVsIso2th2;
        phoSkim.pfnVsIso3th2 = * pho.pfnVsIso3th2;
        phoSkim.pfnVsIso4th2 = * pho.pfnVsIso4th2;
        phoSkim.pfnVsIso5th2 = * pho.pfnVsIso5th2;
        phoSkim.pfpVsIso1 = * pho.pfpVsIso1;
        phoSkim.pfpVsIso2 = * pho.pfpVsIso2;
        phoSkim.pfpVsIso3 = * pho.pfpVsIso3;
        phoSkim.pfpVsIso4 = * pho.pfpVsIso4;
        phoSkim.pfpVsIso5 = * pho.pfpVsIso5;
        phoSkim.pfpVsIso1th1 = * pho.pfpVsIso1th1;
        phoSkim.pfpVsIso2th1 = * pho.pfpVsIso2th1;
        phoSkim.pfpVsIso3th1 = * pho.pfpVsIso3th1;
        phoSkim.pfpVsIso4th1 = * pho.pfpVsIso4th1;
        phoSkim.pfpVsIso5th1 = * pho.pfpVsIso5th1;
        phoSkim.pfpVsIso1th2 = * pho.pfpVsIso1th2;
        phoSkim.pfpVsIso2th2 = * pho.pfpVsIso2th2;
        phoSkim.pfpVsIso3th2 = * pho.pfpVsIso3th2;
        phoSkim.pfpVsIso4th2 = * pho.pfpVsIso4th2;
        phoSkim.pfpVsIso5th2 = * pho.pfpVsIso5th2;

        phoSkim.towerIso1 = * pho.towerIso1;
        phoSkim.towerIso2 = * pho.towerIso2;
        phoSkim.towerIso3 = * pho.towerIso3;
        phoSkim.towerIso4 = * pho.towerIso4;
        phoSkim.towerIso5 = * pho.towerIso5;
        phoSkim.towerVsIso1 = * pho.towerVsIso1;
        phoSkim.towerVsIso2 = * pho.towerVsIso2;
        phoSkim.towerVsIso3 = * pho.towerVsIso3;
        phoSkim.towerVsIso4 = * pho.towerVsIso4;
        phoSkim.towerVsIso5 = * pho.towerVsIso5;
        phoSkim.towerVsSubIso1 = * pho.towerVsSubIso1;
        phoSkim.towerVsSubIso2 = * pho.towerVsSubIso2;
        phoSkim.towerVsSubIso3 = * pho.towerVsSubIso3;
        phoSkim.towerVsSubIso4 = * pho.towerVsSubIso4;
        phoSkim.towerVsSubIso5 = * pho.towerVsSubIso5;

        if(montecarlo){
            phoSkim.nMC = pho.nMC;
            phoSkim.mcPID = * pho.mcPID;
            phoSkim.mcStatus = * pho.mcStatus;
            phoSkim.mcVtx_x = * pho.mcVtx_x;
            phoSkim.mcVtx_y = * pho.mcVtx_y;
            phoSkim.mcVtx_z = * pho.mcVtx_z;
            phoSkim.mcPt = * pho.mcPt;
            phoSkim.mcEta = * pho.mcEta;
            phoSkim.mcPhi = * pho.mcPhi;
            phoSkim.mcE = * pho.mcE;
            phoSkim.mcEt = * pho.mcEt;
            phoSkim.mcMass = * pho.mcMass;
            phoSkim.mcParentage = * pho.mcParentage;
            phoSkim.mcMomPID = * pho.mcMomPID;
            phoSkim.mcMomPt = * pho.mcMomPt;
            phoSkim.mcMomEta = * pho.mcMomEta;
            phoSkim.mcMomPhi = * pho.mcMomPhi;
            phoSkim.mcMomMass = * pho.mcMomMass;
            phoSkim.mcGMomPID = * pho.mcGMomPID;
            phoSkim.mcIndex = * pho.mcIndex;
            phoSkim.mcCalIsoDR03 = * pho.mcCalIsoDR03;
            phoSkim.mcCalIsoDR04 = * pho.mcCalIsoDR04;
            phoSkim.mcTrkIsoDR03 = * pho.mcTrkIsoDR03;
            phoSkim.mcTrkIsoDR04 = * pho.mcTrkIsoDR04;
            phoSkim.pho_genMatchedIndex = * pho.pho_genMatchedIndex;
        }



        outputTreeHLT->Fill();
        outputTreeggHiNtuplizer->Fill();
        outputTreeHiEvt->Fill();
        outputTreeSkim->Fill();

    }

    outFile->cd();
    outputTreeHLT->Write();
    outputTreeggHiNtuplizer->Write();
    outputTreeHiEvt->Write();
    outputTreeSkim->Write();
    configTree->Write();
    //outFile->Close();

    std::cout<< "Loop ENDED : " << treePath.c_str() <<std::endl;
    std::cout << "entries            = " << nentries << std::endl;
    std::cout << "duplicateEntries   = " << duplicateEntries << std::endl;
    std::cout << "entriesAnalyzed    = " << entriesAnalyzed << std::endl;
    std::cout << "outputTreeHLT->GetEntries()           = " << outputTreeHLT->GetEntries() << std::endl;
    std::cout << "outputTreeggHiNtuplizer->GetEntries() = " << outputTreeggHiNtuplizer->GetEntries() << std::endl;
    std::cout << "outputTreeHiEvt->GetEntries() = " << outputTreeHiEvt->GetEntries() << std::endl;

    //! Check
    timer.Stop();
    double rtime  = timer.RealTime();
    double ctime  = timer.CpuTime();
    std::cout << "==== Time Info ====" << std::endl;
    std::cout << "real time : " << rtime << std::endl;
    std::cout << "cpu time : " << ctime << std::endl;

}

int main(int argc, char **argv)
{
    if(argc == 4)
    {
        ggHiNtuplizerSkim(argv[1], argv[2], argv[3]);
        return 0;
    } else if(argc == 3)
    {
        ggHiNtuplizerSkim(argv[1], argv[2]);
        return 0;
    } else {
        std::cout << "Usage : \n" <<
            "./ggHiNtuplizerSkim.exe <configFile> <inputHiForest> <outputSkimFile>"
            << std::endl;
        return 1;
    }
}
