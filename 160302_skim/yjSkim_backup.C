///////////////////////////////////////////////////////////////////                                
// photonMCSkim.C                                                //                                                 
// Creator : Yeonju Go (Korea Univ.), ygo@cern.ch                //                                                 
// Function : Skim photon samples                                //
// (DATA/AllQCDPhoton/EmEnrichedDijet)                           //
// calculate pthatWeight factor and add it to the hiEvtTree      // 
// modified from the xSecCal.C                                   //
// input should be a merged pthat samples                        //
///////////////////////////////////////////////////////////////////         
/////////// There was an critical Segmentation error on this code!!
/////////// It was found out that the vector have to be initialized by 1  
/////////// Don't forget to initialize vector !!!!!!!!!!!

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <TMath.h>
#include "../yjUtility.h"
#include "../phoRaaCuts/phoRaaCuts_v2.h"
#include "../ElectroWeak-Jet-Track-Analyses/CorrelationTuple/EventMatcher.h"
#include <time.h>
using namespace std;
static const long MAXTREESIZE = 10000000000;

float xSecCal(const char* fname_lowestPthat, const char* fname_merged, float pthat_i = 30, float pthat_f = 50);

void photonRaaSkim(/*const TString configFile,*/ sampleType colli=kPPMC) {

    std::cout<<"running photonRaaSkim()"<<std::endl;
    //std::cout<<"configFile  = "<< configFile.Data() <<std::endl;
    //std::cout<<"inputFile   = "<< inputFile.Data() <<std::endl;
    //std::cout<<"outputFile  = "<< outputFile.Data() <<std::endl;
    //std::cout<<"minBiasJetSkimFile = "<< minBiasJetSkimFile.Data() <<std::endl;

#if 0
    InputConfiguration configInput = InputConfigurationParser::Parse(configFile.Data());
    CutConfiguration configCuts = CutConfigurationParser::Parse(configFile.Data());

    // cut configuration
    float cut_vz;
    int cut_pcollisionEventSelection;
    int cut_pPAprimaryVertexFilter;
    int cut_pBeamScrapingFilter;

    std::vector<std::string> jetCollections;
    float cutPhoEt;
    float cutPhoEta;

    int doMix;
    int nMaxEvents_minBiasMixing;
    int nCentralityBins;
    int nVertexBins;
    int nEventsToMix;
    if (configCuts.isValid) {
        cut_vz = configCuts.proc[CUTS::kSKIM].obj[CUTS::kEVENT].f[CUTS::EVT::k_vz];
        cut_pcollisionEventSelection = configCuts.proc[CUTS::kSKIM].obj[CUTS::kEVENT].i[CUTS::EVT::k_pcollisionEventSelection];
        cut_pPAprimaryVertexFilter = configCuts.proc[CUTS::kSKIM].obj[CUTS::kEVENT].i[CUTS::EVT::k_pPAprimaryVertexFilter];
        cut_pBeamScrapingFilter = configCuts.proc[CUTS::kSKIM].obj[CUTS::kEVENT].i[CUTS::EVT::k_pBeamScrapingFilter];

        //jetCollections = ConfigurationParser::ParseList(configCuts.proc[CUTS::kSKIM].obj[CUTS::kJET].s[CUTS::JET::k_jetCollection]);

        cutPhoEt = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].f[CUTS::PHO::k_et];
        cutPhoEta = configCuts.proc[CUTS::kSKIM].obj[CUTS::kPHOTON].f[CUTS::PHO::k_eta];

        //doMix = configCuts.proc[CUTS::kSKIM].obj[CUTS::kGAMMAJET].i[CUTS::GJT::k_doMix];
        //nMaxEvents_minBiasMixing = configCuts.proc[CUTS::kSKIM].obj[CUTS::kGAMMAJET].i[CUTS::GJT::k_nMaxEvents_minBiasMixing];
        //nCentralityBins = configCuts.proc[CUTS::kSKIM].obj[CUTS::kGAMMAJET].i[CUTS::GJT::k_nCentralityBins];
        //nVertexBins = configCuts.proc[CUTS::kSKIM].obj[CUTS::kGAMMAJET].i[CUTS::GJT::k_nVertexBins];
        //nEventsToMix = configCuts.proc[CUTS::kSKIM].obj[CUTS::kGAMMAJET].i[CUTS::GJT::k_nEventsToMix];
    }
    else {
        cut_vz = 15;
        cut_pcollisionEventSelection = 1;
        cut_pPAprimaryVertexFilter = 1;
        cut_pBeamScrapingFilter = 1;

        cutPhoEt = 15;
        cutPhoEta = 3.00;

        // default : no mixing
        doMix = 0;
        nMaxEvents_minBiasMixing = 0;
        nCentralityBins = 0;
        nVertexBins = 0;
        nEventsToMix = 0;
    }
#endif
    bool isMC = true;
    if ((colli==kPPDATA)||(colli==kHIDATA)) isMC=false;
    int nPthat = 5;
    if(colli==kPPMCEmEnr) nPthat=4;
    const char* lowestPthatFileName;
    const char* fileName;
    const char* outFileName; 
    float pthatCut[nPthat];

    if(colli==kPPMC) {
        float temp[] = {15,30,50,80,120,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = temp[j];     
        }
        lowestPthatFileName = "/home/goyeonju/CMS/Files/photon2016/officialMC/HINppWinter16DR/Pythia8_Photon15_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1.root";
        fileName= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINppWinter16DR/Pythia8_Photon15_30_50_80_120_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1.root";
        outFileName= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINppWinter16DR/Pythia8_Photon15_30_50_80_120_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_forest_v1_with_pthatWeight.root";
    } else if(colli==kHIMC) {
        float temp[] = {15,30,50,80,120,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = temp[j];     
        }
        lowestPthatFileName = "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_Photon15_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        fileName= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_Photon15_30_50_80_120_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        outFileName= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_Photon15_30_50_80_120_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1_with_pthatWeight.root";

    } else if(colli==kPPMCEmEnr) {
        float temp[] = {50,80,120,170,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = temp[j];         
        }
        lowestPthatFileName = "Pythia8_EmEnrDijet50_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_forest_v1.root";
        fileName= "Pythia8_EmEnrDijet50_80_120_170_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_forest_v1.root";
        outFileName= "Pythia8_EmEnrDijet50_80_120_170_pp502_TuneCUETP8M1-HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_forest_v1_with_pthatWeight.root";
    } else if(colli==kHIMCEmEnr){
        float tmpPthat[] = {30,50,80,120,170,9999};
        for(int j=0;j<nPthat+1;j++){
            pthatCut[j] = tmpPthat[j];
        }
        lowestPthatFileName = "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_EmEnrDijet30_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        fileName= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_EmEnrDijet30_50_80_120_170_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1.root";
        outFileName= "/home/goyeonju/CMS/Files/photon2016/officialMC/HINPbPbWinter16DR/Pythia8_EmEnrDijet30_50_80_120_170_Hydjet_MB-HINPbPbWinter16DR-75X_mcRun2_HeavyIon_forest_v1_with_pthatWeight.root";
    } else if(colli==kPPDATA) {
        fileName = "/home/goyeonju/CMS/Files/photon2016/2015-Data-promptRECO-photonSkims_pp-photonHLTFilter-v0-HiForest.root";
        outFileName = "/home/goyeonju/CMS/Files/photon2016/forestSkimed_2015-Data-promptRECO-photonSkims_pp-photonHLTFilter-v0-HiForest.root";
    }

    TFile* f = new TFile(fileName);
    TTree* hiforest = (TTree*) f ->Get("HiForest/HiForestInfo");  
    TTree* gen;
    if(isMC) gen = (TTree*) f ->Get("HiGenParticleAna/hi");  

    TTree* evt = (TTree*) f ->Get("hiEvtAnalyzer/HiTree");  
    float pthat, pthatWeight;
    ULong64_t event;
    unsigned run;
    unsigned lumi;
    int hiBin;
    float hiHF;
    if(isMC) evt->Branch("pthatWeight", &pthatWeight,"pthatWeight/F");
    if(isMC) evt->SetBranchAddress("pthat",&pthat);
    evt->SetBranchAddress("evt",&event);
    evt->SetBranchAddress("run",&run);
    evt->SetBranchAddress("lumi",&lumi);
    evt->SetBranchAddress("hiBin",&hiBin);
    evt->SetBranchAddress("hiHF",&hiHF);
    
    TTree* hlt = (TTree*) f ->Get("hltanalysis/HltTree");  
    hlt->SetBranchStatus("*",0);     // disable all branches
    hlt->SetBranchStatus("HLT_HI*SinglePhoton*Eta*v1*",1);     // enable photon branches
    hlt->SetBranchStatus("HLT_HI*DoublePhoton*Eta*v1*",1);     // enable photon branches

    TTree* skim = (TTree*) f ->Get("skimanalysis/HltTree");  
    Int_t HBHENoiseFilterResult;
    Int_t pcollisionEventSelection;
    Int_t pPAprimaryVertexFilter;    // this filter is used for PP.
    Int_t pBeamScrapingFilter;   // this filter is used for PP.
    if(colli==kHIDATA){
        skim->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
        skim->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
    } else if(colli==kPPDATA){
        skim->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
        skim->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
    }

#if 1
    TTree* pho;
    if((colli==kPPDATA) || (colli==kPPMC) || (colli==kPPMCEmEnr)) pho = (TTree*) f ->Get("ggHiNtuplizerGED/EventTree");  
    else if((colli==kHIDATA) || (colli==kHIMC) || (colli==kHIMCEmEnr)) pho = (TTree*) f ->Get("ggHiNtuplizer/EventTree");  
    int nPhos;
    std::vector<float>*  phoEt=0;
    std::vector<float>*  phoEta=0;
    std::vector<float>*  phoHoverE=0;
    std::vector<float>*  phoSigmaIEtaIEta_2012=0;
    std::vector<float>*  pho_swissCrx=0;
    std::vector<float>*  pho_seedTime=0;
    TBranch * b_nPhos=0;
    pho->SetBranchStatus("*",0);
    pho->SetBranchStatus("nPho",1);
    pho->SetBranchStatus("pho*",1);
    pho->SetBranchStatus("pf*",1);
    pho->SetBranchStatus("tower*",1);
    if(isMC) pho->SetBranchStatus("mc*",1);
    pho->SetBranchAddress("nPho", &nPhos,&b_nPhos);
    pho->SetBranchAddress("phoEt", &phoEt);
    pho->SetBranchAddress("phoEta", &phoEta);
    pho->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012);
    pho->SetBranchAddress("phoHoverE", &phoHoverE);
    pho->SetBranchAddress("pho_swissCrx", &pho_swissCrx);
    pho->SetBranchAddress("pho_seedTime", &pho_seedTime);
    cout << "pho tree " << pho << endl;
#endif


    TFile* fout = new TFile(outFileName,"recreate");
    TTree *nhlt = hlt->CloneTree(0);
    TTree *ntpho = pho->CloneTree(0);
    TTree *nskim = skim->CloneTree(0);
    TTree *nevt = evt->CloneTree(0);
    TTree *nhiforest = hiforest->CloneTree(0);
    TTree *ngen;
    if(isMC) ngen = gen->CloneTree(0);
    nskim->SetName("skimTree");

    float tmpWeight[nPthat];
    if(isMC) {
        for(int j=0; j<nPthat ; j++){
            tmpWeight[j] = xSecCal(lowestPthatFileName,fileName, pthatCut[j], pthatCut[j+1]);
            cout << getSampleName(colli) << ", pthatWeight of " << pthatCut[j] << " to " << pthatCut[j+1] << " = " << tmpWeight[j] << endl;
        }
    }

    /////// Event Matching for DATA ///////
    EventMatcher* em = new EventMatcher();
    Long64_t duplicateEntries = 0;
    Long64_t nentries = pho->GetEntries();
    Long64_t entriesAnalyzed = 0;
    std::cout << "entries = " << nentries << std::endl;
    bool eventAdded;


    /////// Event Loop ///////
    //for (Long64_t jentry = 0 ; jentry < 20; jentry++) {
    for (Long64_t jentry = 0 ; jentry < nentries; jentry++) {
        if (jentry% 10000 == 0)  {
            cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
        }
        pho->GetEntry(jentry);
        hlt->GetEntry(jentry);
        skim->GetEntry(jentry);
        evt->GetEntry(jentry);
        hiforest->GetEntry(jentry);
        if(isMC) { 
            //cout << " ### calculate pthat ### " << endl;
            gen->GetEntry(jentry);
            if((pthat>=pthatCut[0]) && (pthat<pthatCut[1])) pthatWeight = tmpWeight[0]; 
            else if((pthat>=pthatCut[1]) && (pthat<pthatCut[2])) pthatWeight = tmpWeight[1]; 
            else if((pthat>=pthatCut[2]) && (pthat<pthatCut[3])) pthatWeight = tmpWeight[2]; 
            else if((pthat>=pthatCut[3]) && (pthat<pthatCut[4])) pthatWeight = tmpWeight[3]; 
            else if((pthat>=pthatCut[4]) && (pthat<pthatCut[5])) pthatWeight = tmpWeight[4]; 
            else continue;
        }
        if(!isMC) {
            eventAdded = em->addEvent(run,lumi,event,jentry);
            if(!eventAdded) // this event is duplicate, skip this one.
            {
                duplicateEntries++;
                continue;
            }
        }
        if(colli==kHIDATA) {
             if ((pcollisionEventSelection < 1))  continue;
        } else if(colli==kPPDATA) {
            if ((pPAprimaryVertexFilter < 1) || (pBeamScrapingFilter < 1))  continue;
        }
#if 1
        int phoIdx = -1;     // index of the leading photon
        double maxPhoEt = -1;
        for(int i=0; i<nPhos;i++){
            bool failedEtCut  = (phoEt->at(i) < 15.0) ;
            bool failedEtaCut = (TMath::Abs(phoEta->at(i)) > 3) ;
            bool failedSpikeRejection;
            if (colli==kPPDATA) {
                failedSpikeRejection = (phoSigmaIEtaIEta_2012->at(i) < 0.002 ||
                        pho_swissCrx->at(i)     > 0.9   ||
                        TMath::Abs(pho_seedTime->at(i)) > 3.);
            }
            else {
                failedSpikeRejection = (phoSigmaIEtaIEta_2012->at(i) < 0.002);
            }

            bool failedHoverE = (phoHoverE->at(i) > 0.2);      // <0.1 cut is applied after corrections
            //               bool failedEnergyRatio = ((float)ggHi.phoSCRawE->at(i)/ggHi.phoE->at(i) < 0.5);

            if (failedEtCut)          continue;
            if (failedEtaCut)         continue;
            if (failedSpikeRejection) continue;
            if (failedHoverE)         continue;
            //               if (failedEnergyRatio)    continue;    // actually applied after corrections
            if (phoEt->at(i) > maxPhoEt)
            {
                maxPhoEt = phoEt->at(i);
                phoIdx = i;
            }
        }
        if (phoIdx == -1) continue;
#endif
        entriesAnalyzed++;
        nhlt->Fill();
        ntpho->Fill();
        nskim->Fill();
        nevt->Fill();
        if(isMC) ngen->Fill();
        nhiforest->Fill();
    }
    cout << entriesAnalyzed << " events are analyzed " << endl;
    fout ->cd();
    nhlt->Write();
    ntpho->Write();
    nskim->Write();
    nevt->Write();
    if(isMC) ngen->Write();
    nhiforest->Write();
}
float xSecCal(const char* fname_lowestPthat, const char* fname_merged, float pthat_i, float pthat_f){
    const int nFile = 2;
    const char *fileName_[nFile];
    fileName_[0] = fname_lowestPthat;
    fileName_[1] = fname_merged;

    TFile *fin[nFile];
    TTree *t[nFile];
    int entries[nFile];
    for(int i=0; i<nFile ; i++){
        entries[i]=0.0;
    }
    for(int ifile=0; ifile<nFile; ifile++){
        fin[ifile] = new TFile(fileName_[ifile]);
        t[ifile] = (TTree*) fin[ifile] -> Get("hiEvtAnalyzer/HiTree");
        Float_t pthat;
        TBranch *b_pthat;
        t[ifile]->SetBranchAddress("pthat",&pthat, &b_pthat);
        entries[ifile] = t[ifile]->GetEntries(Form("pthat>= %.3f && pthat< %.3f", pthat_i, pthat_f));
        cout << "entries " << ifile << " : " << entries[ifile] << endl;
    }
    float weight = (double)entries[0]/(double)entries[1];
    return weight;
}

