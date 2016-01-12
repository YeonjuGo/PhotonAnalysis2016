/*
 * macro to study different photon Reconstruction algorithms
 * modifiedy by Yeonju
 * */
#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif

#include "gedPhotonUtility.h" 
static const long MAXTREESIZE = 10000000000;

void makeGenMatchedPhotonTree_sep(const char* hiForestfileName="/afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_7_5_0/src/PFphoton/nominalForest/crab_AllQCDPhoton30/results/merged_AllQCDPhoton30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV_2nd.root",
        TString outName="AllQCDPhoton30", TString treePath="ggHiNtuplizer", bool isPromptPho=1, float ptThr=20)
{

    TFile* inputFile = new TFile(hiForestfileName, "READ");
    std::cout << "input HiForest : " << inputFile->GetName() << std::endl;
    TString outputFileName = Form("skimFiles/jskim_%s_%s_genMatched_inclusivePho.root",outName.Data(),treePath.Data());
    //TString outputFileName = Form("skimFiles/jskim_%s_%s_isPromptPho%d_genIsoDR3_5GeV.root",outName.Data(),treePath.Data(),(int)isPromptPho);
    TFile* fout = new TFile(outputFileName, "RECREATE");
    fout->cd();


    TTree* hiEvtAnalyzerTree = (TTree*)inputFile->Get("hiEvtAnalyzer/HiTree");
    TTree* ggHiNtuplizerTree = (TTree*)inputFile->Get(Form("%s/EventTree",treePath.Data()));

    Int_t hiBin, hiNevtPlane;
    Float_t hiEvtPlanes[50];
    hiEvtAnalyzerTree->SetBranchAddress("hiBin", &hiBin);
    hiEvtAnalyzerTree->SetBranchAddress("hiNevtPlane", &hiNevtPlane);
    hiEvtAnalyzerTree->SetBranchAddress("hiEvtPlanes", hiEvtPlanes);

    // mc Info.
    Int_t nMC;
    std::vector<float>* mcPt=0;
    std::vector<float>* mcEt=0;
    std::vector<float>* mcEta=0;
    std::vector<float>* mcPhi=0;
    std::vector<int>*   mcPID=0;
    std::vector<int>*   mcMomPID=0;
    std::vector<int>*   mcGMomPID=0;
    std::vector<int>*   mcStatus=0;
    std::vector<int>*   mcCalIsoDR03=0;
    std::vector<int>*   mcCalIsoDR04=0;
    ggHiNtuplizerTree->SetBranchAddress("nMC",&nMC);
    ggHiNtuplizerTree->SetBranchAddress("mcPt",&mcPt);
    ggHiNtuplizerTree->SetBranchAddress("mcEt",&mcPt);
    ggHiNtuplizerTree->SetBranchAddress("mcEta",&mcEta);
    ggHiNtuplizerTree->SetBranchAddress("mcPhi",&mcPhi);
    ggHiNtuplizerTree->SetBranchAddress("mcPID",&mcPID);
    ggHiNtuplizerTree->SetBranchAddress("mcMomPID",&mcMomPID);
    ggHiNtuplizerTree->SetBranchAddress("mcGMomPID",&mcGMomPID);
    ggHiNtuplizerTree->SetBranchAddress("mcStatus",&mcStatus);
    ggHiNtuplizerTree->SetBranchAddress("mcCalIsoDR03",&mcCalIsoDR03);
    ggHiNtuplizerTree->SetBranchAddress("mcCalIsoDR04",&mcCalIsoDR04);

    // RECO photons
    Int_t nPho;
    std::vector<float>* phoEt=0;
    std::vector<float>* phoEta=0;
    std::vector<float>* phoPhi=0;
    std::vector<float>* pho_ecalClusterIsoR2=0;
    std::vector<float>* pho_ecalClusterIsoR3=0;
    std::vector<float>* pho_ecalClusterIsoR4=0;
    std::vector<float>* pho_ecalClusterIsoR5=0;
    std::vector<float>* pho_hcalRechitIsoR2=0;
    std::vector<float>* pho_hcalRechitIsoR3=0;
    std::vector<float>* pho_hcalRechitIsoR4=0;
    std::vector<float>* pho_hcalRechitIsoR5=0;
    std::vector<float>* pho_trackIsoR2PtCut20=0;
    std::vector<float>* pho_trackIsoR3PtCut20=0;
    std::vector<float>* pho_trackIsoR4PtCut20=0;
    std::vector<float>* pho_trackIsoR5PtCut20=0;
    std::vector<float>* phoR9=0;
    std::vector<float>* phoHoverE=0;
    std::vector<float>* phoSigmaIEtaIEta=0;
    std::vector<float>* phoSCE=0;
    std::vector<float>* phoSCEta=0;
    std::vector<float>* phoSCPhi=0;
    std::vector<float>* phoSCEtaWidth=0;
    std::vector<float>* phoSCPhiWidth=0;
    std::vector<float>* phoSCBrem=0;
    std::vector<int>* pho_genMatchedIndex=0;

    ggHiNtuplizerTree->SetBranchAddress("nPho",&nPho);
    ggHiNtuplizerTree->SetBranchAddress("phoEt",&phoEt);
    ggHiNtuplizerTree->SetBranchAddress("phoEta",&phoEta);
    ggHiNtuplizerTree->SetBranchAddress("phoPhi",&phoPhi);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR2",&pho_ecalClusterIsoR2);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR3",&pho_ecalClusterIsoR3);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR4",&pho_ecalClusterIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_ecalClusterIsoR5",&pho_ecalClusterIsoR5);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR2",&pho_hcalRechitIsoR2);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR3",&pho_hcalRechitIsoR3);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR4",&pho_hcalRechitIsoR4);
    ggHiNtuplizerTree->SetBranchAddress("pho_hcalRechitIsoR5",&pho_hcalRechitIsoR5);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR2PtCut20",&pho_trackIsoR2PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR3PtCut20",&pho_trackIsoR3PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR4PtCut20",&pho_trackIsoR4PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("pho_trackIsoR5PtCut20",&pho_trackIsoR5PtCut20);
    ggHiNtuplizerTree->SetBranchAddress("phoR9",&phoR9);
    ggHiNtuplizerTree->SetBranchAddress("phoHoverE",&phoHoverE);
    ggHiNtuplizerTree->SetBranchAddress("phoSigmaIEtaIEta",&phoSigmaIEtaIEta);
    ggHiNtuplizerTree->SetBranchAddress("phoSCE",&phoSCE);
    ggHiNtuplizerTree->SetBranchAddress("phoSCEta",&phoSCEta);
    ggHiNtuplizerTree->SetBranchAddress("phoSCPhi",&phoSCPhi);
    ggHiNtuplizerTree->SetBranchAddress("phoSCEtaWidth",&phoSCEtaWidth);
    ggHiNtuplizerTree->SetBranchAddress("phoSCPhiWidth",&phoSCPhiWidth);
    ggHiNtuplizerTree->SetBranchAddress("phoSCBrem",&phoSCBrem);
    ggHiNtuplizerTree->SetBranchAddress("pho_genMatchedIndex",&pho_genMatchedIndex);

    // RECO photons (pfIso)
    std::vector<float>* pfcIso1=0; 
    std::vector<float>* pfcIso2=0; 
    std::vector<float>* pfcIso3=0; 
    std::vector<float>* pfcIso4=0; 
    std::vector<float>* pfcIso5=0; 
    std::vector<float>* pfpIso1=0; 
    std::vector<float>* pfpIso2=0; 
    std::vector<float>* pfpIso3=0; 
    std::vector<float>* pfpIso4=0; 
    std::vector<float>* pfpIso5=0; 
    std::vector<float>* pfnIso1=0; 
    std::vector<float>* pfnIso2=0; 
    std::vector<float>* pfnIso3=0; 
    std::vector<float>* pfnIso4=0; 
    std::vector<float>* pfnIso5=0; 
    std::vector<float>* pfcVsIso1=0; 
    std::vector<float>* pfcVsIso2=0; 
    std::vector<float>* pfcVsIso3=0; 
    std::vector<float>* pfcVsIso4=0; 
    std::vector<float>* pfcVsIso5=0; 
    std::vector<float>* pfpVsIso1=0; 
    std::vector<float>* pfpVsIso2=0; 
    std::vector<float>* pfpVsIso3=0; 
    std::vector<float>* pfpVsIso4=0; 
    std::vector<float>* pfpVsIso5=0; 
    std::vector<float>* pfnVsIso1=0; 
    std::vector<float>* pfnVsIso2=0; 
    std::vector<float>* pfnVsIso3=0; 
    std::vector<float>* pfnVsIso4=0; 
    std::vector<float>* pfnVsIso5=0;
    ggHiNtuplizerTree->SetBranchAddress("pfcIso1", &pfcIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso2", &pfcIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso3", &pfcIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso4", &pfcIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfcIso5", &pfcIso5); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso1", &pfpIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso2", &pfpIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso3", &pfpIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso4", &pfpIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfpIso5", &pfpIso5); 	
    ggHiNtuplizerTree->SetBranchAddress("pfnIso1", &pfnIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso2", &pfnIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso3", &pfnIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso4", &pfnIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfnIso5", &pfnIso5); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso1", &pfcVsIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso2", &pfcVsIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso3", &pfcVsIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso4", &pfcVsIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfcVsIso5", &pfcVsIso5); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso1", &pfpVsIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso2", &pfpVsIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso3", &pfpVsIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso4", &pfpVsIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfpVsIso5", &pfpVsIso5); 	
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso1", &pfnVsIso1); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso2", &pfnVsIso2); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso3", &pfnVsIso3); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso4", &pfnVsIso4); 
    ggHiNtuplizerTree->SetBranchAddress("pfnVsIso5", &pfnVsIso5); 

    Int_t nPho_;
    Float_t hiEvtPlane1_, hiEvtPlane2_, hiEvtPlane3_, hiEvtPlane4_;
    std::vector<int>* mc_recoMatchedIndex_=0;
    std::vector<float> *mcPt_ = new std::vector<float>();
    std::vector<float> *mcPhi_ = new std::vector<float>();
    std::vector<float> *mcEta_ = new std::vector<float>();
    std::vector<int> *mcPID_ = new std::vector<int>();
    std::vector<int> *mcMomPID_ = new std::vector<int>();
    std::vector<int> *mcGMomPID_ = new std::vector<int>();
    std::vector<int> *mcStatus_= new std::vector<int>();
    std::vector<float> *mcCalIsoDR03_= new std::vector<float>();
    std::vector<float> *mcCalIsoDR04_= new std::vector<float>();
    std::vector<float> *mcDPhi_evtpl1_ = new std::vector<float>();
    std::vector<float> *mcDPhi_evtpl2_ = new std::vector<float>();
    std::vector<float> *mcDPhi_evtpl3_ = new std::vector<float>();
    std::vector<float> *mcDPhi_evtpl4_ = new std::vector<float>();
    std::vector<float> *mc_recoMatchedIndex_= new std::vector<float>();
    std::vector<float> *phoEt_= new std::vector<float>();
    std::vector<float> *phoEta_= new std::vector<float>();
    std::vector<float> *phoPhi_= new std::vector<float>();
    std::vector<float> *phoDPhi_evtpl1_ = new std::vector<float>();
    std::vector<float> *phoDPhi_evtpl2_ = new std::vector<float>();
    std::vector<float> *phoDPhi_evtpl3_ = new std::vector<float>();
    std::vector<float> *phoDPhi_evtpl4_ = new std::vector<float>();
    std::vector<float> *pho_ecalClusterIsoR2_= new std::vector<float>();
    std::vector<float> *pho_ecalClusterIsoR3_= new std::vector<float>();
    std::vector<float> *pho_ecalClusterIsoR4_= new std::vector<float>();
    std::vector<float> *pho_ecalClusterIsoR5_= new std::vector<float>();
    std::vector<float> *pho_hcalRechitIsoR2_= new std::vector<float>();
    std::vector<float> *pho_hcalRechitIsoR3_= new std::vector<float>();
    std::vector<float> *pho_hcalRechitIsoR4_= new std::vector<float>();
    std::vector<float> *pho_hcalRechitIsoR5_= new std::vector<float>();
    std::vector<float> *pho_trackIsoR2PtCut20_= new std::vector<float>();
    std::vector<float> *pho_trackIsoR3PtCut20_= new std::vector<float>();
    std::vector<float> *pho_trackIsoR4PtCut20_= new std::vector<float>();
    std::vector<float> *pho_trackIsoR5PtCut20_= new std::vector<float>();
    std::vector<float> *pho_genMatchedIndex_= new std::vector<float>();
    std::vector<float> *sumIsoR2_= new std::vector<float>();
    std::vector<float> *sumIsoR3_= new std::vector<float>();
    std::vector<float> *sumIsoR4_= new std::vector<float>();
    std::vector<float> *sumIsoR5_= new std::vector<float>();
    std::vector<float> *phoR9_= new std::vector<float>();
    std::vector<float> *phoHoverE_= new std::vector<float>();
    std::vector<float> *phoSigmaIEtaIEta_= new std::vector<float>();
    std::vector<float> *phoSCE_= new std::vector<float>();
    std::vector<float> *phoSCEta_= new std::vector<float>();
    std::vector<float> *phoSCPhi_= new std::vector<float>();
    std::vector<float> *phoSCEtaWidth_= new std::vector<float>();
    std::vector<float> *phoSCPhiWidth_= new std::vector<float>();
    std::vector<float> *phoSCBrem_= new std::vector<float>();
    std::vector<float>* pfcIso1_= new std::vector<float>();
    std::vector<float>* pfcIso2_= new std::vector<float>();
    std::vector<float>* pfcIso3_= new std::vector<float>();
    std::vector<float>* pfcIso4_= new std::vector<float>();
    std::vector<float>* pfcIso5_= new std::vector<float>();
    std::vector<float>* pfpIso1_= new std::vector<float>();
    std::vector<float>* pfpIso2_= new std::vector<float>();
    std::vector<float>* pfpIso3_= new std::vector<float>();
    std::vector<float>* pfpIso4_= new std::vector<float>();
    std::vector<float>* pfpIso5_= new std::vector<float>();
    std::vector<float>* pfnIso1_= new std::vector<float>();
    std::vector<float>* pfnIso2_= new std::vector<float>();
    std::vector<float>* pfnIso3_= new std::vector<float>();
    std::vector<float>* pfnIso4_= new std::vector<float>();
    std::vector<float>* pfnIso5_= new std::vector<float>();
    std::vector<float>* pfSumIso2_= new std::vector<float>();
    std::vector<float>* pfSumIso3_= new std::vector<float>();
    std::vector<float>* pfSumIso4_= new std::vector<float>();
    std::vector<float>* pfSumIso5_= new std::vector<float>();
    std::vector<float>* pfcVsIso1_= new std::vector<float>();
    std::vector<float>* pfcVsIso2_= new std::vector<float>();
    std::vector<float>* pfcVsIso3_= new std::vector<float>();
    std::vector<float>* pfcVsIso4_= new std::vector<float>();
    std::vector<float>* pfcVsIso5_= new std::vector<float>();
    std::vector<float>* pfpVsIso1_= new std::vector<float>();
    std::vector<float>* pfpVsIso2_= new std::vector<float>();
    std::vector<float>* pfpVsIso3_= new std::vector<float>();
    std::vector<float>* pfpVsIso4_= new std::vector<float>();
    std::vector<float>* pfpVsIso5_= new std::vector<float>();
    std::vector<float>* pfnVsIso1_= new std::vector<float>();
    std::vector<float>* pfnVsIso2_= new std::vector<float>();
    std::vector<float>* pfnVsIso3_= new std::vector<float>();
    std::vector<float>* pfnVsIso4_= new std::vector<float>();
    std::vector<float>* pfnVsIso5_= new std::vector<float>();
    std::vector<float>* pfSumVsIso2_= new std::vector<float>();
    std::vector<float>* pfSumVsIso3_= new std::vector<float>();
    std::vector<float>* pfSumVsIso4_= new std::vector<float>();
    std::vector<float>* pfSumVsIso5_= new std::vector<float>();
 
    TTree* t_pho = new TTree("t_pho","gen matched tree distinguished by momId");
    t_pho->SetMaxTreeSize(MAXTREESIZE);
    t_pho->Branch("nPho",&nPho_,"nPho_/I");
    t_pho->Branch("hiBin",&hiBin,"hiBin/I");
    t_pho->Branch("hiEvtPlane1",&hiEvtPlane1_,"hiEvtPlane1_/F");
    t_pho->Branch("hiEvtPlane2",&hiEvtPlane2_,"hiEvtPlane2_/F");
    t_pho->Branch("hiEvtPlane3",&hiEvtPlane3_,"hiEvtPlane3_/F");
    t_pho->Branch("hiEvtPlane4",&hiEvtPlane4_,"hiEvtPlane4_/F");
    t_pho->Branch("mcPt","vector<float>",&mcPt_);
    t_pho->Branch("mcEta","vector<float>",&mcEta_);
    t_pho->Branch("mcPhi","vector<float>",&mcPhi_);
    t_pho->Branch("mcDPhi_evtpl1","vector<float>",&mcDPhi_evtpl1_);
    t_pho->Branch("mcDPhi_evtpl2","vector<float>",&mcDPhi_evtpl2_);
    t_pho->Branch("mcDPhi_evtpl3","vector<float>",&mcDPhi_evtpl3_);
    t_pho->Branch("mcDPhi_evtpl4","vector<float>",&mcDPhi_evtpl4_);
    t_pho->Branch("mcPID","vector<int>",&mcPID_);
    t_pho->Branch("mcMomPID","vector<int>",&mcMomPID_);
    t_pho->Branch("mcGMomPID","vector<int>",&mcGMomPID_);
    t_pho->Branch("mcStatus","vector<int>",&mcStatus_);
    t_pho->Branch("mcCalIsoDR03","vector<float>",&mcCalIsoDR03_);
    t_pho->Branch("mcCalIsoDR04","vector<float>",&mcCalIsoDR04_);
    t_pho->Branch("mc_recoMatchedIndex","vector<int>",&mc_recoMatchedIndex_);
    t_pho->Branch("phoEt","vector<float>",&phoEt_);
    t_pho->Branch("phoEta","vector<float>",&phoEta_);
    t_pho->Branch("phoDPhi_evtpl1","vector<float>",&phoDPhi_evtpl1_);
    t_pho->Branch("phoDPhi_evtpl2","vector<float>",&phoDPhi_evtpl2_);
    t_pho->Branch("phoDPhi_evtpl3","vector<float>",&phoDPhi_evtpl3_);
    t_pho->Branch("phoDPhi_evtpl4","vector<float>",&phoDPhi_evtpl4_);
    t_pho->Branch("phoPhi","vector<float>",&phoPhi_);
    t_pho->Branch("pho_ecalClusterIsoR2","vector<float>",&pho_ecalClusterIsoR2_);
    t_pho->Branch("pho_ecalClusterIsoR3","vector<float>",&pho_ecalClusterIsoR3_);
    t_pho->Branch("pho_ecalClusterIsoR4","vector<float>",&pho_ecalClusterIsoR4_);
    t_pho->Branch("pho_ecalClusterIsoR5","vector<float>",&pho_ecalClusterIsoR5_);
    t_pho->Branch("pho_hcalRechitIsoR2","vector<float>",&pho_hcalRechitIsoR2_);
    t_pho->Branch("pho_hcalRechitIsoR3","vector<float>",&pho_hcalRechitIsoR3_);
    t_pho->Branch("pho_hcalRechitIsoR4","vector<float>",&pho_hcalRechitIsoR4_);
    t_pho->Branch("pho_hcalRechitIsoR5","vector<float>",&pho_hcalRechitIsoR5_);
    t_pho->Branch("pho_trackIsoR2PtCut20","vector<float>",&pho_trackIsoR2PtCut20_);
    t_pho->Branch("pho_trackIsoR3PtCut20","vector<float>",&pho_trackIsoR3PtCut20_);
    t_pho->Branch("pho_trackIsoR4PtCut20","vector<float>",&pho_trackIsoR4PtCut20_);
    t_pho->Branch("pho_trackIsoR5PtCut20","vector<float>",&pho_trackIsoR5PtCut20_);
    t_pho->Branch("pho_genMatchedIndex","vector<int>",&pho_genMatchedIndex_);
    t_pho->Branch("sumIsoR2","vector<float>",&sumIsoR2_);
    t_pho->Branch("sumIsoR3","vector<float>",&sumIsoR3_);
    t_pho->Branch("sumIsoR4","vector<float>",&sumIsoR4_);
    t_pho->Branch("sumIsoR5","vector<float>",&sumIsoR5_);
    t_pho->Branch("phoR9","vector<float>",&phoR9_);
    t_pho->Branch("phoHoverE","vector<float>",&phoHoverE_);
    t_pho->Branch("phoSigmaIEtaIEta","vector<float>",&phoSigmaIEtaIEta_);
    t_pho->Branch("phoSCE","vector<float>",&phoSCE_);
    t_pho->Branch("phoSCEta","vector<float>",&phoSCEta_);
    t_pho->Branch("phoSCPhi","vector<float>",&phoSCPhi_);
    t_pho->Branch("phoSCEtaWidth","vector<float>",&phoSCEtaWidth_);
    t_pho->Branch("phoSCPhiWidth","vector<float>",&phoSCPhiWidth_);
    t_pho->Branch("phoSCBrem","vector<float>",&phoSCBrem_);
    t_pho->Branch("pfcIso1","vector<float>",&pfcIso1_);
    t_pho->Branch("pfcIso2","vector<float>",&pfcIso2_);
    t_pho->Branch("pfcIso3","vector<float>",&pfcIso3_);
    t_pho->Branch("pfcIso4","vector<float>",&pfcIso4_);
    t_pho->Branch("pfcIso5","vector<float>",&pfcIso5_);
    t_pho->Branch("pfnIso1","vector<float>",&pfnIso1_);
    t_pho->Branch("pfnIso2","vector<float>",&pfnIso2_);
    t_pho->Branch("pfnIso3","vector<float>",&pfnIso3_);
    t_pho->Branch("pfnIso4","vector<float>",&pfnIso4_);
    t_pho->Branch("pfnIso5","vector<float>",&pfnIso5_);
    t_pho->Branch("pfpIso1","vector<float>",&pfpIso1_);
    t_pho->Branch("pfpIso2","vector<float>",&pfpIso2_);
    t_pho->Branch("pfpIso3","vector<float>",&pfpIso3_);
    t_pho->Branch("pfpIso4","vector<float>",&pfpIso4_);
    t_pho->Branch("pfpIso5","vector<float>",&pfpIso5_);
    t_pho->Branch("pfSumIso2","vector<float>",&pfSumIso2_);
    t_pho->Branch("pfSumIso3","vector<float>",&pfSumIso3_);
    t_pho->Branch("pfSumIso4","vector<float>",&pfSumIso4_);
    t_pho->Branch("pfSumIso5","vector<float>",&pfSumIso5_);
    t_pho->Branch("pfcVsIso1","vector<float>",&pfcVsIso1_);
    t_pho->Branch("pfcVsIso2","vector<float>",&pfcVsIso2_);
    t_pho->Branch("pfcVsIso3","vector<float>",&pfcVsIso3_);
    t_pho->Branch("pfcVsIso4","vector<float>",&pfcVsIso4_);
    t_pho->Branch("pfcVsIso5","vector<float>",&pfcVsIso5_);
    t_pho->Branch("pfnVsIso1","vector<float>",&pfnVsIso1_);
    t_pho->Branch("pfnVsIso2","vector<float>",&pfnVsIso2_);
    t_pho->Branch("pfnVsIso3","vector<float>",&pfnVsIso3_);
    t_pho->Branch("pfnVsIso4","vector<float>",&pfnVsIso4_);
    t_pho->Branch("pfnVsIso5","vector<float>",&pfnVsIso5_);
    t_pho->Branch("pfpVsIso1","vector<float>",&pfpVsIso1_);
    t_pho->Branch("pfpVsIso2","vector<float>",&pfpVsIso2_);
    t_pho->Branch("pfpVsIso3","vector<float>",&pfpVsIso3_);
    t_pho->Branch("pfpVsIso4","vector<float>",&pfpVsIso4_);
    t_pho->Branch("pfpVsIso5","vector<float>",&pfpVsIso5_);
    t_pho->Branch("pfSumVsIso2","vector<float>",&pfSumVsIso2_);
    t_pho->Branch("pfSumVsIso3","vector<float>",&pfSumVsIso3_);
    t_pho->Branch("pfSumVsIso4","vector<float>",&pfSumVsIso4_);
    t_pho->Branch("pfSumVsIso5","vector<float>",&pfSumVsIso5_);

    TH1::SetDefaultSumw2();
    std::cout << "entering event loop" << std::endl;
    Long64_t entries = ggHiNtuplizerTree->GetEntries();
    std::cout << "number of entries = " << entries << std::endl;

    std::clock_t    start_loop, end_loop;
    start_loop = std::clock();
    for(Long64_t jj = 0; jj < entries; ++jj)
    {
        if (jj % 1000 == 0)  {
            std::cout << "current entry = " <<jj<<" out of "<<entries<<" : "<<std::setprecision(2)<<(double)jj/entries*100<<" %"<<std::endl;
        }

        hiEvtAnalyzerTree->GetEntry(jj);
        ggHiNtuplizerTree->GetEntry(jj);

        nPho_=0;
        hiEvtPlane1_=-99;
        hiEvtPlane2_=-99;
        hiEvtPlane3_=-99;
        hiEvtPlane4_=-99;
        mcPt_->clear(); 
        mcPhi_->clear();  
        mcDPhi_evtpl1_->clear();
        mcDPhi_evtpl2_->clear();
        mcDPhi_evtpl3_->clear();
        mcDPhi_evtpl4_->clear();
        mcEta_->clear(); 
        mcPID_->clear(); 
        mcMomPID_ ->clear(); 
        mcGMomPID_ ->clear(); 
        mcStatus_->clear(); 
        mcCalIsoDR03_->clear(); 
        mcCalIsoDR04_->clear(); 
        mc_recoMatchedIndex_->clear(); 
        phoEt_->clear();
        phoEta_->clear();
        phoPhi_->clear();
        phoDPhi_evtpl1_->clear();
        phoDPhi_evtpl2_->clear();
        phoDPhi_evtpl3_->clear();
        phoDPhi_evtpl4_->clear();
        pho_ecalClusterIsoR2_->clear();
        pho_ecalClusterIsoR3_->clear();
        pho_ecalClusterIsoR4_->clear();
        pho_ecalClusterIsoR5_->clear();
        pho_hcalRechitIsoR2_->clear();
        pho_hcalRechitIsoR3_->clear();
        pho_hcalRechitIsoR4_->clear();
        pho_hcalRechitIsoR5_->clear();
        pho_trackIsoR2PtCut20_->clear();
        pho_trackIsoR3PtCut20_->clear();
        pho_trackIsoR4PtCut20_->clear();
        sumIsoR2_->clear();
        sumIsoR3_->clear();
        sumIsoR4_->clear();
        sumIsoR5_->clear();
        phoR9_->clear();
        phoHoverE_->clear();
        phoSigmaIEtaIEta_->clear();
        phoSCE_->clear();
        phoSCEta_->clear();
        phoSCPhi_->clear();
        phoSCEtaWidth_->clear();
        phoSCPhiWidth_->clear();
        phoSCBrem_->clear();
        pho_genMatchedIndex_->clear(); 
        pfcIso2_->clear();
        pfcIso3_->clear();
        pfcIso4_->clear();
        pfcIso5_->clear();
        pfpIso1_->clear();
        pfpIso2_->clear();
        pfpIso3_->clear();
        pfpIso4_->clear();
        pfpIso5_->clear();
        pfnIso1_->clear();
        pfnIso2_->clear();
        pfnIso3_->clear();
        pfnIso4_->clear();
        pfnIso5_->clear();
        pfSumIso2_->clear();
        pfSumIso3_->clear();
        pfSumIso4_->clear();
        pfSumIso5_->clear();
        pfcVsIso1_->clear();
        pfcVsIso2_->clear();
        pfcVsIso3_->clear();
        pfcVsIso4_->clear();
        pfcVsIso5_->clear();
        pfpVsIso1_->clear();
        pfpVsIso2_->clear();
        pfpVsIso3_->clear();
        pfpVsIso4_->clear();
        pfpVsIso5_->clear();
        pfnVsIso1_->clear();
        pfnVsIso2_->clear();
        pfnVsIso3_->clear();
        pfnVsIso4_->clear();
        pfnVsIso5_->clear();
        pfSumVsIso2_->clear();
        pfSumVsIso3_->clear();
        pfSumVsIso4_->clear();
        pfSumVsIso5_->clear();

        for (int i=0; i < nPho; ++i)
        {
            bool passedPtThr = (phoEt->at(i) > ptThr);
            if(!passedPtThr ) continue;

            bool passedDR;
            bool passedMom; 
            bool passedGMom;
            bool passedPromptPho;      // selections for GEN photon
            bool passedGenIso;      // selections for GEN photon
            double deltaRMin= 999;
            int matchedIndex=-1;           // index of the matched GEN photon in  to this RECO photon

            for (int j=0; j<nMC; ++j){
                if(TMath::Abs(mcPID->at(j))!= PDG_PHOTON) continue;
                if( mcPt->at(j) <= 5 ) continue;
                //if( ! (mcPt->at(j)/phoEt->at(i) >=0.5 && mcPt->at(j)/phoEt->at(i)<=1.5) ) continue;
                double deltaRtmp = getDR(phoEta->at(i), phoPhi->at(i), mcEta->at(j), mcPhi->at(j));
                passedDR = (deltaRtmp < cutdeltaR);
                if (!passedDR) continue;
                //passedGenIso = (mcCalIsoDR04->at(j) <= 5.0);
                //passedGenIso = (mcCalIsoDR03->at(j) <= 5.0);
                //if (!passedGenIso) continue;
                if (deltaRtmp < deltaRMin){
                        deltaRMin = deltaRtmp;
                        matchedIndex = j;
                }
            } // mc loop
            if(matchedIndex==-1) continue;
/*
            passedMom = ( (mcMomPID->at(matchedIndex) == -999) || (abs(mcMomPID->at(matchedIndex))<=22) );
            passedGMom = ( (mcGMomPID->at(matchedIndex) == -999) || (abs(mcGMomPID->at(matchedIndex))<=22) );
            if(isPromptPho) passedPromptPho = passedMom && passedGMom;
            else passedPromptPho = !( passedMom && passedGMom );
*/
            passedPromptPho=1;
            if(passedPromptPho){
                mcPt_->push_back(mcPt->at(matchedIndex));
                mcPhi_->push_back(mcPhi->at(matchedIndex));
                mcEta_->push_back(mcEta->at(matchedIndex));
                mcPID_->push_back(mcPID->at(matchedIndex));
                mcDPhi_evtpl1_->push_back(getDPHI(mcPhi->at(matchedIndex),hiEvtPlanes[2]));
                mcDPhi_evtpl2_->push_back(getDPHI(mcPhi->at(matchedIndex),hiEvtPlanes[8]));
                mcDPhi_evtpl3_->push_back(getDPHI(mcPhi->at(matchedIndex),hiEvtPlanes[15]));
                mcDPhi_evtpl4_->push_back(getDPHI(mcPhi->at(matchedIndex),hiEvtPlanes[21]));
                mcMomPID_->push_back(mcMomPID->at(matchedIndex));
                mcGMomPID_->push_back(mcGMomPID->at(matchedIndex));
                mcStatus_->push_back(mcStatus->at(matchedIndex));
                mcCalIsoDR03_->push_back(mcCalIsoDR03->at(matchedIndex));
                mcCalIsoDR04_->push_back(mcCalIsoDR04->at(matchedIndex));
                phoEt_->push_back(phoEt->at(i));
                phoEta_->push_back(phoEta->at(i));
                phoPhi_->push_back(phoPhi->at(i));
                phoDPhi_evtpl1_->push_back(getDPHI(phoPhi->at(i),hiEvtPlanes[2]));
                phoDPhi_evtpl2_->push_back(getDPHI(phoPhi->at(i),hiEvtPlanes[8]));
                phoDPhi_evtpl3_->push_back(getDPHI(phoPhi->at(i),hiEvtPlanes[15]));
                phoDPhi_evtpl4_->push_back(getDPHI(phoPhi->at(i),hiEvtPlanes[21]));
                pho_ecalClusterIsoR2_->push_back(pho_ecalClusterIsoR2->at(i));
                pho_ecalClusterIsoR3_->push_back(pho_ecalClusterIsoR3->at(i));
                pho_ecalClusterIsoR4_->push_back(pho_ecalClusterIsoR4->at(i));
                pho_ecalClusterIsoR5_->push_back(pho_ecalClusterIsoR5->at(i));
                pho_hcalRechitIsoR2_->push_back(pho_hcalRechitIsoR2->at(i));
                pho_hcalRechitIsoR3_->push_back(pho_hcalRechitIsoR3->at(i));
                pho_hcalRechitIsoR4_->push_back(pho_hcalRechitIsoR4->at(i));
                pho_hcalRechitIsoR5_->push_back(pho_hcalRechitIsoR5->at(i));
                pho_trackIsoR2PtCut20_->push_back(pho_trackIsoR2PtCut20->at(i));
                pho_trackIsoR3PtCut20_->push_back(pho_trackIsoR3PtCut20->at(i));
                pho_trackIsoR4PtCut20_->push_back(pho_trackIsoR4PtCut20->at(i));
                pho_trackIsoR5PtCut20_->push_back(pho_trackIsoR5PtCut20->at(i));
                sumIsoR2_->push_back(pho_ecalClusterIsoR2->at(i)+pho_hcalRechitIsoR2->at(i)+pho_trackIsoR2PtCut20->at(i));
                sumIsoR3_->push_back(pho_ecalClusterIsoR3->at(i)+pho_hcalRechitIsoR3->at(i)+pho_trackIsoR3PtCut20->at(i));
                sumIsoR4_->push_back(pho_ecalClusterIsoR4->at(i)+pho_hcalRechitIsoR4->at(i)+pho_trackIsoR4PtCut20->at(i));
                sumIsoR5_->push_back(pho_ecalClusterIsoR5->at(i)+pho_hcalRechitIsoR5->at(i)+pho_trackIsoR5PtCut20->at(i));
                phoR9_->push_back(phoR9->at(i));
                phoHoverE_->push_back(phoHoverE->at(i));
                phoSigmaIEtaIEta_->push_back(phoSigmaIEtaIEta->at(i));
                phoSCE_->push_back(phoSCE->at(i));
                phoSCEta_->push_back(phoSCEta->at(i));
                phoSCPhi_->push_back(phoSCPhi->at(i));
                phoSCEtaWidth_->push_back(phoSCEtaWidth->at(i));
                phoSCPhiWidth_->push_back(phoSCPhiWidth->at(i));
                phoSCBrem_->push_back(phoSCBrem->at(i));
                pfcIso1_->push_back(pfcIso1->at(i));
                pfcIso2_->push_back(pfcIso2->at(i));
                pfcIso3_->push_back(pfcIso3->at(i));
                pfcIso4_->push_back(pfcIso4->at(i));
                pfcIso5_->push_back(pfcIso5->at(i));
                pfpIso1_->push_back(pfpIso1->at(i));
                pfpIso2_->push_back(pfpIso2->at(i));
                pfpIso3_->push_back(pfpIso3->at(i));
                pfpIso4_->push_back(pfpIso4->at(i));
                pfpIso5_->push_back(pfpIso5->at(i));
                pfnIso1_->push_back(pfnIso1->at(i));
                pfnIso2_->push_back(pfnIso2->at(i));
                pfnIso3_->push_back(pfnIso3->at(i));
                pfnIso4_->push_back(pfnIso4->at(i));
                pfnIso5_->push_back(pfnIso5->at(i));
                pfSumIso2_->push_back(pfpIso2->at(i)+pfcIso2->at(i)+pfnIso2->at(i));
                pfSumIso3_->push_back(pfpIso3->at(i)+pfcIso3->at(i)+pfnIso3->at(i));
                pfSumIso4_->push_back(pfpIso4->at(i)+pfcIso4->at(i)+pfnIso4->at(i));
                pfSumIso5_->push_back(pfpIso5->at(i)+pfcIso5->at(i)+pfnIso5->at(i));
                pfcVsIso1_->push_back(pfcVsIso1->at(i));
                pfcVsIso2_->push_back(pfcVsIso2->at(i));
                pfcVsIso3_->push_back(pfcVsIso3->at(i));
                pfcVsIso4_->push_back(pfcVsIso4->at(i));
                pfcVsIso5_->push_back(pfcVsIso5->at(i));
                pfpVsIso1_->push_back(pfpVsIso1->at(i));
                pfpVsIso2_->push_back(pfpVsIso2->at(i));
                pfpVsIso3_->push_back(pfpVsIso3->at(i));
                pfpVsIso4_->push_back(pfpVsIso4->at(i));
                pfpVsIso5_->push_back(pfpVsIso5->at(i));
                pfnVsIso1_->push_back(pfnVsIso1->at(i));
                pfnVsIso2_->push_back(pfnVsIso2->at(i));
                pfnVsIso3_->push_back(pfnVsIso3->at(i));
                pfnVsIso4_->push_back(pfnVsIso4->at(i));
                pfnVsIso5_->push_back(pfnVsIso5->at(i));
                pfSumVsIso2_->push_back(pfpVsIso2->at(i)+pfcVsIso2->at(i)+pfnVsIso2->at(i));
                pfSumVsIso3_->push_back(pfpVsIso3->at(i)+pfcVsIso3->at(i)+pfnVsIso3->at(i));
                pfSumVsIso4_->push_back(pfpVsIso4->at(i)+pfcVsIso4->at(i)+pfnVsIso4->at(i));
                pfSumVsIso5_->push_back(pfpVsIso5->at(i)+pfcVsIso5->at(i)+pfnVsIso5->at(i));

                nPho_++;
            }
        }// photon loop
        hiEvtPlane1_ = hiEvtPlanes[2];
        hiEvtPlane2_ = hiEvtPlanes[8];
        hiEvtPlane3_ = hiEvtPlanes[15];
        hiEvtPlane4_ = hiEvtPlanes[21];
        if(nPho_!=0) t_pho->Fill();
    } // exited event loop

    end_loop = std::clock();
    std::cout.precision(6);      // get back to default precision
    std::cout << "LOOP finished in             : " << (end_loop - start_loop) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "exited event loop" << std::endl;

    t_pho->Write();
    //=============================================================================================
    //=============================================================================================
    // Save the histograms 
    fout->Close();
    inputFile->Close();
}
