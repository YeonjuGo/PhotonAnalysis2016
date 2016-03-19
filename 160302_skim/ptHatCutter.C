///////////////////////////////////////////////////////////////////                                
// forest2yskim.C                                                //                                                 
// Creator : Yongsun Kim (MIT), jazzitup@mit.edu                 //                                                 
// Function : Transform hiForest files into yskim file           //
// yskims for MinBias1, Minbias2 and photon jet skims            //
///////////////////////////////////////////////////////////////////         

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <TMath.h>
#include "../HiForestAnalysis/hiForest.h"
#include "../yjUtility.h"
#include <time.h>
using namespace std;
static const long MAXTREESIZE = 10000000000;

void ptHatCutter(TString inputFileName, TString outFileName, int ptHatCut1 = 30, int ptHatCut2=50) {

    // start from here
    // path length histogram

    TFile* f = new TFile(inputFileName.Data());
    TTree* hlt = (TTree*) f ->Get("hltanalysis/HltTree");  
    TTree* pho = (TTree*) f ->Get("ggHiNtuplizer/EventTree");  
    TTree* phoGED = (TTree*) f ->Get("ggHiNtuplizerGED/EventTree");  
    TTree* skim = (TTree*) f ->Get("skimanalysis/HltTree");  
    TTree* evt = (TTree*) f ->Get("hiEvtAnalyzer/HiTree");  
    TTree* gen = (TTree*) f ->Get("HiGenParticleAna/hi");  
    TTree* hiforest = (TTree*) f ->Get("HiForest/HiForestInfo");  


    TFile* fout = new TFile(outFileName.Data(),"recreate");
    TTree *nhlt = hlt->CloneTree(0);
    TTree *npho = pho->CloneTree(0);
    TTree *nphoGED = phoGED->CloneTree(0);
    TTree *nskim = skim->CloneTree(0);
    TTree *nevt = evt->CloneTree(0);
    TTree *ngen = gen->CloneTree(0);
    TTree *nhiforest = hiforest->CloneTree(0);
    nphoGED->SetName("EventTreeGED");

    float pthat, pthatWeight;
    nevt->Branch("pthatWeight", &pthatWeight,"pthatWeight/F");
    evt->SetBranchAddress("pthat",&pthat);

    float tmpWeight;
    if(ptHatCut1==15 && ptHatCut2==30) tmpWeight = 74729./74752.; 
    else if(ptHatCut1==30 && ptHatCut2==50) tmpWeight = 63479./180164.; 
    else if(ptHatCut1==50 && ptHatCut2==80) tmpWeight = 32709./234691.; 
    else if(ptHatCut1==80 && ptHatCut2==120) tmpWeight = 10194./222779.; 
    else if(ptHatCut1==120 && ptHatCut2==9999) tmpWeight = 3432./266700.;
    
    //int nentries = 100;
    int nentries = pho->GetEntries();
    cout << "number of entries = " << nentries << endl;
    for (Long64_t jentry = 0 ; jentry < nentries; jentry++) {
        if (jentry% 10000 == 0)  {
            cout <<jentry<<" / "<<nentries<<" "<<setprecision(2)<<(double)jentry/nentries*100<<endl;
        }
        hlt->GetEntry(jentry);
        pho->GetEntry(jentry);
        phoGED->GetEntry(jentry);
        skim->GetEntry(jentry);
        evt->GetEntry(jentry);
        gen->GetEntry(jentry);
        hiforest->GetEntry(jentry);

        //cout << "pthat : " << pthat << endl;
        if ( !((pthat > ptHatCut1) && (pthat < ptHatCut2)) )
            continue;

#if 0    
        bool isHighPtGen = false;
        for (int jpho=0;jpho< t->genp.nPar ;jpho++) {
            if (  ( t->genp.et[jpho] > genPhotonPtCut ) && ( fabs(t->genp.eta[jpho]) < 3 ) && ( abs(t->genp.momId[jpho])<=22 ) && ( t->genp.id[jpho] == 22 ) && ( t->genp.status[jpho]==1 ) ) {
                isHighPtGen = true;
                jpho = t->genp.nPar;
            }
        }

        if (isHighPtGen)
#endif
        //cout << "pthat : " << pthat << endl;
        pthatWeight = tmpWeight;
        nhlt->Fill();
        npho->Fill();
        nphoGED->Fill();
        nskim->Fill();
        nevt->Fill();
        ngen->Fill();
        nhiforest->Fill();
    }
    nhlt->Write();
    npho->Write();
    nphoGED->Write();
    nskim->Write();
    nevt->Write();
    ngen->Write();
    nhiforest->Write();
    //delete t;
}

int main(int argc, char *argv[])
{
    if(argc == 5)
    {
        ptHatCutter(argv[1],
                argv[2],
                atoi(argv[3]),
                atoi(argv[4]));
        return 0;
    }
    else if (argc == 4)
    {
        ptHatCutter(argv[1],
                argv[2],
                atoi(argv[3]));
        return 0;
    }
    else {
        std::cout << "Wrong number of arguments" << std::endl;
        return 1;
    }
}
