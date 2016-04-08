// To get efficiency.
// Apply denominator and numerator cuts as a string in an argument. 
// Using Event, mc, photon For 3 kinds of loops, Fill the *gen* photons one by one. 
// Data : 2016 Apr 07
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
#include "TChain.h"
#include "stdio.h"
#include <vector>
#include <map>
#include <string>
#include "../yjUtility.h"
#include "../phoRaaCuts/phoRaaCuts_v2.h"
#include "../ElectroWeak-Jet-Track-Analyses/TreeHeaders/ggHiNtuplizerTree.h"
#include "../ElectroWeak-Jet-Track-Analyses/Utilities/interface/CutConfigurationParser.h"
#include "../ElectroWeak-Jet-Track-Analyses/Utilities/interface/InputConfigurationParser.h"
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
*/
void stringToMap(string cutst, map<string,int> & cutmap, int mapVal){
    std::string s = cutst;
    std::string delimiter = "_";
    //std::string delimiter = ",";

    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        cutmap[token] = mapVal;
        //std::cout << token << std::endl;
        s.erase(0, pos + delimiter.length());
    }
    //std::cout << s << std::endl;
    cutmap[s] = mapVal;
}

//Cut variation = fragDirPhoCut, fragPhoCut, dirPhoCut, decayPhoCut, hotspotCut, spikeCut, hoeCut(pp/pbpb), sigEtaEtaCut(pp/pbpb), sumIsoCut(pp/pbpb)
void phoEfficiency_v2(const TString inputFile, string den_st="", string num_st="hotspotCut_spikeCut_hoeCut_sigEtaEtaCut_sumIsoCut", const TString coll="pbpb", int cent_i=0, int cent_f=200, const char* ver="v2")
{
    SetHistTitleStyle();
    TH1::SetDefaultSumw2();
    double pur_num_pp[] = {91.0,89.0,88.0,84.0};
    gStyle->SetLabelSize(0.03,"Y");
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleXSize(0.05);
    TChain *treePho;
    if(coll=="pbpb") treePho = new TChain("ggHiNtuplizer/EventTree");
    else if(coll=="pp") treePho = new TChain("ggHiNtuplizerGED/EventTree");
    TChain *treeEvt = new TChain("hiEvtAnalyzer/HiTree");

    std::vector<std::string> inputFiles = InputConfigurationParser::ParseFiles(inputFile.Data());
    std::cout<<"input ROOT files : num = "<<inputFiles.size()<< std::endl;
    for(std::vector<std::string>::iterator it = inputFiles.begin() ; it != inputFiles.end(); ++it) {
        treePho->Add((*it).c_str());
        treeEvt->Add((*it).c_str());
        cout <<"input ROOT file : " << (*it).c_str() << endl;
    }
    TString outputFile;
    if(coll=="pbpb") outputFile = Form("hist_efficiency_pbpb_cent%d-%d_%s_Den_%s_Num_%s.root",(int)cent_i/2,(int)cent_f/2,ver,den_st.data(),num_st.data());
    else if(coll=="pp") outputFile = Form("hist_efficiency_pp_%s_Den_%s_Num_%s.root",ver,den_st.data(),num_st.data());
    cout << "output ROOT file : " << outputFile << endl;

    ggHiNtuplizer pho;
    pho.setupTreeForReading(treePho);
    ULong64_t event=0;
    unsigned run=0;
    unsigned lumi=0;
    int hiBin=0;
    treeEvt->SetBranchAddress("evt",&event);
    treeEvt->SetBranchAddress("run",&run);
    treeEvt->SetBranchAddress("lumi",&lumi);
    treeEvt->SetBranchAddress("hiBin",&hiBin);

    // if map=0,no cut on den&num// map=1, cut on both den and num// map=2, cut on only num but not on den//
    map<string,int> cut_map; 
    stringToMap(den_st,cut_map,1);
    stringToMap(num_st,cut_map,2);
    cout << "############ MAP KEYS ############## " << endl;
    //cout << "0 for non of them, 1 for both of them, 2 for only numberator" << endl;
    for(map<string,int>::iterator it = cut_map.begin() ; it != cut_map.end() ; ++it)
    {
        cout << it->first <<" "  << it->second << endl;
    }
    cout << "#################################### " << endl;

    //reconstructed photons in generated photons (reco pt & reco eta)
    TH2D* h2D_Num = new TH2D("h2D_Num", ";|#eta|;p_{T} (GeV)", nEtaBin, etaBins, nPtBin, ptBins);
    TH2D* h2D_Den = (TH2D*) h2D_Num->Clone("h2D_Den");
    TH2D* h2D_Eff = (TH2D*) h2D_Num->Clone("h2D_Eff");

    TH2D* h2D_Num_recons = new TH2D("h2D_Num_recons", ";|#eta|;p_{T} (GeV)", nEtaBin, etaBins, nPtBin, ptBins);
    TH2D* h2D_Den_recons = (TH2D*) h2D_Num->Clone("h2D_Den_recons");
    TH2D* h2D_Eff_recons = (TH2D*) h2D_Num->Clone("h2D_Eff_recons");
    
    TH2D* h2D_Eff_total;
    if(den_st=="") h2D_Eff_total = (TH2D*) h2D_Num->Clone("h2D_Eff_total");

    Long64_t nentries = treePho->GetEntries();
    for(int ientries = 0; ientries < nentries; ++ientries)
    {
        if (ientries % 10000 == 0)  {
            std::cout << "current entry = " <<ientries<<" out of "<<nentries<<" : "<<std::setprecision(2)<<(double)ientries/nentries*100<<" %"<<std::endl;
        }
        treePho->GetEntry(ientries);
        treeEvt->GetEntry(ientries);

        if(coll=="pbpb" && !(hiBin>=cent_i && hiBin<cent_f)) continue;
        for (int igen = 0; igen < pho.nMC; ++igen){
            if (pho.mcStatus->at(igen) != 1 || (pho.mcPID->at(igen)) != 22 ) continue;
            if (pho.mcPt->at(igen)<30) continue;
            if (pho.mcCalIsoDR04->at(igen)>5) continue;

            //////////////// FLAGS FOR DENOMERATOR and NUMERATOR ////////////
            bool denFlag = 1;
            bool numFlag = 1;
 
            //cout << "passed gen continues " << endl;
            if (cut_map.find("dirPhoCut") != cut_map.end()) {
                bool isPass = (pho.mcMomPID->at(igen) == 22);
                //cout << "momID = " << pho.mcMomPID->at(igen) <<", dirPhoCut isPass = " << isPass << endl; 
                if( (cut_map.at("dirPhoCut") ==1) && (!isPass) ) continue; 
                else if((cut_map.at("dirPhoCut") ==2) && (!isPass) ) numFlag=0; 
            }

            if (cut_map.find("fragPhoCut") != cut_map.end()) {
                bool isPass = (pho.mcMomPID->at(igen) == 22);
                //cout << "momID = " << pho.mcMomPID->at(igen) <<", fragPhoCut isPass = " << isPass << endl; 
                if( (cut_map.at("fragPhoCut") ==1) && (!isPass) ) continue; 
                else if((cut_map.at("fragPhoCut") ==2) && (!isPass) ) numFlag=0; 
            }

            if (cut_map.find("fragDirPhoCut") != cut_map.end()) {
                bool isPass = (pho.mcMomPID->at(igen) == 22);
                //cout << "momID = " << pho.mcMomPID->at(igen) <<", fragDirPhoCut isPass = " << isPass << endl; 
                if( (cut_map.at("fragDirPhoCut") ==1) && (!isPass) ) continue; 
                else if((cut_map.at("fragDirPhoCut") ==2) && (!isPass) ) numFlag=0; 
            }

            if (cut_map.find("decayPhoCut") != cut_map.end()) {
                bool isPass = (pho.mcMomPID->at(igen) == 22);
                //cout << "momID = " << pho.mcMomPID->at(igen) <<", decayPhoCut isPass = " << isPass << endl; 
                if( (cut_map.at("decayPhoCut") ==1) && (!isPass) ) continue; 
                else if((cut_map.at("decayPhoCut") ==2) && (!isPass) ) numFlag=0; 
            }

            ///// MATCHING PART /////
            Float_t currentDiffPt(1000);
            Float_t currentRecoPt(-1);
            int matchedIndex = -1;
            for(int ipho = 0; ipho < pho.nPho; ++ipho){
                //cout << "pho.phoEt = " << pho.phoEt->at(ipho) << endl;
                float tempDR = getDR(pho.phoEta->at(ipho), pho.phoPhi->at(ipho), pho.mcEta->at(igen), pho.mcPhi->at(igen));
                float tempDiffPt = abs( pho.mcPt->at(igen) - pho.phoEt->at(ipho) );
                if( (tempDR < delta) && (tempDiffPt < currentDiffPt) ) {
                    currentDiffPt = tempDiffPt;
                    currentRecoPt = pho.phoEt->at(ipho);
                    matchedIndex = ipho;
                }
            }
           
            ///////////////// RECONSTRUCTION EFFICIENCY ///////////////
            h2D_Den_recons->Fill( abs(pho.mcEta->at(igen)), pho.mcEt->at(igen) );
            if(matchedIndex==-1) continue;
            h2D_Num_recons->Fill( abs(pho.mcEta->at(igen)), pho.mcEt->at(igen) );
            ///////////////////////////////////////////////////////////            
#if 1 
            if(TMath::Abs(pho.phoEta->at(matchedIndex))>3) continue; 
            if(pho.phoEt->at(matchedIndex)<25) continue; 

            ///////////////// START CONDITONS //////////////    
            if (cut_map.find("spikeCut") != cut_map.end()) {
                bool isFail = ( (pho.phoEta->at(matchedIndex)<1.44) &&
                        (pho.phoSigmaIEtaIEta->at(matchedIndex) < 0.002 ||
                         pho.pho_swissCrx->at(matchedIndex)     > 0.9   ||
                         TMath::Abs(pho.pho_seedTime->at(matchedIndex)) > 3) );
                //cout << "is spike(fail) ? = " << isFail << endl; 
                if( (cut_map.at("spikeCut") ==1) && (isFail) ) continue; 
                else if((cut_map.at("spikeCut") ==2) && (isFail) ) numFlag=0; 
            }
   
            if (cut_map.find("hotspotCut") != cut_map.end()) {
                bool isFail = 
                   ( (pho.phoE3x3->at(matchedIndex)/pho.phoE5x5->at(matchedIndex) > 2./3.-0.03 && pho.phoE3x3->at(matchedIndex)/pho.phoE5x5->at(matchedIndex) < 2./3.+0.03) &&
                    (pho.phoE1x5->at(matchedIndex)/pho.phoE5x5->at(matchedIndex) > 1./3.-0.03 && pho.phoE1x5->at(matchedIndex)/pho.phoE5x5->at(matchedIndex) < 1./3.+0.03) &&
                    (pho.phoE2x5->at(matchedIndex)/pho.phoE5x5->at(matchedIndex) > 2./3.-0.03 && pho.phoE2x5->at(matchedIndex)/pho.phoE5x5->at(matchedIndex) < 2./3.+0.03) ); 
                if( (cut_map.at("hotspotCut") ==1) && (isFail) ) continue; 
                else if((cut_map.at("hotspotCut") ==2) && (isFail) ) numFlag=0; 
            }

            if (cut_map.find("hoeCut") != cut_map.end()) {
                if(coll=="pbpb"){
                    bool isPass = (pho.phoHoverE->at(matchedIndex)<0.1);
                    if( (cut_map.at("hoeCut") ==1) && (!isPass) ) continue; 
                    else if((cut_map.at("hoeCut") ==2) && (!isPass) ) numFlag=0; 
                } else if(coll=="pp"){
                    bool isPass = (pho.phoHoverE->at(matchedIndex)<0.05);
                    if( (cut_map.at("hoeCut") ==1) && (!isPass) ) continue; 
                    else if((cut_map.at("hoeCut") ==2) && (!isPass) ) numFlag=0; 
                }
            }

            if (cut_map.find("sigEtaEtaCut") != cut_map.end()) {
                if(coll=="pbpb"){
                    bool isPass = (pho.phoSigmaIEtaIEta_2012->at(matchedIndex)<0.0101);
                    if( (cut_map.at("sigEtaEtaCut") ==1) && (!isPass) ) continue; 
                    else if((cut_map.at("sigEtaEtaCut") ==2) && (!isPass) ) numFlag=0; 
                } else if(coll=="pp"){
                    bool isPass =( 
                        (abs(pho.phoEta->at(matchedIndex))<=1.44 && pho.phoSigmaIEtaIEta_2012->at(matchedIndex)<0.0102) ||
                        (abs(pho.phoEta->at(matchedIndex))>1.44 && pho.phoSigmaIEtaIEta_2012->at(matchedIndex)<0.0268)
                    );
                    if( (cut_map.at("sigEtaEtaCut") ==1) && (!isPass) ) continue; 
                    else if((cut_map.at("sigEtaEtaCut") ==2) && (!isPass) ) numFlag=0; 
                }
            }
            
            if (cut_map.find("sumIsoCut") != cut_map.end()) {
                if(coll=="pbpb"){
                    float sumIso = pho.pho_ecalClusterIsoR4->at(matchedIndex) + pho.pho_hcalRechitIsoR4->at(matchedIndex) + pho.pho_trackIsoR4PtCut20->at(matchedIndex);
                    bool isPass = (sumIso<1);
                    if( (cut_map.at("sumIsoCut") ==1) && (!isPass) ) continue; 
                    else if((cut_map.at("sumIsoCut") ==2) && (!isPass) ) numFlag=0; 
                } else if(coll=="pp"){
                    bool barrelIso = ((pho.pfcIso4->at(matchedIndex)<=1.37) && (pho.pfnIso4->at(matchedIndex)<=1.06+0.014*pho.phoEt->at(matchedIndex)+0.000019*pho.phoEt->at(matchedIndex)*pho.phoEt->at(matchedIndex)) && pho.pfpIso4->at(matchedIndex)<=(0.28+0.0053*pho.phoEt->at(matchedIndex)));                    
                    bool endcapIso = ((pho.pfcIso4->at(matchedIndex)<=1.10) && (pho.pfnIso4->at(matchedIndex)<=2.69+0.0139*pho.phoEt->at(matchedIndex)+0.000025*pho.phoEt->at(matchedIndex)*pho.phoEt->at(matchedIndex)) && pho.pfpIso4->at(matchedIndex)<=(0.39+0.0034*pho.phoEt->at(matchedIndex)));  
                    bool isPass = ( ( (abs(pho.phoEta->at(matchedIndex))<=1.44) && barrelIso ) || 
                            ( (abs(pho.phoEta->at(matchedIndex))>1.44) && endcapIso ) );
                    if( (cut_map.at("sumIsoCut") ==1) && (!isPass) ) continue; 
                    else if((cut_map.at("sumIsoCut") ==2) && (!isPass) ) numFlag=0; 

                    // default photon isolation for PP is "Loose".
//                    cut_phoHOverE_EB = 0.05;            // Barrel
//                    cut_phoSigmaIEtaIEta_EB = 0.0102;
//                    cut_pfcIso4_EB = 3.32;
//                    cut_pfnIso4_c0_EB = 1.92;
//                    cut_pfnIso4_c1_EB = 0.014;
//                    cut_pfnIso4_c2_EB = 0.000019;
//                    cut_pfpIso4_c0_EB = 0.81;
//                    cut_pfpIso4_c1_EB = 0.0053;
//                    cut_phoHOverE_EE = 0.05;            // Endcap
//                    cut_phoSigmaIEtaIEta_EE = 0.0274;
//                    cut_pfcIso4_EE = 1.97;
//                    cut_pfnIso4_c0_EE = 11.86;
//                    cut_pfnIso4_c1_EE = 0.0139;
//                    cut_pfnIso4_c2_EE = 0.000025;
//                    cut_pfpIso4_c0_EE = 0.83;
//                    cut_pfpIso4_c1_EE = 0.0034;
                }
            }
            ///////////////// END OF CONDITON //////////////    
            //cout << "denFlag = " << denFlag << ", numFlag = " << numFlag << endl;

            h2D_Den->Fill( abs(pho.mcEta->at(igen)), pho.mcEt->at(igen) );
            if(numFlag) h2D_Num->Fill( abs(pho.mcEta->at(igen)), pho.mcEt->at(igen) );
#endif
        }//gen loop
    }//event loop

    ///////////////// GET EFFICIENCY /////////////////    

    h2D_Eff->Divide(h2D_Num,h2D_Den,1.,1.,"B");
    h2D_Eff_recons->Divide(h2D_Num_recons,h2D_Den_recons,1.,1.,"B");

    TH1D* h1D_pt_recons= (TH1D*) h2D_Eff_recons->ProjectionY("h1D_Eff_pt_recons",0,1);
    TH1D* h1D_pt = (TH1D*) h2D_Eff->ProjectionY("h1D_Eff_pt",0,1);
 
    TH1D* h1D_pt_recons_pretty = changeToPrettyBin(h1D_pt_recons);
    TH1D* h1D_pt_pretty = changeToPrettyBin(h1D_pt);
    h1D_pt_recons_pretty->SetTitle(";p_{T} (GeV);Reconstruction Efficiency");
    h1D_pt_recons_pretty->SetMarkerStyle(20);
    h1D_pt_recons_pretty->SetAxisRange(0.50,1.0001,"Y");
    h1D_pt_pretty->SetTitle(";p_{T} (GeV);RecoCut Efficiency");
    h1D_pt_pretty->SetMarkerStyle(20);
    h1D_pt_pretty->SetAxisRange(0.50,1.0001,"Y");

    TH1D* h1D_pt_total;
    TH1D* h1D_pt_total_pretty;
    if(den_st=="") { 
        h2D_Eff_total->Divide(h2D_Num,h2D_Den_recons,1.,1.,"B");
        h1D_pt_total = (TH1D*) h2D_Eff_total->ProjectionY("h1D_Eff_pt_total",0,1);
        h1D_pt_total_pretty = changeToPrettyBin(h1D_pt_total);

        h1D_pt_total_pretty->SetTitle(";p_{T} (GeV);Total Efficiency");
        h1D_pt_total_pretty->SetMarkerStyle(20);
        h1D_pt_total_pretty->SetAxisRange(0.50,1.0001,"Y");
    } 
    
//    drawText(cap,0.2,0.2);
//    TEfficiency* hEff = new TEfficiency(h1D_Num_pt,h1D_Den_pt);
//    TCanvas* c = new TCanvas("c1", "", 700,700);
//    TGraphAsymmErrors* g_Eff_pt = new TGraphAsymmErrors(h1D_Num_pt, h1D_Den_pt); 
//    g_Eff_pt->SetMarkerStyle(20);
//    g_Eff_pt->Draw("p");
//    //drawText(cap,0.2,0.2);
//    jumSun(ptBins[0],1,ptBins[nPtBin],1);
//    c->SaveAs(Form("figures/g_Eff_pt_reco_%s_nohotspot.pdf",cap));


    saveHistogramsToPicture(h1D_pt_recons_pretty, "pdf",Form("efficiency_Den_%s_Num_%s_%s_cent%d_%d",den_st.data(),num_st.data(),coll.Data(),(int)cent_i/2,(int)cent_f/2));
    saveHistogramsToPicture(h1D_pt_pretty, "pdf",Form("efficiency_Den_%s_Num_%s_%s_cent%d_%d",den_st.data(),num_st.data(),coll.Data(),(int)cent_i/2,(int)cent_f/2));
    saveHistogramsToPicture(h2D_Eff, "pdf",Form("efficiency_pt_eta_Den_%s_Num_%s_%s_cent%d_%d",den_st.data(),num_st.data(),coll.Data(),(int)cent_i/2,(int)cent_f/2),"figures",0);
    saveHistogramsToPicture(h2D_Eff_recons, "pdf",Form("reconstruction_efficiency_pt_eta_Den_%s_Num_%s_%s_cent%d_%d",den_st.data(),num_st.data(),coll.Data(),(int)cent_i/2,(int)cent_f/2),"figures",0);

    if(den_st=="") saveHistogramsToPicture(h2D_Eff_total, "pdf",Form("total_efficiency_pt_eta_Den_%s_Num_%s_%s_cent%d_%d",den_st.data(),num_st.data(),coll.Data(),(int)cent_i/2,(int)cent_f/2),"figures",0);
    if(den_st=="") saveHistogramsToPicture(h1D_pt_total_pretty, "pdf",Form("efficiency_Den_%s_Num_%s_%s_cent%d_%d",den_st.data(),num_st.data(),coll.Data(),(int)cent_i/2,(int)cent_f/2));

    TFile* outf = new TFile(Form("output/%s",outputFile.Data()), "RECREATE");
    outf->cd();
    h2D_Num->Write();
    h2D_Den->Write();
    h2D_Eff->Write();
    h2D_Num_recons->Write();
    h2D_Den_recons->Write();
    h2D_Eff_recons->Write();
    h1D_pt_recons->Write();
    h1D_pt->Write();
    if(den_st==""){
        h2D_Eff_total->Write();
        h1D_pt_total->Write();
    }
    outf->Close();
}

