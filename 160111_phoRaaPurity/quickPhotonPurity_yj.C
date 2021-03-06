// photonPurity.C
// Author: Alex Barbieri
// Cobbled together from code written by Yongsun Kim

// Calculates and plots the photon purity using a sideband estimation
// of the background and a MC calculation of the signal. This info is
// used to do a template fit in the sigmaIetaIeta variable.

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1D.h"
#include "TCut.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TLatex.h"
#include "stdio.h"
#include <iostream>
#include "../ElectroWeak-Jet-Track-Analyses/Histogramming/PhotonPurity.h"
#include "../ElectroWeak-Jet-Track-Analyses/Utilities/interface/CutConfigurationParser.h"
#include "../ElectroWeak-Jet-Track-Analyses/TreeHeaders/CutConfigurationTree.h"
#include "../ElectroWeak-Jet-Track-Analyses/Plotting/commonUtility.h"
#include "../phoRaaCuts_v1.h"
/*
///// KNU
const TString ppDatafname = "/d3/scratch/goyeonju/files/photons2016/2015-Data-promptRECO-photonSkims_pp-photonHLTFilter-v0-HiForest.root";
const TString ppMCfname = "/d3/scratch/goyeonju/files/photons2016/2015-PP-MC_Pythia8_Photon30_pp502_TuneCUETP8M1.root";
const TString pbpbMCfname = "/d3/scratch/goyeonju/files/photons2016/2015-PbPb-MC_AllQCDPhoton30-v0-HiForest_PYTHIA_HYDJET_160129.root";
const TString pbpbDatafname = "/d3/scratch/goyeonju/files/photons2016/forestSkimed_photonSkim_pbpb_2015data.root";

///// MIT
//const TString ppMCfname = "/mnt/hadoop/cms/store/user/luck/2014-photon-forests/partial_PbPb_gammaJet_MC/HiForest_QCDPhoton30.root";
//const TString ppDatafname = "/mnt/hadoop/cms/store/user/luck/2015-Data-photonSkims/HighPtPhoton30AndZ/pp-photonSkim-wFilter/0.root";
//const TString pbpbMCfname = "/mnt/hadoop/cms/store/user/luck/Pyquen_AllQCDPhoton30_PhotonFilter20GeV_eta24-HiForest.root";
//const TString pbpbDatafname = "/mnt/hadoop/cms/store/user/luck/2015-Data-photonSkims/HIPhoton40AndZ/PbPb-photonSkim-v1/0.root";
//const TString thisConfig = "../ElectroWeak-Jet-Track-Analyses/CutConfigurations/gamma-jet-nominal.conf";
*/
const TString LABEL = "PbPb #sqrt{s}_{_{NN}}=5.02 TeV";
const TCut sampleIsolation_pbpb = "(pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20) < 1.0 && phoHoverE<0.1";
const TCut sampleIsolation_pp = "((pfcIso4<=1.37) && (pfnIso4<=1.06+0.014*phoEt+0.000019*phoEt*phoEt) && pfpIso4<=(0.28+0.0053*phoEt)) && phoHoverE<0.1";


//const Double_t sigShifts[] = {-0.0000989, -0.000131273, -0.00016207, -0.000170555};
const Double_t sigShifts[] = {0, 0, 0, 0};
//const Double_t sigShifts[] = {-0.00015,-0.00015,-0.00015,-0.00015};

// last entry is upper bound on last bin
const Int_t CENTBINS[] = {0,60, 200};//, 40};
//const Int_t CENTBINS[] = {0, 200};//, 40};
//const Int_t CENTBINS[] = {0, 100};
const Int_t nCENTBINS = sizeof(CENTBINS)/sizeof(Int_t) -1;

//const Double_t PTBINS[] = {40, 50, 60, 80, 1000};
//const Double_t PTBINS[] = {40,50,60, 1000};
//const Double_t PTBINS[] = {60,70,80,90,100,220};
const Double_t PTBINS[] = {40,50,60,70,80};
const Int_t nPTBINS = sizeof(PTBINS)/sizeof(Double_t) -1;

const Double_t ETABINS[] = {-1.44, 1.44};
//const Double_t ETABINS[] = {-1.44, -1, -0.5, 0, 0.5, 1, 1.44};
const Int_t nETABINS = sizeof(ETABINS)/sizeof(Double_t) -1;

void quickPhotonPurity(const TString configFile, const TString inputData, const TString inputMC, const TString outputName, const TString coll="pbpb")
{
  TH1::SetDefaultSumw2();
  CutConfiguration config = CutConfigurationParser::Parse(configFile.Data());
  TTree *configTree = setupConfigurationTreeForWriting(config);
   const char* photreeSt="";
   const char* hitreeSt="";
   const char* skimSt="";
  if(coll=="pbpb") {
      photreeSt="EventTree";
      hitreeSt="HiTree";
      skimSt="HltTree";
  } else {
      photreeSt="ggHiNtuplizerGED/EventTree";
      hitreeSt="hiEvtAnalyzer/HiTree";
      skimSt="skimanalysis/HltTree";
  } 
 cout << "s1" << endl;
  TFile *dataFile = TFile::Open(inputData);
  TTree *dataTree = (TTree*)dataFile->Get(photreeSt);
  TTree *dataEvtTree = (TTree*)dataFile->Get(hitreeSt);
  TTree *dataSkimTree= (TTree*)dataFile->Get(skimSt);
//  TTree *dataTree = (TTree*)dataFile->Get("ggHiNtuplizer/EventTree");
//  TTree *dataEvtTree = (TTree*)dataFile->Get("hiEvtAnalyzer/HiTree");
  //TTree *dataTree = (TTree*)dataFile->Get("photonSkimTree");

  dataTree->AddFriend(dataEvtTree);
  dataTree->AddFriend(dataSkimTree);
  TFile *mcFile = TFile::Open(inputMC);
  TTree *mcTree = (TTree*)mcFile->Get(photreeSt);
  TTree *mcEvtTree = (TTree*)mcFile->Get(hitreeSt);
  //TTree *mcTree = (TTree*)mcFile->Get("ggHiNtuplizer/EventTree");
  //TTree *mcEvtTree = (TTree*)mcFile->Get("hiEvtAnalyzer/HiTree");
  //TTree *mcTree = (TTree*)mcFile->Get("photonSkimTree");
    mcTree->AddFriend(mcEvtTree);
 cout << "skk" << endl;
  //TTree *mcTree = (TTree*)mcFile->Get("ggHiNtuplizer/EventTree");
  TFile *outFile = new TFile(outputName,"RECREATE");

  const TCut sidebandIsolation = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)>10) && ((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<20) && phoHoverE<0.1";
  //const TCut sidebandIsolation = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)>5) && ((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<10) && phoHoverE<0.1";
  const TCut mcIsolation = "(pho_genMatchedIndex!= -1) && mcCalIsoDR04[pho_genMatchedIndex]<5 && abs(mcPID[pho_genMatchedIndex])<=22";
  const TCut spikeRejection = "(phoSigmaIEtaIEta_2012>=0.002) && (pho_swissCrx<=0.9) && (abs(pho_seedTime)<=3)";
  TCut evtSel="";
 if(coll=="pbpb") evtSel = "pcollisionEventSelection";
 else if(coll=="pp") evtSel = "(pBeamScrapingFilter && pPAprimaryVertexFilter && pVertexFilterCutEandG)";
    //const TCut spikeRejection = "!((phoSigmaIEtaIEta_2012<0.002) || (pho_swissCrx>0.9) || (abs(pho_seedTime)>3))";

  cout << "JJ" << endl;
  //TCanvas *cPurity[nPTBINS];
  //TCanvas *cPurity = new TCanvas("c1","c1",337*nPTBINS,300*nCENTBINS/**2*/);
  TCanvas *cPurity = new TCanvas("c1","c1",400*nPTBINS,400*nCENTBINS);
  //cPurity->Divide(nPTBINS,2*nCENTBINS,0,0);
  //cPurity->Divide(nPTBINS,nCENTBINS,0,0);
  makeMultiPanelCanvas(cPurity, nPTBINS, nCENTBINS, 0.0, 0.0 , 0.2, 0.15, 0.005);
  cout << "nPTBINS = " << nPTBINS << ", nCENTBINS = " << nCENTBINS << ", nETABINS = " << nETABINS << endl;
  ////////////
  TH1D* h1D_p[nCentBin];
  for(Int_t j = 0; j < nCENTBINS; ++j) {
      h1D_p[j]= new TH1D(Form("h1D_purity_cent%d",j), ";p_{T} (GeV);Purity", nPtBin, ptBins);
  }
  int temp_nCentBin = nCENTBINS;
  if(coll=="pp") temp_nCentBin = 1;
  /////////////////
  for(Int_t i = 0; i < nPTBINS; ++i) {
    cout << "i : " << i << endl;
    //cPurity[i] = new TCanvggas(Form("c1_%d",i),"",1920,1000);
    //cPurity[i]->Divide(nETABINS,2,0,0);
    for(Int_t j = 0; j < temp_nCentBin; ++j) {
    cout << "j : " << j << endl;
      for(Int_t k = 0; k< nETABINS; ++k) {
    cout << "k : " << k << endl;
	TString ptCut = Form("(phoEt >= %f) && (phoEt < %f)",
			     PTBINS[i], PTBINS[i+1]);
	TString centCut = Form("((hiBin) >= %i) && ((hiBin) < %i)",
			     CENTBINS[j], CENTBINS[j+1]);
	TString etaCut = Form("(phoEta >= %f) && (phoEta < %f)",
			      ETABINS[k], ETABINS[k+1]);

    if(coll=="pp") centCut = "(1)"; 
    
	//TString pPbflipetaCut = Form("(eta*((run>211257)*-1+(run<211257)) >=%f) && (eta*((run>211257)*-1+(run<211257)) <%f)",
	//			     ETABINS[k], ETABINS[k+1]);
    TCut sampleIsolation = "";
    if(coll=="pbpb") sampleIsolation = sampleIsolation_pbpb; 
    else if(coll=="pp") sampleIsolation = sampleIsolation_pp; 
	TCut dataCandidateEvtSelCut = sampleIsolation && etaCut && ptCut && centCut && spikeRejection && evtSel;
	TCut dataCandidateCut = sampleIsolation && etaCut && ptCut && centCut && spikeRejection;
	TCut sidebandCut =  sidebandIsolation && etaCut && ptCut && centCut && spikeRejection && evtSel;
	TCut mcSignalCut = dataCandidateCut && mcIsolation;

	// if(nETABINS != 1)
	// {
	//   dataCandidateCut = sampleIsolation && pPbflipetaCut && ptCut && centCut;
	//   sidebandCut =  sidebandIsolation && pPbflipetaCut && ptCut && centCut;
	//   mcSignalCut =  sampleIsolation && etaCut && ptCut && centCut && mcIsolation;
	// }

	PhotonPurity fitr = getPurity(config, dataTree, mcTree,
				      dataCandidateEvtSelCut, sidebandCut,
				      mcSignalCut);

	//cPurity[i*nCENTBINS+j] = new TCanvas(Form("cpurity%d",i*nCENTBINS+j),
	// 					 "",500,500);
    cout << "centBin = " << centCut << ", ptBin : " << ptCut << ",,,, canvas # : "<<2*(k+j)*nPTBINS+i+1 << endl;
    cout << "k = " << k << ", j = " << j << ", i = " << i << endl;
    //cPurity->cd(2*(k+j)*nPTBINS+i+1);
	cPurity->cd((k+j)*nPTBINS+i+1);
	//cPurity[i]->cd(k+1);

	TH1F *hSigPdf = fitr.sigPdf;
	TH1F *hBckPdf = fitr.bckPdf;
	TH1D *hData1  = fitr.data;
	hSigPdf->Add(hBckPdf);

	TString name = "mcfit_total_ptbin";
	name += i;

	// outFile->cd();
	// hSigPdf->SetName(name);
	// hSigPdf->Write();


	// TH1D *err = (TH1D*)hSigPdf->Clone("error");
	// TH1D *tempErr[4];
	// err->Reset();
	// for(int s = 0; s < 4; s++)
	// {
	//   if(s == 0)
	//     tempErr[s] = (TH1D*)TFile::Open("photonPurity_sys_loose.root")->Get(name);
	//   else if(s ==1)
	//     tempErr[s] = (TH1D*)TFile::Open("photonPurity_sys_tight.root")->Get(name);
	//   else if(s ==2)
	//     tempErr[s] = (TH1D*)TFile::Open("photonPurity_sys_sigshift.root")->Get(name);
	//   else if(s ==3)
	//     tempErr[s] = (TH1D*)TFile::Open("photonPurity_sys_bkgshift.root")->Get(name);
	//   tempErr[s]->Divide(hSigPdf);
	//   for (Int_t l=1; l<=tempErr[s]->GetNbinsX();l++)
	//   {
	//     tempErr[s]->SetBinContent(l, TMath::Abs(tempErr[s]->GetBinContent(l))-1);
	//   }
	// }
	// for (Int_t l=1; l<=err->GetNbinsX();l++)
	// {
	//   Double_t errVal = TMath::Sqrt(tempErr[0]->GetBinContent(l)*tempErr[0]->GetBinContent(l) +
	// 				tempErr[1]->GetBinContent(l)*tempErr[1]->GetBinContent(l) +
	// 				tempErr[2]->GetBinContent(l)*tempErr[2]->GetBinContent(l) +
	// 				tempErr[3]->GetBinContent(l)*tempErr[3]->GetBinContent(l)
	//     );
	//   err->SetBinContent(l, errVal);
	// }

	// plot stacked histos
	handsomeTH1(hSigPdf);
	mcStyle(hSigPdf);
	sbStyle(hBckPdf);
	cleverRange(hSigPdf,1.5);
	hSigPdf->SetAxisRange(0.001,0.024,"X");
	hSigPdf->SetNdivisions(505);
	hSigPdf->GetYaxis()->SetTitleOffset(1.75);
	hSigPdf->SetYTitle("Entries");
	hSigPdf->SetXTitle("#sigma_{#eta #eta}");

	hSigPdf->DrawCopy("hist");
	//drawSys(hSigPdf, err, kRed, -1, 0.001);
	hBckPdf->DrawCopy("same hist");
	hData1->DrawCopy("same");

	Float_t xpos = 0.44;
	if(2*(k+j)*nPTBINS+i+1 == 1)
	  xpos = 0.54;

	TLegend *t3=new TLegend(xpos, 0.45, 0.92, 0.71);
	t3->AddEntry(hData1,LABEL,"pl");
	t3->AddEntry(hSigPdf,"Signal","lf");
	t3->AddEntry(hBckPdf,"Background","lf");
	t3->SetFillColor(0);
	t3->SetBorderSize(0);
	t3->SetFillStyle(0);
	t3->SetTextFont(43);
	t3->SetTextSize(20);
	//if(i == 0)
	// TH1D *dummyHist = new TH1D("dummyHist","",10,0,10);
	// dummyHist->Fill(1);
	// dummyHist->SetFillColor(kRed);
	// dummyHist->SetLineColor(kRed);
	// dummyHist->SetFillStyle(1001);
	// t3->AddEntry(dummyHist,"MC Sys. Error","f");
	// if(i == 0)
	//   t3->Draw();

	if(i == 3)
	{
	  drawText("CMS Preliminary", xpos, 0.68,1,20);
	  drawText("PbPb #sqrt{s}_{_{NN}}=5.02 TeV", xpos, 0.60,1,20);
	  drawText("#intL = 404 #ub^{-1}", xpos, 0.50,1,20);
	}



	//drawText("|#eta_{#gamma}| < 1.479",0.5680963,0.9);
	//drawText(Form("%f shift",fitr.sigMeanShift),0.57,0.82);
	//drawText("Background Correction",0.57,0.82);
	//drawText("bkg Tighter",0.57,0.82);
	//if(nPTBINS != 1)
	drawText(Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",
		      PTBINS[i], PTBINS[i+1]),
		 xpos, 0.90,1,20);
	// if(/*nCENTBINS != 1 && */i ==0)
	drawText(Form("%.0f - %.0f%c",
		      CENTBINS[j]/2., CENTBINS[j+1]/2.,'%'),
		 xpos, 0.82,1,20);
	// if(nETABINS != 1)
	//   drawText(Form("%.3f < #eta_{#gamma} < %.3f",
	// 		ETABINS[k], ETABINS[k+1]),
	// 	   xpos, 0.82,1,20);
	drawText(Form("Purity (#sigma_{#eta#eta} < 0.01) : %.2f", (Float_t)fitr.purity),
		 xpos, 0.76,1,20);
	drawText(Form("#chi^{2}/ndf : %.2f", (Float_t)fitr.chisq),
		 xpos, 0.45,1,20);

    /////// MAKE HIST //////////
    //TH2D* h2D_Num = new TH2D("h2D_Num", ";|#eta|;p_{T} (GeV)", nEtaBin, etaBins, nPtBin, ptBins);
    h1D_p[j]->SetBinContent(i,(Float_t)fitr.purity);        



	// //plot ratio
	// cPurity->cd((2*(j+k)+1)*nPTBINS+i+1);
	// //cPurity[i]->cd(nETABINS + k+ 1);
	// TH1D* ratio = (TH1D*)hData1->Clone("ratio");
	// ratio->Divide(hData1, hSigPdf, 1, 1);
	// ratio->SetMinimum(0);
	// ratio->SetMaximum(3);
	// ratio->SetXTitle("#sigma_{#eta #eta}");
	// ratio->GetXaxis()->CenterTitle();
	// ratio->SetYTitle("Data/Fit");
	// ratio->GetYaxis()->CenterTitle();
	// ratio->DrawCopy("E");
	// TLine *line = new TLine(0,1,maxSIGMA,1);
	// line->SetLineStyle(2);
	// line->Draw("same");

	// TString savename = Form("purity_pA_barrel_pt%.0f_hf%.0f_plot",
	// 			PTBINS[i], CENTBINS[j]);
	// cPurity[i*nCENTBINS+j]->SaveAs(savename+".C");
	// cPurity[i*nCENTBINS+j]->SaveAs(savename+".pdf");
	// cPurity[i*nCENTBINS+j]->SaveAs(savename+".png");

      }
    }
    //cPurity[i]->SaveAs(Form("pPb_purity_etadep_wshift_ptbin%.0f.png",PTBINS[i]));
    //cPurity[i]->SaveAs(Form("pPb_purity_etadep_noshift_inclusive.png"));
  }
  outFile->cd();
  configTree->Write();
  cPurity->Write();
    for(int j=0;j<temp_nCentBin;j++){
  h1D_p[j]->Write();
    }
  outFile->Close();
  //cPurity->SaveAs(SAVENAME+".C");
  //cPurity->SaveAs(SAVENAME+".png");
  //cPurity->SaveAs(SAVENAME+".pdf");
}

int main(int argc, char **argv)
{
  if(argc == 6){
    quickPhotonPurity(argv[1], argv[2], argv[3], argv[4], argv[5]);
  return 0;
  } else {
    cout << "wrong argument" << endl;
      return 1;
  }
}
