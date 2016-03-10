// 2016 Mar 3. comparing two datasets (ex. allqcdphoton and emenricheddijet) 
// Author : Yeonju Go

#include "../yjUtility.h"

void compareTwo(TTree* t1=0 ,TTree* t2=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, TCut cut1="",TCut cut2="", const char* cap="", int numbering=1);
void compare_two(TString coll="pbpb"){

    const char* fname_1;
    const char* fname_2;
    if(coll=="pp"){ 
        fname_1="/d3/scratch/goyeonju/files/photons2016/2015-PP-MC_Pythia8_Photon30_pp502_TuneCUETP8M1.root";//AllQCDPhoton
        fname_2="/d3/scratch/goyeonju/files/photons2016/2015-PP-MC_Pythia8_EmEnrDijet30_pp502_TuneCUETP8M1.root";//EmEnrDijet
    } else if(coll=="pbpb"){
        fname_1="/d3/scratch/goyeonju/files/photons2016/2015-PbPb-MC_AllQCDPhoton30_v1/0.root";
        fname_2="/d3/scratch/goyeonju/files/photons2016/2015-PP-MC_Pythia8_EmEnrDijet30_pp502_TuneCUETP8M1.root";
    }
    TFile* f1 = new TFile(fname_1);
    TTree* t1 = (TTree*) f1 -> Get("ggHiNtuplizer/EventTree");
    TFile* f2 = new TFile(fname_2);
    TTree* t2 = (TTree*) f2 -> Get("ggHiNtuplizer/EventTree");
    
    TCut dataTotCut = "";
    TCut mcTotCut = "";
    const char* cap = "";
    compareTwo(t1, t2, "phoEt",50, 60, 250,dataTotCut, mcTotCut,cap);
    compareTwo(t1, t2, "phoR9",50, 0, 1,dataTotCut, mcTotCut,cap);
    compareTwo(t1, t2, "phoSigmaIEtaIEta",50, 0, 0.06,dataTotCut, mcTotCut,cap);
    compareTwo(t1, t2, "phoHoverE",50, 0, 1,dataTotCut, mcTotCut,cap);
    compareTwo(t1, t2, "phoHoverE",50, 0, 0.05,dataTotCut, mcTotCut,cap);
    /* 
    const int Neta = 4;
    TCut etCut;
    TCut trigCut[Neta];
    TCut etaCut[Neta];
    TCut hoeCut[Neta];
    TCut sigmaCut[Neta];
    TCut isoCut[Neta];
    TCut spikeRejCut[Neta];
    etaCut[0] = "abs(phoEta)<=2.4";
    etaCut[1] = "abs(phoEta)<=1.44";
    etaCut[2] = "abs(phoEta)>1.44 && abs(phoEta)<=2";
    etaCut[3] = "abs(phoEta)>2.0 && abs(phoEta)<=2.4";
    spikeRejCut[0] = "(( (abs(phoEta)<=1.44) && phoSigmaIEtaIEta>0.002 && pho_swissCrx<0.9 && abs(pho_seedTime)<3) ||(abs(phoEta)>1.44 && abs(phoEta)<2.4) )";
    spikeRejCut[1] = "((phoSigmaIEtaIEta>0.002 && pho_swissCrx<0.9 && abs(pho_seedTime)<3))";
    spikeRejCut[2] = "";
    spikeRejCut[3] = "";

    if(coll=="pp"){ 
        etCut = "phoEt>40";
        trigCut[0] = "HLT_HISinglePhoton30_Eta3p1_v1";
        trigCut[1] = "HLT_HISinglePhoton30_Eta3p1_v1";
        trigCut[2] = "HLT_HISinglePhoton30_Eta3p1_v1";
        trigCut[3] = "HLT_HISinglePhoton30_Eta3p1_v1";
        //GED
        isoCut[0] = "( (abs(phoEta)<=1.44) && ((pfcIso4<=1.37) && (pfnIso4<=1.06+0.014*phoEt+0.000019*phoEt*phoEt) && pfpIso4<=(0.28+0.0053*phoEt)) ) || ( (abs(phoEta)>1.44 && abs(phoEta)<2.4) && ((pfcIso4<=1.10) && (pfnIso4<=2.69+0.0139*phoEt+0.000025*phoEt*phoEt) && pfpIso4<=(0.39+0.0034*phoEt)) )";  
        isoCut[1] = "((pfcIso4<=1.37) && (pfnIso4<=1.06+0.014*phoEt+0.000019*phoEt*phoEt) && pfpIso4<=(0.28+0.0053*phoEt))";  
        isoCut[2] = "((pfcIso4<=1.10) && (pfnIso4<=2.69+0.0139*phoEt+0.000025*phoEt*phoEt) && pfpIso4<=(0.39+0.0034*phoEt))";  
        isoCut[3] = "((pfcIso4<=1.10) && (pfnIso4<=2.69+0.0139*phoEt+0.000025*phoEt*phoEt) && pfpIso4<=(0.39+0.0034*phoEt))";  
        hoeCut[0] = "phoHoverE<0.05";
        hoeCut[1] = "phoHoverE<0.05";
        hoeCut[2] = "phoHoverE<0.05";
        hoeCut[3] = "phoHoverE<0.05";
        sigmaCut[0] = "( (abs(phoEta)<=1.44) && phoSigmaIEtaIEta<0.0102 ) || ((abs(phoEta)>1.44 && abs(phoEta)<2.4) && phoSigmaIEtaIEta<0.0268)";
        sigmaCut[1] = "phoSigmaIEtaIEta<0.0102";
        sigmaCut[2] = "phoSigmaIEtaIEta<0.0268";
        sigmaCut[3] = "phoSigmaIEtaIEta<0.0268";

    } else if(coll=="pbpb"){

        etCut = "phoEt>60";
        trigCut[0] = "HLT_HISinglePhoton50_Eta3p1_v1";
        trigCut[1] = "HLT_HISinglePhoton50_Eta3p1_v1";
        trigCut[2] = "HLT_HISinglePhoton50_Eta3p1_v1";
        trigCut[3] = "HLT_HISinglePhoton50_Eta3p1_v1";
        //GED
        isoCut[0] = "( (abs(phoEta)<=1.44) && ((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<1.0) ) || ( (abs(phoEta)>1.44 && abs(phoEta)<2.4) && ((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<1.0) )";  
        isoCut[1] = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<1.0)";
        isoCut[2] = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<1.0)";
        isoCut[3] = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<1.0)";
        hoeCut[0] = "phoHoverE<0.10";
        hoeCut[1] = "phoHoverE<0.10";
        hoeCut[2] = "phoHoverE<0.10";
        hoeCut[3] = "phoHoverE<0.10";
        sigmaCut[0] = "( (abs(phoEta)<=1.44) && phoSigmaIEtaIEta<0.0102 ) || ((abs(phoEta)>1.44 && abs(phoEta)<2.4) && phoSigmaIEtaIEta<0.0268)";
        sigmaCut[1] = "phoSigmaIEtaIEta<0.0102";
        sigmaCut[2] = "phoSigmaIEtaIEta<0.0268";
        sigmaCut[3] = "phoSigmaIEtaIEta<0.0268";
    }

    TCut dataTotCut;
    TCut mcTotCut;

#if 1
    for(int j=0;j<Neta;++j){ 
        for(int icut=0;icut<4;++icut){
            const char* cap;
            if(icut==0) {
                // spike rejection 
                dataTotCut = trigCut[j] && etaCut[j] && etCut && spikeRejCut[j];
                mcTotCut = etaCut[j] && etCut && spikeRejCut[j];
                cap = Form("%s_eta%d_spike",coll.Data(),j);
            } else if(icut==1){
                // spike rejection + h/e
                dataTotCut = trigCut[j] && etaCut[j] && etCut && spikeRejCut[j] && hoeCut[j];
                mcTotCut = etaCut[j] && etCut && spikeRejCut[j] && hoeCut[j];
                cap = Form("%s_eta%d_spike_hoe",coll.Data(),j);
            } else if(icut==2){
                // spike rejection + h/e +sigmaIEtaIEta
                dataTotCut = trigCut[j] && etaCut[j] && etCut && spikeRejCut[j] && hoeCut[j] && sigmaCut[j] ;
                mcTotCut = etaCut[j] && etCut && spikeRejCut[j] && hoeCut[j] && sigmaCut[j] ;
                cap = Form("%s_eta%d_spike_hoe_sigma",coll.Data(),j);
            } else if(icut==3){
                // spike rejection + h/e +sigmaIEtaIEta + isolation
                dataTotCut = trigCut[j] && etaCut[j] && etCut && isoCut[j] && spikeRejCut[j] && hoeCut[j] && sigmaCut[j] ;
                mcTotCut = etaCut[j] && etCut && isoCut[j] && spikeRejCut[j] && hoeCut[j] && sigmaCut[j] ;
                cap = Form("%s_eta%d_spike_hoe_sigma_iso",coll.Data(),j);
            }

            if(j==0) compareTwo(t1, t2, "phoEta",50, -2.4, 2.4, dataTotCut, mcTotCut,cap);
            else {
                compareTwo(t1, t2, "phoEt",50, 60, 250,dataTotCut, mcTotCut,cap);
                compareTwo(t1, t2, "phoR9",50, 0, 1,dataTotCut, mcTotCut,cap);
                compareTwo(t1, t2, "phoSigmaIEtaIEta",50, 0, 0.06,dataTotCut, mcTotCut,cap);
                compareTwo(t1, t2, "phoHoverE",50, 0, 1,dataTotCut, mcTotCut,cap);
                compareTwo(t1, t2, "phoHoverE",10, 0, 0.05,dataTotCut, mcTotCut,cap);
            }
        }
    }
#endif
*/
} // main function

void compareTwo(TTree* t1, TTree* t2, TString var, int nBins, double xMin, double xMax, TCut cut1, TCut cut2, const char* cap)  {
    SetHistTitleStyle();
    SetyjPadStyle();
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    static int i = 1;
    TCanvas* c =new TCanvas(Form("c_%s_%d",var.Data(),i),"", 400,800);
    //TCanvas* c =new TCanvas(Form("c_%s_%d",var.Data(),i),"", 400,800);
	c->Divide(1,2);
	c->cd(1);
	gPad->SetLogy();
    double ptbin[] = {60,70,80,90,100,110,120,130,140,150,170,200,250};
    int nptbin = sizeof(ptbin)/sizeof(double) -1;
	TH1D* h1;
    if(var=="phoEt") h1 = new TH1D(Form("h1_%s_%d",var.Data(),i), Form(";%s;",var.Data()),nptbin, ptbin);
    else h1 = new TH1D(Form("h1_%s_%d",var.Data(),i), Form(";%s;",var.Data()), nBins,xMin,xMax);
	TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%d",var.Data(),i));
/*	if(var=="phoEt") h1 = new TH1D(Form("h1_%s_%d",var.Data(),i), Form(";%s;",var.Data()),nptbin, ptbin);
    else h1 = new TH1D(Form("h1_%s_%d",var.Data(),i), Form(";%s;",var.Data()), nBins,xMin,xMax);
	TH1D* h2 = (TH1D*)h1->Clone(Form("h2_%s_%d",var.Data(),i));
*/
    h1->Sumw2();
	h2->Sumw2();
	t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), cut1);
	t2->Draw(Form("%s>>%s",var.Data(),h2->GetName()), cut2);	
	//t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), Form("HLT_HISinglePhoton20_Eta3p1_v1_Prescl*(%s)",cut1.GetTitle()));
	//t1->Draw(Form("%s>>%s",var.Data(),h1->GetName()), Form("HLT_HISinglePhoton20_Eta3p1_v1_Prescl*HLT_HISinglePhoton30_Eta3p1_v1_Prescl*HLT_HISinglePhoton40_Eta3p1_v1_Prescl*HLT_HISinglePhoton50_Eta3p1_v1_Prescl*HLT_HISinglePhoton20_Eta3p1_v1_Prescl*(%s)",cut1.GetTitle()));
	h1->Scale( 1. / h1->Integral(),"width");
	h2->Scale( 1. / h2->Integral(),"width");
	//h1->Scale( 1. / h1->Integral("width"));
	//h2->Scale( 1. / h2->Integral("width"));
	//h1->Scale( 1. / t1->GetEntries(cut1));
	//h2->Scale( 1. / t2->GetEntries(cut2));
	h1->SetMarkerStyle(20);
	h1->SetMarkerSize(0.8);
	h1->SetMarkerColor(2);
	//h2->SetMarkerStyle(20);
	//h2->SetMarkerSize(0.8);
	double range;
    if(var=="phoEt") range = cleverRange(h1,h2,1.5,1.e-4);
    else if(var=="phoHoverE") range = cleverRange(h1,h2,1.5,1.e-4);
    else  range = cleverRange(h1,h2,1.5,1.e-2);
	h1->DrawCopy("L hist e");
	h2->DrawCopy("hist e same");
    TLegend* l1 = new TLegend(0.7,0.7,0.95,0.9);
    legStyle(l1);
    l1->AddEntry(h1,"AllQCDPhoton","pl");
    l1->AddEntry(h1,"EmEnrDijet","pl");
    l1->Draw();

	c->cd(2);
	h1->Divide(h2);
	h1->SetYTitle("AllQCDPhoton / EmEnrDijet Ratio");
	double ratioRange = getCleverRange(h1);
	if(ratioRange > 5) h1->SetAxisRange(0,5,"Y");
	else h1->SetAxisRange(0,2,"Y");
	h1->DrawCopy("le1");
	jumSun(xMin,1,xMax,1);
	c-> SaveAs(Form("pdf/compare_%s_%s.pdf",var.Data(),cap));
    i++;
}
