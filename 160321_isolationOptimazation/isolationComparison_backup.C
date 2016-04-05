// 2016 Mar 3. comparing two datasets (ex. allqcdphoton and emenricheddijet) 
// Author : Yeonju Go

#include "../yjUtility.h"
#include "../phoRaaCuts/phoRaaCuts_v1.h"


class TreeAndHist
{
public:
  TreeAndHist(TTree *t, TString var, int nBins, double xMin, double xMax, TCut cut);
  ~TreeAndHist(){
    delete hist;
  };

  
    TH1D* hist;
  //      int GetBinNumber(double c);
  Double_t evaluate(Double_t *x, Double_t *par);
  TH1D *histSig;
      TH1D *histBck;
      Double_t lowEdge;
      Double_t highEdge;
      Double_t nbin;
};

TreeAndHist::TreeAndHist(TTree *t, TString var, int nBins, double xMin, double xMax, TCut cut)
{
    hist = new TH1D(Form("%s_%s",t->GetName(),var.Data()),"",nBins, xMin, xMax);
    t->Draw(Form("%s>>%s",var.Data(),hist->GetName()),cut);




}
histFunction2::histFunction2(TH1D *h,TH1D *h2)
{
   histSig=(TH1D*)h->Clone();
   histBck=(TH1D*)h2->Clone();

   nbin=h->GetNbinsX();
   lowEdge=h->GetBinLowEdge(1);
   highEdge=h->GetBinLowEdge(nbin+1);

   histSig->SetName("hSig");
   histSig->Scale(1./histSig->Integral(1,histSig->GetNbinsX()+1));
   histBck->SetName("hBck");
   histBck->Scale(1./histBck->Integral(1,histBck->GetNbinsX()+1));
}

Double_t histFunction2::evaluate(Double_t *x, Double_t *par) {

   Double_t xx = x[0];
   int binNum=histSig->FindBin(xx);  //

   return par[0]*(histSig->GetBinContent(binNum)*par[1]+histBck->GetBinContent(binNum)*(1-par[1]));
}




void compareThree(TTree* t1=0 ,TTree* t2=0, TTree* t3=0,TString var="pt", int nBins=10, double xMin=0, double xMax=10, TCut cut1="",TCut cut2="",TCut cut3="", const char* cap="");
void isolationComparison(TString coll="pp"){

    const int nFile =3;
    const char* fname[nFile];
    if(coll=="pp"){ 
        fname[0]="/home/goyeonju/CMS/Files/photon2016/2015-Data-promptRECO-photonSkims_pp-photonHLTFilter-v0-HiForest.root";
        fname[1]="/home/goyeonju/CMS/Files/photon2016/gsfs-Pythia8_Photon_pp_RECO_forest_v28/gsfs-Pythia8_Photon15_30_50_80_120_pp_RECO_forest_v28_with_pthatWeight.root";
        fname[2]="/home/goyeonju/CMS/Files/photon2016/2015-PP-MC_Pythia8_EmEnrDijet30_pp502_TuneCUETP8M1.root";
    } else if(coll=="pbpb"){
        fname[0]="/home/goyeonju/CMS/Files/photon2016/forestSkimed_photonSkim_pbpb_2015data.root";
        fname[1]="/home/goyeonju/CMS/Files/photon2016/Pythia8_Photon_Hydjet_RECO_20160306_forest_v28_2/Pythia8_Photon15_30_50_80_120_Hydjet_RECO_20160306_forest_v28_2_with_pthatWeight.root";
        fname[2]="/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC_Pythia8_EmEnrichedDijet/2015-PbPb-MC_Pythia8_EmEnrichedDijet30_50_80_120_170_with_pthatWeight.root";
    }
    TFile* f[nFile];
    TTree* t[nFile];
    TTree* t_hi[nFile];
    TTree* t_skim[nFile];
    for(int i=0;i<nFile;i++){
        f[i] = new TFile(fname[i]);
        if(coll=="pp" && (i==0 || i==2)) {
            t[i] = (TTree*) f[i] -> Get("ggHiNtuplizerGED/EventTree");
            if(i==1) t_skim[i] = (TTree*) f[i] -> Get("HltTree");
            else t_skim[i] = (TTree*) f[i] -> Get("skimanalysis/HltTree");
            t_hi[i] = (TTree*) f[i] -> Get("hiEvtAnalyzer/HiTree");
        } else{ 
            t[i] = (TTree*) f[i] -> Get("EventTree");
            t_skim[i] = (TTree*) f[i] -> Get("HltTree");
            t_hi[i] = (TTree*) f[i] -> Get("HiTree");
        }
        t[i]->AddFriend(t_hi[i]);
        t[i]->AddFriend(t_skim[i]);
    }

    for(Int_t i = 0; i < 1; ++i) {
        //for(Int_t i = 0; i < nPtBin; ++i) {
        //TCut ptCut = Form("(phoEt >= %f) && (phoEt < %f)", ptBins[i], ptBins[i+1]);
        TCut ptCut = Form("(phoEt >= %f) && (phoEt < %f)",50.0,9999.0);
        TCut etaCut = Form("(abs(phoEta) >= %f) && (abs(phoEta) < %f)", etaBins[0], etaBins[1]);
        TCut filterCut = evtSelFilterCut;
        if(coll=="pp") filterCut = evtSelFilterCut_pp;
        TCut dataTotCut = filterCut && spikeRejection && ptCut && etaCut;
        TCut mcTotCut_bkg = filterCut && mcBkgIsolation && ptCut && etaCut;
        TCut mcTotCut_sig = filterCut && mcIsolation && ptCut && etaCut;

        cout << dataTotCut.GetTitle() << endl;
        cout << mcTotCut_bkg.GetTitle() << endl;
        cout << mcTotCut_sig.GetTitle() << endl;

        const char* cap = Form("%s_pt%dto%d_barrel",coll.Data(),50,9999);
        //const char* cap = Form("%s_pt%dto%d",coll.Data(),(int)ptBins[i],(int)ptBins[i+1]);
        if(coll=="pp") {
            compareThree(t[0], t[1], t[2], "phoSigmaIEtaIEta_2012",50, 0, 0.025, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_ecalClusterIsoR4",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_hcalRechitIsoR4",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_trackIsoR4PtCut20",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_ecalClusterIsoR4+pho_hcalRechitIsoR4+pho_trackIsoR4PtCut20",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(%s)",mcTotCut_sig.GetTitle()), cap);

        } else if(coll=="pbpb") {
            compareThree(t[0], t[1], t[2], "pho_ecalClusterIsoR4",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(pthatWeight)*(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_ecalClusterIsoR4",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(pthatWeight)*(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_hcalRechitIsoR4",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(pthatWeight)*(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_trackIsoR4PtCut20",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(pthatWeight)*(%s)",mcTotCut_sig.GetTitle()), cap);
            compareThree(t[0], t[1], t[2], "pho_ecalClusterIsoR4+pho_hcalRechitIsoR4+pho_trackIsoR4PtCut20",50, -10, 120, dataTotCut, Form("(pthatWeight)*(%s)",mcTotCut_bkg.GetTitle()), Form("(pthatWeight)*(%s)",mcTotCut_sig.GetTitle()), cap);

        }

    }  
    } // main function

    void compareThree(TTree* t1, TTree* t2,  TTree* t3, TString var, int nBins, double xMin, double xMax, TCut cut1, TCut cut2,TCut cut3, const char* cap)  {
        SetHistTitleStyle();
        SetyjPadStyle();
        gStyle->SetOptStat(0);
        TH1::SetDefaultSumw2();
        static int j = 1;
        TCanvas* c =new TCanvas(Form("c_%s_%d",var.Data(),j),"", 400,400);
        //c->Divide(1,2);
        //c->cd(1);
        gPad->SetLogy();
        double ptbin[] = {60,70,80,90,100,110,120,130,140,150,170,200,250};
        int nptbin = sizeof(ptbin)/sizeof(double) -1;
        TH1D* h[3];
        if(var=="phoEt") h[0] = new TH1D(Form("h1_%s_%d",var.Data(),j), Form(";%s;",var.Data()),nptbin, ptbin);
        else h[0] = new TH1D(Form("h1_%s_%d",var.Data(),j), Form(";%s;",var.Data()), nBins,xMin,xMax);
        h[1] = (TH1D*)h[0]->Clone(Form("h2_%s_%d",var.Data(),j));
        h[2] = (TH1D*)h[0]->Clone(Form("h3_%s_%d",var.Data(),j));

        for(int i=0;i<3;i++){
            h[i]->Sumw2();
            h[i]->SetMarkerStyle(20);
            h[i]->SetMarkerSize(0.8);
            h[i]->SetMarkerColor(i+1);
            h[i]->SetLineColor(i+1);
        } 
        t1->Draw(Form("%s>>%s",var.Data(),h[0]->GetName()), cut1);
        t2->Draw(Form("%s>>%s",var.Data(),h[1]->GetName()), cut2);	
        t3->Draw(Form("%s>>%s",var.Data(),h[2]->GetName()), cut3);	

        for(int i=0;i<3;i++){
            h[i]->Scale( 1. / h[i]->Integral());
            //h[i]->GetYaxis()->SetRange(1e-6,1.0);
            h[i]->SetAxisRange(1e-6,1.0,"Y");
            //h[i]->Scale( 1. / h[i]->Integral(),"width");
        }

        h[0]->DrawCopy("hist e");
        h[1]->DrawCopy("hist e same");
        h[2]->DrawCopy("hist e same");
        TLegend* l1 = new TLegend(0.65,0.75,0.95,0.9);
        legStyle(l1);
        l1->AddEntry(h[0],"DATA","p");
        l1->AddEntry(h[1],"AllQCDPhoton","pl");
        l1->AddEntry(h[2],"EmEnrichedDijet","pl");
        l1->Draw();

        drawText(cap,0.2,0.2);
        c-> SaveAs(Form("pdf/isolationComparison_%s_%s.pdf",var.Data(),cap));
        j++;
        /*
         *
         c->cd(2);
         h1->Divide(h2);
         h1->SetYTitle("DATA / MC");
         double ratioRange = getCleverRange(h1);
         if(ratioRange > 5) h1->SetAxisRange(0,5,"Y");
         else h1->SetAxisRange(0,2,"Y");
         h1->DrawCopy("le1");
         jumSun(xMin,1,xMax,1);
         i++;
         */
    }
