#ifndef UTILITY_Yeonju_H
#define UTILITY_Yeonju_H
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLine.h>
#include <TBox.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGaxis.h>
#include <TDatime.h>
#include <iostream>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCut.h>

#include <iostream>     // std::cout
#include <ctime>        // std::clock()
#include <algorithm>    // std::find()
#include <iomanip>      // std::setprecision()
#include <vector>
using namespace std;

const int col[] = {1,2,3,4,6,7,28,46,41};
const int ycol[] = {8,9,28,46,41};
const int marker[] = {24,25,26,27,28,29,31,33,34};

TH1D* changeToPrettyBin(TH1D* h){
    TH1D* hnew = new TH1D(Form("%s_prettyBin",h->GetName()),Form("%s",h->GetTitle()),nPtBin,ptBins_draw);
    for(int i=0;i<h->GetNbinsX();++i){
        hnew->SetBinContent(i+1,h->GetBinContent(i+1));
        float err = h->GetBinError(i+1);
        if(err!=0) hnew->SetBinError(i+1,err);
        else hnew->SetBinError(i+1,0.0001);
    }
    return hnew;
}
void saveHistogramsToPicture(TH1* h, const char* fileType="pdf", const char* caption="", const char* directoryToBeSavedIn="figures", int styleIndex=1, int rebin =1){
    TCanvas* c1=new TCanvas();
    if(rebin!=1)
    {
        h->Rebin(rebin);
    }

    if(styleIndex==1)
    {
        h->Draw("E");
    }
    else
    {
        h->Draw();
        if(h->InheritsFrom("TH2"))
        {
            h->Draw("COLZ TEXT");    // default plot style for TH2 histograms
        }
    }

    if(strcmp(directoryToBeSavedIn, "") == 0)   // save in the current directory if no directory is specified
    {
        c1->SaveAs(Form("%s_%s.%s" ,h->GetName(),caption, fileType));  // name of the file is the name of the histogram
    }
    else
    {
        c1->SaveAs(Form("%s/%s_%s.%s", directoryToBeSavedIn ,h->GetName(),caption, fileType));
    }
    c1->Close();
}

void yjStyleRoot(){
    gStyle -> SetOptStat(0);
    TH1::SetDefaultSumw2();
}
void SetyjPadStyle(){
    gStyle->SetPaperSize(20,26);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.17);
    gStyle->SetPadLeftMargin(0.17);
}
void SetHistTitleStyle(double titlesize=0.06, double labelsize = 0.05){
    gStyle->SetTextFont(42); 
    gStyle->SetTextSize(0.04); 
    gStyle->SetLabelFont( 42, "X" ); 
    gStyle->SetLabelFont( 42, "Y" ); 
    gStyle->SetLabelFont( 42, "Z" ); 
    gStyle->SetTitleSize( titlesize, "X" ); 
    gStyle->SetTitleSize( titlesize, "Y" ); 
    gStyle->SetTitleSize( titlesize, "Z" ); 
    gStyle->SetLabelSize( labelsize, "X" ); 
    gStyle->SetLabelSize( labelsize, "Y" ); 
    gStyle->SetLabelSize( labelsize, "Z" ); 
}
//void SetHistTitleStyle(double titlesize=0.06, double titleoffset=0.04, double labelsize = 0.05, double labeloffset=0.01){
//    gStyle->SetTitleSize( titlesize, "X" ); gStyle->SetTitleOffset(titleoffset, "X");
 //   gStyle->SetTitleSize( titlesize, "Y" ); gStyle->SetTitleOffset(titleoffset, "Y");
  //  gStyle->SetLabelSize( labelsize, "X" ); gStyle->SetLabelOffset(labeloffset, "X");
   // gStyle->SetLabelSize( labelsize, "Y" ); gStyle->SetLabelOffset(labeloffset, "Y");
//}
void SetHistTitleStyle(TH1* h, double titlesize=0.08, double titleoffset=0.01, double labelsize=0.05,double labeloffset=0.01){
    h->SetTitleSize( titlesize, "X" ); h->SetTitleOffset(titleoffset, "X");
    h->SetTitleSize( titlesize, "Y" ); h->SetTitleOffset(titleoffset, "Y");
    h->SetLabelSize( labelsize, "X" ); h->SetLabelOffset(labeloffset, "X");
    h->SetLabelSize( labelsize, "Y" ); h->SetLabelOffset(labeloffset, "Y");
}
void thisPadStyle(){
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.08);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);
}
void SetPadStyle(){
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.08);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTopMargin(0.10);
}
void legStyle( TLegend *a)
{
  a->SetBorderSize(0);
  a->SetFillStyle(0);
//  a->SetHeader(head);
//  a->SetTextFont(62);
//  a->SetTextSize(17);
//  a->SetLineColor(1);
//  a->SetLineStyle(1);
//  a->SetLineWidth(1);
//  a->SetFillColor(0);

}
void graphStyle(TGraph *g1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t mstyle=20, Int_t mcolor=1, Int_t lwidth=1, Double_t msize=1.0)
{
	g1->SetLineStyle(lstyle);
	g1->SetLineColor(lcolor);
	g1->SetLineWidth(lwidth);
	g1->SetMarkerStyle(mstyle);
	g1->SetMarkerColor(mcolor);
	g1->SetMarkerSize(msize);
}

void hLineStyle(TH1 *h1=0, Int_t lstyle=1, Int_t lcolor=1, Int_t lwidth=1, Int_t lfst=0, Int_t lfcolor=0)
{
	h1->SetLineStyle(lstyle);
	h1->SetLineColor(lcolor);
	h1->SetLineWidth(lwidth);
	h1->SetFillStyle(lfst);
	h1->SetFillColor(lfcolor);
}

void hMarkerStyle(TH1 *h1=0, Int_t mstyle=20, Int_t mcolor=1, Double_t msize=1.0)
{
	h1->SetMarkerStyle(mstyle);
	h1->SetMarkerColor(mcolor);
	h1->SetMarkerSize(msize);
}
void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=14){
	TLatex *tex = new TLatex(xp,yp,text);
	tex->SetTextFont(43);
	//   if(bold)tex->SetTextFont(43);
	tex->SetTextSize(textSize);
	tex->SetTextColor(textColor);
	tex->SetLineWidth(1);
	tex->SetNDC();
	tex->Draw();
}
void jumSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
	TLine* t1 = new TLine(x1,y1,x2,y2);
	t1->SetLineWidth(width);
	t1->SetLineStyle(7);
	t1->SetLineColor(color);
	t1->Draw();
}

void onSun(Double_t x1=0,Double_t y1=0,Double_t x2=1,Double_t y2=1,Int_t color=1, Double_t width=1)
{
    TLine* t1 = new TLine(x1,y1,x2,y2);
    t1->SetLineWidth(width);
    t1->SetLineStyle(1);
    t1->SetLineColor(color);
    t1->Draw();
}
double findCross(TH1* h1, TH1* h2, double& frac, double& effi, double& fracErr, double& effiErr){
	Int_t nBins = h1->GetNbinsX();
	double crossVal =0;
	int binAt0 = h1->FindBin(0);
	for(Int_t ix=binAt0; ix<=nBins ;ix++){
		float yy1 = h1->GetBinContent(ix);
		float yy2 = h2->GetBinContent(ix);
		if(yy2>yy1) {
			crossVal= h1->GetBinLowEdge(ix);
			break;
		}
	}
	int crossBin = h1->FindBin(crossVal);
	frac = 1 - (h2->Integral(1,crossBin) / h1->Integral(1,crossBin) );
	effi = ( h1->Integral(1,crossBin) / h1->Integral() );
	fracErr = frac * TMath::Sqrt( (1./h2->Integral(1,crossVal)) + (1./h1->Integral(1,crossVal)) );
	effiErr = ( TMath::Sqrt(h1->Integral(1,crossVal)) / h1->Integral() ) * TMath::Sqrt(1 - (h1->Integral(1,crossVal)/h1->Integral()) );

	return crossVal;
}

void makeMultiPanelCanvas(TCanvas*& canv, const Int_t columns,
		const Int_t rows, const Float_t leftOffset=0.,
		const Float_t bottomOffset=0.,
		const Float_t leftMargin=0.2,
		const Float_t bottomMargin=0.2,
		const Float_t edge=0.05) {
	if (canv==0) {
		//Error("makeMultiPanelCanvas","Got null canvas.");
		return;
	}
	canv->Clear();

	TPad* pad[columns][rows];

	Float_t Xlow[columns];
	Float_t Xup[columns];
	Float_t Ylow[rows];
	Float_t Yup[rows];
	Float_t PadWidth =
		(1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
				(1.0/(1.0-edge))+(Float_t)columns-2.0);
	Float_t PadHeight =
		(1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
				(1.0/(1.0-edge))+(Float_t)rows-2.0);
	Xlow[0] = leftOffset;
	Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
	Xup[columns-1] = 1;
	Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

	Yup[0] = 1;
	Ylow[0] = 1.0-PadHeight/(1.0-edge);
	Ylow[rows-1] = bottomOffset;
	Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

	for(Int_t i=1;i<columns-1;i++) {
		Xlow[i] = Xup[0] + (i-1)*PadWidth;
		Xup[i] = Xup[0] + (i)*PadWidth;
	}
	Int_t ct = 0;
	for(Int_t i=rows-2;i>0;i--) {
		Ylow[i] = Yup[rows-1] + ct*PadHeight;
		Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
		ct++;
	}
	TString padName;
	for(Int_t i=0;i<columns;i++) {
		for(Int_t j=0;j<rows;j++) {
			canv->cd();
			padName = Form("p_%d_%d",i,j);
			pad[i][j] = new TPad(padName.Data(),padName.Data(),
					Xlow[i],Ylow[j],Xup[i],Yup[j]);
			if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
			else pad[i][j]->SetLeftMargin(0);

			if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
			else pad[i][j]->SetRightMargin(0);

			if(j==0) pad[i][j]->SetTopMargin(edge);
			else pad[i][j]->SetTopMargin(0);

			if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
			else pad[i][j]->SetBottomMargin(0);

			pad[i][j]->Draw();
			pad[i][j]->cd();
			pad[i][j]->SetNumber(columns*j+i+1);
		}
	}
}

Double_t getDPHI( Double_t phi1, Double_t phi2) {
        Double_t dphi = phi1 - phi2;

        if ( dphi > 3.141592653589 )
                dphi = dphi - 2. * 3.141592653589;
        if ( dphi <= -3.141592653589 )
                dphi = dphi + 2. * 3.141592653589;

        if ( TMath::Abs(dphi) > 3.141592653589 ) {
                std::cout << " commonUtility::getDPHI error!!! dphi is bigger than 3.141592653589 " << std::endl;
        }

        return dphi;
}

Double_t getDETA(Double_t eta1, Double_t eta2){
        return eta1 - eta2;
}

Double_t getDR( Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2){
        Double_t theDphi = getDPHI( phi1, phi2);
        Double_t theDeta = eta1 - eta2;
        return TMath::Sqrt ( theDphi*theDphi + theDeta*theDeta);
}

Double_t cleverRange(TH1* h,Float_t fac=1.2, Float_t minY=1.e-3)
{
   Float_t maxY =  fac * h->GetBinContent(h->GetMaximumBin());
   //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
   h->SetAxisRange(minY,maxY,"Y");
   return maxY;
}


Double_t getCleverRange(TH1* h)
{
  Double_t maxY = -1000000;
  for ( Int_t ibin = 1 ; ibin <= h->GetNbinsX() ; ibin++) {
    if (maxY < h->GetBinContent(ibin) )
      maxY = h->GetBinContent(ibin);
  }
  return maxY;
}

Double_t cleverRange(TH1* h,TH1* h2, Float_t fac=1.2, Float_t minY=1.e-3)
{
  Float_t maxY1 =  fac * h->GetBinContent(h->GetMaximumBin());
  Float_t maxY2 =  fac * h2->GetBinContent(h2->GetMaximumBin());

  //   cout <<" range will be set as " << minY << " ~ " << maxY << endl;
  h->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  h2->SetAxisRange(minY,max(maxY1,maxY2),"Y");
  return max(maxY1,maxY2);
}

TF1* cleverGaus(TH1* h, const char* title="h", Float_t c = 2.5, bool quietMode=true)
{
    if ( h->GetEntries() == 0 )
    {
        TF1 *fit0 = new TF1(title,"gaus",-1,1);
        fit0->SetParameters(0,0,0);
        return fit0;
    }

    Int_t peakBin  = h->GetMaximumBin();
    Double_t peak =  h->GetBinCenter(peakBin);
    Double_t sigma = h->GetRMS();

    TF1 *fit1 = new TF1(title,"gaus",peak-c*sigma,peak+c*sigma);
    fit1->SetParameter(1, 1.0);
    fit1->SetParameter(1, 0.0005);
    if (quietMode) h->Fit(fit1,"LL M O Q R");
    else    h->Fit(fit1,"LL M O Q R");
    return fit1;
}

void SetHistColor(TH1* h, Int_t color=1)
{
    h->SetMarkerColor(color);
    h->SetLineColor(color);
} 

float mean(float data[], int n)
{
    float mean=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean+=data[i];
    }
    mean=mean/n;
    return mean;
}

float standard_deviation(float data[], int n)
{
    float mean=0.0, sum_deviation=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean+=data[i];
    }
    mean=mean/n;
    for(i=0; i<n;++i)
        sum_deviation+=(data[i]-mean)*(data[i]-mean);
    return sqrt(sum_deviation/n);
}

void normHist(TH1* h=0, TH1* hNominal=0, double cut_i=700, double cut_f=900){
    int cutBinFrom = h->FindBin(cut_i);
    int cutBinTo = h->FindBin(cut_f);
    h->Scale(hNominal->Integral(cutBinFrom,cutBinTo)/h->Integral(cutBinFrom,cutBinTo));
}


#endif
