#ifndef phoRaaCuts_v1_H
#define phoRaaCuts_v1_H

#include <TCut.h>

#define PI 3.141592653589

TString ppDatafname = "/home/goyeonju/CMS/Files/photon2016/2015-Data-promptRECO-photonSkims_pp-photonHLTFilter-v0-HiForest.root";
TString ppMCfname = "/home/goyeonju/CMS/Files/photon2016/2015-PP-MC_Pythia8_Photon30_pp502_TuneCUETP8M1.root";
TString pbpbDatafname = "/home/goyeonju/CMS/Files/photon2016/forestSkimed_photonSkim_pbpb_2015data.root";
TString pbpbMCfname = "/home/goyeonju/CMS/Files/photon2016/2015-PbPb-MC-AllQCDPhoton_v1_ptHat_15_30_50_80_120_with_pthatWeight.root";


//const double ptBins[] = {20,30,40,50,60,70,80,90,100,120,140,180,220};
const double ptBins[] = {40,50,60,70,80};
const int nPtBin = sizeof(ptBins)/sizeof(double) - 1;
//const double etaBins[] = {0.0,1.44,2.0,2.5};
const double etaBins[] = {0.0,1.44};
const int nEtaBin = sizeof(etaBins)/sizeof(double) - 1;
const int centBins[] = {0,60,200};
const int nCentBin = sizeof(centBins)/sizeof(int) -1;


const TCut trigCut = "";
const TCut evtSelFilterCut = "pcollisionEventSelection";
const TCut evtSelFilterCut_pp = "pBeamScrapingFilter && pPAprimaryVertexFilter && pVertexFilterCutEandG";
const TCut spikeRejection = "(phoSigmaIEtaIEta_2012>=0.002) && (pho_swissCrx<=0.9) && (abs(pho_seedTime)<=3)";

const TCut dataCut = trigCut && evtSelFilterCut && spikeRejection;
const TCut dataCut_pp = trigCut && evtSelFilterCut_pp && spikeRejection;

const TCut hoeCut = "phoHoverE<0.1";
const TCut sumIsoCut = "(pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20) < 1.0";
const TCut sigmaCut = "(phoSigmaIEtaIEta_2012)<0.0102";
const TCut hotspotCut = "!(((phoSCEta > -1.194 && phoSCEta < -1.169) && (phoSCPhi > -1.211 && phoSCPhi < -1.162)) || ((phoSCEta > -0.935 && phoSCEta < -0.910) && (phoSCPhi > 2.036 && phoSCPhi < 2.065)) || ((phoSCEta > -0.850 && phoSCEta < -0.819) && (phoSCPhi > -1.723 && phoSCPhi < -1.676)) || ((phoSCEta > -0.764 && phoSCEta < -0.723) && (phoSCPhi > -0.786 && phoSCPhi < -0.720)) || ((phoSCEta > -0.237 && phoSCEta < -0.210) && (phoSCPhi > 0.106 && phoSCPhi < 0.147)) ||((phoSCEta > -0.062 && phoSCEta < -0.036) && (phoSCPhi > 0.860 && phoSCPhi< 0.931)) || ((phoSCEta > 0.973 && phoSCEta < 1.035) && (phoSCPhi > 0.623 && phoSCPhi < 0.697)) || ((phoSCEta > 1.075 && phoSCEta < 1.108) && (phoSCPhi > -3.114 && phoSCPhi < -3.054)) || ((phoSCEta > 1.426 && phoSCEta < 1.447) && (phoSCPhi > 2.744 && phoSCPhi < 2.774)) || ((phoSCEta > 0.915 && phoSCEta < 0.930) && (phoSCPhi > 1.69 && phoSCPhi < 1.72)))";

const TCut phoSignalCut = hoeCut && sumIsoCut && sigmaCut && hotspotCut;


const TCut isoCut_ppGED = "((pfcIso4<=1.37) && (pfnIso4<=1.06+0.014*phoEt+0.000019*phoEt*phoEt) && pfpIso4<=(0.28+0.0053*phoEt))";
const TCut hoeCut_ppGED = "phoHoverE<0.05";
const TCut sigmaCut_ppGED = "phoSigmaIEtaIEta<0.0102";

const TCut phoSignalCut_ppGED = isoCut_ppGED && hoeCut_ppGED && sigmaCut_ppGED;


const TCut sidebandIsolation = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)>10) && ((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<20) && phoHoverE<0.1";
//const TCut sidebandIsolation = "((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)>5) && ((pho_ecalClusterIsoR4 + pho_hcalRechitIsoR4 + pho_trackIsoR4PtCut20)<10) && phoHoverE<0.1";
const TCut mcIsolation = "(pho_genMatchedIndex!= -1) && mcCalIsoDR04[pho_genMatchedIndex]<5 && abs(mcPID[pho_genMatchedIndex])<=22";
const TCut mcBkgIsolation = "(pho_genMatchedIndex!= -1) && (mcCalIsoDR04[pho_genMatchedIndex]>=5 || abs(mcPID[pho_genMatchedIndex])>22)";

TString pbpbEff_fname_cent1 = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pbpb_cent0-60_v1_nohotspot.root";
TString pbpbEff_fname_cent2 = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pbpb_cent60-200_v1_nohotspot.root";
TString ppEff_fname = "/home/goyeonju/CMS/2016/PhotonAnalysis2016/160107_phoRaaEfficiency/output/hist_efficiency_total_pp_v1_nohotspot.root";

//TString pbpbPurity_fname = "../160111_phoRaaPurity/output/";
//TString ppPurity_fname = "../160111_phoRaaPurity/output/";

#endif

