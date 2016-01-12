#void makeGenMatchedPhotonTree_sep(const char* hiForestfileName="/afs/cern.ch/work/y/ygo/private/PRODUCTION/CMSSW_7_5_0/src/PFphoton/nominalForest/crab_AllQCDPhoton30/results/merged_AllQCDPhoton30_PbPb_5TeV_forest_pfisoAdded_genThr500MeV_2nd.root",
#   TString outName="AllQCDPhoton30", TString treePath="ggHiNtuplizer", bool isPromptPho=1, float ptThr=20)

#root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_EmEnrichedDijet30_ParticleFilter20GeV_eta24_TuneZ2_5020GeV.root", "EmEnrichedDijet30", "ggHiNtuplizer",1)'
#root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_EmEnrichedDijet30_ParticleFilter20GeV_eta24_TuneZ2_5020GeV.root", "EmEnrichedDijet30", "ggHiNtuplizerGED",1)'
root -l -b -q 'makeGenMatchedPhotonTree_sep.C++("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV.root", "AllQCDPhoton30", "ggHiNtuplizer",1)'
root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV.root", "AllQCDPhoton30", "ggHiNtuplizerGED",1)'

#root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_EmEnrichedDijet30_ParticleFilter20GeV_eta24_TuneZ2_5020GeV.root", "EmEnrichedDijet30", "ggHiNtuplizer",0)'
#root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_EmEnrichedDijet30_ParticleFilter20GeV_eta24_TuneZ2_5020GeV.root", "EmEnrichedDijet30", "ggHiNtuplizerGED",0)'
#root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV.root", "AllQCDPhoton30", "ggHiNtuplizer",0)'
#root -l -b -q 'makeGenMatchedPhotonTree_sep.C+("/afs/cern.ch/work/y/ygo/public/PFphoton/Pyquen_Unquenched_AllQCDPhoton30_PhotonFilter20GeV_eta24_TuneZ2_PbPb_5020GeV.root", "AllQCDPhoton30", "ggHiNtuplizerGED",0)'


