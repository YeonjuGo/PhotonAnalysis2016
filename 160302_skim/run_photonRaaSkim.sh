#void photonRaaSkim(const TString configFile, const TString inputFile, const TString outputFile, COLL::TYPE)
# COLL::TYPE = COLL::kHI, COLL::kHIMC, COLL::kPP, COLL::kPPMC, COLL::kPA, COLL::kPAMC
root -l -b -q 'photonRaaSkim.C++("../../CutConfigurations/photonRaa.conf","inputFiles/ppAllQCD.txt","/home/goyeonju/CMS/Files/photon2016/photonRaaSkimed_ppAllQCD_5Apr2016.root", COLL::kPPMC)'
root -l -b -q 'photonRaaSkim.C++("../../CutConfigurations/photonRaa.conf","inputFiles/ppEmEnr.txt","/home/goyeonju/CMS/Files/photon2016/photonRaaSkimed_ppEmEnr_5Apr2016.root", COLL::kPPMC)'
root -l -b -q 'photonRaaSkim.C++("../../CutConfigurations/photonRaa.conf","inputFiles/pbpbAllQCD.txt","/home/goyeonju/CMS/Files/photon2016/photonRaaSkimed_pbpbAllQCD_5Apr2016.root", COLL::kHIMC)'
root -l -b -q 'photonRaaSkim.C++("../../CutConfigurations/photonRaa.conf","inputFiles/pbpbEmEnr.txt","/home/goyeonju/CMS/Files/photon2016/photonRaaSkimed_pbpbEmEnr_5Apr2016.root", COLL::kHIMC)'
root -l -b -q 'photonRaaSkim.C++("../../CutConfigurations/photonRaa.conf","inputFiles/ppDATA.txt","/home/goyeonju/CMS/Files/photon2016/photonRaaSkimed_ppDATA_5Apr2016.root", COLL::kPP)'
#root -l -b -q 'photonRaaSkim.C++("../../CutConfigurations/photonRaa.conf","inputFiles/pbpbDATA.txt","/home/goyeonju/CMS/Files/photon2016/photonRaaSkimed_pbpbDATA_5Apr2016.root", COLL::kHI)'
