#!/bin/bash
#conf="../ElectroWeak-Jet-Track-Analyses/CutConfigurations/gamma-jet-nominal.conf"
conf="../ElectroWeak-Jet-Track-Analyses/CutConfigurations/photon.conf"
pbpbData="/d3/scratch/goyeonju/files/photons2016/forestSkimed_photonSkim_pbpb_2015data.root"
#pbpbMC="/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC_AllQCDPhoton30_v1/0.root"
pbpbMC="/u/user/goyeonju/scratch/files/photons2016/2015-PbPb-MC-AllQCDPhoton_v1_ptHat_15_30_50_80_120_with_pthatWeight.root"
#pbpbMC="/d3/scratch/goyeonju/files/photons2016/2015-PbPb-MC_AllQCDPhoton30-v0-HiForest_PYTHIA_HYDJET_160129.root"
ppData="/d3/scratch/goyeonju/files/photons2016/2015-Data-promptRECO-photonSkims_pp-photonHLTFilter-v0-HiForest.root"
ppMC="/d3/scratch/goyeonju/files/photons2016/2015-PP-MC_Pythia8_Photon30_pp502_TuneCUETP8M1.root"
outputname="pbpb_purity_v2.root"
outputname_pp="pp_purity.root"

./quickPhotonPurity_yj.exe $conf $pbpbData $pbpbMC $outputname "pbpb"
#./quickPhotonPurity_yj.exe $conf $ppData $ppMC $outputname_pp "pp"
