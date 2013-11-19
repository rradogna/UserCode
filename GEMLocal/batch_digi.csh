#!/bin/csh
setenv VO_CMS_SW_DIR /opt/exp_soft/cms
source $VO_CMS_SW_DIR/cmsset_default.csh
limit vmem unlim
#mkdir /lustre/cms/store/user/radogna/GEM_reco/
#cd /cmshome/radogna/GEM/CMSSW_5_3_7_patch4/src/VHTauTau/TreeMaker/test/
cd /cmshome/radogna/GEM_Geometry/CMSSW_6_1_2_SLHC6_patch1/src/
setenv SCRAM_ARCH slc5_amd64_gcc472
eval `scramv1 runtime -csh`
cd -
cmsRun /cmshome/radogna/GEM_Geometry/CMSSW_6_1_2_SLHC6_patch1/src/SimMuon/GEMDigitizer/test/runGEMDigiProducer_cfg.py 

