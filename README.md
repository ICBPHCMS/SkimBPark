Original version from  git@github.com:gkaratha/SkimCode.git

imported readme:
Code use for skimming analysis of 1 lepton + 2 tracks. Not the one will run online

Instructions:
cmsrel CMSSW_10_2_5
cd CMSSW_10_2_5/src
cmsenv
git clone git@github.com:ICBPHCMS/SkimBPark.git ./
scram b -j 8

To run:
Use AOD or RECO as inputs. 
All cuts are configured in HLTAnalysis/TriggerAnalyzer/test/runSkim.py
Macros from George are in HLTAnalysis/TriggerAnalyzer/macros



