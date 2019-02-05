Original version from  git@github.com:gkaratha/SkimCode.git

imported readme:
Code use for skimming analysis of 1 lepton + 2 tracks. Not the one will run online
Instructions:
cmsrel CMSSW_10_2_1
cd CMSSW_10_2_1/src
cmsenv
git clone https://github.com/gkaratha/SkimCode
scram b -j 8

To run:
Use AOD or RECO as inputs. All cuts are configured in python/ConfFile_cfg.py. To run simply cmsRun ConfFile_cfg.py.



=> Use in CMSSW_10_2_5 for new MC production and data