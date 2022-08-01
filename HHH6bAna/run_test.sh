#skim
./analyzeHHH testdata.list testdata.root data T F 2017 "((nFatJet>=1 && FatJet_particleNetMD_Xbb[0]/(FatJet_particleNetMD_Xbb[0]+FatJet_particleNetMD_QCD[0]) > 0.5 && nJet>=4) ||  (nJet>=6) ||  (nFatJet>=2 && FatJet_particleNetMD_Xbb[0]/(FatJet_particleNetMD_Xbb[0]+FatJet_particleNetMD_QCD[0]) > 0.5 && FatJet_particleNetMD_Xbb[1]/(FatJet_particleNetMD_Xbb[1]+FatJet_particleNetMD_QCD[1]) > 0.5 && nJet>=2) || (nFatJet>=3 && FatJet_particleNetMD_Xbb[0]/(FatJet_particleNetMD_Xbb[0]+FatJet_particleNetMD_QCD[0]) > 0.5 && FatJet_particleNetMD_Xbb[1]/(FatJet_particleNetMD_Xbb[1]+FatJet_particleNetMD_QCD[1]) > 0.5 && FatJet_particleNetMD_Xbb[2]/(FatJet_particleNetMD_Xbb[2]+FatJet_particleNetMD_QCD[2]) > 0.5))"

./analyzeHHH test.list test.root hhh F F 2017 "((nFatJet>=1 && FatJet_particleNetMD_Xbb[0] > 0.5 && nJet>=4) ||  (nJet>=6) ||  (nFatJet>=2 && FatJet_particleNetMD_Xbb[0] > 0.5 && FatJet_particleNetMD_Xbb[1] > 0.5 && nJet>=2) || (nFatJet>=3 && FatJet_particleNetMD_Xbb[0] > 0.5 && FatJet_particleNetMD_Xbb[1] > 0.5 && FatJet_particleNetMD_Xbb[2] > 0.5))"
#run normal
./analyzeHHH list/testtt.list ana_ttsl_test_2017.root hhh6b_test F F 2017
./analyzeHHH testdata.list testdata.root data T F 2017 > out_data
./analyzeHHH HHH6b_RunIISummer20UL17NANOAODSIM.list ana_hhh_2017.root hhh6b F T 2017 > out
./analyzeHHH test.list ana_hhh_test_2017.root hhh6b_test F T 2017 > out
