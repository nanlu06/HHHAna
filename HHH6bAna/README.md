To set up this code

cmsrel CMSSW_10_6_26

cd CMSSW_10_6_26/src


git clone git@github.com:nanlu06/HHHAna.git

cd HHHAna/HHH6bAna

make

To run the code

./analyzeHHH test.list ana_hhh_test_2017.root hhh6b_test F T 2017 > out
