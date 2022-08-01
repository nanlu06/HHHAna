import uproot
import matplotlib.pyplot as plt

def readbranch(file, name):
    tree = uproot.open('root://cms-xrd-global.cern.ch//'+file+':Runs')
    #print(file.keys())
    branch = tree[name].array(library="np")
    weight = branch[0]
    #print(weight)
    return weight


sumOfw = 0
samplename = 'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8' #'TTToHadronic_TuneCP5_13TeV-powheg-pythia8'#'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'
brn = 'genEventSumw'
# open file in read mode
fin = open('../list/'+samplename+'.list', 'r')

# display content of the file
for x in fin.readlines():
    #print(x)
    sumOfw += readbranch(x, brn)
    #print('sumOfw:', sumOfw)
    
# close the file
fin.close()

# open the file in the write mode
fout = open('output/'+samplename+brn+'.txt', 'w')
fout.write("%.8f" % sumOfw)
# close the file
fout.close()





#file = uproot.open('root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL17NanoAODv9/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_mc2017_realistic_v9-v1/130000/5669AB65-83CF-3D45-9190-6D765216ED67.root:Runs')
#file = uproot.open('/eos/user/m/mstamenk/CxAOD31run/hhh-samples/HHH6b_RunIISummer20UL17/RunIISummer20UL17NANOAODSIM_10.root:Runs')
#print(file.keys())
#branch = file['genEventSumw'].array(library="np")
#print(branch[0])
