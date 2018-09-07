# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 13:43:31 2018

@author: iw596
"""

import pickle
import RBN_v4 as rbn
import bondingProtocol_v4 as bp
import Reactor_v9 as reac
from collections import Counter
string1 = "H:\\YCCSA\\28-08-2018\\ringMolecules\\rings"
string2 =".pickle"
lengthMols = []
mols = []
reactor = reac.Reactor(1,1,10)
for i in range(0,96):
    compString = string1 + str(i) + string2
    pickle_in = open(compString,"rb")
    mol = pickle.load(pickle_in)
    pickle_in.close()
    for i in range(len(mol)):
        lengthMols.append(len(mol[i]))
        mols.append(mol[i])
        #print("The length of the molecule is: " + str(len(mol[i])) + "\n")
lengthMols.sort()
# Get the length of largest mol
largestMol = lengthMols[len(lengthMols) - 1]
numOccurances = []
count = Counter(lengthMols)
# Go through each length of molecule and store how many times a molecule of that legnth appears
for i in range(2,largestMol + 1):
    numOccurances.append(count[i])    
listSizes = list(range(2,largestMol+1))
print ("The list of sizes is \n" + str(listSizes) + "\n")
print ("The num occurances is: " + str(numOccurances) + "\n")
print ("The length of molecules is \n" + str(lengthMols) + "\n")
# Next calculate the probability of finding each molecuke
probability = [i/len(lengthMols) for i in numOccurances]
print ("The probability of each size occuring is:\n" + str(probability) + "\n")
#mols[0][0].spikeArray[0].printNodeProps()

print(Counter(lengthMols))











##reactor.intensityDebug(mols)
#for i in range(len(mols)):
#    #print ("Initial Intenisty \n")
#    #reactor.intensityDebug(mols[i])
#    if reactor.checkStructureDebug(mols[i]) == True:
#        print ("Structure stable \n")
#    #else:
#        #print ("Structure Unstable\n")
#    for j in range(100):
#        reactor.incrementMoleculeState(mols[i])
#    #print ("New Intensity \n")
#    reactor.calculateIntensitySpikes(mols[i])
#    #reactor.intensityDebug(mols[i])
#    if reactor.checkStructureDebug(mols[i]) == True:
#       print ("Structure still stable \n")
#    else:
#        print ("Structure  now Unstable\n")
#    
#    if bp.checkCircularBonding(mols[i]) == False:
#        print ("Molecule is not stable \n")