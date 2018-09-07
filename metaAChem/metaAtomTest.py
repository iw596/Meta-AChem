# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 11:51:22 2018
Test script to try and find errors 
@author: iw596
"""
import pickle
import RBN_v4 as rbn
from numpy import *
import MetaAtom_v2 as metAt
import metaSpikeGenerator_v1 as metSpkGen
import metaBondingProtocol_v1 as mBp


def atomGenerator (molecule,metaAtomNum):
    metSpike1,metSpike2 = metSpkGen.generateMetaSpikes(molecule)
    print ("For atom: " + str(metaAtomNum) + "\n")
    if metSpike1 != -1:
        print ("The intensity of spike 1 is: " + str(metSpike1.intensity) + "\n")
    else:
        print ("No type 1 spike \n")

    if metSpike2 != -1:
        print ("The intensity of spike 2 is: " + str(metSpike2.intensity) + "\n")
    else:
        print ("No type 2 spike \n")
        
    metaAtom = metAt.MetaAtom(metaAtomNum)
    if metSpike1 != -1:
        metaAtom.addMetaSpike(metSpike1)
    if metSpike2 != -1:
        metaAtom.addMetaSpike(metSpike2)
    metaAtom.calculateState()
    return metaAtom
    

    
    
    
    
    
 # First need to import a ring to use
pickle_in = open("C:\\Users\\isaac\\Documents\\Summer Project\\17-08-2018\\new\\ringMolecules\\rings0.pickle","rb")
#Select first ring
mol1 = pickle.load(pickle_in)[4]

pickle_in.close()


 # First need to import a ring to use
pickle_in = open("C:\\Users\\isaac\\Documents\\Summer Project\\17-08-2018\\new\\ringMolecules\\rings0.pickle","rb")
#Select first ring
mol2 = pickle.load(pickle_in)[4]
pickle_in.close()

metaAtom1 = atomGenerator(mol1,0)
#metaAtom2 = atomGenerator(mol2,1)
print ("The length of the molecule is: " + str(len(mol1)) + "\n")
#mBp.bondMetaSpikes(metaAtom1.metaSpikes[0],metaAtom2.metaSpikes[0],metaAtom1,metaAtom2)

#mBp.bondMetaSpikes(metaAtom1.metaSpikes[1],metaAtom2.metaSpikes[1],metaAtom1,metaAtom2)
