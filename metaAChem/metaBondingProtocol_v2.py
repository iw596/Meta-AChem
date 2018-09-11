# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 11:01:16 2018
Need to add function to remove chains of atoms not bonded to meta molecule

This script containins the bonding protocol for metaAtoms, the script consists of a fucntion that takes two metaatoms and attempts
to bond them through the use of their metaspikes. Note that this script is not complete, type 2 metaspike bonding needs to be added
@author: iw596
"""
import type1BondingProtocol_v2 as t1BP
import type2BondingProtocol_v2 as t2BP

# Not happy with this need to make it so random chance of which type of spike will bond, in future want to completely remove meta spikes

def bondMetaAtoms (metaAtom1,metaAtom2,metaMolecule):
    
    origInMol = True # Remains true if metaAtom2 was originally in the molecule
    #If metaAtom2 is not currently in molecule then add it for ease of use
    if metaAtom2 not in metaMolecule.metaAtoms:
        metaMolecule.addMetaAtom(metaAtom2)
        origInMol = False
    for i in range (len(metaAtom1.metaSpikes)):
        for j in range(len(metaAtom2.metaSpikes)):
            # Currently metaspikes need to be same type to bond
            if metaAtom1.metaSpikes[i].typeSpike == 1 and metaAtom2.metaSpikes[j].typeSpike == 1:
                print ("Reaches type 1 bonding algorithm \n")
                t1BP.bondType1MetaSpikes(metaAtom1.metaSpikes[i],metaAtom2.metaSpikes[j],metaAtom1,metaAtom2,metaMolecule)
           
            elif  metaAtom1.metaSpikes[i].typeSpike == 2 and metaAtom2.metaSpikes[j].typeSpike == 2:
                print ("Reaches type 2 bonding algorithm \n")
                t2BP.bondType2MetaSpikes(metaAtom1.metaSpikes[i],metaAtom2.metaSpikes[j],metaAtom1,metaAtom2,metaMolecule)
    
    #Probably should add final bond stability check and check stability of rings
    
    
    
    # If metaAtom2 was not originally in the metaMolecule we need to check that it is now bonded to the metaMolecule
    if origInMol == False:
        bonded = False
        # Go through each dangling node and dangling tail and check to see if it is bonded to an atom in the metamolecule
        for i in range (len(metaAtom2.metaSpikes)):
            if metaAtom2.metaSpikes[i].typeSpike == 1:
                for j in range (len( metaAtom2.metaSpikes[i].danglingNodeList)):
                    if metaAtom2.metaSpikes[i].danglingNodeList[j].bonded == True:
                        bonded = True
            else:
                for j in range (len( metaAtom2.metaSpikes[i].danglingTailList)):
                    if metaAtom2.metaSpikes[i].danglingTailList[j].bonded == True:
                        bonded = True
        
        # If no danglign nodes or tails bonded then remove from molecule
        if bonded == False:
            metaMolecule.removeMetaAtom(metaAtom2) 
               
                
        
