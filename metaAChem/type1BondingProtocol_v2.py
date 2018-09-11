# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 19:07:08 2018
This script contains the implementation of type 1 metaspike bonding, for two type 1 metaspikes to bond
there needs to be at least one connection between dangling nodes in both metaspikes. For two dangling nodes
to combine the sum of the intensity of their spikes (not metaspikes) has to equal 1 if both spikes are
type 2 or equal 0 if one of the spikes is type 1. For the bond to be stable the intensity of the spikes
after bonding still has to meet this condition. When bonding metaspikes the dangling nodes are selected 
at random to be bonded. This is done until all the dangling nodes in the smallest metaspike have had a chance
at bonding.
@author: isaac
"""
import math
import random
import numpy as np
import Reactor_v9 as reactor
from bondingProtocol_v5 import checkCircularBonding

def bondType1MetaSpikes (metaSpike1,metaSpike2,metaAtom1,metaAtom2,metaMolecule):
    """ This function bonds two type 1 metaspikes by bonding their dangling nodes, it then
        checks the stability of the bonds and breaks them if they are unstable
    """
    print ("The ring size metaAtom1 is: " + str(len(metaAtom1.mol)) + "\n")
    print ("The ring size metaAtom2 is: " + str(len(metaAtom2.mol)) + "\n")
    # Find which is the smallest metaspike
    smallest = smallestType1MetaSpike(metaSpike1,metaSpike2)

    # Need to store set of dangling node numbers for both metaspikes 
    set1 = list(range(len(metaSpike1.danglingNodeList)))
    set2 = list(range(len(metaSpike2.danglingNodeList)))
    
    numBonds  = 0 # Stores how many stable bonds have formed
    bondFormed = False
    
    # Randomly shuffle both lists
    random.shuffle(set1)
    random.shuffle(set2)
    # If smallest spike is metaspike1 
    if smallest == 1:
        # Attempt to bond dangling nodes until smallest set is empty
        while len(set1) != 0:
            # Select a dangling node by popping an integer from the set and use this number to select a dangling node from list
            spk1Node = metaSpike1.danglingNodeList[set1.pop()]
            spk2Node = metaSpike2.danglingNodeList[set2.pop()]
            # Attempt to bond these two nodes
            if bondDanglingNodes(spk1Node,spk2Node,metaSpike1,metaSpike2,metaMolecule) == True:
                # Next need to check that rings which atoms are part of are still stable, to do this we combine the molecules
                # into one molecule and recalculate the intensitys and then check to see if the structure is still stable
                # and whether the bonding critieria is still met
                numBonds += 1
                bondFormed = True
                print ("The rbn number of the first dangling node is: " + str(spk1Node.spike.RBN.rbnNumber) + "\n")
                print ("The rbn number of the second dangling node is: " + str(spk2Node.spike.RBN.rbnNumber) + "\n")
                # If the resulting dangling node bond is stable when then check the underlying sturcure of the molecule to see
                # if that is still stable. If this returns false then the structure of the molecule has changed so need to check
                # that the newly formed dangling node bonds are stable
                isMetaMoleculeStable(metaSpike1,metaSpike2,metaMolecule,1,numBonds)
                    
            # Reshuffle Lists
            random.shuffle(set1)
            random.shuffle(set2)
            
    else:
        while len(set2) != 0:
            # Pop a node from both sets, this also removes node from set
            spk1Node = metaSpike1.danglingNodeList[set1.pop()]
            spk2Node = metaSpike2.danglingNodeList[set2.pop()]
            # Attempt to bond these two nodes
            if bondDanglingNodes(spk1Node,spk2Node,metaSpike1,metaSpike2,metaMolecule) == True:
                # Next need to check that rings which atoms are part of are still stable
                numBonds += 1
                bondFormed = True
                print ("The rbn number of the first dangling node is: " + str(spk1Node.spike.RBN.rbnNumber) + "\n")
                print ("The rbn number of the second dangling node is: " + str(spk2Node.spike.RBN.rbnNumber) + "\n")
                
                isMetaMoleculeStable(metaSpike1,metaSpike2,metaMolecule,2,numBonds)
            # Reshuffle Lists
            random.shuffle(set1)
            random.shuffle(set2)
    
    for i in range(len(metaMolecule.metaAtoms)):
        metaMolecule.metaAtoms[i].removeDoubleUnbondedAtoms()
    
    
    # If a bond has been formed we need to double check stability as bonds lead to atom topology change which can 
    # affect stability        


def bondDanglingNodes (danglingNode1,danglingNode2,metaSpike1,metaSpike2,metaMolecule):
    """This function bonds two dangling nodes, for two dangling nodes to bond there spike intensity
        must either equal 1 if both spike types are type 2 or 0  if one spike is type 1. For the bond to be stable this 
        condition must still be met after bonding
    """
    print ("Dangling node 1 intensity: " + str(danglingNode1.spike.intensity) + "\n")
    print ("Dangling node 2 intensity: " + str(danglingNode2.spike.intensity) + "\n")
    print ("Dangling node 1 type: " + str(danglingNode1.spike.type) + "\n")
    print ("Dangling node 2 type: " + str(danglingNode2.spike.type) + "\n")
    print ("Dangling node 1 length of spike: " + str(len(danglingNode1.spike.nodeList)) + "\n")
    print ("Dangling node 2 length of spike: " + str(len(danglingNode2.spike.nodeList)) + "\n")
    # First check that dangling nodes are not already bonded
    if danglingNode1.bonded == True or danglingNode2.bonded == True:
        return False
        
   
    # If type 3 spike bonding criteria is most lenient
    if danglingNode1.spike.type == 3 and danglingNode2.spike.type == 3:
        # Check absolute value of sum is less than or equal to one
        if abs(danglingNode1.spike.intensity + danglingNode2.spike.intensity) <= 2:
            # True so bond both nodes
           danglingNode1.bondFormed(danglingNode2)
           danglingNode2.bondFormed(danglingNode1)
           # After bonding nodes we then need to check stability of the bond, we do this by recalculating the intensity of
           # the spikes and then seeing if the criteria is still met
           metaMolecule.calculateIntensity()
           # If bond is stable return true
           if abs(danglingNode1.spike.intensity + danglingNode2.spike.intensity) <= 2:
               # If this criteria is met we finally need to ensure that the ring the atom is part off is still stable                              
               print ("Type 3 Bond stable \n")
               return True
           else:
               print ("Bond Unstable \n")
               return False
        
        else:
            # Criteria not met so return False
            return False
 
    
    # If both spikes are type 2 bonding crtieria is more lenient
    elif (danglingNode1.spike.type == 2  and danglingNode2.spike.type == 2) or (danglingNode1.spike.type == 3  and danglingNode2.spike.type == 2) or (danglingNode1.spike.type == 2  and danglingNode2.spike.type == 3):
        # Check absolute value of sum is less than or equal to one
        if abs(danglingNode1.spike.intensity + danglingNode2.spike.intensity) <= 1:
            # True so bond both nodes
           danglingNode1.bondFormed(danglingNode2,danglingNode2.spike)
           danglingNode2.bondFormed(danglingNode1,danglingNode1.spike)
           metaMolecule.calculateIntensity()
           if abs(danglingNode1.spike.intensity + danglingNode2.spike.intensity) <= 1:
               print ("Type 2 Bond stable \n")
               return True
           else:
               print ("Bond Unstable \n")
               return False
           
           return True
           # If not stable break bonds and recalculate intensity
        else:

            return False

    
    else:
        # If there is a type 1 spike intensity has to equal zero
        if abs(danglingNode1.spike.intensity + danglingNode2.spike.intensity) == 0:
            # True so bond both nodes
           danglingNode1.bondFormed(danglingNode2,danglingNode2.spike)
           danglingNode2.bondFormed(danglingNode1,danglingNode1.spike)
           metaMolecule.calculateIntensity()
           if abs(danglingNode1.spike.intensity + danglingNode2.spike.intensity) == 0:   
               print ("Type 1 Bond stable \n")
                
               return True
           else:
               print ("Bond unstable \n")
               return False
           
           # If not stable break bonds and recalculate intensity
        else:
#            danglingNode1.bondBroken()
#            danglingNode2.bondBroken()
#            danglingNode1.spike.recalculateIntensity()
#            danglingNode2.spike.recalculateIntensity()
            return False
    
    

def smallestType1MetaSpike (metaSpike1,metaSpike2):
    """ This function returns 1 if metaSpike1 has smallest number of dangling nodes or returns 2 
        if metaspike2 has smallst number of dangling nodes
    """
    if len(metaSpike1.danglingNodeList) <= len(metaSpike2.danglingNodeList):
        return 1
    else:
        return 2    



# Can remove numBonds not needed
def isMetaMoleculeStable (metaSpike1,metaSpike2,metaMolecule,smallest,numBonds):
    """ This function is used to check the internal structure and the dangling node bonds of the metaMolecule, it fisrt
        checks the internal structure and if any bonds have broken then it checks the stability of the dangling node bonds
        this process is repeated until the metaMolecule is stable and no more bonds break
    """
    stable  = True
    stable = metaMolecule.checkMolecularStructure() # Check structure of rings
    metaMolecule.calculateIntensity() # Recaclualte intensity
    # This will keep altering the structure until it becomes stable
    if stable == False:
        # Check dangling tails and spikes
        metaMolecule.checkDanglingNodeStability ()
        metaMolecule.checkDanglingTailStability()
        # Check structure of rings in metaMolecule
        metaMolecule.removedUnbondedAtoms()
        #Recalculate intensity
        metaMolecule.calculateIntensity()
      
        isMetaMoleculeStable (metaSpike1,metaSpike2,metaMolecule,smallest,numBonds)
            