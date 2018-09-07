# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 10:48:47 2018

 NOTE THIS SCRIPT IS NOT CURRENTLY FINISHED
This script is used to bond type 2 metaspikes, type 2 spikes consist of dangling tails, for dangling taisl tob ond the sum of their spike
intensitys must be less than or equal to a number which is dependent on the type of spike
@author: iw596
"""
import math
import random
import numpy as np
import Reactor_v9 as reactor
from bondingProtocol_v5 import checkCircularBonding

def bondType2MetaSpikes (metaSpike1,metaSpike2):
    # Find smallest metaspike
    smallest = smallestMetaSpike(metaSpike1,metaSpike2)
    origInMol = True # Remains true if metaAtom2 was originally in the molecule
    #If metaAtom2 is not currently in molecule then add it for ease of use
    if metaAtom2 not in metaMolecule.metaAtoms:
        metaMolecule.addMetaAtom(metaAtom2)
        origInMol = False
    
    # Generate the set of dangling tails 
    setTails1 = list(range(len(metaSpike1.danglingTailList)))
    setTails2 = list(range(len(metaSpike1.danglingTailList)))
    # Shuffle both sets
    random.shuffle(setTails1)
    random.shuffle(setTails2)
    
    if smallest == 1:
        while len(setTails1) != 0:
            #Select a tail from each set
            tail1 = setTails1.pop()
            tail2 = setTails2.pop()
            # Attemp to react the two tails
            bondDanglingTails (tail1,tail2)
            
            # Shuffle both sets
            random.shuffle(setTails1)
            random.shuffle(setTails2)
            
    else:
        while len(setTails2) != 0 :
            #Select a tail from each set
            tail1 = setTails1.pop()
            tail2 = setTails2.pop()
            # Attemp to react the two tails
            bondDanglingTails (tail1,tail2)
            
            # Shuffle both sets
            random.shuffle(setTails1)
            random.shuffle(setTails2)
            
    
    
    return 1

def bondDanglingTails (tail1,tail2):
    """ This function attempts to bond two danglingTails, to do this the spike intensityies in the tails
        must fufill a critieria which is dependent on the type of spike the tail is part of, if they can
        bond the links and swapped and the spike intensities recalculate to check that the criteria
        is still met, if it is then True is returned, if it is no longer met then the connection is broken,
        the intensity recalcualted and False is returned. If they could not bond in the first place
        False is returned.
    """
    # Check that intensity of both spikes is valid
    if tail1.spike.intensity == 'a':
        return False
    
    if tail2.spike.intensity == 'a':
        return False
    
    
    if tail1.spike.type == 3 and tail2.spike.type == 3:
        # If criteria met then bond
        if abs(tail1.spike.intensity + tail2.spike.intensity) <= 2:
            swapLinks (tail1,tail2)

            # After bonding need to recalculate intensity and check stability
            
    
    elif (tail1.spike.type == 2 and tail2.spike.type == 2) or (tail1.spike.type == 2 and tail2.spike.type == 3) or (tail1.spike.type == 3 and tail2.spike.type == 2):
        # If criteria met then bond
        if abs(tail1.spike.intensity + tail2.spike.intensity) <= 1:
            swapLinks (tail1,tail2)
    
    else:
        # If criteria met then bond
        if tail1.spike.intensity + tail2.spike.intensity == 0:
            swapLinks (tail1,tail2)
        
    return False

def swapLinks (tail1,tail2):
    """ This function is used to swap the connections between the dangling tails, the number of swaps will
        be n where n is the size of the smallest tail, the links the will be swapped starting from
        the highest index in the array of dangling tails. The dangling node at each index will be bonded to the dangling node
        in the other tail at the same index
    """
    # Find smallest tail
    smallest = smallestDanglingTail(tail1,tail2)
    
    # Store the sizes of dangling tails, need to use original values later but also need temps ones
    refMaxT1Index = len(tail1.nodeList) -1
    refMaxT2Index = len(tail2.nodeList) -1
    maxT1Index = refMaxT1Index
    maxT2Index = refMaxT2Index
    
    if smallest == 1:
        numSwaps = refMaxT1Index # stores how many swaps will take place
        for i in range(numSwaps):
            # Connect dangling nodes
            tail1.nodeList[maxT1Index].bondFormed(tail2.nodeList[maxT2Index],tail2.spike)
            tail2.nodeList[maxT2Index].bondFormed(tail1.nodeList[maxT1Index],tail1.spike)
            # Move to next index
            maxT1Index -= 1
            maxT2Index -= 1
        # Update tails with bonding status
        tail1.bondFormed(tail2,numSwaps)
        tail2.bondFormed(tail1,numSwaps)

        
    else:
        numSwaps = refMaxT2Index # stores how many swaps will take place
        for i in range(numSwaps):
            # Connect dangling nodes
            tail1.nodeList[maxT1Index].bondFormed(tail2.nodeList[maxT2Index],tail2.spike)
            tail2.nodeList[maxT2Index].bondFormed(tail1.nodeList[maxT1Index],tail1.spike)
            # Move to the next index
            maxT1Index -= 1
            maxT2Index -= 1
        # Update tails with bonding status
        tail1.bondFormed(tail2,numSwaps)
        tail2.bondFormed(tail1,numSwaps)
    
        
        
def smallestMetaSpike (metaSpike1,metaSpike2):
    """ This function returns 1 if metaspike2 has the fewest danglinng tails or returns 2 if metaspike 2
        has the fewest dangling tails
    """

    if len(metaSpike1.danglingTailList) <= len(metaSpike2.danglingTailList):
        return 1
    else:
        return 2

def smallestDanglingTail(danglingTail1,danglingTail2):
    """ This function returns 1 if danglingTail1 has less (or equall) dangling nodes than danglingTail2. The
        function returns 2 if the opposite is true
    """
    
    if len(danglingTail1.nodeList) <= len(danglingTail2.nodeList):
        return 1
    else:
        return 2


def checkStability ():
    """ Need to complete this function, used to check meta molecule until it is stable """
    return 1
    