# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 11:20:51 2018

@author: iw596
"""

import MetaSpike_v2 as metaSpk
import DanglingNode_v1 as dn
import danglingTail_v1 as dt
from numpy import *
def generateMetaSpikes(molecule):
   # print ("The length of the molecule passed is: " + str(len(molecule)) + "\n")
    # Initally generate a type 1 and a type 2 metaspike
    type1Spike = metaSpk.MetaSpike(1,0)
    type2Spike = metaSpk.MetaSpike(2,1)
   
    fixDanglingNodes(molecule)
    
    # First step is to add dangling nodes to the metaSpikes
    # Go through each atom in molecule
    for i in range(size(molecule)):
        #print ("Reaches here \n")
        # Go throughn each spike in atom
        for j in range(size(molecule[i].spikeArray)):
            #print ("Inside spikearray \n")
            spk = molecule[i].spikeArray[j]
            #print ("The number of dangling bonds is: " + str(molecule[i].spikeArray[j].numDanglingBonds) + "\n")
            if spk.numDanglingBonds == 1:
                #print ("One dangling bond \n")
                # Convert to dangling node
                dangNode = dn.DanglingNode(spk.danglingBonds[0],spk)
                # Add dangling node to type 1 spike
                type1Spike.addDanglingNode(dangNode)
            elif spk.numDanglingBonds > 1:
                #print ("The number of dangling bonds is: " + str(spk.numDanglingBonds) + "\n")
                #print ("The length of the dangling bond array is: " + str(len(spk.danglingBonds)) + "\n")
                # Need below if statement to fix initial error where array was not reset after bonding
                # can remove in future generations of ring molecules
                if spk.bonded == True:
                    numDanNodes = len(spk.nodeList) - len(spk.bondedRBN.spikeArray[spk.bondedSpikeNum].nodeList)
                    spk.addDanglingBonds(numDanNodes)
                    #print ("The number of dangling bonds is: " + str(spk.numDanglingBonds) + "\n")
                    #print ("The length of dangling bond array is: " + str(len(spk.danglingBonds)) + "\n")
                newTail = dt.danglingTail(molecule[i].spikeArray[j])   
                for k in range(spk.numDanglingBonds):
                    # If more than two dangling nodes we need to add these to the type 2 spike
                    # Need to convert each node to dangling node, then append to list, then add list to type 2 spike
                    #print ("The value of k is: " + str(k) + "\n")
                    dangNode = dn.DanglingNode(spk.danglingBonds[k],spk)
                    newTail.addNode(dangNode)
                type2Spike.addTailDanglingBonds(newTail)
    
    # Next step is to calculate the intensity of each spike, note if spike has no nodes then its
    # intensity is not calculated
    
    if size(type1Spike.danglingNodeList) == 0:
       # print ("Reaches here \n")
        type1Spike = -1

    if size(type2Spike.danglingTailList) == 0:
        #print ("Reaches here \n")
        type2Spike = -1
    

    
    return type1Spike,type2Spike
        
    
def fixDanglingNodes(molecule):
    """ Dangling bonds less than zero due to programming error so need to make positive and
        select right number of dangling bonds
    """
    #print (" Function called \n")
    for i in range(size(molecule)):
        for j in range(size(molecule[i].spikeArray)):
            spk = molecule[i].spikeArray[j]
            if spk.numDanglingBonds < 0:
                spk.numDanglingBonds = spk.numDanglingBonds * -1 # Make positive
                # Call function to recalculate danglung nodes
                spk.addDanglingBonds(spk.numDanglingBonds)
                
                
                    
                