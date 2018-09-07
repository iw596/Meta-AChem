# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 09:16:32 2018

@author: iw596
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 10:37:15 2018
This scipt is used to bond two spikes as it it is important and used often
I have put it in its own individual script to highlight its complexity and 
importance to the overall project
@author: iw596
"""
from numpy import *
import Reactor_v9

def reactRBNS (spike1,spike2,rbn1,rbn2,mol,reactor):
    """ This function is used as a base to call the other functions
        associated with bonding two random boolean networks
    """
    if checkIntensityError(rbn1) == False or checkIntensityError(rbn2) == False:
        return False

    #rbn1.rbnBonded(spike1.returnSpkNum())
    #rbn2.rbnBonded(spike2.returnSpkNum())
    # Change spike bonding status
    spike1.hasBonded(rbn2,spike2.spikeNumber)
  #  print ("Spike 1 test after initial bond func: " + str(spike1.intensity) + "\n")
    spike2.hasBonded(rbn1,spike1.spikeNumber)

 #   print ("Spike 2 test after initial bond func: " + str(spike2.intensity) + "\n")

   # Determine which of the spikes is smallest and store size of smallest
    smallest,numSwaps = findSmallestSpike(spike1,spike2,rbn1,rbn2)
#
#    # Maximum number of swaps is size of smallest spike-1
    numSwaps = numSwaps - 1
    #print ("For spike 1 test intensity is: " + str(currentIntensitySpike1) + "\n")
    #print ("For spike 2 test intensity is: " + str(currentIntensitySpike2) + "\n")
##    
##    
##    print ("-------------------------------- \n")
##    print ("Spike 1 Before Bonding: \n ")
##    spike1.printNodeProps()
##    print ("Spike 2 Before Bonding: \n ")
##    spike2.printNodeProps()
#    print ("-------------------------------- \n")
    swapLinks(spike1,spike2,numSwaps,smallest,rbn1,rbn2)
##    
#
#    print ("RBN 2 number is now: " + str(spike2.RBN.rbnNumber) + "\n")
   # print ("old spk 1 int is: " + str(spike1.intensity) + "\n")
    #print ("old spk 2 int is: " + str(spike2.intensity) + "\n")
    reactor.calculateIntensitySpikes(mol)
    newIntensitySpk1 = spike1.intensity
    newIntensitySpk2 = spike2.intensity
    
    # This used to handle situation in which cyclength could not be found
    if newIntensitySpk1 == 'a' or newIntensitySpk2 == 'a' or checkCircularBonding(mol) == False:
        #print ("Error thrown \n")
        stable = False
        spike1.bondBreak()
        spike2.bondBreak()
        reactor.calculateIntensitySpikes(mol)

        return stable
    
    #print ("new spk 1 int is: " + str(spike1.intensity) + "\n")
    #print ("new spk 2 int is: " + str(spike2.intensity) + "\n")
    stable = 0
##    
   # print ("For spike 1 post bond intensity is: " + str(newIntensitySpk1) + "\n")
    #print ("For spike 2 post bond intensity is: " + str(newIntensitySpk2) + "\n")
    if spike1.type == 3 and spike2.type == 3:
        
        if newIntensitySpk1 + newIntensitySpk2 <= 2 and newIntensitySpk1 + newIntensitySpk2 >= -2:
            #print ("Stable bond \n")
            stable = True
        else: 
            spike1.bondBreak ()
            spike2.bondBreak ()
            reactor.calculateIntensitySpikes(mol)
            #print("TYPE2 Unstable bond \n")
            stable = False 

    elif (spike1.type == 2  and spike2.type == 2) or (spike1.type == 2  and spike2.type == 3) or (spike1.type == 3  and spike2.type == 2):
        
        if newIntensitySpk1 + newIntensitySpk2 <= 1 and newIntensitySpk1 + newIntensitySpk2 >= -1:
            #print ("Stable bond \n")
            stable = True
        else: 
            spike1.bondBreak ()
            spike2.bondBreak ()
            reactor.calculateIntensitySpikes(mol)
            #print("TYPE2 Unstable bond \n")
            stable = False
        
    else:
        
        if newIntensitySpk1 + newIntensitySpk2 == 0:
            #print ("Stable bond \n")
            stable = True
        else: 
            # If unstable break both bonds after doing this we need to recalculate the intensity
            spike1.bondBreak ()
            spike2.bondBreak ()
            reactor.calculateIntensitySpikes(mol)
            #print("Unstable bond \n")
            stable = False
        
    # Check that intensity of spikes is valid
    if checkIntensityError(rbn1) == False or checkIntensityError(rbn2) == False:
        spike1.bondBreak ()
        spike2.bondBreak ()
        reactor.calculateIntensitySpikes(mol)
        print ("Return error thrown \n")
        return False

#    print ("Spike 1 AFTER Bonding: \n ")
#    spike1.printNodeProps()
#    print ("Spike 2 AFTER Bonding: \n ")
#    spike2.printNodeProps()
#    
   # print ("For spike 1 new intensity is: " + str(newIntensitySpk1) + "\n")
   # print ("For spike 2 new intensity is: " + str(newIntensitySpk2) + "\n") 
    
    return stable

def findSmallestSpike (spike1,spike2,rbn1,rbn2):
    """ This function is used to determine which spike is the smalllest
        and returns the smallest spike and its size
    """
    sizeSpike1 = size(spike1.returnNodeArray()) 
    #spike1.printNodeArray()
    sizeSpike2 = size(spike2.returnNodeArray())
    #spike2.printNodeArray()
    #print ("The size of spike 1 is: " + str(sizeSpike1) + "\n")
    #print("The size of spike 2 is: " + str(sizeSpike2) + "\n")
    if (sizeSpike1<= sizeSpike2):
        return 1,sizeSpike1
    else:
        return 2,sizeSpike2

def swapLinks (spike1,spike2,numSwaps,smallestSpike,rbn1,rbn2):
    spike1NodeArray = spike1.returnNodeArray()
    spike2NodeArray = spike2.returnNodeArray()
    # Store the sizes of spikes, need to use original values later but also need temps ones
    refMaxS1Index = size(spike1NodeArray) -1
    refMaxS2Index = size(spike2NodeArray) -1
    maxS1Index = refMaxS1Index
    maxS2Index = refMaxS2Index
    
    if smallestSpike == 1:
        #print ("If statement called \n")
        spike2.addDanglingBonds(size(spike2.nodeList) - size(spike1.nodeList))
        for i in range(numSwaps):
            spike2NodeArray[maxS2Index-1].involvedInBond(spike2NodeArray[maxS2Index],spike1NodeArray[maxS1Index])  
            spike1NodeArray[maxS1Index-1].involvedInBond(spike1NodeArray[maxS1Index],spike2NodeArray[maxS2Index]) 

            maxS1Index = maxS1Index-1
            maxS2Index = maxS2Index-1
    
    else:
        #print ("Else statement called \n")
        spike1.addDanglingBonds(size(spike1.nodeList) - size(spike2.nodeList))
        for i in range(numSwaps):
            
            spike1NodeArray[maxS1Index-1].involvedInBond(spike1NodeArray[maxS1Index],spike2NodeArray[maxS2Index])
            spike2NodeArray[maxS2Index-1].involvedInBond(spike2NodeArray[maxS2Index],spike1NodeArray[maxS1Index]) 
            
            maxS1Index = maxS1Index-1
            maxS2Index = maxS2Index-1
        
    spike1.changeNodeArray(spike1NodeArray)
    spike2.changeNodeArray(spike2NodeArray)
    
#    print ("Node info for spike 1 \n")
#    for i in range(size(spike1.nodeList)):
#        spike1.nodeList[i].printProps()
#    
#    
#    print ("Node info for spike 2 \n")
#    for i in range(size(spike2.nodeList)):
#        spike2.nodeList[i].printProps()
#        
       
    

def reanalyseRBN (RBN):
    """ This function is used to re calculate the intensity of all the spikes in a RBN, it is called after a
        RBN has bonded or unbonded as the topology of the new work has changed which could have lead to
        a change in the intensity of the spikes
    """
    
    for i in range(size(RBN.spikeArray)):
     #print ("The RBN is: " + str(RBN.rbnNumber) + "\n")
     #print ("The spike number is: " + str(RBN.spikeArray[i].spikeNumber) + "\n")
     #print ("The old Intensity is: "  + str(RBN.spikeArray[i].intensity) + "\n")
     if RBN.spikeArray[i].checked ==  False:
        RBN.spikeArray[i].checked = True
        if RBN.spikeArray[i].bonded == True:
            reanalyseRBN(RBN.spikeArray[i].bondedRBN)
        
       # print ("The RBN is: " + str(RBN.rbnNumber) + "\n")
        #print ("The spike number is: " + str(RBN.spikeArray[i].spikeNumber) + "\n")  
        #print ("Being recalculated \n")    
        RBN.spikeArray[i].recalculateIntensity()
     
    # print ("The RBN is: " + str(RBN.rbnNumber) + "\n")
    # print ("The spike number is: " + str(RBN.spikeArray[i].spikeNumber) + "\n")   
    # print ("The new Intensity is: "  + str(RBN.spikeArray[i].intensity) + "\n")   
        



def checkCircularBonding (molecule):
    """ This function takes an array of RBNs which are bonded to each other  in a linear fashion, it iterates through the RBNs and checks
        that the bond is still stable, if it is not false is returned otherwise true is returned
    """
    stable = True
    
    for i in range(size(molecule)):
        for j in range(size(molecule[i].spikeArray)):
            if molecule[i].spikeArray[j].bonded == True:
                
                spike1 = molecule[i].spikeArray[j]
                spike2 = molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum]
                
                if (spike1.type == 2  and spike2.type == 2) or (spike1.type == 2  and spike2.type == 3) or (spike1.type == 3  and spike2.type == 2):
                    atomIntensity = spike1.intensity
                    bondedAtomIntensity = spike2.intensity
                    #Check if both intensitys are valid
                    if atomIntensity == 'a' or bondedAtomIntensity == 'a':
                        stable = False
                     #   print ("Strucure unstable \n")
                        return stable
                    
                    if abs(atomIntensity + bondedAtomIntensity) > 1:
#                        print  ("Spike 1 type is: " + str(molecule[i].spikeArray[j].type) + "\n")
#                        print ("Spike 1 intensity is: " + str(molecule[i].spikeArray[j].intensity) + "\n")
#                        print ("Spike 2 type is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].type) + "\n")
#                        print ("Spike2 intensity is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity) + "\n")
#                        print ("Absolute intensity is: " + str( abs(molecule[i].spikeArray[j].intensity + molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity)) + "\n")
#                        print ("Bond unstable \n")
                        stable = False
                   #     print ("Strucure unstable \n")
                        return stable
                
                elif spike1.type == 3 and spike2.type == 3:
                    atomIntensity = spike1.intensity
                    bondedAtomIntensity = spike2.intensity
                    #Check if both intensitys are valid
                    if atomIntensity == 'a' or bondedAtomIntensity == 'a':
                        stable = False
                      #  print ("Strucure unstable \n")
                        return stable
                    
                    if abs(atomIntensity + bondedAtomIntensity) > 1:
#                        print  ("Spike 1 type is: " + str(molecule[i].spikeArray[j].type) + "\n")
#                        print ("Spike 1 intensity is: " + str(molecule[i].spikeArray[j].intensity) + "\n")
#                        print ("Spike 2 type is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].type) + "\n")
#                        print ("Spike2 intensity is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity) + "\n")
#                        print ("Absolute intensity is: " + str( abs(molecule[i].spikeArray[j].intensity + molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity)) + "\n")
#                        print ("Bond unstable \n")
                        stable = False
                     #   print ("Strucure unstable \n")
                        return stable
                
                
                
                else:
                    atomIntensity = spike1.intensity
                    bondedAtomIntensity = spike2.intensity
                    if atomIntensity == 'a' or bondedAtomIntensity == 'a':
                        stable = False
                     #   print ("Strucure unstable \n")
                        return stable
                    if atomIntensity + bondedAtomIntensity != 0:
#                        print ("Bond unstable \n")
                        stable = False
                        #print ("Strucure unstable \n")
                        return stable
                    
                
                
    return stable


def checkCircularBondingDebug (molecule,reactor):
    """ This function takes an array of RBNs which are bonded to each other  in a linear fashion, it iterates through the RBNs and checks
        that the bond is still stable, if it is not false is returned otherwise true is returned. Note this is a debugging
        version of the function and prints off information about the bonds 
    """
    stable = True
    
    for i in range(size(molecule)):
        for j in range(size(molecule[i].spikeArray)):
            if molecule[i].spikeArray[j].bonded == True:
                
                spike1 = molecule[i].spikeArray[j]
                spike2 = molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum]
                
                if (spike1.type == 2  and spike2.type == 2) or (spike1.type == 2  and spike2.type == 3) or (spike1.type == 3  and spike2.type == 2):
                    atomIntensity = spike1.intensity
                    bondedAtomIntensity = spike2.intensity
                    #Check if both intensitys are valid
                    if atomIntensity == 'a' or bondedAtomIntensity == 'a':
                        stable = False
                        return stable
                    
                    if abs(atomIntensity + bondedAtomIntensity) > 1:
                        print  ("Spike 1 type is: " + str(molecule[i].spikeArray[j].type) + "\n")
                        print ("Spike 1 intensity is: " + str(molecule[i].spikeArray[j].intensity) + "\n")
                        print ("Spike 2 type is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].type) + "\n")
                        print ("Spike2 intensity is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity) + "\n")
                        print ("Absolute intensity is: " + str( abs(molecule[i].spikeArray[j].intensity + molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity)) + "\n")
                        print ("Bond unstable \n")
                        stable = False
                        return stable
                
                elif spike1.type == 3 and spike2.type == 3:
                    atomIntensity = spike1.intensity
                    bondedAtomIntensity = spike2.intensity
                    #Check if both intensitys are valid
                    if atomIntensity == 'a' or bondedAtomIntensity == 'a':
                        stable = False
                        return stable
                    
                    if abs(atomIntensity + bondedAtomIntensity) > 1:
                        print  ("Spike 1 type is: " + str(molecule[i].spikeArray[j].type) + "\n")
                        print ("Spike 1 intensity is: " + str(molecule[i].spikeArray[j].intensity) + "\n")
                        print ("Spike 2 type is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].type) + "\n")
                        print ("Spike2 intensity is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity) + "\n")
                        print ("Absolute intensity is: " + str( abs(molecule[i].spikeArray[j].intensity + molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity)) + "\n")
                        print ("Bond unstable \n")
                        stable = False
                        return stable
                
                
                
                else:
                    atomIntensity = spike1.intensity
                    bondedAtomIntensity = spike2.intensity
                    if atomIntensity == 'a' or bondedAtomIntensity == 'a':
                        stable = False
                        return stable
                    if atomIntensity + bondedAtomIntensity != 0:
#                        print ("Bond unstable \n")
                        stable = False
                        return stable
                    
                
                
    return stable


def checkIntensityError (atom):
    """ This function checks every spike in a atom to see if its intensity is valid if it is not valid
        then false is returned
    """
    
    for i in range(size(atom.spikeArray)):
        if atom.spikeArray[i].intensity == 'a':
            return False
    return True