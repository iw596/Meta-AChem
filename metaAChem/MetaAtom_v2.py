# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 14:00:17 2018
This class desribes a meta node, a meta node is made up of bonded RBNs and is
used to form a higher level tree structure with other meta nodes in order to 
form a meta network
@author: iw596
"""
from numpy import *
import MetaSpike_v2 as metSpk
 

class MetaAtom:
    """ Meta node class has data and functions which allow it to be generated, store 
        connectoons and states, very similar to a normal RBN node
    """
    
    def __init__(self,atomNumber,mol):
        self.metaSpikes = []
        self.bonded  = False
        self.state = 0 # The state of the metaatom is the sum of the states of the node, with 1 adding plus 1 and zero state subtracting 1
        self.atomNumber = atomNumber
        self.mol = mol # Stores the molecule the metaatom is made out of
        self.metaMolecule = 0 # If the metaAtom is bonded then it is part of a metaMolecule

    
    def addMetaSpike(self,metaSpike):
        
        metaSpike.addMetaAtom (self)
        self.metaSpikes.append(metaSpike)
        #self.updateMetaAtom()
    
    def addMetaMolecule (self,metaMolecule):
        """ This function adds the metaMolecule the metaAtom is in """
        self.metaMolecule = metaMolecule
    
    
    def calculateState (self):
        """ This function is used to calculate the state of the metaatom, it does this by first resetting
            the state of the metaatom to zero and then going through each dangling node
            and if the dangling node has a state of one the state of the metaatom is incremented and if it is zero
            then the state is decremented by one
        """
        newState = 0
       # print ("Inside state function the states DNs are: \n")
       # print ("Before starting \n")
        self.stateDanglingNodes()
        #for i in range(len(self.metaSpikes)):
        #    if self.metaSpikes[i].typeSpike == 1:
         #       print ("Meta atom number is: " + str(self.atomNumber) + "\n")
        
        insideMetState = []
        # To calculate the state we need to update every atom the metaatom consistrs off then see
        # the states of every dangling node in the metaspikes
        for i in range(len(self.metaSpikes)):
            if self.metaSpikes[i].typeSpike == 1:
                #print ("Inside type 1 \n")
                #print ("Number of type 1 nodes: " + str(len(self.metaSpikes[i].danglingNodeList)) + "\n")
                for j in range(len(self.metaSpikes[i].danglingNodeList)):
                    insideMetState.append(self.metaSpikes[i].danglingNodeList[j].state)
                    if self.metaSpikes[i].danglingNodeList[j].state == 1:
          #              print ("Adding one \n" )
                        newState += 1
                    else:
            #            print ("Subracting one \n")
                        newState -= 1
            else:
                
               # print ("Inside type 2 \n")
              #  print ("Number od type 1 tales: " + str(len(self.metaSpikes[i].danglingTailList)) + "\n")
                for j in range(len(self.metaSpikes[i].danglingTailList)):
                    #print ("Size of tail: " + str(len(self.metaSpikes[i].danglingTailList[j].nodeList)) + "\n")
                    for k in range(len(self.metaSpikes[i].danglingTailList[j].nodeList)):
                        insideMetState.append(self.metaSpikes[i].danglingTailList[j].nodeList[k].state)
                        if self.metaSpikes[i].danglingTailList[j].nodeList[k].state == 1:
                            newState += 1
                        else:
                            newState -= 1 
       
       # print ("The state of analysed nodes: \n" + str(insideMetState) + "\n")
      #  print ("The length of analysed nodes: \n" + str(len(insideMetState)) + "\n")
      #  print ("The new state is: " + str(newState) + "\n")                    
        self.state = newState
    
    def updateMetaAtom (self):
        """ To update a metaAtom every node in the metaspike must be updated simultaneouslys, to do this
            we generate a list of nodes using the molecule stored in the metaAtom, next we update each node
            simultaneously before passing the new values back to each node in the spike
        """
       # print ("Old state DNS: \n")
       # self.stateDanglingNodes()
        synchList = []
        synchListState = []
        for i in range(len(self.mol)):
            for j in range(len(self.mol[i].nodeArray)):
                synchList.append(self.mol[i].nodeArray[j])
                synchListState.append(synchList[i].state)
        #print ("The original state is: \n" + str(synchListState) + "\n")
        # Find new state for every node
        newStates = []
        for i in range(len(synchList)):
            oldState = synchList[i].state
            synchList[i].calculateState()
            newStates.append(synchList[i].state)
            synchList[i].state = oldState
        
        for i in range(len(synchList)):
            synchList[i].state = newStates[i]
            synchListState[i] = synchList[i].state
        
        offSet = 0 
        for i in range(len(self.mol)):
            for j in range(len(self.mol[i].nodeArray)):
                self.mol[i].nodeArray[j].state = synchListState[offSet]
                offSet += 1
        stateMol = []
        
        for i in range(len(self.mol)):
            for j in range(len(self.mol[i].nodeArray)):
               stateMol.append(self.mol[i].nodeArray[j].state)
        
     #   print ("The new state is: \n" + str(synchListState) + "\n")
      #  print ("The state of the mol array is: " + str(stateMol) + "\n")
        #print ("Post update \n")
        self.stateDanglingNodes()
        offSet = 0 
        oldStateNodes = [] # Store the old state of nodes in molecule
        newStateNodes = [] # Stores the new state
        # The code below goes through each metaspike and ensures that the dangling nodes have been updated with the correct
        # new state
        for i in range(len(self.metaSpikes)):
            if self.metaSpikes[i].typeSpike == 1:
                #print ("Inside type 1 \n")
                #print ("The number of DNs is: " + str(len(self.metaSpikes[i].danglingNodeList)) + "\n")
                for j in range(len(self.metaSpikes[i].danglingNodeList)):
                    # Find the location of the dangling node in the synch list and change the dangling nodes state to match
                    # state locted in the synch list
                    if self.metaSpikes[i].danglingNodeList[j].node in synchList:
                        oldStateNodes.append(self.metaSpikes[i].danglingNodeList[j].node.state)
                        indexNode = synchList.index(self.metaSpikes[i].danglingNodeList[j].node)
                      #  print ("The current value is: " + str(self.metaSpikes[i].danglingNodeList[j].node.state) + "\n")
                  #      print ("The index of node is: " + str(indexNode) + "\n")
                  #      print ("The new value should be: " + str(synchListState[indexNode]) + "\n")
                        self.metaSpikes[i].danglingNodeList[j].changeState(synchListState[indexNode]) 
                        newStateNodes.append(self.metaSpikes[i].danglingNodeList[j].node.state)
       #                 print ("Node in list \n")
            else:
               # print ("The number of DTs is: " + str(len(self.metaSpikes[i].danglingTailList)) + "\n")
                #print ("Inside type 2 \n")
                # With dangling tails we need an extra for loop to iterate across each nodelist of the tail
                for j in range(len(self.metaSpikes[i].danglingTailList)):
                    for k in range(len(self.metaSpikes[i].danglingTailList[j].nodeList)):
                            if self.metaSpikes[i].danglingTailList[j].nodeList[k].node in synchList:
                                oldStateNodes.append(self.metaSpikes[i].danglingTailList[j].nodeList[k].state)
                                indexNode = synchList.index(self.metaSpikes[i].danglingTailList[j].nodeList[k].node)
                                self.metaSpikes[i].danglingTailList[j].nodeList[k].changeState(synchListState[indexNode]) 
                                newStateNodes.append(self.metaSpikes[i].danglingTailList[j].nodeList[k].state)
                                #print ("Node in list \n")
                                
       # print ("After running update code \n")
        self.stateDanglingNodes()
        # Recalculate the state of the metaatom
        self.calculateState()
        #print ("The old state is:\n" + str(oldStateNodes) + "\n")
        #print ("The new state is:\n" + str(newStateNodes) + "\n")
            
        
        # Next need to give each node in mol its state 
        
            #print ("Intensity before update: " + str(self.metaSpikes[i].intensity) + "\n")
            #print ("Intensity after update: " + str(self.metaSpikes[i].intensity) + "\n")
        # Now need to recalculate state
        
    
    def bondMetaAtoms(self):
        """ This function is called then the metaAtoms is involved in a bond, it sets the bonding status of the metaAtom to true
        """
        self.bonded = True
     
    def stateDanglingNodes(self):
        """ This function appends the state of each dangling node to a an array then returns the array, this was 
            used as a debugging tool but kept in as could be usefull later for more debugging
        
        """
        state = []
        for i in range(len(self.metaSpikes)):
            if self.metaSpikes[i].typeSpike == 1:
                #print ("Inside type 1 \n")
                #print ("The number of DNs is: " + str(len(self.metaSpikes[i].danglingNodeList)) + "\n")
                for j in range(len(self.metaSpikes[i].danglingNodeList)):
                        state.append(self.metaSpikes[i].danglingNodeList[j].node.state)

       #                 print ("Node in list \n")
            else:
                #print ("The number of DTs is: " + str(len(self.metaSpikes[i].danglingTailList)) + "\n")
                #print ("Inside type 2 \n")
                for j in range(len(self.metaSpikes[i].danglingTailList)):
                    for k in range(len(self.metaSpikes[i].danglingTailList[j].nodeList)):
                        state.append(self.metaSpikes[i].danglingTailList[j].nodeList[k].state)

        #print ("The state of dangling nodes is: \n" + str(state) + "\n")
        #print ("The length of dangling nodes is: \n" + str(len(state)) + "\n")   
        
        return state
    
    def updateIntensities (self,listAtoms):
        """ This function is called when the intensities of the spikes located atoms in the metaAtoms may have changed so there
            values need to be updated, this is done by going through every atom the molecule and giving each atoms spikes their
            new intenity values
        """
        
        for i in range(len(listAtoms)):
            for j in range(len(listAtoms[i].spikeArray)):
                self.mol[i].spikeArray[j].intensity = listAtoms[i].spikeArray[j].intensity
        
    
    def updateNodeStates (self,listAtoms):
        """ This function is called when the states of the nodes in the Atoms which make up the molecule have changed so 
            their values need to be updated this is done by going through every node and changing their value
        """
        
        for i in range(len(listAtoms)):
            for j in range(len(listAtoms[i].nodeArray)):
                self.mol[i].nodeArray[j].state = listAtoms[i].nodeArray[j].state
    
    
    def checkSpikeBonding (self):
        """ This function goes through every atom in the every molecule of the meta atom to see if the bonds are still stable
            if the bonds are not stable they broken, after breaking unstable bonds this process is run again until no more bonds
            break and the resulting structure is now stable. Note that at each stage the change in intensity and spike bonding status
            is passed to the underlying strucure
        """
        stable = True # If any bonds break this will be set to false
        stabilityChecker = True # Checks the result of each function call, if set to false then stable will be set to false
        # Go through each atom
        for i in range(len(self.mol)):
            # Go through each spike
            for j in range(len(self.mol[i].spikeArray)):
                if self.mol[i].spikeArray[j].bonded == True:
                    stabilityChecker = self.stabilitySpike(self.mol[i].spikeArray[j])
                    if stabilityChecker == False:
                        stable = False
                    #print (stable)
        if stable == True:
            print("No Bonds have broken \n")
        else:
            print ("Bonds have broken \n")
        return stable
    
    def stabilitySpike(self,spike):
        """ This function checks if the sum of the intensities of the spike and the spike it is bonded too meet the requred value
            in order to remain stable, if they dont then False is returned otherwise True is returned
        """
        stable = True
       # print ("The bonded RBN is: " + str(spike.bondedRBN) +  "\n")
        bondedSpike = spike.bondedRBN.spikeArray[spike.bondedSpikeNum]
        # If type 3 spikes then sum of an intensity can be plus or minus two for bond to be stable
        if spike.type == 3 and bondedSpike.type == 3:
            if abs(spike.intensity + bondedSpike.intensity) >= 2:
                # If unstable break bonds and set stable to false
                bondedSpike.bondBreak()
                spike.bondBreak()
                print ("Bond Broken in metaAtom: " + str(self.atomNumber) + "\n")
                print ("The RBN number broken is: " + str(spike.RBN.rbnNumber) + "\n")
                
                stable = False
        # If one or boths spikes is type 2 and neither is type 1 then the sum of intensitys can be plus or minus one to be stable
        elif (spike.type == 2 and bondedSpike.type == 2) or  (spike.type == 3 and bondedSpike.type == 2) or  (spike.type == 2 and bondedSpike.type == 3):
            if abs(spike.intensity + bondedSpike.intensity) >= 1:
                # If unstable break bonds and set stable to false
                bondedSpike.bondBreak()
                spike.bondBreak()
                print ("Bond Broken in metaAtom: " + str(self.atomNumber) + "\n")
                print ("The RBN number broken is: " + str(spike.RBN.rbnNumber) + "\n")
                stable = False
        # If one of the spikes or both is a type 1 spike then the intesnity needs to sum to zero to be stable
        else:
            if spike.intensity + bondedSpike.intensity != 0:
                # If unstable break bonds and set stable to false
                bondedSpike.bondBreak()
                spike.bondBreak()
                print ("Bond Broken in metaAtom: " + str(self.atomNumber) + "\n")
                print ("The RBN number broken is: " + str(spike.RBN.rbnNumber) + "\n")
                stable = False
        
#        # If any bonds have broken False will be returned otherwise True will be returned
#        if stable == True:
#            print ("Stability has not changed \n")
#        else:
#            print ("Stability has changed \n")
        return stable        
            
    def removeDoubleUnbondedAtoms (self):
        """ This function removes any atoms in the molecule which are not bonded to anything (i.e.Both spikes unbonded), this is done by going through every
            atom in the ring then going through the atoms spikes, if both spikes are not bonded and the dangling nodes from the spikes are not connected
            to a metaAtom then the atom is removed from the molecule. If the dangling nodes/ tails are bonded then the atom is removed
            from the ring and converted into a metaAtom and then is added to the metaMolecule
        """
        atomsToRemove = [] # Stores index of atoms we will need to remove
        
        # Go through each mol
        for i in range(len(self.mol)):
            # Atom is disconnected if number of unbonded spikes is equal to the number of spikes in the atom
            numUnbondedSpikes = 0
            for j in range(len(self.mol[i].spikeArray)):
                if self.mol[i].spikeArray[j].bonded == False:
                    # Spike not bonded so increment counter
                    numUnbondedSpikes += 1
            # If atom disconnected then need to check to see if dangling nodes or tails are bonded
            if numUnbondedSpikes == len(self.mol[i].spikeArray):
                print ("Atom: " + str(self.mol[i].rbnNumber) + " is being removed \n")
                anyBondedDanglingNodes  = False
                for j in range(len(self.mol[i].spikeArray)):
                    if self.isUnbondedAtomConnected(self.mol[i].spikeArray[j]) == True:
                        anyBondedDanglingNodes =  True
                # If atom has connected dangling nodes then need to convert atom to metaAtom, add metaAtom to metaMolecule and
                # remove atom from ring
                if anyBondedDanglingNodes == True:
                    print ("A new metaAtom is being created \n")
                    newMetaAtom = self.convertUnbondedAtomToMetaAtom(self.mol[i])
                    self.metaMolecule.addMetaAtom(newMetaAtom)
                    atomsToRemove.append(i)
        
        # Now need to remove atoms
        for i in range(len(atomsToRemove)):
            self.mol.pop(atomsToRemove[i])
        
        # Finally need to update metaMolecule with new mol    
        self.metaMolecule.updateListMols(self)            
    
    def isUnbondedAtomConnected (self,spike):
        """ This function is used by the function above, it takes a spike from an unbonded atom then checks to see if it has any
            dangling nodes or tails, if it doesnt then False is returned if it does then the bonding status of these dangling nodes
            or tails is checked, if they are unbonded False is returned. If they at least one is bonded then True is returned
        """   
        for i in range(len(self.metaSpikes)):
            if self.metaSpikes[i].typeSpike == 1:
                for j in range(len(self.metaSpikes[i].danglingNodeList)):
                    if self.metaSpikes[i].danglingNodeList[j].spike == spike and self.metaSpikes[i].danglingNodeList[j].bonded == True:
                        # If spike has a dangling node and the this is bonded then return True
                        return True
            else:
                for j in range(len(self.metaSpikes[i].danglingTailList)):
                    if self.metaSpikes[i].danglingTailList[j].spike == spike and self.metaSpikes[i].danglingTailList[j].bonded == True:
                        # If spike has a dangling tail and the this is bonded then return True
                        return True
        # If True has not been returned by this point either spike has no dangling tails or nodes or if it does have them then none
        # are bonded
        return False
    
    def convertUnbondedAtomToMetaAtom (self,atom):
        """ This function takes an unbonded atom and converts it to a metaAtom, it does this by generating a new metaAtom and 
            new metaSpikes, the dangling nodes and dangling tails of the bonds are added to the metaSpikes 
        """
        
        
        # Create a molecule consisting of just the atom
        mol = []
        mol.append(atom)
        # metaAtomNumber will be the length of the current metaMolecule
        metaAtomNumber =  len(self.metaMolecule.metaAtoms)
        newMetaAtom = MetaAtom(metaAtomNumber,atom)
        # Next need to generate metaspikes we can do this by generating a metaspike object then finding the dangling nodes
        # and tails already created for this spike
        type1MetaSpike = metSpk.MetaSpike(1,0)
        type2MetaSpike = metSpk.MetaSpike(2,0)
        numDanglingNodes = 0
        numDanglingTails = 0
        danglingNodesToRemove = []
        danglingTailsToRemove = []
        indexType1Spike = 0
        indexType2Spike = 0
        for i in range(len(self.metaSpikes)):
            if self.metaSpikes[i].typeSpike == 1:
                indexType1Spike = i
                for j in range(len(self.metaSpikes[i].danglingNodeList)):
                    print ("The length of the dangling node list is: " + str(len(self.metaSpikes[i].danglingNodeList)) + "\n")
                    print ("The value of j is: " + str(j) + "\n")
                    # See if dangling node belongs to atom being removed
                    if self.metaSpikes[i].danglingNodeList[j].spike.RBN == atom:
                        # If so then add to list of nodes being removed
                        type1MetaSpike.addDanglingNode(self.metaSpikes[i].danglingNodeList[j])# Add to metaspike
                        danglingNodesToRemove.append(j)
                        numDanglingNodes += 1
            else:
                indexType2Spike = i
                for j in range(len(self.metaSpikes[i].danglingTailList)):
                    print ("The length of the dangling tail list is: " + str(len(self.metaSpikes[i].danglingTailList)) + "\n")
                    print ("The value of j is: " + str(j) + "\n")
                    # See if dangling tail belongs to atom being removed
                    if self.metaSpikes[i].danglingTailList[j].spike.RBN == atom:
                        # If so then add to list of rails being removed
                        type2MetaSpike.addTailDanglingBonds(self.metaSpikes[i].danglingTailList[j]) # Add to metaspike
                        danglingTailsToRemove.append(j)
                        numDanglingTails += 1 
                
        # If any dangling nodes have been found then the spike can be added to the metaAtom
        if numDanglingNodes > 0:             
            newMetaAtom.addMetaSpike(type1MetaSpike)
            # Next need to remove danglingNodes no longer located in this metaAtom
            danglingNodesToRemove.sort(reverse = True) # Sort so popping values doesnt affect other indexes
            for i in range (len(danglingNodesToRemove)):
                self.metaSpikes[indexType1Spike].danglingNodeList.pop(danglingNodesToRemove[i])
                
       # If any dangling tails have been found then the spike can be added to the metaAtom
        if numDanglingTails > 0:             
            newMetaAtom.addMetaSpike(type2MetaSpike)
            # Next need to remove danglingtails no longer located in this metaAtom
            danglingTailsToRemove.sort(reverse = True) # Sort so popping values doesnt affect other indexes
            for i in range (len(danglingTailsToRemove)):
                self.metaSpikes[indexType2Spike].danglingNodeList.pop(danglingTailsToRemove[i])

        
        
        
        # We repeated this process for dangling tails
        
        type2MetaSpike = metSpk.MetaSpike(2,1)
        
        numDanglingTails = 0
        for i in range(len(self.metaSpikes)):
            if self.metaSpikes[i].typeSpike == 2:
                for j in range(len(self.metaSpikes[i].danglingTailList)):
                    if self.metaSpikes[i].danglingTailList[j].spike.RBN == atom:
                        type2MetaSpike.addTailDanglingBonds(self.metaSpikes[i].danglingTailList[j])
                        danglingTailsToRemove.append(j)
                        numDanglingTails += 1
        # If any dangling nodes have been found then the spike can be added to the metaAtom
        if numDanglingTails > 0:             
            newMetaAtom.addMetaSpike(type2MetaSpike)
            # Next need to remove danglingNodes no longer located in this metaAtom
            danglingTailsToRemove.sort(reverse = True) # Sort so popping values doesnt affect other indexes
            for i in range (len(danglingTailsToRemove)):
                self.metaSpikes[i].danglingTailList.pop(danglingTailsToRemove[i])
            
        return newMetaAtom
        
        
        
    