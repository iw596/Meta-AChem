# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 13:07:48 2018
This class defines the metamolecule, a metamolecule is a collection of metaatoms which are bonded together. The class stores the bonded
metaatoms and a list of all the atoms in the metaAtoms
@author: iw596
"""
import MetaAtom_v2 as metAt
import Reactor_v9 as reac
import numpy as np
class MetaMolecule:
    
    def __init__ (self,largestAtomSize):
        self.metaAtoms = []
        self.listMols = [] # stores the molecules the which make up the metaAtom
        self.largestAtomSize = largestAtomSize
        
    def addMetaAtom (self,metaAtom):
        """ This function adds a metaAtom to the metaMolecule, it does this by appending the metaAtom to the metaMolecule metaAtom list
            and appending the metaAtoms mol to this the list of mols
        """
        self.listMols.append(metaAtom.mol)
        self.metaAtoms.append(metaAtom)
        # Add MetaMolecule to metaAtom
        metaAtom.addMetaMolecule(self)
        print ("The length of listMols is now: " + str(len(self.listMols)) + "\n")
    
    def updateMetaMolecule (self):
        """ This function updates every metaAtom in the metaMolecule synchronously this is done by creating a list of every atom
            which makes up the metaAtoms molecules and updating all those synchrnously before passing the new states back to the 
            metaAtoms and recaulating the state of the metaAtoms
        """
        conMol = []
        for i in range(len(self.listMols)):
            for j in range(len(self.listMols[i])):
                conMol.append(self.listMols[i][j])
        # We now have a single large molecule which we can use to calculate the intensity of the spikes
        self.updateMolecule(conMol)
        
        # Now need to give the list of mols their new states
        offset = 0
        #print ("The length of the list of molecule is: " + str(len(self.listMols)) + "\n")
        for i in range(len(self.listMols)):
            for j in range(len(self.listMols[i])):
                #print ("The number of spikes is: " + str(len(self.listMols[i][j].spikeArray)) + "\n")
                for k in range(len(self.listMols[i][j].nodeArray)):
                 #   print ("For rbn: " + str(conMol[offset].rbnNumber) + " the intensity is now: " + str(conMol[offset].spikeArray[k].intensity) + "\n")
                    
                    self.listMols[i][j].nodeArray[k].state = conMol[offset].nodeArray[k].state
                offset += 1    
    
        # Update molecules in metaatoms with new states
        for i in range(len(self.metaAtoms)):
            self.metaAtoms[i].updateNodeStates(self.listMols[i])
            # Then calculate their new state
            self.metaAtoms[i].calculateState()
    
    def updateListMols (self,metaAtom):
        """ This updates the list of mols which is associated with the metaAtom passed as an argument """
        
        # First need to find index
        index  = -1
        for i in range(len(self.metaAtoms)):
            if self.metaAtoms[i] == metaAtom:
                index = i
                break
        # Update list mols
        self.listMols[index] = metaAtom.mol
   
    
    
    def calculateIntensity (self):
        """ This function calculates the intensity of every spike in the list of molecules, this is done by converting the individual
            list into a single large list before calling a function which calculates the intensity, this new intensity is then
            passed onto the smallest lists
        """
        conMol = []
        for i in range(len(self.listMols)):
            for j in range(len(self.listMols[i])):
                conMol.append(self.listMols[i][j])
        # We now have a single large molecule which we can use to calculate the intensity of the spikes
        origStates = []
        
        for i in range(len(conMol)):
            origStates.append(conMol[i].states)
            conMol[i].zeroRBN()

            #print ("The state matrix for atom " + str(i) + " is: \n" + str(mol[i].states) + "\n")
        for i in range(self.largestAtomSize + 30):
            self.updateMolecule(conMol)
        #for i in range (len(conMol)):
           # for j in range(len(conMol[i].spikeArray)):
               # print ("For rbn: " + str(conMol[i].rbnNumber) + " the intensity is: " + str(conMol[i].spikeArray[j].intensity) + "\n")
        
        # Recalcaulate the intensity of every spike inside molecule
        self.analyseMolecule(conMol)   
       
        for i in range(len(conMol)):
            conMol[i].setState(origStates[i],np.size(origStates[i],0)-1)
        offset = 0
        #print ("The length of the list of molecule is: " + str(len(self.listMols)) + "\n")
        for i in range(len(self.listMols)):
            for j in range(len(self.listMols[i])):
                #print ("The number of spikes is: " + str(len(self.listMols[i][j].spikeArray)) + "\n")
                for k in range(len(self.listMols[i][j].spikeArray)):
                 #   print ("For rbn: " + str(conMol[offset].rbnNumber) + " the intensity is now: " + str(conMol[offset].spikeArray[k].intensity) + "\n")
                    
                    self.listMols[i][j].spikeArray[k].intensity = conMol[offset].spikeArray[k].intensity
                offset += 1    
        
        
        for i in range(len(self.metaAtoms)):
            for j in range(len(self.metaAtoms[i].metaSpikes)):
                self.metaAtoms[i].metaSpikes[j].debugIntensity()
            
        # Update molecules in metaatoms with new intensity
        for i in range(len(self.metaAtoms)):
            self.metaAtoms[i].updateIntensities(self.listMols[i])
    
    def updateMolecule (self,molecule):
        poppedStates = [] # Stores the popped states 
        for i in range(len(molecule)):
            molecule[i].updateRBN()
            poppedStates.append(molecule[i].popState())
        
        # Afer all updated states found add these states back to atoms
        for i in range(len(molecule)):
            molecule[i].appendState(poppedStates[i])
        
    def analyseMolecule (self,molecule):
        """ This function goes through every Atpm in a molecule and calculates the intensity of the spikes, this is needed
            as bonding and unbonding causes spike intensity values to change
        """
        for i in range(len(molecule)):
            self.analyseAtom(molecule[i])
    
    def analyseAtom (self,rbn):
        """ Recalculculates intensity of a atom """
        for i in range(np.size(rbn.spikeArray)):
           # print ("The value of i is: " + str(i) + "\n")
            rbn.spikeArray[i].calcMolIntenisty()
    
    def checkMolecularStructure (self):
        """ This function goes through every metaatom and checks the stability of the underlying molecular structure if any bonds
            break then the list of mols is updated and the function is run again until no more bonds break and the resulting structure
            is stable
        """
        stable = True # If any bonds break this will be set to false
        stabilityChecker = True # Checks the result of each function call, if set to false then stable will be set to false
        
        # First calculate the intensity
        self.calculateIntensity()
    
        # Next go through every metaatom and check the stability of the metaatoms molecular structure
        for i in range(len(self.metaAtoms)):
            stabilityChecker = self.metaAtoms[i].checkSpikeBonding()
            if stabilityChecker == False:
                stable = False
        
        if stable == True:
            print ("No molecular changes \n")
        else:
            print ("Some molecular changes \n")
    
        if stable == False:
            self.checkMolecularStructure()
        # After obtaining final strucrure need to update the listMols
        for i in range(len(self.metaAtoms)):
            self.listMols[i] = self.metaAtoms[i].mol
        
        return stable
     
    def checkDanglingNodeStability (self):        
        """ This function checks the stability of all the dangling node bonds in the metaMolecule, it does this by going through every
            dangling node bond in each meta atom and checking to see if the bonding criteria is still satisfied, if not then 
            the dangling node bond is broken
        """
        
        stable = True
        # Go through each meta atom
        for i in range(len(self.metaAtoms)):
            for j in range(len(self.metaAtoms[i].metaSpikes)):
                # If a type 2 spike then go through each dangling tail
                if self.metaAtoms[i].metaSpikes[j].typeSpike == 1:
                    for k in range(len(self.metaAtoms[i].metaSpikes[j].danglingNodeList)):
                        danglingNode = self.metaAtoms[i].metaSpikes[j].danglingNodeList[k]
                        if danglingNode.bonded == True:
                            bondedDanglingNode = danglingNode.connectedDanglingNode
                            # Check to see if spike bonding criteria still met
                            if danglingNode.spike.type == 3 and bondedDanglingNode.spike.type == 3:
                                # If criteria not met then break bond
                                if abs(danglingNode.spike.intensity + bondedDanglingNode.spike.intensity) > 2:
                                    # Note need to finish functions below
                                    bondedDanglingNode.bondBroken()
                                    danglingNode.bondBroken()
                                    print ("dangling tail bond is unstable \n")
                                    stable = False
                            
                            elif (danglingNode.spike.type == 2 and bondedDanglingNode.spike.type == 2) or (danglingNode.spike.type == 2 and bondedDanglingNode.spike.type == 3) or (danglingNode.spike.type == 3 and bondedDanglingNode.spike.type == 2):
                                if abs(danglingNode.spike.intensity + bondedDanglingNode.spike.intensity) > 1:
                                    # Note need to finish functions below
                                    bondedDanglingNode.bondBroken()
                                    danglingNode.bondBroken()
                                    print ("dangling tail bond is unstable \n")
    
                                    stable = False
                            else:
                                if abs(danglingNode.spike.intensity + bondedDanglingNode.spike.intensity) != 0:
                                    # Note need to finish functions below
                                    bondedDanglingNode.bondBroken()
                                    danglingNode.bondBroken()
                                    print ("dangling tail bond is unstable \n")
                                    stable = False 
            
        
        return stable                        
                            
    def checkDanglingTailStability (self):
        """ This function checks the stability of all the dangling tail bonds in the metaMolecule, it does this by going through every
            dangling tail bond in each meta atom and checking to see if the bonding criteria is still satisfied, if not then 
            the dangling tail bond is broken
        """
        
        stable = True
        # Go through each meta atom
        for i in range(len(self.metaAtoms)):
            for j in range(len(self.metaAtoms[i].metaSpikes)):
                # If a type 2 spike then go through each dangling tail
                if self.metaAtoms[i].metaSpikes[j].typeSpike == 2:
                    
                    for k in range(len(self.metaAtoms[i].metaSpikes[j].danglingTailList)):
                        danglingTail = self.metaAtoms[i].metaSpikes[j].danglingTailList[k]
                        if danglingTail.bonded == True:
                            bondedDanglingTail = danglingTail.bondedTail
                            # Check to see if spike bonding criteria still met
                            if danglingTail.spike.type == 3 and bondedDanglingTail.spike.type == 3:
                                # If criteria not met then break bond
                                if abs(danglingTail.spike.intensity + bondedDanglingTail.spike.intensity) > 2:
                                    # Note need to finish functions below
                                    bondedDanglingTail.bondBroken()
                                    danglingTail.bondBroken()
                                    print ("dangling node bond is unstable \n")
                                    stable = False
                            
                            elif (danglingTail.spike.type == 2 and bondedDanglingTail.spike.type == 2) or (danglingTail.spike.type == 2 and bondedDanglingTail.spike.type == 3) or (danglingTail.spike.type == 3 and bondedDanglingTail.spike.type == 2):
                                if abs(danglingTail.spike.intensity + bondedDanglingTail.spike.intensity) > 1:
                                    # Note need to finish functions below
                                    bondedDanglingTail.bondBroken()
                                    danglingTail.bondBroken()
                                    print ("dangling node bond is unstable \n")
                                    stable = False
                            else:
                                if abs(danglingTail.spike.intensity + bondedDanglingTail.spike.intensity) != 0:
                                    # Note need to finish functions below
                                    bondedDanglingTail.bondBroken()
                                    danglingTail.bondBroken()
                                    print ("dangling node bond is unstable \n")
                                    stable = False 
            
        
        return stable     
                
    
    
    
    
    def intensityOFEverySpikeDebug (self):
        """ This is a debugging function which returns a list containing the intensity of every spike in the metamolecule"""
        intensitys = []
        for i in range(len(self.listMols)):
            for j in range (len(self.listMols[i])):
                for k in range(len(self.listMols[i][j].spikeArray)):
                    intensitys.append(self.listMols[i][j].spikeArray[k].intensity)
        print ("The intensitys is: " + str(intensitys) + "\n")
        return intensitys
    
    
    def removeMetaAtom (self,metaAtom):
        """ This function removes the metaAtom passed as an argument from the metaMolecule, if the metaAtom is not in the
            meta molecule then an error message is printed
        """
        
        if metaAtom not in self.metaAtoms:
            print ("Meta atom not in meta molecule")
        else:
            # Find index and remove metaAtom and molecule associated with meta atom
            indexMetaAtom = self.metaAtoms.index(metaAtom)
            self.metaAtoms.pop(indexMetaAtom)
            self.listMols.pop(indexMetaAtom)
    
    def removedUnbondedAtoms (self):
        """This function goes through every metaatom and removes atoms which are no longer atttached to
           the ring
        """
        
        for i in range(len(self.metaAtoms)):
            # For every meta atom call function which removes atoms with both spikes unbonded
            self.metaAtoms[i].removeDoubleUnbondedAtoms()
    
    
    
    def stabilityDebugTest (self):
        """ This function tests that the intensty of the metamolecule is stable over time by calculating the initiail intensity, then
            the state of the metaatom before updating the metaatom again and recalculating the stability if the test is succesfull
            then the intensity will not have changed
        """
        print ("The initial intensitys is: \n" + str(self.intensityOFEverySpikeDebug()) + "\n")
        self.calculateIntensity()
        print ("The initial intensitys after calling recalc function is: \n" + str(self.intensityOFEverySpikeDebug()) + "\n")
        for i in range(100):
            self.updateMetaMolecule()
        self.calculateIntensity()
        print ("The intensity after update is: \n" + str(self.intensityOFEverySpikeDebug()) + "\n")