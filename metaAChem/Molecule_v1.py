# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 18:36:22 2018
This class defines a molecule, a molecule is a collection of atoms, a molecule can either be a signle ended or
a double ended chain, a single ended chain has only one unbonded spike, a double ended chain has two unbonded
spike. A double ended chain has a chance that the two unbonded spikes will attempt to react and form a ring
@author: isaac
"""
import numpy as np
import Reactor_v10

class Molecule:
    
    def __init__ (self):
        self.atoms = []
        self.sizeMol = 0
        self.type = 2
        self.head = 0
        self.tail  = 0
   
    def addInitialPair (self,atom1,atom2,reactor):
        """ This function adds an initial pair of bonded atoms to the molecule, it then checks the
         stability of the bond by recalculating the intensity, if the bond is stable true is returned, if
         the bond is not stable then 
        """
        self.atoms.extend((atom1,atom2))
        self.reactAtoms(self.atoms[0],self.atoms[1],reactor)
        self.sizeMol += 2
        # Next we need to react the two atoms
        
        return 1
    
    
    
    
    def addAtom (self,atom):
        self.atoms.append(atom)
        self.sizeMol += 1
        # If atom only has one spike then the molecule will become a type 1 molecule
        if np.size(atom.spikeArray) == 1:
            self.type = 1
    
    def calculateIntensitySpikes (self,reactor):
        """ This function calculates the intensity of every spike in the molecule, it does this by updating the molecule
            and then calculating the intensity of evey spike
        """
        
        origStates = []
        for i in range(self.sizeMol):
            origStates.append(self.atoms[i].states)
            self.atoms[i].zeroRBN()
            #print ("The state matrix for atom " + str(i) + " is: \n" + str(mol[i].states) + "\n")
        for i in range(reactor.maxSizeAtom + 30):
            self.updateMolecule()

        
        self.analyseMolecule()   
       
        for i in range(self.sizeMol):
            self.atoms[i].setState(origStates[i],np.size(origStates[i],0)-1)
    

    def updateMolecule(self):
        """ This function updates every Atom in a molecule synchronously, it does this by updating an Atom then popping
            its new state, before moving on to the next Atom, this ensures that every Atom is updated synchronously with
            the other Atoms. After this update and pop loop the popped states are appeneded to the Atoms 
        """
        poppedStates = [] # Stores the popped states 
        for i in range(self.sizeMol):
            self.atoms[i].updateRBN()
            poppedStates.append(self.atoms[i].popState())
        
        # Afer all updated states found add these states back to atoms
        for i in range(self.sizeMol):
            self.atoms[i].appendState(poppedStates[i])            
            
    def analyseMolecule (self):
        """ This function goes through every Atpm in a molecule and calculates the intensity of the spikes, this is needed
            as bonding and unbonding causes spike intensity values to change
        """
        for i in range(self.sizeMol):
            self.analyseAtom(self.atoms[i])
    
    def analyseAtom (self,atom):
        """ Recalculculates intensity of a atom """
        for i in range(np.size(atom.spikeArray)):
           # print ("The value of i is: " + str(i) + "\n")
            atom.spikeArray[i].calcMolIntenisty()
    
    def checkFullyBonded (self):
        """ This function goes through every spike in the head and tail of the molecule and checks its bonding status, if
            there are any unbonded spikes False is returned as further bonding can take place
            other wise True is returned as no further bonds can take place
        """
        # Check every spike in head
        for i in range(len(self.head.spikeArray)):
            if self.head.spikeArray[i].bonded == False:
                return False
        
        # Check every spike in tail
        for i in range(len(self.tail.spikeArray)):
            if self.tail.spikeArray[i].bonded == False:
                return False
            
        return True
    
    
    def reactAtoms (self,atom1,atom2,reactor):
        """ This function is used to try and react two RBNs, it works by going across
            the spike list of one RBN and comparing it agaisnt the spike list of another RBN it the sum
            is zero the RBNs are reacted together, if the product is stable then True is returned, if the 
            product is not stable we attempt to find another spike paiur we can use, if there are no spike
            pairs left and the RBNs still havent reacted then we return False to indicate no succesfull reaction
            has taken place
        """
       # Check to see that both atoms have valid intensities
        for i in range(np.size(atom1.spikeArray)):
            if atom1.spikeArray[i].intensity == 'a':
                if atom1 in self.soup:
                    self.soup.remove(atom1)
                    self.soup.append(self.rbnGenerator())
                return False
        
        for i in range(np.size(atom2.spikeArray)):
            if atom2.spikeArray[i].intensity == 'a':
                if atom2 in self.soup:
                    self.soup.remove(atom2)
                    self.soup.append(self.rbnGenerator())
                return False
        
        
        # Next check that the intensities of the spikes can sum to meet a critieria which
        # depends on the type of spike
        for i in range(np.size(atom1.spikeArray)):
            for j in range(np.size(atom2.spikeArray)):
                if atom1.spikeArray[i].type  +  atom2.spikeArray[j].type == 6 :
                    #print ("Reacting type 2 spikes \n")
                    if abs(atom1.spikeArray[i].intensity + atom2.spikeArray[j].intensity) <= 1 and atom1.spikeArray[i].bonded == False and atom2.spikeArray[j].bonded == False:
                        #print ("Reacting RBNs \n")
                        #print ("Intensity of spike 1: " + str(RBN1.spikeArray[i].intensity) + "\n")
                        ##print ("Intensity of spike 2: " + str(RBN2.spikeArray[j].intensity) + "\n")
                        sucessReacted = bp.reactRBNS(atom1.spikeArray[i],atom2.spikeArray[j],atom1,atom2,mol,self)
                        if sucessReacted == True:
                            return True # Stable bond so return true

                
                if atom1.spikeArray[i].type + atom2.spikeArray[j].type >= 4 + atom1.spikeArray[i].type + atom2.spikeArray[j] .type <6:
                    #print ("Reacting type 2 spikes \n")
                    if abs(atom1.spikeArray[i].intensity + atom2.spikeArray[j].intensity) <= 1 and atom1.spikeArray[i].bonded == False and atom2.spikeArray[j].bonded == False:
                        #print ("Reacting RBNs \n")
                        #print ("Intensity of spike 1: " + str(RBN1.spikeArray[i].intensity) + "\n")
                        ##print ("Intensity of spike 2: " + str(RBN2.spikeArray[j].intensity) + "\n")
                        sucessReacted = bp.reactRBNS(atom1.spikeArray[i],atom2.spikeArray[j],atom1,atom2,mol,self)
                        if sucessReacted == True:
                            return True # Stable bond so return true
                else:
                   # print ("Reacting type 1 spikes \n")
                    if atom1.spikeArray[i].intensity + atom2.spikeArray[j].intensity == 0 and atom1.spikeArray[i].bonded == False and atom2.spikeArray[j].bonded == False:
                        #print ("Reacting RBNs \n")
                        sucessReacted = bp.reactRBNS(atom1.spikeArray[i],atom2.spikeArray[j],atom1,atom2,mol,self)
                        if sucessReacted == True:
                            return True # Stable bond so return true
                    
        return False # If no stable bond found return False