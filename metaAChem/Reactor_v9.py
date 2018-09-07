# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 11:21:32 2018

@author: iw596
"""
import random
import math
import pickle
import bondingProtocol_v5 as bp
import RBN_v4 as rbn
import Spike_v4 as spk
import os.path

import numpy as np
#from numpy import *
class Reactor:
    """ The reactor is used to generate rings of atoms, it does this by creating a soup of atoms and then randomly selecting
        atoms from this soup and joining to form a molecule, this is done until either a ring forms or the soup is empty
        Note each atom in the soup has two spikes and a cyclength of three
    """
    def __init__ (self,numProducts,sizeSoup,maxSizeAtom):
        print ("Building Reactor \n")
        self.wantedNumProducts  = numProducts
        self.sizeSoup = sizeSoup
        self.maxSizeAtom = maxSizeAtom
        self.soup = []
        self.genSoup()
        print ("Reactor Built \n")
    
    def runReactor (self):
        """ This function runs the reactor. It does this by calling a function which attempts to generatre rings, it does
            this until the wanted number of products has been reached. The products are also pickled by a function called
            by the reactor so they can be used later
        """
        print ("Reactor started \n")
        setRings = []
        currentNumProducts = 0
        while currentNumProducts != self.wantedNumProducts:
            ring = self.generateRing ()
            if ring != -1:
                print ("The final size of the molecule is: " + str(len(ring)) + "\n")
                #self.checkStability (ring)
               
                setRings.append(ring)
                currentNumProducts += 1
                print ("The current number of rings is: " + str(currentNumProducts) + "\n")
                if len(setRings)%5 == 0 and len(setRings) > 0 :
                    self.pickleSet(setRings)
                    setRings = []
            else:
                print ("Reactor failed \n")
            self.genSoup()    

    def pickleSet (self,setMolecules):
        # Path where all pickled files stored need to change depending on computer
        path = 'H:\\YCCSA\\28-08-2018\\ringMolecules' 
        # First need to write number of files to text file
        numSetsFileName = "numSets" # File which stores the number of steps
        numSetsFileName = os.path.join(path,numSetsFileName+".txt")
        # Open file so we can read and write to it
        numSetsFile = open(numSetsFileName,"a+") # If file does not exist it will be created
        numSetsFile.seek(0)
        numFiles = numSetsFile.read()
        print ("The num of files is: " + str(numFiles) + "\n")
        # Can now close this file
        numSetsFile.close()
        # If file did not exist before opening
        if numFiles == "":
            # Reopen file and write 1 to it
            numSetsFile = open(numSetsFileName,"w+") 
            numFiles = 1
            numSetsFile.write(str(numFiles) + "\n")
        # Otherwise reset file by opening in w+ mode and write new value    
        else:
            numSetsFile = open(numSetsFileName,"w+") 
            numFiles = int(numFiles)  + 1
            numSetsFile.write(str(numFiles) + "\n")
            
        numSetsFile.close()    
            
        # Now need to pickle set of molecule
        pickleFileName = "rings" + str(numFiles - 1) + ".pickle" # Name of of the fail,minus one as zero indexing
        pickleFileName = os.path.join(path,pickleFileName) # Join the path with the name
        pickle_in = open(pickleFileName,"wb") # Open new pickle file
        pickle.dump(setMolecules,pickle_in) # Dump data in file
        pickle_in.close() # Close the file            
                
    
    
    def generateRing(self):
        """ This function is used to attempt to generate a ring of atoms it does this by selecting atoms from the soup and bonding
            them to form a molecule, it keeps doing this until either a ring forms or the soup is empty or the max number of 
            attempts has been reached. If a ring is formedthen the ring is returned, otherwise False is returned to indicate 
            no ring was produced
        """
        print ( "Reactor run started \n")
        mol = []
        mol= self.createInitialPair(mol)
        # If initial pair has not been generated then return false so new soup can be made
        if mol == -1:
            return False
        headAtom = mol[0]
        tailAtom = mol[1]
        
        if self.completeRing(headAtom,tailAtom,mol) == True:
            return mol
        
        
        self.soup.extend((self.rbnGenerator(),self.rbnGenerator()))
        maxTimeSteps = self.sizeSoup * 10 + 500
        curTimeSteps = 0
        while len(self.soup) != 0 and curTimeSteps != maxTimeSteps:
            testAtom = self.selectAtom() 
            #print ("New atom selected \n")
            headOrTail = random.random()
            # If less than 0.5 react atom with head of molecule
            if headOrTail <= 0.5:
                if self.reactAtomWithMol (testAtom,headAtom,mol) == True:
                    #print ("Head Reacted \n")
                  
                    print ("Size of molecule is now: " + str(len(mol)) + "\n")
                    #self.intensityDebug(mol)
                    self.calculateIntensitySpikes(mol)
                 #   print ("RBN number of head atom was: " + str(headAtom.rbnNumber) + "\n")
                    headAtom = testAtom
                   # print ("RBN number of head atom is now: " + str(headAtom.rbnNumber) + "\n")
                    self.soup.append(self.rbnGenerator())
                    if self.completeRing(headAtom,tailAtom,mol) == True:
                        self.calculateIntensitySpikes(mol)
                        return mol
                    else:
                         self.calculateIntensitySpikes(mol)
                else:
                    reReact = random.random()
                    if reReact >= 0.3:
                        if self.reactAtomWithMol (testAtom,tailAtom,mol) == True:
                            #print ("Tail Re-Reacted \n")
                            print ("Size of molecule is now: " + str(len(mol)) + "\n")
                            #self.intensityDebug(mol)
                            self.calculateIntensitySpikes(mol)
                            self.soup.append(self.rbnGenerator())
#                            print ("RBN number of tail atom was: " + str(tailAtom.rbnNumber) + "\n")
                            tailAtom = testAtom
                          #  print ("RBN number of tail atom is: " + str(tailAtom.rbnNumber) + "\n")
                            if self.completeRing(headAtom,tailAtom,mol) == True:
                                return mol
                            else:
                                self.calculateIntensitySpikes(mol)
                         #reReact = random.random()
                        else:
                            self.calculateIntensitySpikes(mol)
                            self.soup.append(testAtom) 

                    else:
                        self.calculateIntensitySpikes(mol)
                        self.soup.append(testAtom)
            else:
                if self.reactAtomWithMol (testAtom,tailAtom,mol) == True:
                    #print ("Tail Reacted \n")
                    print ("Size of molecule is now: " + str(len(mol)) + "\n")
                    #self.intensityDebug(mol)
                    self.calculateIntensitySpikes(mol)
                   # print ("RBN number of tail atom was: " + str(tailAtom.rbnNumber) + "\n")
                    tailAtom = testAtom
                    #print ("RBN number of tail atom is now: " + str(tailAtom.rbnNumber) + "\n")
                    self.soup.append(self.rbnGenerator())
                    if self.completeRing(headAtom,tailAtom,mol) == True:
                         return mol
                    else:
                         self.calculateIntensitySpikes(mol)
                     
                else:
                    reReact = random.random()
                    if reReact >= 0.3:
                        if self.reactAtomWithMol (testAtom,headAtom,mol) == True:
                           # print ("Head Re-Reacted \n")
                            print ("Size of molecule is now: " + str(len(mol)) + "\n")
                            #self.intensityDebug(mol)
                            #self.calculateIntensitySpikes(mol)
                            self.soup.append(self.rbnGenerator())
                            #print ("RBN number of head atom was: " + str(headAtom.rbnNumber) + "\n")
                            headAtom = testAtom
                            #print ("RBN number of head atom is now: " + str(headAtom.rbnNumber) + "\n")
                            if self.completeRing(headAtom,tailAtom,mol) == True:
                                return mol
                            else:
                                self.calculateIntensitySpikes(mol)
                         #reReact = random.random()
                        else:
                            self.calculateIntensitySpikes(mol)
                            self.soup.append(testAtom) 

                    else:
                        self.calculateIntensitySpikes(mol)
                        self.soup.append(testAtom)
            
            self.incrementMoleculeState(mol)
            curTimeSteps += 1
            if len(self.soup) < self.sizeSoup:
                print ("Soup is not the max size \n")
                print ("Size of soup is: " + str(len(self.soup)) + "\n")
            
            
            
            if curTimeSteps % 200 == 0:
                print ("Program progressing \n")
                #self.intensityDebug(mol)
                #self.printNumberUnbondedSpikes(mol)
                #print ("The head RBN number is: " + str(headAtom.rbnNumber) + "\n")
                #print ("The tail RBN number is: " + str(tailAtom.rbnNumber) + "\n")
        
        #print ("Length of failed molecule is: " + str(len(mol)) + "\n")
        return -1    
    

    
    
    def intensityDebug (self,mol):
        intensitys = []
        lengths = []
        spkNumbers = []
        
        for i in range(len(mol)):
            for j in range(len(mol[i].spikeArray)):
              #  print ("The number of spikes is: " + str(len(mol[i].spikeArray)) + "\n")
                intensitys.append(mol[i].spikeArray[j].intensity)
                lengths.append(len(mol[i].spikeArray[j].nodeList))
                spkNumbers.append(mol[i].spikeArray[j].spikeNumber)
        print ("The intensity of every spike is: \n" + str(intensitys) + "\n")
       # print ("The length of each spike is \n" + str(lengths) + "\n")
      #  print ("The number of every spike is: \n" + str(spkNumbers) + "\n")
    
    
    
    
    def createInitialPair (self,mol):
    
        """ This function is used to generate the initial molecule, it does this by bonding two atoms from the soup and 
            returning them. If this cannot be achieved after a certain number of time steps then -1 is returned to 
            indicate that no initial pair could be formed
        """
        
        # Formula found here https://stackoverflow.com/questions/18859430/how-do-i-get-the-total-number-of-unique-pairs-of-a-set-in-the-database
        # Answer 1 wher k = 2
        maxNumTimeSteps = math.factorial(len(self.soup))/(4*math.factorial(len(self.soup) - 2))
        curNumSteps = 0
        
        while curNumSteps != maxNumTimeSteps:
            atom1 = self.selectAtom()
            atom2 = self.selectAtom()
            mol.extend((atom1,atom2))
            check = self.reactAtoms(atom1,atom2,mol)
            
            while check != True and curNumSteps != maxNumTimeSteps:
                mol = []
                self.soup.append(atom1)
                self.soup.append(atom2)
                atom1 = self.selectAtom()
                atom2 = self.selectAtom()
                mol.extend((atom1,atom2))
                check = self.reactAtoms(atom1,atom2,mol)
                curNumSteps += 1
            
            if check == True:
                return mol
            else:
                return -1,-1
    
    
    def reactAtomWithMol(self,atom,molAtom,molecule):
        """ This function attempts to react an atom with a molecule, if the reaction is sucesfull and the
            molecule stable then true is returned otherwise false is returned. 
        """
        molecule.append(atom)
        self.calculateIntensitySpikes(molecule)
        if self.reactAtoms(molAtom,atom,molecule) == True:
            #molecule.append(atom)
            self.calculateIntensitySpikes(molecule)
            if self.checkStructure(molecule) == False:
                #print ("Stabiity test failed bottom bonding \n")
                molecule.pop()
                for i in range(np.size(molAtom.spikeArray)):
                    if molAtom.spikeArray[i].bondedRBN == atom:
                        molAtom.spikeArray[i].bondBreak()
                self.calculateIntensitySpikes(molecule)
                self.analyseAtom(atom)
                return False
            else:
                self.calculateIntensitySpikes(molecule)
                return True
        else:
            self.calculateIntensitySpikes(molecule)
            molecule.pop()
            
            return False        
        
    
    def completeRing(self,headAtom,tailAtom,molecule):
        #self.debugHeadTail(headAtom,tailAtom)
        if self.reactAtoms(headAtom,tailAtom,molecule) == True:
            self.calculateIntensitySpikes(molecule)
            if self.checkStructure (molecule) == True:
                print ("Ring Created \n")
                return True
            else:
                spikeToBreakIndex = headAtom.activeSpikes[np.size(headAtom.activeSpikes) - 1]
                # Breaks the tail atoms bond to head atom
                headAtom.spikeArray[spikeToBreakIndex].bondedRBN.spikeArray[headAtom.spikeArray[spikeToBreakIndex].bondedSpikeNum].bondBreak()
                # Breaks the head atoms bond to tail atom
                headAtom.spikeArray[spikeToBreakIndex].bondBreak()
                # Recalculates the intensity
                self.calculateIntensitySpikes(molecule)
                print ("Ring unstable \n")
                self.printNumberUnbondedSpikes(molecule)
                return False
        else:
            #print ("Head and Tail cannot bond \n")
            self.calculateIntensitySpikes(molecule)
            #self.debugHeadTail(headAtom,tailAtom)
#            self.printNumberUnbondedSpikes(molecule)
            return False
     
   
    
    
    
    
    def rbnGenerator (self):
        """ This funcion is used to generate an RBN with a random  number of nodes, the max possible given
            by maxSizeRBN, max number of connections is numberNodes -1, the function generates RBN and checks
            that its cycle length is equal to three before returning this RBN, note all RBNs currently have
            initial state of 0 (all nodes state 0)
        """
    
        
        while True:
            numNodes = random.randint (3,self.maxSizeAtom + 1) # Generates random number for node, note need to add 1 so maxSizeRBN is possible
            numConnections = 2 # At moment fix k = 2 as otherwise RBN too unstable
            testRBN = rbn.RBN(numNodes,numConnections,len(self.soup) - 1) # Generate the RBN 
            cycleLength = testRBN.findCycleLength()
            
            if np.size(testRBN.spikeArray) == 2 and cycleLength == 3 :
                #print ("The number of spikes is: " + str(size(testRBN.spikeArray)) + "\n")
                return testRBN # If cyclelength is correct size then return RBN
    
    def genSoup (self):
        """ This function generates a soup of atoms by calling the RBN generator function until the wanted
            number of atoms in the soup has been reached
        """
        self.soup = []
        #print ("Generating soup \n")
        for i in range (self.sizeSoup):
            self.soup.append(self.rbnGenerator())
            #print ("Size of soup is now: " + str(len(self.soup)) + "\n")
        #print ("Soup of size: " + str(len(self.soup)) + " has been produced \n")
    
    def selectAtom (self):
        """ This function is used to select an atom from the soup, this is achieved by randomly sorting the list
            and then popping an atom from the soup
        """
        random.shuffle (self.soup)
        return self.soup.pop()
    
    
    def reactAtoms (self,atom1,atom2,mol):
        """ This function is used to try and react two RBNs, it works by going across
            the spike list of one RBN and comparing it agaisnt the spike list of another RBN it the sum
            is zero the RBNs are reacted together, if the product is stable then True is returned, if the 
            product is not stable we attempt to find another spike paiur we can use, if there are no spike
            pairs left and the RBNs still havent reacted then we return False to indicate no succesfull reaction
            has taken place
        """
       
        for i in range(np.size(atom1.spikeArray)):
            if atom1.spikeArray[i].intensity == 'a':
                return False
        
        for i in range(np.size(atom2.spikeArray)):
            if atom2.spikeArray[i].intensity == 'a':
                return False
        
        
        
        for i in range(np.size(atom1.spikeArray)):
            for j in range(np.size(atom2.spikeArray)):
                spike1 = atom1.spikeArray[i]
                spike2 = atom2.spikeArray[j]
                if isinstance(spike1.intensity, str) or isinstance(spike2.intensity, str):
                    return False
                if spike1.type == 3  and spike2.type == 3 :
                    #print ("Reacting type 2 spikes \n")
                    if abs(atom1.spikeArray[i].intensity + atom2.spikeArray[j].intensity) <= 1 and spike1.bonded == False and spike2.bonded == False:
                        #print ("Reacting RBNs \n")
                        #print ("Intensity of spike 1: " + str(RBN1.spikeArray[i].intensity) + "\n")
                        ##print ("Intensity of spike 2: " + str(RBN2.spikeArray[j].intensity) + "\n")
                        sucessReacted = bp.reactRBNS(spike1,spike2,atom1,atom2,mol,self)
                        if sucessReacted == True:
                            return True # Stable bond so return true

                
                elif (spike1.type == 2  and spike2.type == 2) or (spike1.type == 2  and spike2.type == 3) or (spike1.type == 3  and spike2.type == 2):
                    #print ("Reacting type 2 spikes \n")
                    if abs(spike1.intensity + spike2.intensity) <= 1 and spike1.bonded == False and spike2.bonded == False:
                        #print ("Reacting RBNs \n")
                        #print ("Intensity of spike 1: " + str(RBN1.spikeArray[i].intensity) + "\n")
                        ##print ("Intensity of spike 2: " + str(RBN2.spikeArray[j].intensity) + "\n")
                        sucessReacted = bp.reactRBNS(spike1,spike2,atom1,atom2,mol,self)
                        if sucessReacted == True:
                            return True # Stable bond so return true
                else:
                   # print ("Reacting type 1 spikes \n")
                    if spike1.intensity + spike2.intensity == 0 and spike1.bonded == False and spike2.bonded == False:
                        #print ("Reacting RBNs \n")
                        sucessReacted = bp.reactRBNS(spike1,spike2,atom1,atom2,mol,self)
                        if sucessReacted == True:
                            return True # Stable bond so return true
                    
        return False # If no stable bond found return False


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
    
    
    def updateMolecule(self,molecule):
        """ This function updates every RBN in a molecule synchronously, it does this by updating an RBN then popping
            its new state, before moving on to the next RBN, this ensures that every RBN is updated synchronously with
            the other RBNs. After this update and pop loop the popped states are appeneded to the RBNs 
        """
        poppedStates = [] # Stores the popped states 
        for i in range(len(molecule)):
            molecule[i].updateRBN()
            poppedStates.append(molecule[i].popState())
        
        # Afer all updated states found add these states back to atoms
        for i in range(len(molecule)):
            molecule[i].appendState(poppedStates[i])
    
    def checkStructure (self,mol):
        """ To check the structure of the molecule we update it many times, then calculate the
            intensity of each spike before checking that the mol is stable, if it is not then False
            is returned if it is then true is returned
        """
        self.calculateIntensitySpikes(mol)
        if self.checkCircularBonding(mol) == False:
            print ("Stucture unstable \n")
            return False
        return True
    
    def checkCircularBonding (self,molecule):
        """ This function takes an array of RBNs which are bonded to each other  in a linear fashion, it iterates through the RBNs and checks
            that the bond is still stable, if it is not false is returned otherwise true is returned
        """
        stable = True
        
        for i in range(np.size(molecule)):
            for j in range(np.size(molecule[i].spikeArray)):
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
    #                        print  ("Spike 1 type is: " + str(molecule[i].spikeArray[j].type) + "\n")
    #                        print ("Spike 1 intensity is: " + str(molecule[i].spikeArray[j].intensity) + "\n")
    #                        print ("Spike 2 type is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].type) + "\n")
    #                        print ("Spike2 intensity is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity) + "\n")
    #                        print ("Absolute intensity is: " + str( abs(molecule[i].spikeArray[j].intensity + molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity)) + "\n")
    #                        print ("Bond unstable \n")
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
    #                        print  ("Spike 1 type is: " + str(molecule[i].spikeArray[j].type) + "\n")
    #                        print ("Spike 1 intensity is: " + str(molecule[i].spikeArray[j].intensity) + "\n")
    #                        print ("Spike 2 type is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].type) + "\n")
    #                        print ("Spike2 intensity is: " + str(molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity) + "\n")
    #                        print ("Absolute intensity is: " + str( abs(molecule[i].spikeArray[j].intensity + molecule[i].spikeArray[j].bondedRBN.spikeArray[molecule[i].spikeArray[j].bondedSpikeNum].intensity)) + "\n")
    #                        print ("Bond unstable \n")
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
    
    
    
    def calculateIntensitySpikes (self,mol):
        """ This function calculates the intensity of every spike in the molecule, it does this by updating the molecule
            and then calculating the intensity of evey spike
        """
        
        origStates = []
        for i in range(len(mol)):
            origStates.append(mol[i].states)
            mol[i].zeroRBN()

            #print ("The state matrix for atom " + str(i) + " is: \n" + str(mol[i].states) + "\n")
        for i in range(self.maxSizeAtom + 30):
            self.updateMolecule(mol)

        
        self.analyseMolecule(mol)   
       
        for i in range(len(mol)):
            mol[i].setState(origStates[i],np.size(origStates[i],0)-1)
        
    def calculateIntensitySpikesDebug (self,mol):
        """ This function calculates the intensity of every spike in the molecule, it does this by updating the molecule
            and then calculating the intensity of evey spike
        """
        
        origStates = []
        for i in range(len(mol)):
            origStates.append(mol[i].states)
            mol[i].zeroRBN()
            #print ("The state matrix for atom " + str(i) + " is: \n" + str(mol[i].states) + "\n")
        for i in range(self.maxSizeAtom + 30):
            self.updateMolecule(mol)
            
        #for i in range(len(mol)):    
            #print ("The state matrix for atom " + str(i) + " is now: \n" + str(mol[i].states) + "\n")
        self.analyseMolecule(mol)   
       
        for i in range(len(mol)):
            mol[i].setState(origStates[i],np.size(origStates[i],0)-1) 
    
    
    def analyseMoleculeDebug (self,molecule):
        """ This function goes through every Atpm in a molecule and calculates the intensity of the spikes, this is needed
            as bonding and unbonding causes spike intensity values to change
        """
        for i in range(len(molecule)):
            self.analyseAtom(molecule[i])
    
    def analyseAtomDebug (self,rbn):
        """ Recalculculates intensity of a atom """
        for i in range(np.size(rbn.spikeArray)):
            rbn.spikeArray[i].calcMolIntenistyDebug()
    
    def incrementMoleculeState (self,mol):
        """ This function is called after ever time step, it updates each molecule and then sets the state matrix to be the new state
            this is done to prevent huge state matrices from being generted in every RBN which uses too much memory
        """
        self.updateMolecule(mol)
        for i in range(len(mol)):
            mol[i].selectMostRecentState()
    
    
    def debugHeadTail (self,head,tail):
        headNumUnbonded = np.size(head.spikeArray) - np.size(head.activeSpikes)
        print ("For head number of unbonded spikes: " + str(headNumUnbonded) + "\n")
        for i in range(np.size(head.spikeArray)):
            if head.spikeArray[i].bonded == False:
                print ("The intensity of unbonded spike is: " + str(head.spikeArray[i].intensity) + "\n")
            if headNumUnbonded == 0:
                print ("Head has no unbonded spikes \n")
                print ("Length of head spike array is: " + str(np.size(head.spikeArray)) + "\n")
        tailNumUnbonded = np.size(tail.spikeArray) - np.size(tail.activeSpikes)
        print ("For tail number of unbonded spikes: " + str(tailNumUnbonded) + "\n")
        for i in range(np.size(tail.spikeArray)):
            if tail.spikeArray[i].bonded == False:
                print ("The intensity of unbonded spike is: " + str(tail.spikeArray[i].intensity) + "\n")
            if tailNumUnbonded == 0:
                print ("Tail has no unbonded spikes \n")
                print ("Length of tail spike array is: " + str(np.size(tail.spikeArray)) + "\n")
                
    def printNumberUnbondedSpikes(self,molecule):
        """ This is a debugging function used to the print the number of unbonded spikes
            in each rbn in a molecule
        """
        for i in range(len(molecule)):
            print ("For RBN: " + str(molecule[i].rbnNumber) + " the number of unbonded spikes is: " + str(np.size(molecule[i].spikeArray) - np.size(molecule[i].activeSpikes)) + "\n")
    
    
    def checkStability (self,mol):
        if self.checkStructure(mol) == True:
            print ("Molecule is initially stable \n")
        else:
            print ("Molecule is inititally unstable \n")
        self.updateMolecule(mol)
        if self.checkStructure(mol) == True:
            print ("Molecule is still stable \n")
            
        self.updateMolecule(mol)
        if self.checkStructure(mol) == True:
            print ("Molecule is still stable after 2 updates \n")
        
        self.updateMolecule(mol)
        if self.checkStructure(mol) == True:
            print ("Molecule is still stable after 3 updates \n")
            
    def checkStructureDebug (self,mol):
        """ To check the structure of the molecule we update it many times, then calculate the
            intensity of each spike before checking that the mol is stable, if it is not then False
            is returned if it is then true is returned
        """
        self.calculateIntensitySpikes(mol)
        if bp.checkCircularBondingDebug(mol,self) == False:
            print ("Stucture unstable \n")
            return False
        return True