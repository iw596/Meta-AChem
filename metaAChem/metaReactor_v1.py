# -*- coding: utf-8 -*-
"""
Created on Thu Aug 23 16:29:44 2018
NOTE THIS REACTOR IS NOT CLOSE TO BEING FULLY COMPLETED, ALL IT CAN DO IS REACT TWO META ATOMS TOGETHER

This class defines the metaReactor, this reactor is used to bond metaatoms to generate metarings. To do this a soup of metaatoms
is generated and reacted for a certain number of time steps or until a metaring is formed.
@author: iw596
"""
import os.path
import random
import pickle
import metaAtomGenerator_v2 as metGen
import MetaMolecule_v1 as MetaMolecule
import metaBondingProtocol_v2 as mBP
class MetaReactor:
    
    def __init__ (self,numProducts,soupSize):
        self.wantedNumProducts = numProducts
        self.sizeSoup = soupSize
        self.generateSoup()
        
    def generateSoup (self):
        """ This function generates a soup of metaAtoms, it does this by openeing a random file containing a set
            of metaAtoms then selecting a random metAtom from this set and adding it to the soup, this is done until
            the size of the soup is the wanted size
        """
        self.soup = [] # Array to store soup
        # Store path where rings and supporting files are located
        path = 'H:\\YCCSA\\28-08-2018\\ringMolecules' # NEED TO CHANGE DEPENDING ON LOCATION
        # Need to find how many files there are, we do this by opening a supporting file containing this info
        numSetsFileName = "numSets" # File which stores the number of sets
        numSetsFileName = os.path.join(path,numSetsFileName+".txt")
        numSetsFile = open(numSetsFileName,"r+") # Open the file in read mode
        numSetsFile.seek(0) # Go to start of file
        numFiles = int(numSetsFile.read()) # Read the number and convert to an integer
        numSetsFile.close()
        curSize = 0 
        while len(self.soup) != self.sizeSoup:
            # Need to randomly select which set to open
            setNumber = random.randint(0,numFiles-1)
            # Next to load this file and depickle its contents
            setRingsFileName = "rings" + str(setNumber) + ".pickle"
            setRingsFileName = os.path.join(path,setRingsFileName)
            setRingsFile = open(setRingsFileName,"rb") # Open File
            setRings = pickle.load(setRingsFile) # Load the set of rings
            setRingsFile.close()
            ringNumber = random.randint(0,len(setRings)-1) # See how many rings there are and select which ring to pick
            ring = setRings[ringNumber] # Select the ring
            setRings = [] # Delete the set of rings
            
            # Next need to convert ring into metaatom
            metaAtom = metGen.atomGenerator(curSize,ring)
            if len(metaAtom.metaSpikes) == 2:
            #    print ("MetaAtom state is " + str(metaAtom.state) + "\n")
                self.soup.append(metaAtom)
                curSize += 1
        
    def runReactor (self):
        currentNumProducts = 0
        setMetaRings = []
        while currentNumProducts != self.wantedNumProducts:
            metaRing = []
            metaRing = self.generateMetaRing
            
        return 1
    def bondMetaAtomsTest (self,metaAtom1,metaAtom2):
        """ This is a test function which takes two metaAtoms and attempts to bond them by calling functions located in
            the bonding protocol script
        """
        metaMol = MetaMolecule.MetaMolecule(20)
        metaMol.addMetaAtom(metaAtom1)
        mBP.bondMetaAtoms(metaAtom1,metaAtom2,metaMol)
        
    
    def updateMetaAtoms (self):
        stateTransitions = []
        stateTransitions.append(self.soup[0].state)
        for i in range(1):
            self.soup[0].updateMetaAtom()
            stateTransitions.append(self.soup[0].state)
        
        print ("The states are:\n " + str(stateTransitions) + "\n")
 
    
    def generateMetaRing (self):
        """ This function generates a ring of metaAtoms, it does this by generating an intial pair of bonded metaAtoms and then
            reacting it with metaAtoms in the soup in order to generate a ring of metaAtoms where all metaSpikes are bonded
        """
        
        return 1 
            
    def generateInitialPair (self):
        """ This function is used to generate an initial meta molecule from two meta atoms it does this by selecting meta atoms
            until it finds two which can react together sucessfully
        """
        # This function is not complete, need to fully complete bonding code + complete a way of reforming meta atom if
        # they don't initally bond
        
        # Select two metaAtoms at random to bond
        metaAtom1 = self.selectMetaAtom()
        metaAtom2 = self.selectMetaAtom()
        # Generate a metaMolecule
        metaMol = MetaMolecule.MetaMolecule(20)
        # Add first metaAtom to metaMolecule
        metaMol.addMetaAtom(metaAtom1)

        # Attempt to react second MetaAtom with first MetaAtom
                
        return 1    
        #self.soup.  
    
    def selectMetaAtom (self):
        """ This function returns a metaAtom from the soup, it does this by randomly shuffling the soup and popping a metaAtom from
            it, note this function removes the metaAtom from the soup
        """
        random.shuffle(self.soup)
        return self.soup.pop()

    def debugMetaMolecule (self):
        """ This a debugging function which creates a metamolecule and runs one of its debugging functions """
        metaAtom = self.selectMetaAtom()
        metaMol = MetaMolecule.MetaMolecule(20)
        metaMol.addMetaAtom(metaAtom)
        metaMol.stabilityDebugTest()

reactor = MetaReactor (5,2)
#reactor.debugMetaMolecule()
reactor.bondMetaAtomsTest(reactor.selectMetaAtom(),reactor.selectMetaAtom())
#reactor.updateMetaAtoms()
            