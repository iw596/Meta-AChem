# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:35:10 2018
--- Need to review this and stip out unused code --- 


This class defines the metaspikes, these are made up of dangling nodes which are nodes which are not involved in the bonding
between two RBNs. MetaSpikes come in two types, the first is made up of danling ndoes obtained from bonding between two spikes which results
in only a singal danglind node being present. The second type is made of chains of dangling nodes obtained from bonding between two spikes which
results in multiple dangling nodes being present
@author: iw596
"""
import DanglingNode_v1 as dn
import metaSpikeInt_v1 as spkint
class MetaSpike:
    
    def __init__ (self,typeSpike,spikeNumber):
        """ This function initalises the class with the type of spike the MetaSpike will be, it also creates a empty lsit to store
            dangling nodes and sets the bonding status of metaspike to be false
        """
        
        self.bonded = False
        self.typeSpike  = typeSpike # What type of metaspike this instance will be
        # If type 1 spike we will have an array of individual dnagling nodes
        if self.typeSpike == 1:
            self.danglingNodeList = [] # List of dangling nodes which make up metaspike
            
        # If type 2 spike will store a list of arrays of dangling bonds
        else:
            self.danglingTailList = []
            
        self.intensity = 0
        self.spikeNumber  = spikeNumber
        self.metaAtom = 0
    
    def addMetaAtom (self,metaAtom):
        self.metaAtom = metaAtom
   
    def addDanglingNode (self,danglingNode):
        """ This function appends dangling node to the list """
        
        self.danglingNodeList.append(danglingNode)
        if danglingNode.bonded == True:
            self.bonded == True
    def addTailDanglingBonds (self,danglingNodeTail):
        """ This function appends list of dangling nodes to danglingTailList """
        self.danglingTailList.append(danglingNodeTail)
    
    def calculateIntensity (self):
        self.intensity = spkint.calculateIntensity(self)
        
    def bonded (self,bondedMetaSpike):
        self.bonded = True
        self.bondedMetaSpike = bondedMetaSpike
    
    def updateMetaSpike (self):
        if self.typeSpike == 1:
            # Need to updat each node simultaneously
            newStateDanglingNodes = []
            for i in range(len(self.danglingNodeList)):
                origState = self.danglingNodeList[i].state
                newStateDanglingNodes.append(self.danglingNodeList[i].calculateState())
                self.danglingNodeList[i].state = origState
            for i in range(len(self.danglingNodeList)):
                self.danglingNodeList[i].state = newStateDanglingNodes[i]        
        else:
            newStateDanglingTales = []
            for i in range(len(self.danglingTailList)):
                newStateDanglingNodes = []
                for j in range(len(self.danglingTailList[i])):
                    origState = self.danglingTailList[i][j].state
                    newStateDanglingNodes.append(self.danglingTailList[i][j].calculateState())
                    self.danglingTailList[i][j].state = origState
                newStateDanglingTales.append(newStateDanglingNodes)
            
            for i in range(len(newStateDanglingTales)):
                for j in range(len(newStateDanglingTales[i])):
                    self.danglingTailList[i][j].state = newStateDanglingTales[i][j]

    def changeState (self,states):
        """ This function changes the state of each node in the metaspike for a type 1 spike
            we simply iterate across the list of nodes and change the states. For a type 2 spike
            we decompose the dangling tail list into a list of nodes before iterating across each node
            and changing the states again.
        """
        
        if self.typeSpike == 1:
            for i in range(len(states)):
                self.danglingNodeList[i].state = states[i]
        else:
            return 2
        
    def returnStates (self):
        states = []
        
        if self.typeSpike == 1:
            for i in range(len(self.danglingNodeList)):
                states.append(self.danglingNodeList[i].state)
        
        else:
            for i in range(len(self.danglingTailList)):
                for j in range(len(self.danglingTailList[i])):
                    states.append(self.danglingTailList[i][j].state)
        
        return states
    
    def debugIntensity (self):
        """ This debugging function goes through every dangling node in the metaspike and prints off the intesnity of the spike
            it is associated with and the rbn Number of the RBN it is associated with
        """
        intens = []        
        rNum = []
        if self.typeSpike == 1:
            for i in  range(len(self.danglingNodeList)):
              #  print ("Single dangling nodes \n")
                intens.append(self.danglingNodeList[i].spike.intensity)
                rNum.append(self.danglingNodeList[i].spike.RBN.rbnNumber)
        else:
            for i in range(len(self.danglingTailList)):
             #   print ("Dangling tail nodes \n")
           #     print ("The size of the dangling tail is: " + str(len(self.danglingTailList[i].nodeList)) + "\n")
                intens.append(self.danglingTailList[i].spike.intensity)
                rNum.append(self.danglingTailList[i].rbnNumber)
            
    #    print ("The RBN numbers are: \n" + str(rNum) + "\n")
    #    print ("The intensities are: \n" + str(intens) + "\n")

            