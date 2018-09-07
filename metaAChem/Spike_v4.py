# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 11:00:48 2018
This script contains the code to implement the spike class
@author: iw596
"""
from numpy import *
import Node_v4 as node
import spikeInt_v7 as spkInt # Module used to calculate intensity of spike and dangling bonds

class Spike:
    """ The spike class is used by the RBN to determine whether the RBN should bond to another RBN
        bonding occurs when two spikes interact and share links between nodes. Spike class contains data such as
        nodes in spikes, rbn spike is part of and whether the spike is involved in a bond or not
    """
   
    def __init__ (self,spikeNumber,RBN):
        self.nodeList = array([],dtype = node.Node) # Stores nodes in the spike
        self.bonded = False # Boolean to indicate if spike is involved in a bond or not initally all spikes are not bonded
        self.RBN = RBN # Stores the RBN the spike is part of
        
        self.spikeNumber = spikeNumber
        self.bondedRBN = 0
        self.checked = False # This boolean is used to indicate whether the spike has been checked when recalcualting intensity
        self.danglingBonds = [] # This stores dangling bonds in the spike, dangling bonds are bonds not involved in bonding
        self.numDanglingBonds = 0 # Stores number of dangling bonds
    
        # The type of node is dependent on the size of the spike
        self.type = 1
    
    
    
    def addNode (self,node):
        """ This function adds a node to the spike and adds the nodes
            flashiness to the flashiness array
        """
        self.nodeList = append(self.nodeList,node)
        # Check to see if type of spike can change
        if len(self.nodeList) >= 5 and len(self.nodeList) < 10:
            self.type = 2
        elif len(self.nodeList) >= 10:
            self.type = 3
            
        if self.bonded == True:
            print ("called when bonded \n")
    def setSpikeNum (self,spikeNum):
        """ This function sets the spike number """
        self.spikeNumber = spikeNum
    
    def returnSpkNum (self):
        return self.spikeNumber

    def hasBonded (self,RBN,bondedSpikeNum):
        """ This function is called when spike is involved in a bond, it sets the boolean variable
            indicating whether spike has bonded to true
        """
        self.bonded = True # Update bonding status of spike
        self.bondedRBN = RBN # Store RBN  spike is bonded too
        self.bondedSpikeNum = bondedSpikeNum # Stores the spike num the spike is bonded to
        # Update RBN that spike is part of with new bonding infor
        self.RBN.rbnBonded(self.spikeNumber,RBN)

    def addDanglingBonds (self,numDangleNodes):
        """ This function is caleld after a spike has bonded, it takes the number of dangling bonds the spike will have,dangling
            bonds are bonds  which has not been involved in the linking with another RBN, the function will then place the 
            first nodes in the node array into the dangling bond array as these nodes will be the ones which are not involved in 
            the bonding, the varaible containg number of bonds will also be updated
        """
        # Add nodes to dangle bond array
        #print ("The number of dangling nodes is: " + str(numDangleNodes) + "\n")
        #print ("The length of the spike is: " + str(size(self.nodeList)) + "\n")
        #print ("The length of the bonded spike is: " + str(size(self.bondedRBN.spikeArray[self.bondedSpikeNum].nodeList)) + "\n")
        self.danglingBonds= []
        for i in range(numDangleNodes):
            self.danglingBonds = append(self.danglingBonds,self.nodeList[i])
        
        self.numDanglingBonds = numDangleNodes # update number

    
    
    
    
    def bondBreak (self):
        """ This function is called when the bond between two bonds breaks, this function works
            by starting at the bottom index of the list of nodes and connecting it to the node
            below. Finally the bonding state of the node is set to false and the RBN is updated to
            reflect the fact that the bond has been broken
        """
        self.bonded = False # Set the bonding status of the bond to false as spike is now unbonded
        
        # Next go through each node reconnecting it to the node after it in order to 
        # reform the spike to the conenction list it had before the bond formed
        for i in range(size(self.nodeList)-1):
            self.nodeList[i].bondBroken (self.nodeList[i+1])
       # print ("Before In the function intensity is: " + str (self.intensity) + "\n")
        #print ("The state of the the RBN is: \n")
        #self.RBN.printMostRecentState()
        # Need to recalculate intensity
        #self.recalculateIntensity()
        
        #print ("after In the function intensity is: " + str (self.intensity) + "\n")
        
        # Need to update RBN to inform it that bond has broken
        self.RBN.rbnUnbonded(self.spikeNumber,self.bondedRBN)
        # Need to update dangling bonds
        self.danglingBonds = []
        self.numDanglingBonds = 0
        
        # Can remove reference to bonded RBN as no longer needed
        self.bondedRBN = 0
        self.bondedSpikeNum = -1 # Set to -1 as invalid spike num
            
    
    
    
    def returnNodeArray (self):
        return self.nodeList
    
    def calculateInitialIntensity(self):
        """ This function calculates the initial intensity of the spike
           it works by calling a function which works by 
           finding the cycle length of RBN (only stores states of nodes in spikes)
           the intensity is sum of weighted transistions, 0-1 transistion is +1, 0-1 transition is -1
        """
        
        spkInt.calculateIntensity(self) # See function script for more detail
        # print ("The intensity is: " + str(self.intensity) + "\n")
        self.intensity = spkInt.calculateIntensity(self) 
        return self.intensity       
    
    
    def recalculateIntensity (self):
        """ This function is used to recalculate the intensity of a spike after the spike has formed a bond 
            it works in the same way as the function above but updates the other RBN when findinging the cycle
        """
        self.intensity = spkInt.calculateIntensity(self) # See function script for more detail
        return self.intensity

    def calcMolIntenisty (self):
        self.intensity = spkInt.findMolecularIntensity(self)
    
    def calcMolIntenistyDebug (self):
        self.intensity = spkInt.findMolecularIntensityDebug(self)
    
    def returnIntensity (self):
         """ This function returns the intenisty of the spike """
         return self.intensity
    
    def changeNodeArray (self,newNodeArray):
        """ This function changes the nodes (or properties of the nodes), this is called when spike is involved in a bond """ 
        self.nodeList = newNodeArray
    
    def printNodeArray (self):
        print ("The node list is: \n")
        for i in range(size(self.nodeList)):
            print (str(self.nodeList[i].nodeNumber))
            print ("  ")
    
    def returnSize(self):
        """ Returns the number of nodes in the spike """
        return size(self.nodeList)
    
    def returnRBNNumber(self):
    
        return self.rbnNumber
    
    def returnBondedRBN (self):
        """ Returns the RBN the spike is bonded to """
        return self.bondedRBN
    
    def setBondedRBN (self,RBN):
        """ Sets the RBN the spike is bonded to """
        self.bondedRBN = RBN
    
    def checkedPrint (self):
        if self.checked == True:
            print ("Spike has been checked \n")
        else:
            print ("Spike has not been checked \n")
     
             
    def printNodeProps (self):
        """ This debugging function goes through every node in spikes and prints out the nodes state,function ,node number and 
            connection list (including which RBN connected node is from)
        """
        for i in range (size(self.nodeList)):
            self.nodeList[i].printProps()
    
