# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 14:46:34 2018
This script contains the dangling tail class, a dangling tail is a collection of dangling nodes which belong to the same
spike, a dangling tail contains all of these nodes along with extra information like the bonding status of the tail and 
the spike it belongs too
@author: isaac
"""

class danglingTail:
    
    def __init__ (self,spike):
        self.nodeList = []
        #self.spike = spike
        self.bonded = False # Indicates if the tail is bonded
        self.fullyBonded = False # Indicates if all the nodes in the tail are fully bonded
        self.bondedTail = 0
        self.spike = 0
        self.rbnNumber = 0
    
    def addNode (self,node):
        """ This function adds a dangling node to the tail """
        self.nodeList.append(node)
        self.spike = node.spike
        self.rbnNumber = node.spike.RBN.rbnNumber
    
    def bondFormed (self,bondedTail,numBonds):
        """ This function is called when a dangling tail is bonded, it sets the bonding status of the tail to true and stores
            the dangling tail that it is bonded to. It also takes the number of bonds which is then used to determine
            if the tail is fully bonded
        """
        self.bonded = True
        self.bondedTail = bondedTail
        
        if numBonds == len(self.nodeList):
            self.fullyBonded = True
    
    def bondBroken (self):
        """ This function is used to unbond a dangling tail, it does this by going through each bonded dangling node and
            unbonding it.
        """
        self.bonded = False # Set bonding status to false
        
        

        