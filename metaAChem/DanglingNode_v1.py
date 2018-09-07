# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 10:19:23 2018

@author: iw596
"""
import numpy as np
class DanglingNode:
    
    def __init__ (self,node,spike):
        self.bonded = False
       
        self.node= node
        self.spike = spike
        self.connectedNode = 0
        for i in range (np.size(self.spike.nodeList)):
            if self.spike.nodeList[i] == self.node:
                self.connectedNode = self.spike.nodeList[i+1]
                break
        
        self.bondedSpike = 0 
        # Stores the connection of the node in spike nodelist which this node takes input from this is the node which will be swapped
        # if dangling node is involved in bonding
        
        self.state = node.state
       # print ("The boolean function is: " + str(self.boolFunc) + "\n")
        
    def findState (self):
        self.state = self.node.state
    
    def bondFormed (self,bondedNode,bondedSpike):
        """ This function bonds the dangling node to another node, this changes the connection list of the
            node and sets the bonded status of the dangling node to be true.
        """
        self.connectedNode = bondedNode.node
        self.bondedSpike = bondedSpike
        self.bonded = True
     #   print ("Before bonding connected node rbn is: " + str(self.connectedNode.rbn.rbnNumber) + "\n")
        # Go through connection list in node and find the current connected node and replace with
        # new connection
        for i in range(len(self.node.connections)):
            if self.node.connections[i] == self.connectedNode:
                # Change connection
                self.node.connections[i] = bondedNode.node
                # Change bonded status of node
                self.node.bonded = True
        self.connectedDanglingNode = bondedNode
        self.spike.RBN.addBondedRBN (bondedSpike.RBN)
        
  #      print ("After bonding connected node rbn is: " + str(self.connectedNode.node.rbn.rbnNumber) + "\n")
        
    def bondBroken (self):
        """ Function called if bond between two metaspikes is unstable and needs to broken,this function works by
            removing the dangling nodes connection the other dangling node in the bonded metaspike and replacing it
            with a connection to the node which comes after the dangling node in the danglings nodes spike in atom
        """
        self.bonded = False
        tempNode = self.connectedNode # Temporaly stores a node so a swap can take place
        
        # First need to find node which node the dangling node should bond to once bond between metaspikes is broken
        for i in range (np.size(self.spike.nodeList)):
            if self.spike.nodeList[i] == self.node:
                self.connectedNode = self.spike.nodeList[i+1]
                break
       
        # Next need to update the connection list of the node
        for i in range (len(self.node.connections)):
            if self.node.connections[i] == tempNode:
                self.node.connections[i] = self.connectedNode
                
        # Finnaly we need to remove bonded node from RBN and remove bonded spike
 #       print ("The spike is: " + str(self.spike) + "\n")
 #       print ("The spike bonded RBN is: " + str(self.spike) + "\n")
        self.spike.RBN.removeBondedRBN(self.bondedSpike.RBN)
        self.bondedSpike = 0
    
    def changeState (self,newState):
#        print ("Old state is: " + str(self.node.state) + "\n")
        self.node.state = newState
#        print ("New state is: " + str(self.node.state) + "\n")
        self.state = self.node.state
        