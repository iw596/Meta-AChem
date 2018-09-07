""" Class which describes Node object which is used in the sructure of RBN
    node consists of node number, RBN  it is part of, current state and connections
"""
from numpy import *
import pickle

#from random import randint

class Node:
    """ A node is a building block of an RBN it takes k number of connections from other nodes
        and has a boolean function which determines how the state of the node changes
        in response to the state of its inputs
    """
    
    

    def __init__ (self,nodeNumber,rbn,boolFunc,numConnections):
        """ This method initalises the node object with its node number, rbnNumber
            and function
        """
        self.nodeNumber = nodeNumber
        self.rbn = rbn # RBN node is part of
        self.state = random.randint(0,2)
        self.boolFunc = boolFunc
        self.connections = array([],dtype = Node)
        
        self.numConnections = numConnections
        self.bonded = False # This is triffered if the RBN is involved in a bond
        

    def addConnection(self,inputNode):
        """ This method adds a input to the node it takes the
            input node number and the rbn number the node is part of
        """
        self.connections = append(self.connections,inputNode)  

    def calculateState (self):
        """ This function is used to determine how to calculate the state of the node
            if the node is involved with a bond a special function has to be called in order to calculate
            the state, if it is not involved in a bond then a simpler function can be called
        """
        sumOfStates = 0
        
        for i in range(self.numConnections):
            #print ("The value of i is: " + str(i) + "\n")
            connectedNode = self.connections[i]
            power = 2**i
            #print ("The value of power is: " + str(power) + "\n")
            mostRecentState = connectedNode.state
            if mostRecentState != 0 and mostRecentState > 1:
                mostRecentState = 1
                connectedNode.rbn.fixStateMatrix()
            sumOfStates = (power * mostRecentState ) + sumOfStates
            
           # print ("The value of sum of states is: " + str(sumOfStates) + "\n")
       
        if sumOfStates >= size(self.boolFunc):
            print ("Error occuring \n")
            print ("The number of connections is: " + str(self.numConnections) + "\n")
            print ("The node RBN number is: " + str(self.rbn.rbnNumber) + "\n")
            print ("The power value is: " + str(power) + "\n")
            print ("Connections and their rbn numbers \n")
            for i in range(size(self.connections)):
                print ("The connected number is: " + str(self.connections[i].nodeNumber) + " the rbn is "  + str(self.connections[i].rbn.rbnNumber) + "\n")
            print ("The boolean function is:\n" + str(self.boolFunc) + "\n")
            print ("Actual number of connections is: " + str(size(self.connections)) + "\n")
            print ("RBN State matrix is: \n")
            self.rbn.printStateMatrix()
            
            pickle_out= open("ErroBond","wb")
            pickle.dump(self.rbn,pickle_out)
            pickle_out.close() 
        
        self.state = self.boolFunc[sumOfStates]
        return self.state
    
    def returnConnections (self):
        """ This method returns the matrix which shows the connections of the node and the rbns those connections
            come from
        """
        return self.connections

   
    def changeState (self,newState):

        if newState != 1 and newState != 0:
            print ("Error invalid state: " + str(newState) + "\n")
            newState  = 1
        
        self.state = newState
    
    def changeConnection (self,changedConnection,newConnection):
        """ This function changes the connectoon list of the new node, this means that the node will now
            be recieving an input from somewhere else 
        """
        # First have to find index of connection we are going to change
        for i in range(self.numConnections):
           if self.connections[i] == changedConnection:
               self.connections[i] = newConnection
               break


    def involvedInBond (self, changedConnection,newConnection):
        self.bonded = True
        #self.rbnArray = append(self.rbnArray,RBNOrigin) # Appends RBN node is part of to list
        #self.rbnArray = append(self.rbnArray,RBNNewConnection) # Appends RBN the node is bonded too to the list
        self.changeConnection (changedConnection,newConnection)
    
    def bondBroken (self,expectedNode):
        """ This function is called when the spike the node is in has a bond which becomes unstable, this function links the 
            the node up to another node passed in as an argument, this other node is located in the same RBN as this node. The nodes
            connection list is then changed
        """
        self.bonded = False

        # Search through bottom row of connection list until bonded RBN is found, as bond is breaking
        # we want to remove this connection to the other RBN

        for i in range(self.numConnections):
            if self.connections[i].rbn.rbnNumber != self.rbn.rbnNumber: # Search until connection to other RBN is found
                # Once found replace the connection with the new node
                self.connections[i] = expectedNode


    def printProps (self):
        """ This is a debugging functions and prints the connection list and node number of the node """
        print ("The node number is: " + str(self.nodeNumber) + "\n")
        conList = []
        conListStates = []
        conListRBN = []
        for i in range(self.numConnections):
            conList.append(self.connections[i].nodeNumber)
            conListRBN.append(self.connections[i].rbn.rbnNumber)
            conListStates.append(self.connections[i].state)        
        #print ("The connection is list: \n" + str(conList) + "\n")
        
        print ("The con list is: \n" + str(conList) + "\n")
        print ("The con list RBN: \n" + str(conListRBN) + "\n")
        print ("The con list states is: \n" + str(conListStates) + "\n")
        print ("The state of the node is: " + str(self.state) + "\n")
        print ("The boolean function is: " + str(self.boolFunc) + "\n")
        

