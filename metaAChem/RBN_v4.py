# -*- coding: utf-8 -*-
""" Class which describeds RBN object which consists of nodes and connections
    between nodes, class has methods which desrcibed the dynamics of the
    RBN
"""

import Node_v4 as node
import Spike_v4 as spike
from numpy import *
from scipy.stats import mode

class RBN:

    def __init__(self,numNodes,numConnections,rbnNumber):
        """ Method which initalises the RBN by assinging internal varaibles
            there value and calling a method which creates the internal
            structure of the RBN
        """
        self.n = numNodes
        self.k = numConnections
        self.rbnNumber = rbnNumber
        self.nodeArray = array([],dtype = node.Node) # Create array which can be filled with nodes
        self.spikeArray = array([],dtype = spike.Spike)
        self.bonded = False # Boolean used to indicate if RBN is bonded to another RBN
        self.bondedRBNs = []
        self.states = array([],dtype = int) # Matrix to hold states of each node
        self.numIterations = 0 # Number of times RBN has been run
        self.activeSpikes = array([],dtype = int) # Stores the spike numbers of the spikes which are currently invovled in bonds
        self.type = 0 # The type of RBN is determined by the number of spikes it has
        self.createRBN()
        self.generateSpikes()
        
    def test(self):
       """ A quick test function to check on what RBN looks like """
       print ("The RBN number is " + str(self.rbnNumber) + "\n " + "The num of nodes is: " + str(self.n) + "\n")
       for i in range(self.n):
           self.nodeArray[i].printProps()

    def createRBN (self):
        """ This method creates an RBN by generating an array of nodes and
            assigning each node its connectons to other nodes and its
            internal function
        """

        # First generate connection matrix
        con = apply_along_axis(random.permutation, 1, tile(range(self.n), (self.n,1)))[:, 0:self.k]
        #print ("The original connection matrix is: " + str(con) + "\n")

        #Next generate boolean function matrix which maps how node reacts to inputs
        booleanFuncs = random.randint(0,2,(self.n,2**self.k))

        # Fill array with approprate number of nodes and give each node its
        # connections
        for i in range(self.n):
        
            self.nodeArray = append(self.nodeArray,node.Node(i,self,booleanFuncs[i,],self.k))
        
        for i in range(self.n):
            # Generate a new node object
            #self.nodeArray = append(self.nodeArray,node.Node(i,self,booleanFuncs[i,],self.k))
            for j in range(self.k):
                # Add connections to newly generated node
                #print ("The node number to add is: " + str(con[i,j]) + "\n")
                self.nodeArray[i].addConnection(self.nodeArray[con[i,j]])
            self.states = append(self.states,self.nodeArray[i].state) # Add state of node to create initial state
    

    def resetRBN (self):
        """ This function resets an RBN by giving the state matrix of the RBN
            the initial value of the State matrix and then resetting number of iterations
            after this the state of the nodes is reset
        """
        
        if ndim(self.states) > 1:
            #print ("Original states is: \n" + str(self.states) + "\n" )
            self.states = self.states[0,]
            #print (" New states is: \n" + str(self.states) + "\n" )
            self.numIterations = 0
            
            for i in range(size(self.n)):
                self.nodeArray[i].changeState(self.states[i])
        else:
            #print ("Original states is: \n" + str(self.states) + "\n" )
            self.states = self.states
            #print (" New states is: \n" + str(self.states) + "\n" )
            self.numIterations = 0
            
            for i in range(size(self.n)):
                self.nodeArray[i].changeState(self.states[i])
    
     
    def findCycleLength (self):
        """ This calculates the cyclength of a RBN with starting conditions given by current state and connectons
            of each node in the RBN, this function is called after a bond occurs as cycle length could
            change
        """
#
        originalStateMatrix = self.states # Store the original matrix
        originalNumIteration = self.numIterations
        numAttempts = self.n
        numRBNUpdate = self.n + 30
        #print ("The number of attemps is: " + str(numAttempts) + "\n")
        for i in range(numAttempts):
            for j in range(i):
                self.updateRBN()
            self.selectMostRecentState()
            stateOfRBN = self.states
            # Run RBN for this number of time
            for k in range(1,numRBNUpdate):
                self.updateRBN()
                stateOfRBN = self.states
                #print ("The state of the RBN is now: " + str(spike.RBN.states) + "\n")
                #print ("The first row is: " + str(stateOfRBN[0,:]) + "\n")
                #print ("The last row is: " + str(stateOfRBN[k,:]) + "\n")
                if array_equal(stateOfRBN[0,:],stateOfRBN[k,:]) == True:
                    self.popState()
             #       print ("The cyclength is: " + str(spike.RBN.numIterations) + "\n")
                    # Need to pop last state
                   # for i in range(size(spike.nodeList)):
                       # print ("Node in spike: " + str(spike.nodeList[i].nodeNumber) + "\n")
                    #print ("The state matrix is: " + str(self.states) + "\n")
                    count = k                    
                    self.setState(originalStateMatrix,originalNumIteration)
                    return count
            self.zeroRBN()
        
        return -1
        



    def runRBN (self,numTimeSteps,isBonded):
        """ This function uses a for loop to run the RBN for a given number
            of time steps passed to it as an argument
        """
        numSpikes = size(self.activeSpikes)
        originalRBNs  = []
        if isBonded == True:
            for i in range(numSpikes):
                originalRBNs = append(originalRBNs,self.activeSpikes[i].returnBondedRBN())
            runList = originalRBNs
        
        for i in range(numTimeSteps):
            
            if isBonded == True:
                # Need to have copies of RBN at each state
                for i in range(numSpikes):
                    self.activeSpike[i].setBondedRBN(runList[i])
                    runList[i].updateRBN()
                
                self.updateRBN()
                
                for i in range(numSpikes):
                    self.activeSpike[i].setBondedRBN(runList[i])
            else:
                self.updateRBN()
    
        # If the RBN is bonded after updating we need to reset the rbns of all the spikes
        if isBonded == True:
            for i in range(numSpikes):
                self.activeSpike[i].setBondedRBN(originalRBNs[i]) 
            
        
    
    def updateRBN(self):
        """ This function works by calling the update state function for each
            node then appending thew new state of each node to a new
            row in the state matrix
        """
        # Generate new array to store the new state of node
        newStates  = empty([self.n],dtype = int)
        rowValue = []
        
        for i in range (self.n):
            origState = self.nodeArray[i].state
            self.nodeArray[i].calculateState()
            #print ("Return state is: " + str(self.nodeArray[i].returnState()) + "\n")
            newStates[i] = self.nodeArray[i].state
            self.nodeArray[i].state = origState
        
        self.numIterations += 1
        # Append new states to state matrix
        self.states = vstack((self.states,newStates))
        
        # Update nodes with new state
        for i in range(self.n):
            self.nodeArray[i].state = newStates[i]
      

    def generateSpikes(self):
        """ This function generates the spikes that the RBN uses to bond with
            other RBNS, the function is split into two smaller functions
            one calculates the nodes in the spike and the other
            the properties of the spike
        """
        # Generates the spikes 
        self.findSpikeNodes()
        # Finally spikes with zero intensity are removed as they are considered 'inert' and give each remaing spike a fixed spike number 
        self.intensityOfSpikes()
        self.removeZeroIntensitySpikes ()
        for i in range(size(self.spikeArray)):
            self.type += 1
            self.spikeArray[i].setSpikeNum(i)
    
        #self.resetRBN()
    
    
    def findSpikeNodes(self):
        """ This function finds the number of spikes and the nodes in them
            the algorithm to do this is described throughout this function 
        """
        # First generate set of nodes by randomly sorting the list of node numbers
        setOfNodes = arange(self.n)
    
        random.shuffle(setOfNodes)
        
        # While the set of nodes is not empty 
        while size(setOfNodes) != 0:
            # Next we select a node from this list
            node = setOfNodes[0]
            # Remove node from set
            setOfNodes = delete(setOfNodes,[0])    
            # Generate a new spike
            self.spikeArray = append(self.spikeArray,spike.Spike(0,self))
            # Append node to spike
            self.spikeArray[size(self.spikeArray)-1].addNode(self.nodeArray[node])
            
            # Get list of node connections, note we can discard second row as all nodes will be bonded to this RBN
            inputList = self.nodeArray[node].returnConnections()
            # While there are still nodes in the input list
            while (size(inputList) != 0):
                # The next node is randomly selected from the input lis
                nextNodeIndex = random.randint(0,size(inputList))
                nextNode = inputList[nextNodeIndex]
                # The next node is then removed from the input list
                inputList = delete(inputList,nextNodeIndex)
          
                
                if nextNode.nodeNumber in setOfNodes:
                    indexToDelete = int(where(setOfNodes  == nextNode.nodeNumber)[0])
                    setOfNodes = delete(setOfNodes,indexToDelete)
                    self.spikeArray[size(self.spikeArray)-1].addNode(self.nodeArray[nextNode.nodeNumber])
                    node = nextNode
                    inputList = self.nodeArray[node.nodeNumber].returnConnections()
        
        
            
        
        
    def intensityOfSpikes(self):
        """ This function is used to calculate the intensity of each spike it works by going through each spike
            and calling a funcion in the spike class which calculates the intensity. Note the flashiness of each node is stored
            is the spike class so the only data that needs to be passed to it is the cycle length
        """
        for i in range(size(self.spikeArray)):
            self.spikeArray[i].calculateInitialIntensity()
    
    def removeZeroIntensitySpikes (self):
        """ This function scans through the list of spikes and removes the ones with zero intensity, this is done
            as zero intensity is considered an 'inert' spike which can't bond with anything
            also removes spikes with only one node as they cannot swap links with other spikes
        """
        spikesToBeDeleted = array([],dtype = int)
        for i in range(size(self.spikeArray)):
            if self.spikeArray[i].returnSize() == 1 or self.spikeArray[i].intensity == "a":
                spikesToBeDeleted = append(spikesToBeDeleted,i)
        for i in range(size(spikesToBeDeleted)-1,-1,-1):
            self.spikeArray = delete(self.spikeArray,spikesToBeDeleted[i])
   
    
    def printStateMatrix (self):
        """ This method is for debugging, prints state matrix """
        print("The state matrix is: \n" + str(self.states) + "\n")
    
    def printMostRecentState (self):
        """ This method prints the most recent state of RBN """
        if self.numIterations > 0:
            print("The most recent state is: \n" + str(self.states[self.numIterations,]) + "\n")
        else:
            print ("The most recent state is: \n" + str(self.states) + "\n")
    
    
    def returnMostRecentNodeState(self,nodeNumber):

        if self.numIterations == 0:
            #print ("In the IF Statement \n")
            #print ("The node numeber is: " + str(nodeNumber) + "\n")
            #print ("The state matrix is: \n" + str(self.states) + "\n")

            return self.states[nodeNumber]
        else:
            #print ("In else statemtn \n")
            return self.states[self.numIterations,nodeNumber]
    
    def returnNodCons(self):
        for i in range(self.n):
            self.nodeArray[i].returnConnections()
    
    def returnNumSpikes(self):
        """ Returns the number of bonding spikes the RBN has """
        return size(self.spikeArray)
    
    def addBondedSpike (self,spikeNum):
        """This function adds an active spike to the array """
        self.activeSpikes = append(self.activeSpikes,spikeNum)
    
    def returnSpikeArray (self):
        """ This function returns all the spikes asssocaited with this RBN """
        return self.spikeArray
    
    def returnRBNNumber(self):
        """ This returns the RBN number associated with this RBN """
        return self.rbnNumber

    def rbnBonded (self,spikeNum,bondedRBN):
        """ This function adds a spike to the list of spikes involved in a bond """
        self.bonded = True
        self.activeSpikes = append(self.activeSpikes,spikeNum)
        self.bondedRBNs.append(bondedRBN)
    
    def returnNumberNodes(self):
        """ This function returns the number of nodes the RBN is made up of """
        return self.n

    def rbnUnbonded (self,spikeNum,bondedRBN):
        """ This function removes a spike from the list of spikes involved in a bond
            and if there are no more spikes in the active spike array then the state of
            the rbn is set to unbonded
        """
        # Need to remove spike num from list of active spikes
        #print ("Before spike removal: \n" + str(self.activeSpikes) + "\n")
        self.activeSpikes = setdiff1d(self.activeSpikes,spikeNum) # This numpy function will remove spikeNUm
        #print ("After spike removal: \n" + str(self.activeSpikes) + "\n")
        
        if size(self.activeSpikes) == 0:
            self.bonded = False
        
        for i in range(len(self.bondedRBNs)):
            if self.bondedRBNs[i] == bondedRBN:
                self.bondedRBNs.pop(i)
                break
    def returnNumIters (self):
        """ Returns the number of times the RBN has beeen run """
        return self.numIterations
    
    def popState (self):
        """ This function removes state given by iterNumber and returns the state value
            and decrements the number of iterations
        """
        #print ("The number of rows is: " + str(size(self.states,0)) + "\n")
        #print ("The number of iterations is: " + str(self.numIterations) + "\n")
        state = self.states[size(self.states,0) -1,]
        self.states = delete(self.states,size(self.states,0) -1 ,0)
        self.numIterations -= 1
        #print ("Number of iterations is: " + str(self.numIterations) + "\n")
        if self.numIterations == 0:
            self.states = self.states.flatten()
            for i in range (self.n):
                self.nodeArray[i].changeState(self.states[i])
        else:
            for i in range (self.n):
                self.nodeArray[i].changeState(self.states[size(self.states,0) -1,i])
        #print ("The state to be returned is: " + str(state) + "\n")
        return state 
    
    def appendState (self,state):
        """ This function appends a state passed in as an argument and 
            increments number of states and updates node values
        """
        self.states = vstack((self.states,state))
        
        #print ("The number of iterations is: " + str(self.numIterations) + "\n")
        for i in range (self.n):
            if self.states[size(self.states,0) -1,i] != 0 and self.states[size(self.states,0) -1,i] != 1:
                print ("Error the most recent state is \n"  + str(self.states) + "\n")
            self.nodeArray[i].changeState(self.states[size(self.states,0) -1,i])
        self.numIterations += size(self.states,0) -1
    def returnMostRecentRow (self):
        if self.numIterations == 0:
            return self.states
        else:
            return self.states[size(self.states,0) -1,]
    
    
    def returnStateMatrix (self):
        """ This returns the state matrix associated with this RBN """
        return self.states

    def setState (self,stateMatrix,numIterations):
        #print ("New state is: " + str(stateMatrix) + "\n")
        #print ("New num iters is: " + str(numIterations) + "\n")
        self.numIterations = numIterations
        self.states = stateMatrix
        if self.states.ndim == 1:
            for i in range (self.n):
                self.nodeArray[i].changeState(stateMatrix[i])
        else:
            for i in range (self.n):
                #print ("States is: " + str(self.states) + "\n")
                #print ("Rows: " + str(size(self.states,0)) + "\n")
                self.nodeArray[i].changeState(stateMatrix[(size(self.states,0)-1),i])
        self.numIterations = size(self.states,0)-1
    def selectMostRecentState (self):
        """ This function resets the RBN, with the starting state being equal
            to the final state before the reset
        """
        
        
        #print ("Number of iteratiosn before reset is: " + str(self.numIterations) + "\n") 
        #print ("The rbn number is: " + str(self.rbnNumber) + "\n")
        self.numIterations = 0
        if self.states.ndim == 1:
            for i in range(self.n):
             self.nodeArray[i].state = self.states[i]
           # print ("In this area \n")
            return
        else:
            self.states = self.states[size(self.states,0)-1,]
        #print ("The state matrix is: "  + str(self.states) + "\n")
        
        for i in range(size(self.states)):
            self.nodeArray[i].state = self.states[i]

    def fixStateMatrix(self):
        """ This is a quick fix for a bug where states become corrupted with integers besides zero and one
            need to find root cause of this bug at some point
        """
        if self.numIterations == 0:
            for i in range(size(self.states)):
                if self.states[i] > 1:
                    self.states[i] = 1
        else:
            for i in range(self.numIterations):
                for j in range(self.n):
                    if self.states[i,j] > 1:
                        self.states[i,j] = 1

    def zeroRBN (self):
        self.states = array([],dtype = int)
        self.numIterations = 0
        for i in range(self.n):
            self.states = append(self.states,0)
        for i in range(self.n):
            self.nodeArray[i].state = 0
    def addBondedRBN (self,bondedRBN):
        """ This function adds an RBN to the bonded RBN list this is used when metaspikes bond which causes two RBNs in separate
            metaAtoms to bond together
        """
        self.bondedRBNs.append(bondedRBN)
    
    def removeBondedRBN (self,bondedRBN):
        """ This function removes an RBN from the bonded RBN list this is used when two metaspikes break which causes two RBNs
            in separate metaAtoms to break appart
        """
        
        for i in range(len(self.bondedRBNs)):
            if self.bondedRBNs[i] == bondedRBN:
                self.bondedRBNs.pop(i)
                break
    






#testRBN = RBN(5,2,0)
#for i in range(15):
#    testRBN.updateRBN()
#print (testRBN.states)


##testRBN.test()
#print ("The number of spikes is: " + str(size(testRBN.spikeArray)) + "\n")
#for i in range(size(testRBN.spikeArray)):
#    print ("The intensity is: " + str(testRBN.spikeArray[i].intensity))
#
#for i in range(size(testRBN.spikeArray)):
#    testRBN.spikeArray[i].recalculateIntensity()
#    print ("The intensity is now: " + str(testRBN.spikeArray[i].intensity))
##
##testRBN.printStateMatrix()
#
#
#
#print ("The number of spikes is: " +str(size(testRBN.spikeArray)) + "\n")
#for i in range(size(testRBN.spikeArray)):
#    print ("For spike: " + str(testRBN.spikeArray[i].spikeNumber) + "\n")
#    print ("The intensity is: " + str(testRBN.spikeArray[i].intensity) + "\n")
#    for j in range(size(testRBN.spikeArray[i].nodeList)):
#        print ("The node number is: " + str(testRBN.spikeArray[i].nodeList[j].nodeNumber) + "\n")
#
#
#

#testRBN.test()
#testRBN.printStateMatrix()
#testRBN.updateRBN()
#testRBN.printStateMatrix()
#
#
#
#print("Cyle length is: " + str(testRBN.findCycleLength()[0]) + "\n")  
#        
