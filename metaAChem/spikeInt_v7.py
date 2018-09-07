# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 18:28:31 2018

@author: isaac
"""

from numpy import *
from scipy.stats import mode


def calculateIntensity (spike):
    # First need to store the original state of RBN Nodes
   # print ("The original states: " + str(spike.RBN.states) + "\n")

    
     intensity = 0
     rbnList = 0
     originalStateRBNNodes = spike.RBN.states
     origNumIters = spike.RBN.numIterations
     spike.RBN.zeroRBN()
      # print ("State RBN after reset is: " + str(spike.RBN.states) + "\n")
     intensity = findIntensity(spike)
     spike.RBN.setState(originalStateRBNNodes,origNumIters) 

        
        #spike.bondedRBN.setState(originalStateBRBNNodes,origNumItersBonded)
    #print ("before resettting states: " + str(originalStateRBNNodes) + "\n")
       
    #print ("After resettting states: " + str(originalStateRBNNodes) + "\n")
    #print ("Returns back to here \n")
     return intensity
        



def newFinalStage (states,spike):
    
    # Fist store list node number
    nodeNumbers = []
    for i in range(size(spike.nodeList)):
        nodeNumbers.append(spike.nodeList[i].nodeNumber)
    #print ("The node numbers are: " + str(nodeNumbers) + "\n")
    #print ("Number of columns: " + str(size(states,1)) + "\n")
    #print ("Number of rows: " + str(size(states,0)) + "\n")
   # print ("The states are \n " + str(states) + "\n")
    #print ("Number of dimenions is: " + str(states.ndim) + "\n")
    
    #First check that a cycle has been found if not then return a string to indicate failure
    if states is None:
        return 'a'
    
    if states.ndim == 1:
        states = atleast_2d(states)
    intensity = 0
    # Go through each column
    for i in range(size(states,1)):
        #Go through each row
        #print ("Column number " + str(i) +  "\n" )
        currentSum = 0 
        if i in nodeNumbers:
            #print ("Starting value for sum: " + str(currentSum) + "\n")
            for j in range(size(states,0)):
                #print ("State matrix is: " + str(states) + "\n")
                #print ("Row number " + str(j) +  "\n" )
                #print (" Element value is: " + str(states[j,i]) + "\n")
                if states[j,i] == 1:
                    currentSum += 1
                    #print ("Current sum value: " + str(currentSum) + "\n")
                else:
                    currentSum -= 1
                    #print ("Current sum value: " + str(currentSum) + "\n")
            
           # print ("Current sum final value: " + str(currentSum) + "\n")
            if currentSum == size(states,0):
                intensity += 1
            elif currentSum == (size(states,0) *-1):
                intensity -= 1
        
    #print ("The intensity is: " + str(intensity) + "\n")
    return intensity








def findMolecularAttractorCycle (spike):
    for i in range(size(spike.RBN.states,0)):
        for j in range(i+1,size(spike.RBN.states,0)):
            #print ("The state matrix is currently: " + str(spike.RBN.states) + "\n")
            #print ("First row: " + str(spike.RBN.states[i,:]) + "\n")
            #print ("Second row: " + str(spike.RBN.states[j,:]) + "\n")
            if array_equal(spike.RBN.states[i,:],spike.RBN.states[j,:]) == True:
                #print ("The spike number is: " + str(spike.spikeNumber) + "\n")
                #print ("The RBN number is: " + str(spike.RBN.rbnNumber) + "\n")
                #print ("Returning:\n" + str(spike.RBN.states[i:j,]) + "\n")
                return spike.RBN.states[i:j,]
        
 
def findMolecularIntensity (spike):
    stateNodes = findMolecularAttractorCycle (spike)
    intensity =  newFinalStage(stateNodes,spike)
    #print ("The intensity is:\n " + str(intensity) + "\n")
    return intensity
    
 
 

def findUnbondedAttractorCycle (spike):
    
        # We will run RBN for n + 30 times where n is number of nodes in RBM
   
    indexOfStart = 0  # Index when mode is first found
    indexOfEnd = -1 # Index when  mode is found again
    numAttempts = spike.RBN.n
    numRBNUpdate = spike.RBN.n + 30
    stateOfRBN = spike.RBN.states
    #print ("The number of attemps is: " + str(numAttempts) + "\n")
    for i in range(numAttempts):
        #print ("Initial state is: " + str(stateOfRBN) + "\n")
        for j in range(i):
            #print ("Inside for loop \n")
            spike.RBN.updateRBN()
        spike.RBN.selectMostRecentState()
        #print ("The states are now: " + str(spike.RBN.states) + "\n") 
        stateOfRBN = spike.RBN.states
        #print ("The other matrix has the value: " + str(stateOfRBN) + "\n")
        
        # Run RBN for this number of time
        for k in range(1,numRBNUpdate):
            spike.RBN.updateRBN()
            stateOfRBN = spike.RBN.states
            #print ("The state of the RBN is now: " + str(spike.RBN.states) + "\n")
            #print ("The first row is: " + str(stateOfRBN[0,:]) + "\n")
            #print ("The last row is: " + str(stateOfRBN[k,:]) + "\n")
            if array_equal(stateOfRBN[0,:],stateOfRBN[k,:]) == True:
                spike.RBN.popState()
                #print ("The cyclength is: " + str(spike.RBN.numIterations) + "\n")

            # Need to pop last state
               # for i in range(size(spike.nodeList)):
                   # print ("Node in spike: " + str(spike.nodeList[i].nodeNumber) + "\n")
               # print ("The states are: \n" + str(spike.RBN.states) + "\n")
                return spike.RBN.states
        
        spike.RBN.resetRBN()
        #print ("After fail states: " + str(spike.RBN.states) + " \n")
    

def findIntensity (spike):
        # First need to store the original state of RBN Nodes
    originalStateRBNNodes = []
    for i in range(spike.RBN.n):
        originalStateRBNNodes.append(spike.RBN.nodeArray[i].state)
    
    
    # Stores the state of nodes in each update
    stateNodes = zeros(spike.RBN.n,dtype = int)
    
    # Reset the RBN
    spike.RBN.resetRBN()
    
    #Add initial states to array
    for i in range(spike.RBN.n):
        stateNodes[i] = spike.RBN.states[i]

    stateNodes = findUnbondedAttractorCycle(spike)

        
    
    
 #   print ("After function state nodes is \n" + str(stateNodes) + "\n")



    #print ("The transposed array is: \n"  + str(trans) + "\n")

    intensity = newFinalStage(stateNodes,spike)
    #for i in range (len(spike.nodeList)):
     #   print ("Nodes in spike: " + str(spike.nodeList[i].nodeNumber) + "\n")
    #print ("The intensity is: " + str(intensity) + "\n")
    return intensity                
        
        
    
    
    
    
    trans = trans[indexOfStart:indexOfEnd]
    # Next need to convert integers to binary numbers to see individual nodes
    # Line below taken from
    # https://stackoverflow.com/questions/22227595/convert-integer-to-binary-array-with-suitable-padding
    trans = (((trans[:,None] & (1 << arange(size((danglingNodeList))))) > 0).astype(int))
    print ("After converting back to binary: " + str(trans) + "\n")
    
    intensity = 0
    sumNodes = array([],dtype = int)
    sumIndNode = 0  # sum of an individual node
    # If cyclength is one then all nodes are frozen by definition
       
    if size(trans,0) == 1:
        #print ("In cycle length one if statement \n")
        for i in range(size(trans,1)):
            #print ("The size is: " + str(size(trans,1)) + "\n")
            if trans[0,i] == 1:
                intensity += 1
            else:
                intensity -= 1
                    
    else:
        for i in range(size(trans,1)):
            for j in range(size(trans,0)):   
                   if trans[j,i] == 1:
                       sumIndNode += 1
            sumNodes = append(sumNodes,sumIndNode)
            sumIndNode = 0

        for i in range(size(sumNodes)):
            if sumNodes[i] == size(trans,0):
                intensity += 1
            
            elif sumNodes[i] == 0:
                intensity -= 1
    #print ("The rbn num is " + str(spike.RBN.rbnNumber) + "\n")
    #print ("The cycle length is: " + str(count) + "\n")
    #spike.intensity = intensity
    #print ("The intensity is: " + str(intensity) + "\n")
    return intensity        










def findMolecularAttractorCycleDebug (spike):
    print ("The starting matrix is: " + str(spike.RBN.states))
    for i in range(size(spike.RBN.states,0)):
        for j in range(i+1,size(spike.RBN.states,0)):
            #print ("The state matrix is currently: " + str(spike.RBN.states) + "\n")
            #print ("First row: " + str(spike.RBN.states[i,:]) + "\n")
            #print ("Second row: " + str(spike.RBN.states[j,:]) + "\n")
            if array_equal(spike.RBN.states[i,:],spike.RBN.states[j,:]) == True:
                print ("Returning: " + str(spike.RBN.states[i:j,]) + "\n")
                return spike.RBN.states[i:j,]
        
 
def findMolecularIntensityDebug (spike):
    stateNodes = findMolecularAttractorCycleDebug (spike)
    intensity = newFinalStage(stateNodes,spike) 
    print ("The intensity being returned is: " + str(intensity) + "\n")
    return intensity 