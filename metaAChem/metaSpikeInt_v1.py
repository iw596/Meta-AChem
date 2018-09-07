# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 13:24:39 2018

@author: iw596
"""
import numpy as np

def calculateIntensity (metaSpike):
    originalStateNodes = returnStateNodes(metaSpike)
    numNodes = len(originalStateNodes)
    
    if metaSpike.bonded == True:
        return -1
    else:
        intensity = calculateUnbondedIntensity(metaSpike,numNodes)
    
    return intensity
    


def returnStateNodes (metaSpike):
    stateNodes = np.array([],dtype = int)
    if metaSpike.typeSpike == 1:
        for i in range(len(metaSpike.danglingNodeList)):
            stateNodes = np.append(stateNodes,metaSpike.danglingNodeList[i].state)
    else:
        for i in range(len(metaSpike.danglingTailList)):
            for j in range(len(metaSpike.danglingTailList[i])):
                stateNodes = np.append(stateNodes,metaSpike.danglingTailList[i][j].state)
    
    return stateNodes



def unbondedAttractorCycle (numNodes,metaSpike):
    numAttempts = numNodes
    numUpdates = numNodes + 30
    stateNodes = np.array([],dtype = int)
    for i in range(numAttempts):
        for j in range(i):
            metaSpike.updateMetaSpike()
        stateNodes = returnStateNodes(metaSpike)
        print ("The initial state is: " + str(stateNodes) + "\n")
        states = []
        for k in range(numUpdates):
            metaSpike.updateMetaSpike()
            states = returnStateNodes(metaSpike)
            stateNodes = np.vstack((stateNodes,states))
            states = []
            if np.array_equal(stateNodes[0,:],stateNodes[k,:]) == True:
                stateNodes = np.delete(stateNodes,np.size(stateNodes,0)-1,0)
                return stateNodes
            
    return -1    
    
def calculateUnbondedIntensity (metaSpike,numNodes):
    
    stateNodes = unbondedAttractorCycle (numNodes,metaSpike)
    
    intensity = newFinalStage(stateNodes,metaSpike)
    #for i in range (len(spike.nodeList)):
     #   print ("Nodes in spike: " + str(spike.nodeList[i].nodeNumber) + "\n")
    #print ("The intensity is: " + str(intensity) + "\n")
    return intensity         

def newFinalStage (states,metaSpike):
    
    #print ("Number of dimenions is: " + str(states.ndim) + "\n")
    print ("The attractor cycle is: " + str(states) + "\n")
    
    #First check that a cycle has been found if not then return a string to indicate failure
    if states is None:
        return 'a'
    
    if states.ndim == 1:
        states = np.atleast_2d(states)
    intensity = 0
    # Go through each column
    for i in range(np.size(states,1)):
        #Go through each row
        #print ("Column number " + str(i) +  "\n" )
        currentSum = 0 
        #print ("Starting value for sum: " + str(currentSum) + "\n")
        for j in range(np.size(states,0)):
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
        if currentSum == np.size(states,0):
            intensity += 1
        elif currentSum == (np.size(states,0) *-1):
            intensity -= 1
        
    #print ("The intensity is: " + str(intensity) + "\n")
    return intensity

