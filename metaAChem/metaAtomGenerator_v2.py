# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 09:55:27 2018
This script is used to generate 
@author: iw596
"""


from numpy import *
import pickle
import metaSpikeGenerator_v1 as metSpkGen
import MetaAtom_v2 as metAt

def atomGenerator (metaAtomNum,molecule):
   # To generate metaatom we first need to analyse the dangling nodes of the molecule in order to calculate the metaspikes 
   # to do this we call a function located in the metaSpikeGeneration script
   metSpike1,metSpike2 = metSpkGen.generateMetaSpikes(molecule)
   metaAtom = metAt.MetaAtom(metaAtomNum,molecule)
   if metSpike1 != -1:
       metaAtom.addMetaSpike(metSpike1)
   if metSpike2 != -1:
        metaAtom.addMetaSpike(metSpike2)
#        print ("Initial state \n")
   metaAtom.calculateState() 
   
   return metaAtom

