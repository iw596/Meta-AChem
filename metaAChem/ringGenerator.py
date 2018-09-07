# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 14:19:54 2018
This script is used to test the v8 soup base reactor
@author: iw596
"""
import Reactor_v9 as reac

print ("Welcome to the reactor which generates the rings used by meta atoms.\n" + " You have three parameter to select, the number of wanted products, the size of the reactor and the max size of atom. \n The number of wanted products is the number of rings which will be produced in the reactor run. The size of the reactor is the number of atoms in the reactor that are avaialle to react and the max size of atom is the maximum number of nodes in the RBN which makes up an atom.  I would go for 100,20,15 as this matches the settings of the rings already generated and included in this program ")


wantedNumProducts = int(input("Please enter how many rings you want to generate \n"))
sizeSoup = int(input("Please enter how many atoms you want the reactor to contain \n"))
sizeRBn  = 100int(input("Please enter the max number of nodes an atom can have \n"))
reactor = reac.Reactor(150,20,15)
reactor.runReactor()