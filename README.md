# Meta-AChem
A proof of concept for a higher level artificial chemistry

--- This is a work in progress not all feature are complete and there will be bugs ----



Requirements
Python 3.6
Numpy

Included in this program:

- Set of already generated rings
- Code to generate new rings
- Code to generate meta atoms
- Code to react meta atom dangling nodes together


Running the low level reactor:
This part of the program generates rings of atoms which are then pickled. To run this program download it and run it from the script called "ringGenerator.py".
Upon running a it will present a paragraph of text informing you of what parameters you can changes. Do not alter the numSets.txt file, doing this will
mean the script cannot correctly store the ring pickle file.

Note - The file paths in Reactor_v9.py will need to be changed, this file path is located on lines 57

Running the high level reactor:
The high level reactor is currently a work in progress and can only react a single pair of meta atoms at a time and is currently not able to save the resulting product. To run
the high level reactor download everything and run the script called "metaReactor_v1.py"

Note- The file path in metaReactor_v1.py will need to be changed, this file path is located on line 30 and will need to be changed to X:\\XXXX\\XXXX\\ringMolecules where
	the X's are the path the ringMolecules file is located
