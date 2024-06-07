# BHV Extension Spaces

BHVExtMinDistance is a Java package for the construction of Extension Spaces in the BHV tree space and the computation of distances between these extension spaces. 

The main executor ExtDistExecutor.jar performs our algorithm, based on reduced gradient methods, to find the minimum-length path between two extension spaces. The details of the algorithm can be found in... (future paper). It takes as input a text file with a list of Trees in Newick Format and for each pair (two trees listed one after the other) it: 

- Finds the union of the leaf sets and builds the extension spaces for each tree in the BHV-space for the complete leaf set. 
- Performs the optimization method to find the shortest path in between the extension spaces.
- Prints out a report and (if asked for it) a table detailing  results obtained by the algorithm.
    
## Installation 

The executable .jar is free to download of use. It was build on Java 20. 

## Input

The only required input to run ExtDistExecutor.jar is a text file containing a list of trees in Newick format. Each tree must occupy its own line, paired up trees must be one after the other without any empty line between them and an empty line must separate different pairs (see InputList.txt in Example1). Since the input file can contain several pairs of trees, and for each pair a separate report will be written, the user has the option to add an extra line on topo of each tree pair to give them a useful name (see InputList2.txt in Example2). The default name for the pairs in the list are Pair\#.   

## Use

To run ExtDistExecutor.jar from the command line (at the folder containing the executable) use the command

```sh
java -jar ExtDistExecutor.jar <InputFileName.txt> [options]


Options are: 

- -d: Allows you to pick the name of the new folder results are going into. Default is ExtDistResults
- -o: Allows you to pick the name of the output documents (previous to the specific name for each tree pair). Default is ExtensionDistance. 

