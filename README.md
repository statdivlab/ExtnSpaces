# BHV Extension Spaces

BHVExtMinDistance is a Java package for the construction of Extension Spaces in the BHV tree space and the computation of distances between these extension spaces. 

The main executor ExtDistExecutor.jar performs our algorithm, based on reduced gradient methods, to find the minimum-length path between two extension spaces. The details of the algorithm can be found in (link to preprint TBD). It takes as input a text file with a list of Trees in Newick Format and for each pair (two trees listed one after the other) it: 

- Finds the union of the leaf sets and builds the extension spaces for each tree in the BHV-space for the complete leaf set. 
- Performs the optimization method to find the shortest path in between the extension spaces.
- Prints out a report and (if asked for it) a table detailing  results obtained by the algorithm.
    
## Installation 

The executable .jar is free to download. It was built on Java 20. 

## Input

The only required input to run ExtDistExecutor.jar is a text file containing a list of trees in Newick format. Each tree must occupy its own line, paired up trees must be one after the other without any empty line between them and an empty line must separate different pairs (see InputList.txt in Example1). Since the input file can contain several pairs of trees, and for each pair a separate report will be written, the user has the option to add an extra line on topo of each tree pair to give them a useful name (see InputList2.txt in Example2). The default name for the pairs in the list are Pair\#.   

## Use

To run ExtDistExecutor.jar from the command line (at the folder containing the executable) use the command

```sh
java -jar ExtDistExecutor.jar <InputFileName.txt> [options]
```

Options are: 

- -d: Allows you to pick the name of the new folder results are going into. Default is ExtDistResults
- -o: Allows you to pick the name of the output documents (previous to the specific name for each tree pair). Default is ExtensionDistance. 
- -k: Number of Threads to use in case of multi-threading. It performs sequentially if not specified. 
- -T: Indicator. If present, it not only prints a Report on a text file, but also creates a csv table with information in each orthant pair in our algorithm. 
- -l: When -T is present, but csv should not contain all information of the whole set of orthant pairs, indicate with this option how many rows (number of orthant pairs) should the table contain maximum. 
- -Tol1: Use to change the tolerance threshold for the gradient to be considered equal to zero. Default is 1E-8.
- -Tol2: Use to change the tolerance threshold for the partial derivative of line search to be considered equal to zero. Default is 1E-16. 

## Example

To obtain the results shown in Example1 the default code can be runned: 

```sh
java -jar ExtDistExecutor.jar InputList.txt
```

To obtain the results in Example2, we run the command,

```sh
java -jar ExtDistExecutor.jar InputList2.txt -d MyDirectoryName -o Distances -T -l 15 -Tol1 0.00001
```

in which we indicate results should be stored in the folder "MyDirectoryName", each reports name should start with the name "Distances", we indicate we want the algorithm to print a table with a total of 15 lines and the tolerance threshold is changed to 0.00001. 
