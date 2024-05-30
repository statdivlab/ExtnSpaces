/** This is intended as the representation of the subset of the extension space of a tree retricted to a particular orthant. 

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;

public class OrthExt{
    private PhyloTree originalTree; //The phylogenetic tree to which leaves are being added. 
    private Vector<String> completeLeafSet; //All leaves in the 'maximal' BHV space.
    private BitSet originalLeaves; //Vector indicating which leaves in the complete Leaf Set are part of the Original Tree
    private int[] orgLeaves2compLeaves; //Vector that maps the leaves in the original tree to the leaves in the complete leaf set. The entry in position k being of value r indicates that the k-th leaf in the original tree is equivalent to the r-th leaf in the complete leaf set.
    private int[] compLeaves2orgLeaves; //Vector that indicates which position corresponds to the leaf in the Original leaf set
    
    private Vector<Bipartition> orthantAxis; //Bipartitions representing the axes of the orthant (only internal splits). 
    private int Dim; //Number of free degrees. It coincides with the number of leaves being added. 
    
    private int[] Axis2Edges;
    private Vector<Integer> Edges2Axis;
    
    private double[] fixedLengths; //Vector to which the matrix maps every vector in the orthant to. In the restricted version, the size of this vector is the number of internal branches in the original tree; in the unrestricted version, the external edges in the original tree is also added.
    private extMatrix mapMatrix; //Matrix that describes how edges in extension trees are merged into the edges of the original tree (mapping matrix). In the restricted case, we are mapping interior to interior edges. In the unrestricted case, some exterior edges will also map into potentially exterior edges. 
    
    private Map<Integer, Vector<Integer>> mapList; // A hash map representing the mapMatrix; Keys are the rows and the vector is the entries in the row that should be 1. 
    private int[] backMap; //List of integer entries that indicate to which final edge (from original tree) the edge in the extension space is adding to. In the unrestricted case, we are adding external edges to leaves in the original tree. So length coincides with the number of axes in the orthant in the restricted case, or to the number of axes plus number of original leaves in the unrestricted case. 
    
    private PhyloTree startingTree; //Phylogenetic tree right in the "middle" of the extension space. The gradient of descent algorithm to find distances between Extension spaces will start in this tree. 
    
    private int oID;
    
    // Extra functions to internally do the construction of the Orthant Extension Space
    
    //Grabs a BitSet representing an axis in the orthant and returns its mapping to the edges of the original tree, i.e. the edge to which the bipartition goes to when performing leaf pruning. 
    private BitSet reducedBitSet(BitSet High){
        BitSet rBit = new BitSet();
        
        for (int i = 0; i<orgLeaves2compLeaves.length; i++){ //For each leaf in the Original Tree we verify if this leaf is part of the bipartition of the axis. 
            if(High.get(orgLeaves2compLeaves[i])){
                rBit.set(i);
            }
        }
        
        return(rBit);
    }
    
    //Quick function to determine if two BitSets represent the same bipartition.
    private boolean equivalentBip(BitSet Bit1, BitSet Bit2, int numberLeaves){
        return (Bit2.equals(Bit1) || ((!Bit2.intersects(Bit1)) && (Bit1.cardinality() + Bit2.cardinality() == numberLeaves)));
    }
    
    
    //Constructor
    public OrthExt(PhyloTree t, Vector<Bipartition> axis, Vector<String> cLeafSet){
        originalTree = t;
        completeLeafSet = cLeafSet;
        orthantAxis = axis;
        oID = 0;
        
        Axis2Edges = new int[axis.size()];
        for (int i = 0; i < axis.size(); i++){
            Axis2Edges[i] = -1;
        }
        Edges2Axis = new Vector<Integer>();
        
        EdgeAttribute[] startingleafEdgeLengths = new EdgeAttribute[cLeafSet.size()];//Leaf edges attributes in the trees in the extension space
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(t.getLeaf2NumMap());
        
        orgLeaves2compLeaves = new int[oLeafSet.size()];
        compLeaves2orgLeaves = new int[cLeafSet.size()];
        originalLeaves = new BitSet(cLeafSet.size());
        
        for (int i = 0; i< compLeaves2orgLeaves.length; i++){
            compLeaves2orgLeaves[i] = -1;
        }
        
        //For the i-th leaf in the original tree, we find its position r in the complete leaf set, and declare the i-th entry in orgLeaves2compLeaves to be r, and the r-th entry in compLeaves2orgLeaves to be i. We also set the r-th entry in the bitset originalLeaves to be true, since this leaf is effectively part of the original tree.
        for (int i = 0; i < oLeafSet.size(); i++){
            int temp = cLeafSet.indexOf(oLeafSet.get(i));
            if (temp == -1){
                System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
			    System.exit(1);
            }//end of if
            orgLeaves2compLeaves[i] = temp;
            compLeaves2orgLeaves[temp] = i;
            originalLeaves.set(temp);
            
        }//end of for loop
        
        //Loop to find all the attributes of the external edges in the "new" tree (the tree after attaching the extra leaves). These attributes coincide with the attributes in the original tree if the leaf was already part of that tree, and it is zero otherwise. 
        for (int i=0; i < cLeafSet.size(); i++){
            if(originalLeaves.get(i)){
                startingleafEdgeLengths[i] = t.getLeafEdgeAttribs()[compLeaves2orgLeaves[i]];
            } else {
                EdgeAttribute tempEdgeAttribute = new EdgeAttribute("[0]");
                startingleafEdgeLengths[i] = tempEdgeAttribute;
            }
        }
        
        fixedLengths = new double[t.getEdges().size()];
        
        //Iteration through the interior edges of the original tree to retrieve the values of the edges to which the coordinates in the orthant inside the extension space have to map to (get added up to be equal to). 
        Iterator<PhyloTreeEdge> edgesIter = t.getEdges().iterator();
        int count = 0;
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
			fixedLengths[count] = e.getAttribute().get(0);
            count++;
		}
        
        mapMatrix = new extMatrix(fixedLengths.length, orthantAxis.size());
        
        mapList = new HashMap<Integer, Vector<Integer>>();
        
        backMap = new int[axis.size()];
        
        for (int i = 0; i< backMap.length; i++){
            backMap[i] = -1;
        }
        
        Dim = orthantAxis.size() - fixedLengths.length;
        
        int[] numEdgesCombined = new int[fixedLengths.length]; //Vector that will record how many edges (coordinates of Axes of the orhans) are being merged to produce the edge in the original tree.
        
        //Loop through the interior edges in the original tree and the axes in the current orthant to verify which axis map to each edge, defining the mapping matrix and the number of edges being merged to form the interior edge in the original tree. 
        for (int i = 0; i < fixedLengths.length; i++){
            Vector<Integer> tempVect = new Vector<Integer>();
            for (int j = 0; j < orthantAxis.size(); j++){
                if (equivalentBip(t.getEdge(i).getOriginalEdge().getPartition(), reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                    mapMatrix.setItem(i,j,1);
                    numEdgesCombined[i]++;//We increment the respective vector whenever a new 1 appears in the mapMatrix
                    tempVect.add(j);
                    backMap[j] = i;
                }
            }
            mapList.put(i, tempVect);
        }
        
        
        Vector<PhyloTreeEdge> startingEdges = new Vector<PhyloTreeEdge>(); //Vector of edges for the starting tree.
        int CountEdges = 0;
        
        //We define the attributes of the internal edges in the starting tree as the attribute of the edge in the original tree divided by the number of edges being merged. 
        for (int j = 0; j < orthantAxis.size(); j++){
            for (int i = 0; i < fixedLengths.length; i++){
                if (mapMatrix.element(i,j) == 1){
                    double[] tempVecEA = {fixedLengths[i]/numEdgesCombined[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    PhyloTreeEdge tempEdge = new PhyloTreeEdge(axis.get(j), tempEA, j);
                    startingEdges.add(tempEdge);
                    Axis2Edges[j] = CountEdges;
                    Edges2Axis.add(Integer.valueOf(j));
                    CountEdges++;
                }
            }
        }
        
        
        startingTree = new PhyloTree(startingEdges, cLeafSet, startingleafEdgeLengths, false);
        
    }//end of constructor 1
    
    //Constructor 1 with ID: 
    public OrthExt(PhyloTree t, Vector<Bipartition> axis, Vector<String> cLeafSet, int oid){
        originalTree = t;
        completeLeafSet = cLeafSet;
        orthantAxis = axis;
        oID = oid;
        
        Axis2Edges = new int[axis.size()];
        for (int i = 0; i < axis.size(); i++){
            Axis2Edges[i] = -1;
        }
        Edges2Axis = new Vector<Integer>();
        
        EdgeAttribute[] startingleafEdgeLengths = new EdgeAttribute[cLeafSet.size()];//Leaf edges attributes in the trees in the extension space
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(t.getLeaf2NumMap());
        
        orgLeaves2compLeaves = new int[oLeafSet.size()];
        compLeaves2orgLeaves = new int[cLeafSet.size()];
        originalLeaves = new BitSet(cLeafSet.size());
        
        for (int i = 0; i< compLeaves2orgLeaves.length; i++){
            compLeaves2orgLeaves[i] = -1;
        }
        
        //For the i-th leaf in the original tree, we find its position r in the complete leaf set, and declare the i-th entry in orgLeaves2compLeaves to be r, and the r-th entry in compLeaves2orgLeaves to be i. We also set the r-th entry in the bitset originalLeaves to be true, since this leaf is effectively part of the original tree.
        for (int i = 0; i < oLeafSet.size(); i++){
            int temp = cLeafSet.indexOf(oLeafSet.get(i));
            if (temp == -1){
                System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
			    System.exit(1);
            }//end of if
            orgLeaves2compLeaves[i] = temp;
            compLeaves2orgLeaves[temp] = i;
            originalLeaves.set(temp);
            
        }//end of for loop
        
        //Loop to find all the attributes of the external edges in the "new" tree (the tree after attaching the extra leaves). These attributes coincide with the attributes in the original tree if the leaf was already part of that tree, and it is zero otherwise. 
        for (int i=0; i < cLeafSet.size(); i++){
            if(originalLeaves.get(i)){
                startingleafEdgeLengths[i] = t.getLeafEdgeAttribs()[compLeaves2orgLeaves[i]];
            } else {
                EdgeAttribute tempEdgeAttribute = new EdgeAttribute("[0]");
                startingleafEdgeLengths[i] = tempEdgeAttribute;
            }
        }
        
        fixedLengths = new double[t.getEdges().size()];
        
        //Iteration through the interior edges of the original tree to retrieve the values of the edges to which the coordinates in the orthant inside the extension space have to map to (get added up to be equal to). 
        Iterator<PhyloTreeEdge> edgesIter = t.getEdges().iterator();
        int count = 0;
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
			fixedLengths[count] = e.getAttribute().get(0);
            count++;
		}
        
        mapMatrix = new extMatrix(fixedLengths.length, orthantAxis.size());
        
        mapList = new HashMap<Integer, Vector<Integer>>();
        
        backMap = new int[axis.size()];
        
        for (int i = 0; i< backMap.length; i++){
            backMap[i] = -1;
        }
        
        Dim = orthantAxis.size() - fixedLengths.length;
        
        int[] numEdgesCombined = new int[fixedLengths.length]; //Vector that will record how many edges (coordinates of Axes of the orhans) are being merged to produce the edge in the original tree.
        
        //Loop through the interior edges in the original tree and the axes in the current orthant to verify which axis map to each edge, defining the mapping matrix and the number of edges being merged to form the interior edge in the original tree. 
        for (int i = 0; i < fixedLengths.length; i++){
            Vector<Integer> tempVect = new Vector<Integer>();
            for (int j = 0; j < orthantAxis.size(); j++){
                if (equivalentBip(t.getEdge(i).getOriginalEdge().getPartition(), reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                    mapMatrix.setItem(i,j,1);
                    numEdgesCombined[i]++;//We increment the respective vector whenever a new 1 appears in the mapMatrix
                    tempVect.add(j);
                    backMap[j] = i;
                }
            }
            mapList.put(i, tempVect);
        }
        
        
        Vector<PhyloTreeEdge> startingEdges = new Vector<PhyloTreeEdge>(); //Vector of edges for the starting tree.
        int CountEdges = 0;
        
        //We define the attributes of the internal edges in the starting tree as the attribute of the edge in the original tree divided by the number of edges being merged. 
        for (int j = 0; j < orthantAxis.size(); j++){
            for (int i = 0; i < fixedLengths.length; i++){
                if (mapMatrix.element(i,j) == 1){
                    double[] tempVecEA = {fixedLengths[i]/numEdgesCombined[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    PhyloTreeEdge tempEdge = new PhyloTreeEdge(axis.get(j), tempEA, j);
                    startingEdges.add(tempEdge);
                    Axis2Edges[j] = CountEdges;
                    Edges2Axis.add(Integer.valueOf(j));
                    CountEdges++;
                }
            }
        }
        
        
        startingTree = new PhyloTree(startingEdges, cLeafSet, startingleafEdgeLengths, false);
        
    }//end of constructor 1
    
    
    //Constructor 2: Allows for new leaves to be attached to external edges. 
    public OrthExt(PhyloTree t, Vector<Bipartition> axis, Vector<String> cLeafSet, boolean restricted){
        originalTree = t;
        completeLeafSet = cLeafSet;
        orthantAxis = axis;
        oID = 0;
        
        Axis2Edges = new int[axis.size()];
        for (int i = 0; i < axis.size(); i++){
            Axis2Edges[i] = -1;
        }
        Edges2Axis = new Vector<Integer>();
        
        EdgeAttribute[] startingleafEdgeLengths = new EdgeAttribute[cLeafSet.size()];//Leaf edges attributes in the trees in the extension space
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(t.getLeaf2NumMap());
        
        orgLeaves2compLeaves = new int[oLeafSet.size()];
        compLeaves2orgLeaves = new int[cLeafSet.size()];
        originalLeaves = new BitSet(cLeafSet.size());
        
        for (int i = 0; i< compLeaves2orgLeaves.length; i++){
            compLeaves2orgLeaves[i] = -1;
        }
        
        
        //For the i-th leaf in the original tree, we find its position r in the complete leaf set, and declare the i-th entry in orgLeaves2compLeaves to be r, and the r-th entry in compLeaves2orgLeaves to be i. We also set the r-th entry in the bitset originalLeaves to be true, since this leaf is effectively part of the original tree.
        for (int i = 0; i < oLeafSet.size(); i++){
            int temp = cLeafSet.indexOf(oLeafSet.get(i));
            if (temp == -1){
                System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
			    System.exit(1);
            }//end of if
            orgLeaves2compLeaves[i] = temp;
            compLeaves2orgLeaves[temp] = i;
            originalLeaves.set(temp);
            
        }//end of for loop
        
        
        //If restricted == TRUE; then the fixed lenghts focus only on internal branches. 
        //If restricted == FALSE; then fixed lenghts include the attributes of the original leaves. 
        if (restricted){
            fixedLengths = new double[t.getEdges().size()];
        } else {
            fixedLengths = new double[oLeafSet.size() + t.getEdges().size()];
        }
        
        
        //Iteration through the exterior (in the case restricted == false) and interior edges of the original tree to retrieve the values of the edges to which the coordinates in the orthant inside the extension space have to map to (get added up to be equal to).
        EdgeAttribute[] oLeafAtt = t.getLeafEdgeAttribs();
        Iterator<PhyloTreeEdge> edgesIter = t.getEdges().iterator();
        int count = 0;
        if (restricted == false){
            for (int i = 0; i < oLeafAtt.length; i++){
                fixedLengths[count] = oLeafAtt[i].get(0);
                count++;
            }
        }
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
			fixedLengths[count] = e.getAttribute().get(0);
            count++;
		}
        
        if (restricted){
            mapMatrix = new extMatrix(fixedLengths.length, orthantAxis.size());
        } else {
            mapMatrix = new extMatrix(fixedLengths.length, oLeafSet.size() + orthantAxis.size());
        }
        
        mapList = new HashMap<Integer, Vector<Integer>>();
        
        //backMap indicates to which final edge (from original tree) the edge in the extension space is adding to.
        if (restricted){
            backMap = new int[axis.size()];
        } else {
            backMap = new int[oLeafSet.size() + axis.size()];
        }
        
        for (int i = 0; i< backMap.length; i++){
            backMap[i] = -1;
        }
        
        if(restricted){
            Dim = orthantAxis.size() - t.getEdges().size();   
        } else {
            Dim = 0;
        }
        
        int[] numEdgesCombined = new int[fixedLengths.length]; //Vector that will record how many edges (coordinates of Axes of the orthans and possibly external edges to original leaves) are being merged to produce the edge in the original tree.
        
        //Loop through the interior (and possible exterior) edges in the original tree and the axes (and possibly some of the exterior edges) in the current orthant to verify which axis map to each edge, defining the mapping matrix and the number of edges being merged to form the interior edge in the original tree. 
        if (restricted){
            for (int i = 0; i < fixedLengths.length; i++){
                Vector<Integer> tempVect = new Vector<Integer>();
                for (int j = 0; j < orthantAxis.size(); j++){
                    if (equivalentBip(t.getEdge(i).getOriginalEdge().getPartition(), reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                        mapMatrix.setItem(i,j,1);
                        numEdgesCombined[i]++;//We increment the respective vector whenever a new 1 appears in the mapMatrix
                        tempVect.add(j);
                        backMap[j] = i;
                    }
                }
                mapList.put(i, tempVect);
            }
        } else {
            for (int i = 0; i < oLeafSet.size(); i++){
                Vector<Integer> tempVect = new Vector<Integer>();
                mapMatrix.setItem(i,i,1);
                numEdgesCombined[i]++; 
                tempVect.add(i);
                backMap[i] = i;
                Dim++;
                BitSet tempExtBitSet = new BitSet(oLeafSet.size());
                tempExtBitSet.set(i);
                for (int j = 0; j < orthantAxis.size(); j++){
                    if (equivalentBip(tempExtBitSet, reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                        mapMatrix.setItem(i,j+oLeafSet.size(), 1);
                        numEdgesCombined[i]++;
                        tempVect.add(j+oLeafSet.size());
                        backMap[j+oLeafSet.size()] = i;
                        Dim++;
                    }
                }
                mapList.put(i, tempVect);
            }
            for (int i = 0; i < t.getEdges().size(); i++){
                Vector<Integer> tempVect = new Vector<Integer>();
                for (int j = 0; j < orthantAxis.size(); j++){
                    if (equivalentBip(t.getEdge(i).getOriginalEdge().getPartition(), reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                        mapMatrix.setItem(i+oLeafSet.size(),j+oLeafSet.size(),1);
                        numEdgesCombined[i+oLeafSet.size()]++;//We increment the respective vector whenever a new 1 appears in the mapMatrix
                        tempVect.add(j+oLeafSet.size());
                        backMap[j+oLeafSet.size()] = i+oLeafSet.size();
                        Dim++;
                    }
                }
                mapList.put(i+oLeafSet.size(), tempVect);
            }
            Dim = Dim - (oLeafSet.size() + t.getEdges().size());
        }
        
        
        Vector<PhyloTreeEdge> startingEdges = new Vector<PhyloTreeEdge>(); //Vector of edges for the starting tree.
        int CountEdges = 0;
        //Loop to find all the attributes of the external edges in the "new" tree (the tree after attaching the extra leaves). If restricted == TRUE, these attributes coincide with the attributes in the original tree if the leaf was already part of that tree, and it is zero otherwise. 
        //If restricted == FALSE, some of the attributes will be shared with external leaves.
        if (restricted){
            for (int i=0; i < cLeafSet.size(); i++){
                if(originalLeaves.get(i)){
                    startingleafEdgeLengths[i] = t.getLeafEdgeAttribs()[compLeaves2orgLeaves[i]];
                } else {
                    EdgeAttribute tempEdgeAttribute = new EdgeAttribute("[0]");
                    startingleafEdgeLengths[i] = tempEdgeAttribute;
                }
            }
        } else {
            for (int i=0; i < cLeafSet.size(); i++){
                if(originalLeaves.get(i)){
                    double[] tempAtt = new double[]{t.getLeafEdgeAttribs()[compLeaves2orgLeaves[i]].get(0)/numEdgesCombined[compLeaves2orgLeaves[i]]};
                    startingleafEdgeLengths[i] = new EdgeAttribute(tempAtt);
                } else {
                    EdgeAttribute tempEdgeAttribute = new EdgeAttribute("[0]");
                    startingleafEdgeLengths[i] = tempEdgeAttribute;
                }
            }
        }
        
        //We define the attributes of the internal edges in the starting tree as the attribute of the edge in the original tree divided by the number of edges being merged. 
        if (restricted){
            for (int j = 0; j < orthantAxis.size(); j++){
                for (int i = 0; i < fixedLengths.length; i++){
                    if (mapMatrix.element(i,j) == 1){
                        double[] tempVecEA = {fixedLengths[i]/numEdgesCombined[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        PhyloTreeEdge tempEdge = new PhyloTreeEdge(axis.get(j), tempEA, j);
                        startingEdges.add(tempEdge);
                        Axis2Edges[j] = CountEdges;
                        Edges2Axis.add(Integer.valueOf(j));
                        CountEdges++;
                    }
                }
            }
        } else {
            for (int j = 0; j < orthantAxis.size(); j++){
                for (int i = 0; i < fixedLengths.length; i++){
                    if (mapMatrix.element(i,j+oLeafSet.size()) == 1){
                        double[] tempVecEA = {fixedLengths[i]/numEdgesCombined[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        PhyloTreeEdge tempEdge = new PhyloTreeEdge(axis.get(j), tempEA, j);
                        startingEdges.add(tempEdge);
                        Axis2Edges[j] = CountEdges;
                        Edges2Axis.add(Integer.valueOf(j));
                        CountEdges++;
                    }
                }
            }   
        }
        
        startingTree = new PhyloTree(startingEdges, cLeafSet, startingleafEdgeLengths, false);
        
    }//end of constructor 2
    
    //Constructor 2 wiht oID: Allows for new leaves to be attached to external edges. 
    public OrthExt(PhyloTree t, Vector<Bipartition> axis, Vector<String> cLeafSet, boolean restricted, int oid){
        originalTree = t;
        completeLeafSet = cLeafSet;
        orthantAxis = axis;
        oID = oid;
        
        Axis2Edges = new int[axis.size()];
        for (int i = 0; i < axis.size(); i++){
            Axis2Edges[i] = -1;
        }
        Edges2Axis = new Vector<Integer>();
        
        EdgeAttribute[] startingleafEdgeLengths = new EdgeAttribute[cLeafSet.size()];//Leaf edges attributes in the trees in the extension space
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(t.getLeaf2NumMap());
        
        orgLeaves2compLeaves = new int[oLeafSet.size()];
        compLeaves2orgLeaves = new int[cLeafSet.size()];
        originalLeaves = new BitSet(cLeafSet.size());
        
        for (int i = 0; i< compLeaves2orgLeaves.length; i++){
            compLeaves2orgLeaves[i] = -1;
        }
        
        
        //For the i-th leaf in the original tree, we find its position r in the complete leaf set, and declare the i-th entry in orgLeaves2compLeaves to be r, and the r-th entry in compLeaves2orgLeaves to be i. We also set the r-th entry in the bitset originalLeaves to be true, since this leaf is effectively part of the original tree.
        for (int i = 0; i < oLeafSet.size(); i++){
            int temp = cLeafSet.indexOf(oLeafSet.get(i));
            if (temp == -1){
                System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
			    System.exit(1);
            }//end of if
            orgLeaves2compLeaves[i] = temp;
            compLeaves2orgLeaves[temp] = i;
            originalLeaves.set(temp);
            
        }//end of for loop
        
        
        //If restricted == TRUE; then the fixed lenghts focus only on internal branches. 
        //If restricted == FALSE; then fixed lenghts include the attributes of the original leaves. 
        if (restricted){
            fixedLengths = new double[t.getEdges().size()];
        } else {
            fixedLengths = new double[oLeafSet.size() + t.getEdges().size()];
        }
        
        
        //Iteration through the exterior (in the case restricted == false) and interior edges of the original tree to retrieve the values of the edges to which the coordinates in the orthant inside the extension space have to map to (get added up to be equal to).
        EdgeAttribute[] oLeafAtt = t.getLeafEdgeAttribs();
        Iterator<PhyloTreeEdge> edgesIter = t.getEdges().iterator();
        int count = 0;
        if (restricted == false){
            for (int i = 0; i < oLeafAtt.length; i++){
                fixedLengths[count] = oLeafAtt[i].get(0);
                count++;
            }
        }
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
			fixedLengths[count] = e.getAttribute().get(0);
            count++;
		}
        
        if (restricted){
            mapMatrix = new extMatrix(fixedLengths.length, orthantAxis.size());
        } else {
            mapMatrix = new extMatrix(fixedLengths.length, oLeafSet.size() + orthantAxis.size());
        }
        
        mapList = new HashMap<Integer, Vector<Integer>>();
        
        //backMap indicates to which final edge (from original tree) the edge in the extension space is adding to.
        if (restricted){
            backMap = new int[axis.size()];
        } else {
            backMap = new int[oLeafSet.size() + axis.size()];
        }
        
        for (int i = 0; i< backMap.length; i++){
            backMap[i] = -1;
        }
        
        if(restricted){
            Dim = orthantAxis.size() - t.getEdges().size();   
        } else {
            Dim = 0;
        }
        
        int[] numEdgesCombined = new int[fixedLengths.length]; //Vector that will record how many edges (coordinates of Axes of the orthans and possibly external edges to original leaves) are being merged to produce the edge in the original tree.
        
        //Loop through the interior (and possible exterior) edges in the original tree and the axes (and possibly some of the exterior edges) in the current orthant to verify which axis map to each edge, defining the mapping matrix and the number of edges being merged to form the interior edge in the original tree. 
        if (restricted){
            for (int i = 0; i < fixedLengths.length; i++){
                Vector<Integer> tempVect = new Vector<Integer>();
                for (int j = 0; j < orthantAxis.size(); j++){
                    if (equivalentBip(t.getEdge(i).getOriginalEdge().getPartition(), reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                        mapMatrix.setItem(i,j,1);
                        numEdgesCombined[i]++;//We increment the respective vector whenever a new 1 appears in the mapMatrix
                        tempVect.add(j);
                        backMap[j] = i;
                    }
                }
                mapList.put(i, tempVect);
            }
        } else {
            for (int i = 0; i < oLeafSet.size(); i++){
                Vector<Integer> tempVect = new Vector<Integer>();
                mapMatrix.setItem(i,i,1);
                numEdgesCombined[i]++; 
                tempVect.add(i);
                backMap[i] = i;
                Dim++;
                BitSet tempExtBitSet = new BitSet(oLeafSet.size());
                tempExtBitSet.set(i);
                for (int j = 0; j < orthantAxis.size(); j++){
                    if (equivalentBip(tempExtBitSet, reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                        mapMatrix.setItem(i,j+oLeafSet.size(), 1);
                        numEdgesCombined[i]++;
                        tempVect.add(j+oLeafSet.size());
                        backMap[j+oLeafSet.size()] = i;
                        Dim++;
                    }
                }
                mapList.put(i, tempVect);
            }
            for (int i = 0; i < t.getEdges().size(); i++){
                Vector<Integer> tempVect = new Vector<Integer>();
                for (int j = 0; j < orthantAxis.size(); j++){
                    if (equivalentBip(t.getEdge(i).getOriginalEdge().getPartition(), reducedBitSet(axis.get(j).getPartition()),oLeafSet.size())){
                        mapMatrix.setItem(i+oLeafSet.size(),j+oLeafSet.size(),1);
                        numEdgesCombined[i+oLeafSet.size()]++;//We increment the respective vector whenever a new 1 appears in the mapMatrix
                        tempVect.add(j+oLeafSet.size());
                        backMap[j+oLeafSet.size()] = i+oLeafSet.size();
                        Dim++;
                    }
                }
                mapList.put(i+oLeafSet.size(), tempVect);
            }
            Dim = Dim - (oLeafSet.size() + t.getEdges().size());
        }
        
        
        Vector<PhyloTreeEdge> startingEdges = new Vector<PhyloTreeEdge>(); //Vector of edges for the starting tree.
        int CountEdges = 0;
        //Loop to find all the attributes of the external edges in the "new" tree (the tree after attaching the extra leaves). If restricted == TRUE, these attributes coincide with the attributes in the original tree if the leaf was already part of that tree, and it is zero otherwise. 
        //If restricted == FALSE, some of the attributes will be shared with external leaves.
        if (restricted){
            for (int i=0; i < cLeafSet.size(); i++){
                if(originalLeaves.get(i)){
                    startingleafEdgeLengths[i] = t.getLeafEdgeAttribs()[compLeaves2orgLeaves[i]];
                } else {
                    EdgeAttribute tempEdgeAttribute = new EdgeAttribute("[0]");
                    startingleafEdgeLengths[i] = tempEdgeAttribute;
                }
            }
        } else {
            for (int i=0; i < cLeafSet.size(); i++){
                if(originalLeaves.get(i)){
                    double[] tempAtt = new double[]{t.getLeafEdgeAttribs()[compLeaves2orgLeaves[i]].get(0)/numEdgesCombined[compLeaves2orgLeaves[i]]};
                    startingleafEdgeLengths[i] = new EdgeAttribute(tempAtt);
                } else {
                    EdgeAttribute tempEdgeAttribute = new EdgeAttribute("[0]");
                    startingleafEdgeLengths[i] = tempEdgeAttribute;
                }
            }
        }
        
        //We define the attributes of the internal edges in the starting tree as the attribute of the edge in the original tree divided by the number of edges being merged. 
        if (restricted){
            for (int j = 0; j < orthantAxis.size(); j++){
                for (int i = 0; i < fixedLengths.length; i++){
                    if (mapMatrix.element(i,j) == 1){
                        double[] tempVecEA = {fixedLengths[i]/numEdgesCombined[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        PhyloTreeEdge tempEdge = new PhyloTreeEdge(axis.get(j), tempEA, j);
                        startingEdges.add(tempEdge);
                        Axis2Edges[j] = CountEdges;
                        Edges2Axis.add(Integer.valueOf(j));
                        CountEdges++;
                    }
                }
            }
        } else {
            for (int j = 0; j < orthantAxis.size(); j++){
                for (int i = 0; i < fixedLengths.length; i++){
                    if (mapMatrix.element(i,j+oLeafSet.size()) == 1){
                        double[] tempVecEA = {fixedLengths[i]/numEdgesCombined[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        PhyloTreeEdge tempEdge = new PhyloTreeEdge(axis.get(j), tempEA, j);
                        startingEdges.add(tempEdge);
                        Axis2Edges[j] = CountEdges;
                        Edges2Axis.add(Integer.valueOf(j));
                        CountEdges++;
                    }
                }
            }   
        }
        
        startingTree = new PhyloTree(startingEdges, cLeafSet, startingleafEdgeLengths, false);
        
    }//end of constructor 2
    
    //Clone Constructor: 
    public OrthExt(OrthExt OE){
        this.originalTree = new PhyloTree(OE.getOriginalTree());
        this.completeLeafSet = polyAlg.Tools.myVectorCloneString(OE.getCompleteLeafSet());
        this.originalLeaves = (BitSet) OE.getOriginalLeaves().clone();
        
        this.orgLeaves2compLeaves = OE.getOrgLeaves2compLeaves().clone();
        this.compLeaves2orgLeaves = OE.getCompLeaves2orgLeaves().clone();
        
        this.orthantAxis = new Vector<Bipartition>();    
        for (Bipartition bip : OE.getOrthantAxis()){
            this.orthantAxis.add(bip.clone());
        }
        this.Dim = OE.getDim();
        
        this.Axis2Edges = OE.getCloneAxis2Edges();
        this.Edges2Axis = OE.getCloneEdges2Axis();
        
        this.fixedLengths = OE.getFixedLengths().clone();
        this.mapMatrix = new extMatrix(OE.getMapMatrix());
        
        this.mapList = OE.getMapListClone();
        this.backMap = OE.getCloneBackMap();
        
        this.startingTree = new PhyloTree(OE.getStartTree());
        this.oID = OE.getOID();
        
    }
    
    
    //Getters and Printers
    public PhyloTree getOriginalTree(){
        return originalTree;
    }
    
    public BitSet getOriginalLeaves(){
        return originalLeaves;
    }
    
    public int[] getOriginalLeavesAsVector(){
        int[] vectAns = new int[originalLeaves.length()];
        for(int i = 0; i< originalLeaves.length(); i++){
            if (originalLeaves.get(i)){
                vectAns[i] = 1;
            } else{
                vectAns[i] = 0;
            }
        }
        return vectAns;
    }
    
    public int[] getOrgLeaves2compLeaves(){
        return orgLeaves2compLeaves;
    }
    
    public int getOrgLeaves2compLeaves(int i){
        return orgLeaves2compLeaves[i];
    }
    
    public int[] getCompLeaves2orgLeaves(){
        return compLeaves2orgLeaves;
    }
    
    public int getCompLeaves2orgLeaves(int i){
        return compLeaves2orgLeaves[i];
    }
    
    
    public Vector<Bipartition> getOrthantAxis(){
        return orthantAxis;
    }
    
    public Bipartition getOrthantAxis(int i){
        return orthantAxis.get(i);
    }
    
    public Vector<String> getCompleteLeafSet(){
        return completeLeafSet;
    }
    
    public double[] getFixedLengths(){
        return fixedLengths;
    }
    
    public double getFixedLengths(int i){
        return fixedLengths[i];
    }
    
    public extMatrix getMapMatrix(){
        return mapMatrix;
    }
    
    public Map<Integer, Vector<Integer>> getMapList(){
        return mapList;
    }
    
    public Map<Integer, Vector<Integer>> getMapListClone(){
        Map<Integer, Vector<Integer>> mapCopy = new HashMap<Integer, Vector<Integer>>();
        
        for (Integer key : this.mapList.keySet()){
            Vector<Integer> intObj = new Vector<Integer>();
            for (int i  = 0; i < this.mapList.get(key).size(); i++){
                intObj.add(Integer.valueOf(this.mapList.get(key).get(i)));
            }
            mapCopy.put(Integer.valueOf(key),intObj);   
		}
        
        return mapCopy;
    }
    
    public int[] getBackMap(){
        return backMap;
    }
    
    public int getBackMap(int i){
        return backMap[i];
    }
    
    public int[] getCloneBackMap(){
        return this.backMap.clone();
    }
    
    public int getDim(){
        return Dim;
    }
    
    public int getOID(){
        return oID;
    }
    
    public PhyloTree getStartTree(){
        return startingTree;
    }
    
    public int[] getCloneAxis2Edges(){
        return this.Axis2Edges.clone();
    }
    
    public Vector<Integer> getCloneEdges2Axis(){
        return new Vector<Integer>(this.Edges2Axis);
    }
    
    //Function to print the axes of the orthants with a nice format. 
    public void PrintOrthantAxes(){
         System.out.print("Axes in the extension orthant: ");
         for (int i=0; i < this.orthantAxis.size(); i++){
             Bipartition eClone = this.orthantAxis.get(i).clone();
             eClone.complement(this.completeLeafSet.size());
             System.out.print("  {" + this.orthantAxis.get(i).toStringVerbose(this.completeLeafSet) + "|" + eClone.toStringVerbose(this.completeLeafSet) + "}");
            /**if (this.orthantAxis.get(i).getPartition().cardinality() > (this.completeLeafSet.size()/2)){
                Bipartition eClone = this.orthantAxis.get(i).clone();
                eClone.complement(this.completeLeafSet.size());
                System.out.print("  {" + eClone.toStringVerbose(this.completeLeafSet) + "}");
            } else {
                System.out.print("  {" + this.orthantAxis.get(i).toStringVerbose(this.completeLeafSet) + "}");
            }*/
            if (i < this.orthantAxis.size()-1){
                System.out.print(", ");
            }
        }
        
     }
    
    //Function to print a big Summary of the final product of the constructor. 
    public void printLN(){
        PhyloNicePrinter treePrint = new PhyloNicePrinter();
        System.out.println("EXTENSION SPACE CREATED of dimension " + this.Dim);
        System.out.println("Original tree: ");
        System.out.println(treePrint.toString(originalTree));
        System.out.println("");
        
        System.out.println("Map of the Original Leaves:");
        System.out.println("  " + Arrays.toString(this.getOriginalLeavesAsVector()));
        System.out.println("");
        
        
        System.out.println("Original leaves to Complete leaves:");
        System.out.println("  " + Arrays.toString(orgLeaves2compLeaves));
        System.out.println("");
        
        System.out.println("Complete leaves to Original leaves:");
        System.out.println("  " + Arrays.toString(compLeaves2orgLeaves));
        System.out.println("");
        
        System.out.println("The vector of final lengths:");
        System.out.println("  " + Arrays.toString(fixedLengths));
        System.out.println("");
        
        this.PrintOrthantAxes();
        System.out.println("");
        
        System.out.println("The mapping matrix:");
        System.out.println("Size of the matrix: (" + mapMatrix.getNRow() + ", " + mapMatrix.getNCol() +")");
        mapMatrix.PrintMat();
        System.out.println("");
        
        System.out.println("Mapping matrix, other representation:");
        for (int i = 0; i < mapMatrix.getNRow(); i++){
            if (mapList.containsKey(i)){
                System.out.println("  "+i+" --> "+ mapList.get(i));
            }
        }
        System.out.println("With the back Map:" + Arrays.toString(backMap));
        System.out.println("");
        
        System.out.println("From Axis to Edges:" + Arrays.toString(Axis2Edges));
        System.out.println("");
        
        
        System.out.println("From Edges to Axis:" + Edges2Axis);
        System.out.println("");
        
        System.out.println("New starting tree: ");
        System.out.println(treePrint.toString(startingTree));
        System.out.println("");
    }
    
    //Printing a smaller summary thant the previous function. 
    public void PrintReduced(){
        System.out.println("------ Orthant Extension Description ------");
        System.out.println("The dimension is: " + this.Dim);
        this.PrintOrthantAxes();
        System.out.println("");
        System.out.println("The vector of final lengths:" + Arrays.toString(fixedLengths));
        System.out.println("");
        System.out.println("The mapping matrix: ( "+ this.mapMatrix.getNRow() + ", " + this.mapMatrix.getNCol() +")");
        this.mapMatrix.PrintMat();
        System.out.println("-------------------------------------------");
    }
    
    //Printing a smaller summary thant the previous function, with the option to include starting tree; 
    public void PrintReduced(boolean withTree){
        System.out.println("------ Orthant Extension Description ------");
        System.out.println("The dimension is: " + this.Dim);
        this.PrintOrthantAxes();
        System.out.println("");
        System.out.println("The vector of final lengths:" + Arrays.toString(fixedLengths));
        System.out.println("");
        System.out.println("And the Mapping to return is: " + Arrays.toString(this.backMap));
        System.out.println("");
        System.out.println("The mapping matrix: ( "+ this.mapMatrix.getNRow() + ", " + this.mapMatrix.getNCol() +")");
        this.mapMatrix.PrintMat();
        if(withTree){
            System.out.println("Starting tree: ");
            PhyloNicePrinter treePrint = new PhyloNicePrinter();
            System.out.println(treePrint.toString(this.startingTree));
        }
        System.out.println("-------------------------------------------");
    }
    
//End of getters
    
}