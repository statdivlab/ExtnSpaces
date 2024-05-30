/** This is intended as the representation of the of the extension space of a tree.

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;
import polyAlg.*;


public class ExtensionSpace{
    //The whole tree to which leaves are been added
    private PhyloTree originalTree;
    //Vector indicating which leaves in the completeLeafSet are part of the Original Tree
    private BitSet originalLeaves; 
    //Vector that maps the leaves in the original tree to the leaves in the complete leaf set.
    //private int[] orgLeaves2compLeaves; 
    //All leaves in the 'maximal' BHV space.
    private Vector<String> completeLeafSet; 
    
    //All orthants to which the Extension Space of the tree belongs to, coonected by Rotation.
    private orthantGraph connectCluster;
    //List of all the maximal orthants to which the Extension Space belongs to.
    private Vector<OrthExt> listOrthants;
    //Number of maximal orthants to which the Extension Space belongs to.
    private int numOrthants;
    
    //constructors
    
    public ExtensionSpace(PhyloTree t, Vector<String> cLeafSet){
        this.originalTree = new PhyloTree(t);
        this.completeLeafSet = Tools.myVectorCloneString(cLeafSet);
        
        originalLeaves = new BitSet(this.completeLeafSet.size());
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(t.getLeaf2NumMap());
        
        //REALIZED I DO NOT USE THIS HERE, AND CALCULATE AGAIN INSIDE edgeCrossGraph. COMMENT --> DELETE WHEN SURE ALL GOOD
        //orgLeaves2compLeaves = new int[oLeafSet.size()]; 
        
        //For each leaf in the original tree's leaf set, we find its position in the complete set, and modify the arrays and BitSet accordingly. 
        //for (int i = 0; i < oLeafSet.size(); i++){
        //    int temp = cLeafSet.indexOf(oLeafSet.get(i));
        //    if (temp == -1){
        //        System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
		//	    System.exit(1);
        //    }
        //    orgLeaves2compLeaves[i] = temp;
        //    originalLeaves.set(temp);
        //}
        
        //Create the edgeCrossGraph calling its constructor with this tree and complete leaf set.
        edgeCrossGraph GraphAllEdges = new edgeCrossGraph(this.originalTree, this.completeLeafSet);
        
        //Compute the Maximal Independen Sets of the previous graph.
        GraphAllEdges.MIScalculator();
        
        //From the edgeCrossGraph, we find all maximal orthants in this extension space.
        connectCluster = new orthantGraph(GraphAllEdges);
        numOrthants = connectCluster.getVertexNum();
        
        listOrthants = new Vector<OrthExt>();
        
        for (int i = 0; i < numOrthants; i++){
            Vector<Bipartition> tempAxes = connectCluster.getAxesClone(i);
            OrthExt tempOrthExt = new OrthExt(new PhyloTree(t), tempAxes, Tools.myVectorCloneString(cLeafSet), i);
            listOrthants.add(tempOrthExt);
        }
    } //End constructor 1
    
    //constructor 2: admitting unrestricted version
    public ExtensionSpace(PhyloTree t, Vector<String> cLeafSet, boolean restricted){
        this.originalTree = new PhyloTree(t);
        this.completeLeafSet =  Tools.myVectorCloneString(cLeafSet);
        
        originalLeaves = new BitSet(this.completeLeafSet.size());
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(t.getLeaf2NumMap());
        //orgLeaves2compLeaves = new int[oLeafSet.size()]; 
        
        //REALIZED I DO NOT USE THIS HERE, AND CALCULATE AGAIN INSIDE edgeCrossGraph. COMMENT --> DELETE WHEN SURE ALL GOOD
        //For each leaf in the original tree's leaf set, we find its position in the complete set, and modify the arrays and BitSet accordingly. 
        //for (int i = 0; i < oLeafSet.size(); i++){
        //    int temp = cLeafSet.indexOf(oLeafSet.get(i));
        //    if (temp == -1){
        //        System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
        //            System.exit(1);
        //    }
        //    orgLeaves2compLeaves[i] = temp;
        //    originalLeaves.set(temp);
        //}
        
        //Create the edgeCrossGraph calling its constructor with this tree and complete leaf set.
        edgeCrossGraph GraphAllEdges = new edgeCrossGraph(this.originalTree, this.completeLeafSet, restricted);
        
        //Compute the Maximal Independen Sets of the previous graph.
        GraphAllEdges.MIScalculator();
        
        //From the edgeCrossGraph, we find all maximal orthants in this extension space.
        connectCluster = new orthantGraph(GraphAllEdges);
        numOrthants = connectCluster.getVertexNum();
        
        listOrthants = new Vector<OrthExt>();
        
        for (int i = 0; i < numOrthants; i++){
            Vector<Bipartition> tempAxes = connectCluster.getAxesClone(i);
            OrthExt tempOrthExt = new OrthExt(new PhyloTree(t), tempAxes, Tools.myVectorCloneString(cLeafSet), restricted, i);
            listOrthants.add(tempOrthExt);
        }
    }// end constructor 2
    
    //Printers and Getters. 
    
    public void PrintSummary(){
        PhyloNicePrinter treePrint = new PhyloNicePrinter();
        System.out.println("EXTENSION SPACE FOR: ");
        System.out.println(treePrint.toString(this.originalTree));
        System.out.println("");
        System.out.println("To complete Leaf Set: "+ this.completeLeafSet);
        System.out.println("");
        System.out.println("");
        System.out.println("The extension space is composed by "+this.numOrthants+" orthant extensions:");
        
        System.out.println("");
        for (int i = 0; i < numOrthants; i++){
            System.out.println("Orthant number "+i+ " adjacent to: " + Arrays.toString(connectCluster.getAdjIDs(i)));
            listOrthants.get(i).PrintReduced();
            System.out.println("");
        }
    }
    
    //Same as before, but with the option to include the starting trees of the orthants.
    public void PrintSummary(boolean withTrees){
        PhyloNicePrinter treePrint = new PhyloNicePrinter();
        System.out.println("EXTENSION SPACE FOR: ");
        System.out.println(treePrint.toString(this.originalTree));
        System.out.println("");
        System.out.println("To complete Leaf Set: "+ this.completeLeafSet);
        System.out.println("");
        System.out.println("");
        System.out.println("The extension space is composed by "+this.numOrthants+" orthant extensions:");
        
        System.out.println("");
        for (int i = 0; i < numOrthants; i++){
            System.out.println("Orthant number "+i+ " adjacent to: " + Arrays.toString(connectCluster.getAdjIDs(i)));
            listOrthants.get(i).PrintReduced(withTrees);
            System.out.println("");
        }
    }
    
    public void PrintSummaryLN(){
        PhyloNicePrinter treePrint = new PhyloNicePrinter();
        System.out.println("EXTENSION SPACE FOR: ");
        System.out.println(treePrint.toString(this.originalTree));
        System.out.println("");
        System.out.println("To complete Leaf Set: "+ this.completeLeafSet);
        System.out.println("");
        System.out.println("");
        System.out.println("The extension space is composed by "+this.numOrthants+" orthant extensions:");
        
        System.out.println("");
        for (int i = 0; i < numOrthants; i++){
            System.out.println("Orthant number "+i+ " adjacent to: " + Arrays.toString(connectCluster.getAdjIDs(i)));
            listOrthants.get(i).printLN();
            System.out.println("");
        }
    }
    
    public void PrintOrthants(){
        for (int i = 0; i < numOrthants; i++){
            System.out.println("Orthant number "+i+ " adjacent to: " + Arrays.toString(connectCluster.getAdjIDs(i)));
            System.out.print("  ");
            listOrthants.get(i).PrintOrthantAxes();
            System.out.println("");
        }
    }
    
    public Vector<OrthExt> getOrthExts(){
        return this.listOrthants;
    }
    
    public OrthExt getOrthExts(int i){
        return this.listOrthants.get(i);
    }
    
    public orthantGraph getConnectCluster(){
        return this.connectCluster;
    }
    
    public int getNumOrthants(){
        return this.numOrthants;
    }
    
}