/** This is intended as the class defining the vertex of a graph, which represents an interior edge in a binary tree, and it will be connected to other vertices representing other edges which are not compatible to the respective one. 

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;

public class edgeVertex{
    private int ID;
    private Bipartition edge;
    
    //Constructor
    public edgeVertex(int newID, Bipartition newEdge){
        this.ID = newID;
        this.edge = newEdge;
    }
    
    //Getters 
    
    public Bipartition getEdge(){
        return this.edge;
    }
    
    public int getID(){
        return this.ID;
    }
    
    //Function to print the bipartition in the edge in a nice format
    //TO DO: Use nice printer here in the future
    public String toStringVerbose(Vector<String> cLeafSet){
        return this.edge.toStringVerbose(cLeafSet);
    }
}