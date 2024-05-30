/** This is intended as the class defining the vertex of a graph, which represents a complete maximal orthant that is part of the extension spaces of a tree. It will be connected to other vertices representing other orthants that are neighbouring by rotation. 

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;

public class orthantVertex{
    private int ID;
    private Vector<Bipartition> orthantAxis;//All the axes defining the orthant in this vertex.
    
    //Constructor
    public orthantVertex(int newID, Vector<Bipartition> newOrthantAxis){
        this.ID = newID;
        this.orthantAxis = newOrthantAxis;
    }
    
    //Getters
    
    public int getID(){
        return this.ID;
    }
    
    public Vector<Bipartition> getOrthantAxis(){
        return orthantAxis;
    }
    
    public Vector<Bipartition> getOrthantAxisClone(){
        Vector<Bipartition> copyOrthantAxis = new Vector<Bipartition>();
        for (Bipartition bip : this.orthantAxis){
            copyOrthantAxis.add(bip.clone());
        }
        return copyOrthantAxis;
    }
    
}