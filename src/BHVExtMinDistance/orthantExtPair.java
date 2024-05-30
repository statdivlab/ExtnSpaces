/** This is intended as the class defining the vertex of a graph, which represents a complete maximal orthant that is part of the extension spaces of a tree. It will be connected to other vertices representing other orthants that are neighbouring by rotation. 

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;

public class orthantExtPair{
    private OrthExt OrthE1;
    private OrthExt OrthE2;
    
    //Constructor
    public orthantExtPair(OrthExt OE1, OrthExt OE2){
        this.OrthE1 = OE1;
        this.OrthE2 = OE2;
    }
    
    //Getters
    
    public OrthExt getOrthE1(){
        return this.OrthE1;
    }
    
    public OrthExt getOrthE2(){
        return this.OrthE2;
    }
    
}