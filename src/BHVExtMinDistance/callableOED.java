/** This is intended as the class defining the distance between two extension spaces. It will include methods to compute this distance, as well as a list of pairs of trees that obtain smaller distances per each maximal orthant covered by these extension spaces.

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import java.util.concurrent.*; 
import distanceAlg1.*;

public class callableOED implements Callable<OrthExtDistance> {
    private orthantExtPair Opair;
    
    public callableOED(orthantExtPair OEpair){
        this.Opair = OEpair;
    }
    
    public OrthExtDistance call() {
        return new OrthExtDistance(Opair, false);
    }
}