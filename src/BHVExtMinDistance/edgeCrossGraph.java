/** This is intended as the class defining the a graph connecting incompatible edges. The edges represented by the vertices in this graph are all those that may be obtain by attaching new leaves to a fixed tree. This graph is used to find all the possible maximal orthants where the extension space for this fixed tree may have trees on. 

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;
import polyAlg.*;

public class edgeCrossGraph{
    //Adjacency lists that define the graph
    private Map<edgeVertex, List<edgeVertex>> adjVertices;
    //Number of vertices in the graph;
    private int vertexNum;
    //List of the vertices in the graph, so we always have the i-th vertex with the ID equal to i.
    private Vector<edgeVertex> orderedVertices; 
    
    //Variables necessary for the computation of the Maximal Independent Sets, which will correspond to the axes in all maximal orthants included in the extension spaces. 
    //The algorithm to compute these maximal independent sets is that described in Tsukiyama, S., Ide, M., Ariyoshi, I., Shirakawa, I. (1977).
    private int[] IS; 
    private List<Set<Integer>> Bucket;
    private List<Vector<Integer>> MIS; //List of all maximal independent sets (only the IDs of the edge vertices).
    
    //My own function to determine if splits crosses
    private boolean Crosses(Bipartition Split1, Bipartition Split2, int NumberLeaves){
        if ((Split1.disjointFrom(Split2) || Split1.contains(Split2) || Split2.contains(Split1))){
            return false;
        } else {
            Bipartition temp1 = Split1.clone();
            temp1.complement(NumberLeaves);
            Bipartition temp2 = Split2.clone();
            temp2.complement(NumberLeaves);
            
            return(!temp1.disjointFrom(temp2));
        }
    }
    
    //Constructor
    public edgeCrossGraph(PhyloTree T, Vector<String> cLeafSet){
        
        //Vector indicating which leaves in the completeLeafSet are part of the Original Tree
        BitSet originalLeaves = new BitSet(cLeafSet.size());
        
        //Vector of the original tree's leaf set.
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(T.getLeaf2NumMap()); 
        
        //Vector that maps the leaves in the original tree to the leaves in the complete leaf set. 
        int[] orgLeaves2compLeaves = new int[oLeafSet.size()]; 
        
        //Number of internal edges in the original tree. 
        int m = T.getEdges().size(); 
        
        //Number of leaves to be added. 
        int l = cLeafSet.size() - oLeafSet.size(); 
        
        //The number of potentially new edges is the number of internal edges in the original tree times the different ways the l extra leaves can be added to this edges.
        this.vertexNum = (int) (m*Math.pow(2,l));   
        
        //Array of the leaves that are being added
        int[] listAddedLeaves = new int[l];
        
        //For each leaf in the original tree's leaf set, we find its position in the complete set, and modify the arrays and BitSet accordingly. 
        for (int i = 0; i < oLeafSet.size(); i++){
            int temp = cLeafSet.indexOf(oLeafSet.get(i));
            if (temp == -1){
                System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
			    System.exit(1);
            }
            orgLeaves2compLeaves[i] = temp;
            originalLeaves.set(temp);
            
        }
        
        
        //System.out.println("Org 2 Comp is: " + Arrays.toString(orgLeaves2compLeaves));
        
        //System.out.println("original Leaves" + originalLeaves);
        
        //We iterate over the set of complete leaves to find which leafs must be added to the original tree.
        int count = 0;
        for (int i = 0; i < cLeafSet.size(); i++){
            if (!originalLeaves.get(i)){
                listAddedLeaves[count] = i;
                count++;
            }
        }
        
        //System.out.println("listAddedLeaves: " + Arrays.toString(listAddedLeaves));
        //System.out.println("Number of leaves added = " + l);
        
        //For each edge in the original tree, we loop through all the ways the extra leaves can be added to both parts of the Bipartition in that edge, and add the resulting edges after doing that to the vertices of the graph. 
        Vector<edgeVertex> VertexList = new Vector<edgeVertex>();
        Iterator<PhyloTreeEdge> edgesIter = Tools.myVectorClonePhyloTreeEdge(T.getEdges()).iterator();
        
        //counter to assign IDs from 0 to vertexNum - 1. 
        int VertexCount = 0;
        
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
            //Bitset that indicates how the original leaves appear in the partition of e, but with the positions in the complete leaf set
            BitSet moldPartition = new BitSet(); 
            for (int i = 0; i < oLeafSet.size(); i++){
                if (e.getOriginalEdge().getPartition().get(i)){
                    moldPartition.set(orgLeaves2compLeaves[i]);
                }
            }
            //System.out.println("Mold Partition for " + VertexCount + " is : " + moldPartition);
            //Loop on all the ways to add the new leaves to moldPartition
			for (int bitc = 0; bitc < Math.pow(2,l); bitc++){
                BitSet tempPartition = (BitSet) moldPartition.clone();
                int rest = bitc;
                //Adding the leaves in the bipartition in concordance to bitc (bitc to binary to decide which side)
                for (int j = 0; j < l; j++){
                    if(rest % 2 != 0){
                        tempPartition.set(listAddedLeaves[j]);
                    }
                    rest = rest/2; 
                }
                //Adding the new edge created by adding those leaves to the Vertex list, with proper ID.
                if (tempPartition.get(cLeafSet.size()-1)){//We always use the split not including the last leaf (considered the root as per Megan's Owen code) as the representative
                    tempPartition.flip(0,cLeafSet.size());
                }
                
                //System.out.println("One temp Partition for " + VertexCount + " is : " + tempPartition);
                edgeVertex newV = new edgeVertex(VertexCount, new Bipartition(tempPartition));
                VertexList.add(newV);
                VertexCount++;
            }
		}
        
        adjVertices = new HashMap<>();
        orderedVertices = new Vector<edgeVertex>();
        
        //We iterate through all edges in the tree to find which are non-compatible, and connect the respective vertices in the graph.
        Iterator<edgeVertex> vertexIter = VertexList.iterator();
        while (vertexIter.hasNext()){
            edgeVertex eKey = (edgeVertex) vertexIter.next();
            Iterator<edgeVertex> vertexIter2 = VertexList.iterator();
            List<edgeVertex> tempAdjList = new ArrayList<edgeVertex>();
            while(vertexIter2.hasNext()){
                edgeVertex potAdj = (edgeVertex) vertexIter2.next();
                //System.out.println("For the eKey = " + eKey.getEdge().toString());
                //System.out.println("And the potAdj = " + potAdj.getEdge().toString());
                //System.out.println("Disjoint part = " + eKey.getEdge().disjointFrom(potAdj.getEdge()));
                //System.out.println("Crosses according to my own function: " + Crosses(eKey.getEdge(), potAdj.getEdge(), cLeafSet.size()));
                if (Crosses(eKey.getEdge(), potAdj.getEdge(), cLeafSet.size())){
                    //System.out.println("It crosses");
                    tempAdjList.add(potAdj);
                }
            }
            
            adjVertices.put(eKey, tempAdjList);
            orderedVertices.add(eKey);
        }
        
        
        //We initialize the vector of all Maximal Independent Sets, but it will be filled up in an independent function.
        MIS = new ArrayList<Vector<Integer>>();
        
    } //End of first constructor
    
    //Constructor 2: allows for unrestricted case
    public edgeCrossGraph(PhyloTree T, Vector<String> cLeafSet, boolean restricted){
        
        //Vector indicating which leaves in the completeLeafSet are part of the Original Tree
        BitSet originalLeaves = new BitSet(cLeafSet.size());
        
        //Vector of the original tree's leaf set.
        Vector<String> oLeafSet = polyAlg.Tools.myVectorCloneString(T.getLeaf2NumMap()); 
        
        //Vector that maps the leaves in the original tree to the leaves in the complete leaf set. 
        int[] orgLeaves2compLeaves = new int[oLeafSet.size()]; 
        
        //Number of internal edges in the original tree. 
        int m = T.getEdges().size(); 
        
        //Number of leaves to be added. 
        int l = cLeafSet.size() - oLeafSet.size(); 
        
        //In the restricted case, the number of potentially new internal edges is the number of internal edges in the original tree times the different ways the l extra leaves can be added to this edges.
        //in the unrestricted case, we also obtain those edges obtained by adding to the external edges in the original tree, minus those that end up being the external edge again. 
        this.vertexNum = (int) (m*Math.pow(2,l));  
        if (!restricted){
            this.vertexNum += (int) (oLeafSet.size()*(Math.pow(2,l) - 1)); 
            this.vertexNum += (int) Math.pow(2,l) - l - 1;
        }
        
        //Array of the leaves that are being added
        int[] listAddedLeaves = new int[l];
        
        //For each leaf in the original tree's leaf set, we find its position in the complete set, and modify the arrays and BitSet accordingly. 
        for (int i = 0; i < oLeafSet.size(); i++){
            int temp = cLeafSet.indexOf(oLeafSet.get(i));
            if (temp == -1){
                System.err.println("Error: The original tree has a leaf that is not part of the complete leaf set");
			    System.exit(1);
            }
            orgLeaves2compLeaves[i] = temp;
            originalLeaves.set(temp);
            
        }
        
        //We iterate over the set of complete leaves to find which leafs must be added to the original tree.
        int count = 0;
        for (int i = 0; i < cLeafSet.size(); i++){
            if (!originalLeaves.get(i)){
                listAddedLeaves[count] = i;
                count++;
            }
        }
        
        //For each edge in the original tree, we loop through all the ways the extra leaves can be added to both parts of the Bipartition in that edge, and add the resulting edges after doing that to the vertices of the graph. 
        Vector<edgeVertex> VertexList = new Vector<edgeVertex>();
        Iterator<PhyloTreeEdge> edgesIter = Tools.myVectorClonePhyloTreeEdge(T.getEdges()).iterator();
        
        //counter to assign IDs from 0 to vertexNum - 1. 
        int VertexCount = 0;
        
		while (edgesIter.hasNext()){
			PhyloTreeEdge e = (PhyloTreeEdge) edgesIter.next();
            //Bitset that indicates how the original leaves appear in the partition of e, but with the positions in the complete leaf set
            BitSet moldPartition = new BitSet(); 
            for (int i = 0; i < oLeafSet.size(); i++){
                if (e.getOriginalEdge().getPartition().get(i)){
                    moldPartition.set(orgLeaves2compLeaves[i]);
                }
            }
            //Loop on all the ways to add the new leaves to moldPartition
			for (int bitc = 0; bitc < Math.pow(2,l); bitc++){
                BitSet tempPartition = (BitSet) moldPartition.clone();
                int rest = bitc;
                //Adding the leaves in the bipartition in concordance to bitc (bitc to binary to decide which side)
                for (int j = 0; j < l; j++){
                    if(rest % 2 != 0){
                        tempPartition.set(listAddedLeaves[j]);
                    }
                    rest = rest/2; 
                }
                //Adding the new edge created by adding those leaves to the Vertex list, with proper ID.
                if (tempPartition.get(cLeafSet.size()-1)){//We always use the split not including the last leaf (considered the root as per Megan's Owen code) as the representative
                    tempPartition.flip(0,cLeafSet.size());
                }
                edgeVertex newV = new edgeVertex(VertexCount, new Bipartition(tempPartition));
                VertexList.add(newV);
                VertexCount++;
            }
		}
        
        //We are adding extra edges in the unrestricted case: 
        
        if (!restricted){
            for(int i = 0; i < oLeafSet.size(); i++){
                BitSet moldPartition = new BitSet(); 
                moldPartition.set(orgLeaves2compLeaves[i]);
                for(int bitc = 1; bitc < Math.pow(2,l); bitc++){//We start in 1 to ensure at least one leaf is added to the side of the leaf in the external branch
                    BitSet tempPartition = (BitSet) moldPartition.clone();
                    int rest = bitc;
                    //Adding the leaves in the bipartition in concordance to bitc (bitc to binary to decide which side)
                    for (int j = 0; j < l; j++){
                        if(rest % 2 != 0){
                            tempPartition.set(listAddedLeaves[j]);
                        }
                        rest = rest/2; 
                    }
                    //Adding the new edge created by adding those leaves to the Vertex list, with proper ID.
                    if (tempPartition.get(cLeafSet.size()-1)){//We always use the split not including the last leaf (considered the root as per Megan's Owen code) as the representative
                        tempPartition.flip(0,cLeafSet.size());
                    }
                    
                    edgeVertex newV = new edgeVertex(VertexCount, new Bipartition(tempPartition));
                    VertexList.add(newV);
                    VertexCount++;
                }
            }
            
            BitSet moldEmptyPartition = new BitSet();
            for (int pow2 = 1; pow2 < l; pow2++){//We are using the formula 2^pow2 + bitext to generate all bitc numbers that have at least 2 ones in the binary representation
                //System.out.println("pow 2 = "+ pow2);
                int vPow2 = (int) Math.pow(2,pow2);
                //System.out.println("V pow 2 = "+ vPow2);
                for (int bitext = 1; bitext < vPow2; bitext++){
                    int rest = vPow2 + bitext;
                    //System.out.println("Rest = "+ rest);
                    BitSet tempPartition = (BitSet) moldEmptyPartition.clone();
                    for (int j = 0; j < l; j++){
                        if(rest % 2 != 0){
                            tempPartition.set(listAddedLeaves[j]);
                        }
                        rest = rest/2; 
                    }
                    //Adding the new edge created by adding those leaves to the Vertex list, with proper ID.
                    if (tempPartition.get(cLeafSet.size()-1)){//We always use the split not including the last leaf (considered the root as per Megan's Owen code) as the representative
                        tempPartition.flip(0,cLeafSet.size());
                    }
                    
                    edgeVertex newV = new edgeVertex(VertexCount, new Bipartition(tempPartition));
                    VertexList.add(newV);
                    VertexCount++;
                }
            }
        }
        
        adjVertices = new HashMap<>();
        orderedVertices = new Vector<edgeVertex>();
        
        //We iterate through all edges in the tree to find which are non-compatible, and connect the respective vertices in the graph.
        Iterator<edgeVertex> vertexIter = VertexList.iterator();
        while (vertexIter.hasNext()){
            edgeVertex eKey = (edgeVertex) vertexIter.next();
            Iterator<edgeVertex> vertexIter2 = VertexList.iterator();
            List<edgeVertex> tempAdjList = new ArrayList<edgeVertex>();
            while(vertexIter2.hasNext()){
                edgeVertex potAdj = (edgeVertex) vertexIter2.next();
                if (Crosses(eKey.getEdge(), potAdj.getEdge(), cLeafSet.size())){
                    tempAdjList.add(potAdj);
                }
            }
            
            adjVertices.put(eKey, tempAdjList);
            orderedVertices.add(eKey);
        }
        
        //We initialize the vector of all Maximal Independent Sets, but it will be filled up in an independent function.
        MIS = new ArrayList<Vector<Integer>>();
        
    }
    
    //COMPUTING ALL MAXIMAL INDEPENDENT SETS
    
    //Function that is used to build the MIS, described in  Tsukiyama, S., Ide, M., Ariyoshi, I., Shirakawa, I. (1977).
    private void Backtrack(int i){
        if (i < vertexNum-1){
            int x = i+1;
            int c = 0;
            
            List<edgeVertex> adjIter = this.adjVertices.get(this.orderedVertices.get(x));
            for(edgeVertex e_y : adjIter){
                int y = e_y.getID();
                if(y<=i){
                    if(IS[y] == 0){c = c + 1;}
                }
            }
            
            if (c==0){
                List<edgeVertex> adjIter2 = this.adjVertices.get(this.orderedVertices.get(x));
                for(edgeVertex e_y : adjIter2){
                    int y = e_y.getID();
                    if(y <= i){IS[y] = IS[y] + 1;}
                    
                }
                
                Backtrack(x);
                
                List<edgeVertex> adjIter3 = this.adjVertices.get(this.orderedVertices.get(x));
                for(edgeVertex e_y : adjIter3){
                    int y = e_y.getID();
                    if(y <= i){IS[y] = IS[y] - 1;}   
                }
            }else{
                IS[x] = c;
                Backtrack(x);
                IS[x] = 0;
                
                boolean fo = true;
                List<edgeVertex> adjIter2 = this.adjVertices.get(this.orderedVertices.get(x));
                for(edgeVertex e_y : adjIter2){
                    int y = e_y.getID();
                    if(y <= i){
                        if (IS[y] == 0){
                            Bucket.get(x).add(y);
                            List<edgeVertex> adjIter3 = this.adjVertices.get(this.orderedVertices.get(y));
                            for (edgeVertex e_z: adjIter3){
                                int z = e_z.getID(); 
                                if (z <= i){
                                    IS[z] = IS[z] - 1;
                                    if (IS[z] == 0){fo = false;}
                                }
                            }
                        }
                        IS[y] = IS[y] + 1;
                    }
                }
                if (fo) {Backtrack(x);}
                
                List<edgeVertex> adjIter4 = this.adjVertices.get(this.orderedVertices.get(x));
                for(edgeVertex e_y : adjIter4){
                    int y = e_y.getID();
                    if (y <= i){IS[y] = IS[y] - 1;}
                }
                
                Integer[] bucketIter = new Integer[this.Bucket.get(x).size()];
                bucketIter = this.Bucket.get(x).toArray(bucketIter);
                for (int j = 0; j < bucketIter.length; j++){//(bucketIter.hasNext()){
                    int y = bucketIter[j];
                    List<edgeVertex> adjIter5 = this.adjVertices.get(this.orderedVertices.get(y));
                    for(edgeVertex e_z : adjIter5){
                        int z = e_z.getID(); 
                        if (z<=i) {IS[z] = IS[z] + 1;}
                    }
                    this.Bucket.get(x).remove(y);
                }
                
            }
            
        } else {
            Vector<Integer> tempMIS = new Vector<Integer>();
            for (int j = 0; j < vertexNum; j++){
                if (IS[j] == 0){
                    tempMIS.add(j);
                }
            }
            MIS.add(tempMIS);
        }
    }
    
    
    //Important construction: Finding the maximal independent sets, which is the possible orthants. 
    //Function that starts the recursive function described in described in Tsukiyama, S., Ide, M., Ariyoshi, I., Shirakawa, I. (1977).
    public void MIScalculator(){
        IS = new int[vertexNum];
        Bucket = new ArrayList<Set<Integer>>();
        
        for (int j = 0; j < vertexNum; j++){
            IS[j] = 0;
            Bucket.add(new HashSet<Integer>());
        }
        
        Backtrack(0);
        
    }
    
    //Printers and Getters
    
    public void PrintMap(Vector<String> cLeafSet){
        System.out.println("This graph has " + this.vertexNum + " vertices: ");
        Iterator<edgeVertex> keyIter = this.orderedVertices.iterator();
        while (keyIter.hasNext()){
            edgeVertex eKey = (edgeVertex) keyIter.next();
            System.out.println(eKey.getID() + ": "+eKey.getEdge().toStringVerbose(cLeafSet));
            System.out.print("  Adjacent to: ");
            Iterator<edgeVertex> adjIter = this.adjVertices.get(eKey).iterator();
            while (adjIter.hasNext()){
                edgeVertex adjE = (edgeVertex) adjIter.next();
                System.out.print(" " + adjE.getID() + " ");
            }
            System.out.println("");
        }
    }
    
    public List<Vector<Integer>> getMIS(){
        return MIS;
    }
    
    public Vector<edgeVertex> getOV(){
        return orderedVertices;
    }
    
}