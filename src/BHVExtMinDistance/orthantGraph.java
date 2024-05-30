/** This is intended as the class defining the a graph connecting NNI orthants, all belonging to extension space of a tree. The orthants represented by the vertices in this graph are all those that may be obtained by attaching new leaves to a tree. 

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;

public class orthantGraph{
    //Adjacency lists that define the graph
    private Map<orthantVertex, List<orthantVertex>> adjVertices;
    //List of the vertices in the graph, so we always have the i-th vertex with the ID equal to i.
    private Vector<orthantVertex> orderedVertices;
    //Number of vertices in the graph;
    private int vertexNum;
    
    //Constructor (We use the Maximal Independen Sets created for a edgeCrossGraph)
    public orthantGraph(edgeCrossGraph edgeCG){
        //Get the edgeVertex and the Maximal Independent Sets from the edgeCrossGraph passed as parameter.
        Vector<edgeVertex> orderedEdges = edgeCG.getOV();
        List<Vector<Integer>> MIS = edgeCG.getMIS();
        
        adjVertices = new HashMap<>();
        orderedVertices = new Vector<orthantVertex>();
        //The number of orthants are the same as the number of Maximal Independent Sets.
        vertexNum = MIS.size();
        
        
        //For every Maximal Independent Set, we find the bipartitions of the edges inside that set, add to the axis of the orthant represented by these edges, and add this as a new vertex of our graph.
        for(int i=0; i < vertexNum; i++){
            Vector<Bipartition> tempOAxes = new Vector<Bipartition>();
            for (int j : MIS.get(i)){
                tempOAxes.add(orderedEdges.get(j).getEdge().clone());
            }
            orthantVertex newVertex = new orthantVertex(i, tempOAxes);
            orderedVertices.add(newVertex);
        }
        
        //We iterate through all pairs of vertices, and find whichones share all but one edges, to find which are neighbours by rotation
        for(int i=0; i < vertexNum; i++){
            List<orthantVertex> tempAdjList = new ArrayList<orthantVertex>();
            for(int j = 0; j < vertexNum; j++){
                Vector<Integer> tempMISline = (Vector<Integer>) MIS.get(j).clone();
                tempMISline.retainAll(MIS.get(i));
                if (tempMISline.size() == MIS.get(i).size()-1){
                    tempAdjList.add(orderedVertices.get(j));
                }
            }
            adjVertices.put(orderedVertices.get(i), tempAdjList);
        }
    }
    
    //Printers and Getters
    public void Print(Vector<String> cLeafSet){
        System.out.println("This graph has " + this.vertexNum + " vertices: ");
        Iterator<orthantVertex> keyIter = this.orderedVertices.iterator();
        while (keyIter.hasNext()){
            orthantVertex eKey = (orthantVertex) keyIter.next();
            System.out.print(eKey.getID() + ": {");
            for (Bipartition bipe : eKey.getOrthantAxis()){
                System.out.print(" {"+bipe.toStringVerbose(cLeafSet)+"}");
            }
            System.out.println("};");
            System.out.print("   Adjacent to: ");
            Iterator<orthantVertex> adjIter = this.adjVertices.get(eKey).iterator();
            while (adjIter.hasNext()){
                orthantVertex adjE = (orthantVertex) adjIter.next();
                System.out.print(" " + adjE.getID() + " ");
            }
            System.out.println("");
        }
    }
    
    public Vector<orthantVertex> getOV(){
        return this.orderedVertices;
    }
    
    public Vector<Bipartition> getAxes(int k){
        return this.orderedVertices.get(k).getOrthantAxis();
    }
    
    public Vector<Bipartition> getAxesClone(int k){
        return this.orderedVertices.get(k).getOrthantAxisClone();
    }
    
    public int getVertexNum(){
        return this.vertexNum;
    }
    
    //Get the IDs of the neighbour orthants to the orthant in the k-th vertex of the graph as an Array.
    public int[] getAdjIDs(int k){
        orthantVertex tempKey = this.orderedVertices.get(k);
        List<orthantVertex> tempAdj = this.adjVertices.get(tempKey);
        int[] ReturnValue = new int[tempAdj.size()]; 
        for(int i = 0; i< tempAdj.size(); i++){
            ReturnValue[i] = tempAdj.get(i).getID();
        }
        
        return(ReturnValue);
    }
    
    
    //Get the IDs of the neighbour orthants to the orthant in the k-th vertex of the graph as a List.
    public List<Integer> getAdjIDsList(int k){
        orthantVertex tempKey = this.orderedVertices.get(k);
        List<orthantVertex> tempAdj = this.adjVertices.get(tempKey);
        List<Integer> ReturnValue = new ArrayList<Integer>(); 
        for(orthantVertex e : tempAdj){
            ReturnValue.add(e.getID());
        }
        
        return(ReturnValue);
    }
    
}