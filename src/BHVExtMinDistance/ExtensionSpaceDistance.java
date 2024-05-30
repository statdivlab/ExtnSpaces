/** This is intended as the class defining the distance between two extension spaces. It will include methods to compute this distance, as well as a list of pairs of trees that obtain smaller distances per each maximal orthant covered by these extension spaces.

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/

package BHVExtMinDistance;

import java.util.*;
import java.util.concurrent.*; 
import java.util.concurrent.AbstractExecutorService;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.stream.Collectors;
import distanceAlg1.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ExtensionSpaceDistance{
    //List of all the distances between the subset of the extension spaces restricted to particular orthants.
    //This list will be ordered from the shorter distance to the higher.
    private List<OrthExtDistance> orderedOrthExtDistances;
    //List of the ID's of the orthants in the first and second extension space that are involved in the respective orthant extension distance
    //private OrthExtDistance[] DistancesTemp;
    //The two trees belonging to each of the extension spaces that produces the smaller distance possible.
    private PhyloTree bestTree1;
    private PhyloTree bestTree2;
    //The distance (smaller of the orthant extension distances) between the extension spaces.
    private double Distance;
    //The geodesic between the extension spaces that produces the smaller distance.
    //private Geodesic bestGeode;
    
    //Constructor
    public ExtensionSpaceDistance(ExtensionSpace ES1, ExtensionSpace ES2){
        orderedOrthExtDistances = new ArrayList<OrthExtDistance>();
        
        Vector<OrthExt> OEs1 = ES1.getOrthExts();
        Vector<OrthExt> OEs2 = ES2.getOrthExts();
        
        int oNum1 = OEs1.size();
        int oNum2 = OEs2.size();
        
        //For each pair of orthant extensions in the extension spaces we compute the Orthant Extension Distances in between them, find how it compares to the other distances already added to the list, and we add it to the correct position, also adding the orthants that produced this distance to each the list of orthants. 
        for (int k1 = 0; k1 < oNum1; k1++){
            OrthExt OE1 = OEs1.get(k1);
            for (int k2 = 0; k2 < oNum2; k2++){
                //System.out.println("************");
                //System.out.println("STARTING O pair ("+k1+", "+k2+")");
                OrthExt OE2 = OEs2.get(k2);
                OrthExtDistance tempOED = new OrthExtDistance(OE1, OE2);
                //System.out.println("THE DISTANCE WAS "+ tempOED.getDistance());
                //System.out.println("************");
                //System.out.println("");
                
                orderedOrthExtDistances.add(tempOED);
                
                
                /*boolean Added = false;
                for(int i = 0; i < orderedOrthExtDistances.size(); i++){
                    if (tempOED.getDistance() < orderedOrthExtDistances.get(i).getDistance()){
                        Added = true;
                        orderedOrthExtDistances.add(i, tempOED);
                        break;
                    }
                }
                if (!Added){
                    orderedOrthExtDistances.add(orderedOrthExtDistances.size(), tempOED);
                }*/
                
                /*if ((orderedOrthExtDistances.size()>0) && (tempOED.getDistance() <= orderedOrthExtDistances.get(0).getDistance())){
                    orderedOrthExtDistances.add(0, tempOED);
                } else {
                    orderedOrthExtDistances.add(tempOED);
                }*/
            }
        }
        
        Collections.sort(orderedOrthExtDistances, Comparator.comparingDouble(OrthExtDistance::getDistance));
        //The best trees, distance and geodesic will be those at the beginning of our list. 
        bestTree1 = orderedOrthExtDistances.get(0).getFirstTree();
        bestTree2 = orderedOrthExtDistances.get(0).getSecondTree();
        Distance = orderedOrthExtDistances.get(0).getDistance();
        //bestGeode = orderedOrthExtDistances.get(0).getFinalGeode();
    }//end of constructor 
    
    //Constructor 2: allowing for unrestricted version
    public ExtensionSpaceDistance(ExtensionSpace ES1, ExtensionSpace ES2, boolean restricted){
        orderedOrthExtDistances = new ArrayList<OrthExtDistance>();
        
        Vector<OrthExt> OEs1 = ES1.getOrthExts();
        Vector<OrthExt> OEs2 = ES2.getOrthExts();
        
        int oNum1 = OEs1.size();
        int oNum2 = OEs2.size();
        
        //For each pair of orthant extensions in the extension spaces we compute the Orthant Extension Distances in between them, find how it compares to the other distances already added to the list, and we add it to the correct position, also adding the orthants that produced this distance to each the list of orthants. 
        for (int k1 = 0; k1 < oNum1; k1++){// k1 < oNum1; k1++){
            OrthExt OE1 = OEs1.get(k1);
            for (int k2 = 0; k2 < oNum2; k2++){//k2 < oNum2; k2++){
                OrthExt OE2 = OEs2.get(k2);
                //System.out.println("************");
                //System.out.println("STARTING O pair ("+k1+", "+k2+")");
                //OE1.printLN();
                //OE2.printLN();
                OrthExtDistance tempOED = new OrthExtDistance(OE1, OE2, restricted);
                //System.out.println("************");
                //System.out.println("");
                //if (tempOED.getWarning() != null){
                    //System.out.println("************");
                    //System.out.println("STARTING O pair ("+k1+", "+k2+")");
                    //System.out.println(tempOED.getWarning());
                    //System.out.println("************");
                    //System.out.println("");
                //}
                //long End = System.currentTimeMillis();
                //double TimeSeconds = ((double)(End - Start))/1000;
                //System.out.println("THE DISTANCE WAS "+ tempOED.getDistance());
                
                //System.out.println("Time needed: " + TimeSeconds);
                
                /*boolean Added = false;
                for(int i = 0; i < orderedOrthExtDistances.size(); i++){
                    if (tempOED.getDistance() < orderedOrthExtDistances.get(i).getDistance()){
                        Added = true;
                        orderedOrthExtDistances.add(i, tempOED);
                        break;
                    }
                }
                if (!Added){
                    orderedOrthExtDistances.add(orderedOrthExtDistances.size(), tempOED);
                }*/
                
                /*if ((orderedOrthExtDistances.size()>0) && (tempOED.getDistance() <= orderedOrthExtDistances.get(0).getDistance())){
                    orderedOrthExtDistances.add(0, tempOED);
                } else {
                    orderedOrthExtDistances.add(orderedOrthExtDistances.size(),tempOED);
                }*/
                orderedOrthExtDistances.add(tempOED);
            }
        }
        
        Collections.sort(orderedOrthExtDistances, Comparator.comparingDouble(OrthExtDistance::getDistance));
        //The best trees, distance and geodesic will be those at the beginning of our list. 
        bestTree1 = orderedOrthExtDistances.get(0).getFirstTree();
        bestTree2 = orderedOrthExtDistances.get(0).getSecondTree();
        Distance = orderedOrthExtDistances.get(0).getDistance();
        //bestGeode = orderedOrthExtDistances.get(0).getFinalGeode();
    } //end of constructor 2
    
    //Constructor 3: allowing for unrestricted version and parallelizing 
    /*public List<OrthExtDistance> ParallelComputation(Vector<OrthExt> OESS1, Vector<OrthExt> OESS2, boolean restricted, int numT) throws InterruptedException, ExecutionException {
        ExecutorService service = Executors.newFixedThreadPool(numT);
        
        List<Future<OrthExtDistance>> futures = new ArrayList<Future<OrthExtDistance>>();
        for (final OrthExt OE1 : OESS1){
            for (final OrthExt OE2 : OESS2){
                Callable<OrthExtDistance> callable = new Callable<OrthExtDistance>(){
                    public OrthExtDistance call() throws Exception{
                        OrthExtDistance output = new OrthExtDistance(OE1, OE2, restricted);
                        return output;
                    }
                };
                futures.add(service.submit(callable));
            }
        }
        
        service.shutdown();

        List<OrthExtDistance> outputs = new ArrayList<OrthExtDistance>();
        for (Future<OrthExtDistance> future : futures) {
            outputs.add(future.get());
        }
        return outputs;
    }*/
    
    private void processParallelyWithExecutorService(List<orthantExtPair> Input, int numT){// throws InterruptedException {
    //ExecutorService executorService = Executors.newFixedThreadPool(numThreads);
    //List<CompletableFuture<Void>> futures = new ArrayList<>();
    
    final ExecutorService executor = Executors.newFixedThreadPool(numT);
    List<Future<OrthExtDistance>> futures = new ArrayList<>();
    
    for (int i = 0; i < Input.size(); i++) {
        final orthantExtPair inpT = Input.get(i);
        try{
            Future<OrthExtDistance> future = executor.submit(new callableOED(inpT));
            futures.add(future);
        } catch (Exception e){
            System.out.println("This happened " + i);
            System.out.println(e.getCause());
        }
        
            
            /*CompletableFuture.runAsync(() -> {
            //try {
                DistancesTemp[i] = new OrthExtDistance(inpT, false);
            /*} catch (InterruptedException e) {
                e.printStackTrace();
            }
        }, executorService);*/
    
    }
    //System.out.println("The length of futures is: " + futures.size());
    int tempCount = 0;
    for (Future<OrthExtDistance> future : futures) {
        tempCount++;
        try{
            orderedOrthExtDistances.add(future.get());   
        } catch (Exception e){
            System.out.println("Catched " + tempCount);
            System.out.println(e);
        }
    }
    /*try {
        for (Future<OrthExtDistance> future : futures) {
            orderedOrthExtDistances.add(future.get()); // do anything you need, e.g. isDone(), ...
        }
    } catch (Exception e) {
        System.out.println("And this");
        System.out.println(e);
    }*/
    executor.shutdown();
}
    
    public ExtensionSpaceDistance(ExtensionSpace ES1, ExtensionSpace ES2, boolean restricted, int numThreads){
        orderedOrthExtDistances = new ArrayList<OrthExtDistance>();
        
        Vector<OrthExt> OEs1 = ES1.getOrthExts();
        Vector<OrthExt> OEs2 = ES2.getOrthExts();
        
        
        List<orthantExtPair> OEpairs = new ArrayList<orthantExtPair>();
        
        //System.out.println("About to create list of pairs to run in 'Parallel'");
        
        for (OrthExt OE1 : OEs1){
            for (OrthExt OE2 : OEs2){
                OEpairs.add(new orthantExtPair(OE1, OE2));
            }
        }
        
        //DistancesTemp = new OrthExtDistance[OEpairs.size()];
        
        //System.out.println("About to run it in 'Parallel'");
        processParallelyWithExecutorService(OEpairs, numThreads);
        //orderedOrthExtDistances = OEs1.parallelStream().flatMap(OE1 -> OEs2.parallelStream().map(OE2 -> new OrthExtDistance(OE1, OE2, restricted))).collect(Collectors.toList());
        
        //orderedOrthExtDistances = OEpairs.parallelStream().map(OEpair -> new OrthExtDistance(OEpair, false)).collect(Collectors.toList());
        
        //for (orthantExtPair OEpar : OEpairs){
        //    orderedOrthExtDistances.add(new OrthExtDistance(OEpar, restricted));
        //}

        
        Collections.sort(orderedOrthExtDistances, Comparator.comparingDouble(OrthExtDistance::getDistance));
        //The best trees, distance and geodesic will be those at the beginning of our list. 
        bestTree1 = orderedOrthExtDistances.get(0).getFirstTree();
        bestTree2 = orderedOrthExtDistances.get(0).getSecondTree();
        Distance = orderedOrthExtDistances.get(0).getDistance();
        //bestGeode = orderedOrthExtDistances.get(0).getFinalGeode();
    } //end of constructor 3*/
    
    //Printers and Getters
    
    public List<OrthExtDistance> getOOED(){
        return orderedOrthExtDistances;
    }
    
    public PhyloTree getBestTree1(){
        return bestTree1;
    }
    
    public PhyloTree getBestTree2(){
        return bestTree2;
    }
    
    public double getDistance(){
        return Distance;
    }
    
    public void PrintSummary(boolean withTrees, boolean withIterCount){
        System.out.println("There are a total of "+ orderedOrthExtDistances.size()+" orthant pairs");
        for(int i = 0; i < orderedOrthExtDistances.size(); i++){
            System.out.println("Pair "+ (i+1) +": (" + orderedOrthExtDistances.get(i).getO1ID() + ", " + orderedOrthExtDistances.get(i).getO2ID() + ")\n The distance is " + orderedOrthExtDistances.get(i).getDistance());
            if (withTrees){
                PhyloNicePrinter treePrinter = new PhyloNicePrinter();
                System.out.println("  Best Tree 1: ");
                System.out.println(treePrinter.toString(orderedOrthExtDistances.get(i).getFirstTree()));
                System.out.println("");
                System.out.println("  Best Tree 2: ");
                System.out.println(treePrinter.toString(orderedOrthExtDistances.get(i).getSecondTree()));
            }
            if (withIterCount){
                System.out.println("   Number of Iterations: " + orderedOrthExtDistances.get(i).getIterCount());
            }
            System.out.println("---------------------------------------------------------------");
        }
    }
    
    public void PrintSummary(boolean withTrees, boolean withIterCount, int NumOrthants){
        if (NumOrthants < 1){
            System.out.println("The number of orthant pairs to be printed should be at least 1");
        }
        System.out.println("There are a total of "+ orderedOrthExtDistances.size()+" orthant pairs");
        for(int i = 0; i < NumOrthants; i++){
            System.out.println("Pair "+ (i+1) +": (" + orderedOrthExtDistances.get(i).getO1ID() + ", " + orderedOrthExtDistances.get(i).getO2ID() + ")\n The distance is " + orderedOrthExtDistances.get(i).getDistance());
            if (withTrees){
                PhyloNicePrinter treePrinter = new PhyloNicePrinter();
                System.out.println("  Best Tree 1: ");
                System.out.println(treePrinter.toString(orderedOrthExtDistances.get(i).getFirstTree()));
                System.out.println("");
                System.out.println("  Best Tree 2: ");
                System.out.println(treePrinter.toString(orderedOrthExtDistances.get(i).getSecondTree()));
            }
            if (withIterCount){
                System.out.println("   Number of Iterations: " + orderedOrthExtDistances.get(i).getIterCount());
            }
            System.out.println("---------------------------------------------------------------");
        }
    }
    
    public void PrintSummary(boolean withTrees, boolean withIterCount, String fileName){
        try {
            File myObj = new File(fileName);
            if (myObj.createNewFile()) {
                System.out.println("Report created in: " + myObj.getName());
            } else {
                System.out.println("Report already exists in: " + myObj.getName());
            }
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        
        try {
            FileWriter myWriter = new FileWriter(fileName);
            myWriter.write("There are a total of "+ orderedOrthExtDistances.size()+" orthant pairs \n");
            for(int i = 0; i < orderedOrthExtDistances.size(); i++){
                myWriter.write("Pair "+ (i+1) +": (" + orderedOrthExtDistances.get(i).getO1ID() + ", " + orderedOrthExtDistances.get(i).getO2ID() + ")\n The distance is " + orderedOrthExtDistances.get(i).getDistance() + "\n");
                if(withTrees){
                    PhyloNicePrinter treePrinter = new PhyloNicePrinter();
                    myWriter.write("  Best Tree 1: \n" + treePrinter.toString(orderedOrthExtDistances.get(i).getFirstTree()) + "\n \n");
                    myWriter.write("  Best Tree 2: \n" + treePrinter.toString(orderedOrthExtDistances.get(i).getSecondTree()) + "\n");   
                }
                if (withIterCount){
                    myWriter.write("   Number of Iterations: " + orderedOrthExtDistances.get(i).getIterCount());
                }
                myWriter.write("---------------------------------------------------------------\n");
                myWriter.write(" \n");
            }
            myWriter.close();
            System.out.println("Successfully wrote the report.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        
    }
    
    public void PrintSummary(boolean withTrees, boolean withIterCount, String fileName, int NumOrthants){
        try {
            File myObj = new File(fileName);
            if (myObj.createNewFile()) {
                System.out.println("Report created in: " + myObj.getName());
            } else {
                System.out.println("Report already exists in: " + myObj.getName());
            }
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        
        try {
            FileWriter myWriter = new FileWriter(fileName);
            if (NumOrthants < 1){
                myWriter.write("The number of orthant pairs to be printed should be at least 1");
            }
            myWriter.write("There are a total of "+ orderedOrthExtDistances.size()+" orthant pairs \n");
            for(int i = 0; i < NumOrthants; i++){
                myWriter.write("Pair "+ (i+1) +": (" + orderedOrthExtDistances.get(i).getO1ID() + ", " + orderedOrthExtDistances.get(i).getO2ID() + ")\n The distance is " + orderedOrthExtDistances.get(i).getDistance() + "\n");
                if(withTrees){
                    PhyloNicePrinter treePrinter = new PhyloNicePrinter();
                    myWriter.write("  Best Tree 1: \n" + treePrinter.toString(orderedOrthExtDistances.get(i).getFirstTree()) + "\n \n");
                    myWriter.write("  Best Tree 2: \n" + treePrinter.toString(orderedOrthExtDistances.get(i).getSecondTree()) + "\n");   
                }
                if (withIterCount){
                    myWriter.write("   Number of Iterations: " + orderedOrthExtDistances.get(i).getIterCount());
                }
                myWriter.write("---------------------------------------------------------------\n");
                myWriter.write(" \n");
            }
            myWriter.close();
            System.out.println("Successfully wrote the report.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        
    }
    
    // Graphical representations to explore the relationship between NNI and shorter distances. 
    
    //We are creating a graph that shows, for each pair of orthants, which are "adjacent", which means: one of the orthants is identical to the other for one of the extension spaces, and the other ones in the pair are neighbours by rotation.
    public Map<Integer, List<Integer>> JointNNI(ExtensionSpace ES1, ExtensionSpace ES2){
        orthantGraph connectCluster1 = ES1.getConnectCluster();
        orthantGraph connectCluster2 = ES2.getConnectCluster();
        
        Map<Integer, List<Integer>> ReturnValue = new HashMap<>();
        
        for(int i = 0; i < orderedOrthExtDistances.size(); i++){
            List<Integer> adjTemp = new ArrayList<Integer>();
            for(int j = 0; j < orderedOrthExtDistances.size(); j++){
                if(((orderedOrthExtDistances.get(i).getO1ID() == orderedOrthExtDistances.get(j).getO1ID()) && (connectCluster2.getAdjIDsList(orderedOrthExtDistances.get(i).getO2ID()).contains(orderedOrthExtDistances.get(j).getO2ID()))) || ((orderedOrthExtDistances.get(i).getO2ID() == orderedOrthExtDistances.get(j).getO2ID()) && (connectCluster1.getAdjIDsList(orderedOrthExtDistances.get(i).getO1ID()).contains(orderedOrthExtDistances.get(j).getO1ID())))){
                    adjTemp.add(j);
                }
            }
            ReturnValue.put(i, adjTemp);
            
        }
        
        return(ReturnValue);
    }
    
    public Map<Integer, List<Integer>> StartTreeNNI(ExtensionSpace ES1){
        orthantGraph connectCluster1 = ES1.getConnectCluster();
        
        Map<Integer, List<Integer>> ReturnValue = new HashMap<>();
        
        for(int i = 0; i < orderedOrthExtDistances.size(); i++){
            List<Integer> adjTemp = new ArrayList<Integer>();
            for(int j = 0; j < orderedOrthExtDistances.size(); j++){
                if((orderedOrthExtDistances.get(i).getO1ID() == orderedOrthExtDistances.get(j).getO1ID()) || (connectCluster1.getAdjIDsList(orderedOrthExtDistances.get(i).getO1ID()).contains(orderedOrthExtDistances.get(j).getO1ID()))){
                    adjTemp.add(j);
                }
            }
            ReturnValue.put(i, adjTemp);
            
        }
        
        return(ReturnValue);
    }
    
    public Map<Integer, List<Integer>> EndTreeNNI(ExtensionSpace ES2){
        orthantGraph connectCluster2 = ES2.getConnectCluster();
        
        Map<Integer, List<Integer>> ReturnValue = new HashMap<>();
        
        for(int i = 0; i < orderedOrthExtDistances.size(); i++){
            List<Integer> adjTemp = new ArrayList<Integer>();
            for(int j = 0; j < orderedOrthExtDistances.size(); j++){
                if((connectCluster2.getAdjIDsList(orderedOrthExtDistances.get(i).getO2ID()).contains(orderedOrthExtDistances.get(j).getO2ID())) || (orderedOrthExtDistances.get(i).getO2ID() == orderedOrthExtDistances.get(j).getO2ID())){
                    adjTemp.add(j);
                }
            }
            ReturnValue.put(i, adjTemp);
            
        }
        
        return(ReturnValue);
    }
    
}