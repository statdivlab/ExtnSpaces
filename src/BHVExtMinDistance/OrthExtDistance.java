/** This is intended as the class defining the distance between two orthant extension spaces. It will include methods to compute this distance, as well was the pair of trees that obtain this smaller distance

Part of the package BHVExtMinDistance and it is constructed using tools from the packages: 
 * distanceAlg1; PolyAlg; constructed by Megan Owen

Part of the package that computes distances between Extension Spaces.
*/


package BHVExtMinDistance;

import java.util.*;
import distanceAlg1.*;
import polyAlg.*;
import java.util.Scanner;

public class OrthExtDistance{
    //Trees in each Orthant extension space from which the shorter geodesic is obtained
    private PhyloTree Tree1;
    private PhyloTree Tree2;
    
    //Shorter distance in between the orthant extension space.
    private double Distance;
    private int IterCount; //Number of iterations used to compute the distance. 
    
    private int O1ID;// ID's of the Orthants in their respective Orthant Extensions
    private int O2ID;
    
    private String Message;
    
    
    
    //The following function is no longer necessary, every since we standarized all edges to be represented by the split in which the root (the last leaf) is 0. We keep it here in case I need to re-utilize it in future versions of the code. 
    /*//Funcion to determine the edge ID on a tree, dealing with the fact that sometimes the edge will be listed as the complement of the edge in the tree. 
    private int edgeIDonT(PhyloTreeEdge e, PhyloTree T, int numberLeaves){
        int reID = T.getSplits().indexOf(e.asSplit());
        
        if (reID == -1){
            Bipartition eClone = e.getOriginalEdge().clone();
            eClone.complement(numberLeaves);
            reID = T.getSplits().indexOf(eClone);
            if (T.getSplits().indexOf(eClone) != -1){
                System.out.println("It entered the non-problem");
            }
        }
        
        return(reID);
    }*/
    
    //Constructor
    private void Constructor1(OrthExt OE1, OrthExt OE2){
        double TolLimit = 0.00000001; 
        PhyloNicePrinter treePrinter = new PhyloNicePrinter();
        //We start by the starting trees in each orthant extension.
        PhyloTree T1 = new PhyloTree(OE1.getStartTree());
        PhyloTree T2 = new PhyloTree(OE2.getStartTree());
        
        //Find the the geodesic in between these trees. 
        Geodesic tempGeode = parPolyMain.getGeodesic(T1, T2);
        
        //System.out.println("");
        //System.out.println("PROCESSING THE DISTANCE ALGORITHM... Starting at distance "+ tempGeode.getDist());
        
        //Some useful counters
        int k1 = OE1.getDim();//Dimension of the first Orthant Extension Space
        int k2 = OE2.getDim();//Dimension of the second Orthant Extension Space
        int n = OE1.getOrthantAxis().size(); //Dimension of the rows in the orthogonal matrices for both extension orthants. This should be the number of interior edges in binary trees with the complete leaf set for both, and should coincide. 
        int m1 = OE1.getFixedLengths().length;
        int m2 = OE2.getFixedLengths().length;
        
        //TO DO: add code to verify both orthant extensions are in fact inside the same BHV tree space. For now, I just assume every user will be careful about this. 
        
        
        //The following while will perform reduced gradient method algorithm, with a conjugate gradient method in each classification of variables. In each iteration the gradient of the "active variables" (those clasified into S1 and S2) function from the current trees is computed, the optimal descent direction is selected following the conjugate gradient method, and the minimum in that direction is computed. If we hit a boundary, we reclasify variables in order to increment those forced to be zero. We continue until finding a semi-stationary point, and corroborate this is the optimum or add new non-basic variables otherwise.
        
        //Initializing the indexes sets B, S and N, with some extra structures to easy change. 
        //THIS COULD POTENTIALLY BE A PART OF OrthExt class TO AVOID IT BEING COMPUTED EVERY TIME A DISTANCE IS COMPUTED
        
        Vector<Integer> B1 = new Vector<Integer>();
        Vector<Integer> B2 = new Vector<Integer>();
        
        Vector<Integer> S1 = new Vector<Integer>();
        Vector<Integer> S2 = new Vector<Integer>();
        
        Vector<Integer> N1 = new Vector<Integer>();
        Vector<Integer> N2 = new Vector<Integer>();
        
        //We will keep a vector of indexes that have already been non-basic variables, to give priority to new potential non-basic variables with possible, trying to prevent cycling. 
        Vector<Integer> alreadyN1 = new Vector<Integer>();
        Vector<Integer> alreadyN2 = new Vector<Integer>();
        
        for (int i = 0; i < m1; i++){//For each row in the map matrix
            Vector<Integer> tempVect = OE1.getMapList().get(i); //Get the edges that merge into the final edge in the original tree
            B1.add(tempVect.get(0)); //Add the first entry of this list of edges into B1
            S1.addAll(tempVect.subList(1,tempVect.size())); //The rest is added to S1
        }
        
        for (int i = 0; i < m2; i++){//For each row in the map matrix for the second extension
            Vector<Integer> tempVect = OE2.getMapList().get(i); //Get the edges that merge into the final edge in the original tree
            B2.add(tempVect.get(0)); //Add the first entry of this list of edges into B1
            S2.addAll(tempVect.subList(1,tempVect.size())); //The rest is added to S1
        }
        
        //Some values before the iterations start
        
        int iterCount = 0; //Counter of the number of iterations performed.  
        
        boolean optimNotReached = true;//We will stop the loop when the gradient is small enough to guarantee we have reach the minimum. 
        
        int conjugate_initial_counter = 0; //Counter for re-initialization of the conjugate gradient method
        
        //We need to keep track on gradients and change directions
        
        double[] gradientxs1 = new double[S1.size()];
        double[] gradientxs2 = new double[S2.size()];
        
        double[] dDirectionxs1 = new double[S1.size()];
        double[] dDirectionxs2 = new double[S2.size()];
        
        //System.out.println("ABOUT TO ENTER THE MAIN LOOP");
        //System.out.println("");
        
        while ((optimNotReached)){ // && (iterCount<50) 
            iterCount++;
            if (iterCount == 10000){
                this.Message = "Warning: Number of iterations reached second threshold (10000) and the process was stopped. The last gradient for free variables was: " + Arrays.toString(gradientxs1) + " " + Arrays.toString(gradientxs2);
                break;
            }
            /**
            System.out.println(":::: ITERATION "+iterCount+"::::");
            System.out.println("   T1: \n" + treePrinter.toString(T1)+"\n \n");
            System.out.println("   T2: \n" + treePrinter.toString(T2)+"\n \n");
            System.out.println("   B1 = " + B1);
            System.out.println("   S1 = " + S1);
            System.out.println("   N1 = " + N1);
            System.out.println("   B2 = " + B2);
            System.out.println("   S2 = " + S2);
            System.out.println("   N2 = " + N2);
            System.out.println("");*/
            
            if (conjugate_initial_counter > 15){
                conjugate_initial_counter = 0;
            }
            //System.out.println("Iteration number " + iterCount);
            double[] gradient1 = new double[n];
            double[] gradient2 = new double[n];
            
            RatioSequence currentRSeq = tempGeode.getRS();//The derivaties will depend on the ratio sequence in the geodesic of the geodesic between current trees T1 and T2. 
            
            //And it also depends on which common edges they have
            Vector<PhyloTreeEdge> currentECEs = tempGeode.geteCommonEdges(); 
            Vector<PhyloTreeEdge> currentFCEs = tempGeode.getfCommonEdges();
            
            //For each ratio we compute the contribution of the expression relating to the ratio in the final geodesic length in the derivative with respect to the edge. 
            Iterator<Ratio> rsIter = currentRSeq.iterator();
            while(rsIter.hasNext()){
                Ratio rat = (Ratio) rsIter.next();
                for (PhyloTreeEdge e : rat.getEEdges()){
                    //int eID = e.getOriginalID();  
                    int eID = T1.getEdges().indexOf(e);
                    if (rat.getELength() == 0){
                        //Subgradient case, we replace the non-defined gradient term with 0
                        gradient1[eID] += rat.getFLength()/Math.sqrt((double) rat.getEEdges().size());
                    } else {
                        gradient1[eID] += e.getNorm()*(1 + (rat.getFLength()/rat.getELength()));
                    }       
                }
                for (PhyloTreeEdge e : rat.getFEdges()){
                    int eID = T2.getEdges().indexOf(e);
                    if (rat.getFLength() == 0){
                        //Subgradient case, we replace the non-defined gradient term with 0
                        gradient2[eID] += rat.getELength()/Math.sqrt((double) rat.getFEdges().size());
                    } else {
                        gradient2[eID] += e.getNorm()*(1 + (rat.getELength()/rat.getFLength()));
                    }   
                }
            }
            
            //For each common edge, we compute the contribution to the derivatives in the gradient.
            
            for(PhyloTreeEdge e : currentECEs){
                int eID = T1.getEdges().indexOf(e);
                if (eID == -1){
                    continue;
                }
                EdgeAttribute T2EAtt = T2.getAttribOfSplit(e.asSplit());
                if (T2EAtt == null){
                    Bipartition eClone = e.getOriginalEdge().clone();
                    eClone.complement(OE2.getCompleteLeafSet().size());
                    T2EAtt = T2.getAttribOfSplit(eClone);
                }
                
                gradient1[eID] += (e.getNorm() - T2EAtt.norm());
            } 
            for(PhyloTreeEdge e : currentFCEs){
                int eID = T2.getEdges().indexOf(e);
                if (eID == -1){
                    continue;
                }
                EdgeAttribute T1EAtt = T1.getAttribOfSplit(e.asSplit());
                if (T1EAtt == null){
                    Bipartition eClone = e.getOriginalEdge().clone();
                    eClone.complement(OE1.getCompleteLeafSet().size());
                    T1EAtt = T1.getAttribOfSplit(eClone);
                }
                gradient2[eID] += (e.getNorm() - T1EAtt.norm());
            } 
            
            
            //Using the gradients for each "variable" (the values of the edges for each current tree) we compute the gradients of the free variables in the reduced gradient method. But first, we need to save the previous values if we are not in the first iteration of a re-initialization of the conjugate gradient method. 
            
            double[] gradientxs1Prev = gradientxs1.clone();
            double[] gradientxs2Prev = gradientxs2.clone();
            
            double akDenom = 0;
            
            if (conjugate_initial_counter > 0){
                for (int i = 0; i < gradientxs1Prev.length; i++){
                    akDenom += gradientxs1Prev[i]*gradientxs1Prev[i];
                }
                for (int i = 0; i < gradientxs2Prev.length; i++){
                    akDenom += gradientxs2Prev[i]*gradientxs2Prev[i];
                }
            }
            
            boolean gradient_small = true; // as we compute the new gradient, we assess if the size is big enough to justify another loop or we have arrive to an stationary point. 
            
            gradientxs1 = new double[S1.size()];
            gradientxs2 = new double[S2.size()];
        
            
            if (iterCount == 2000){
                this.Message = "Warning: Number of iterations reached first threshold (2000) and the tolerance for the gradient was increased to 0.001";
                TolLimit = 0.001;
            }
            
            for (int i = 0; i < S1.size(); i++){
                gradientxs1[i] = gradient1[S1.get(i)] - gradient1[B1.get(OE1.getBackMap(S1.get(i)))];
                if((gradientxs1[i] < -TolLimit) || (gradientxs1[i] > TolLimit)){
                    gradient_small = false;
                }
            }
            
            for (int i = 0; i < S2.size(); i++){
                gradientxs2[i] = gradient2[S2.get(i)] - gradient2[B2.get(OE2.getBackMap(S2.get(i)))];
                if((gradientxs2[i] < -TolLimit) || (gradientxs2[i] > TolLimit)){
                    gradient_small = false;
                }
            }
            
            //System.out.println("The gradient for xs1 is: " + Arrays.toString(gradientxs1));
            //System.out.println("The gradient for xs2 is: " + Arrays.toString(gradientxs2));
            
            //We use continue; in case we have arrived to an stationary point in the current face being considered.
            
            if(gradient_small){
                //If the gradient is small, we have arrived to an semi-stationary point. We will check if it holds the condition to be the optimum or we need to shuffle things around to find the potential one. 
                //System.out.println("   It entered the gradient small if...");
                Vector<Integer> promisingEN1 = new Vector<Integer>();
                Vector<Integer> promisingEN2 = new Vector<Integer>();
                
                optimNotReached = false; //Assume at first that the current semi-stationary point is in fact the optimum. 
                
                for (int i = 0; i < N1.size(); i++){
                    if ((gradient1[N1.get(i)] - gradient1[B1.get(OE1.getBackMap(N1.get(i)))]) < 0){
                        promisingEN1.add(N1.get(i));
                    }
                }
                
                for (int i = 0; i < N2.size(); i++){
                    if ((gradient2[N2.get(i)] - gradient2[B2.get(OE2.getBackMap(N2.get(i)))]) < 0){
                        promisingEN2.add(N2.get(i));
                    }
                }
                
                if ((promisingEN1.size()>0) || (promisingEN2.size()>0)){
                    //System.out.println("But promising N1 is: " + promisingEN1);
                    //System.out.println("But promising N2 is: " + promisingEN2);
                    N1.removeAll(promisingEN1);
                    S1.addAll(promisingEN1);
                    
                    N2.removeAll(promisingEN2);
                    S2.addAll(promisingEN2);
                    
                    conjugate_initial_counter = 0;
                    optimNotReached = true;
                }
                
                continue;//We go back to the beginning of the loop. 
            }
            
            //We know need to determine the best direction of change depending on whether we are in the first iteration of a re=initialization of the conjutage gradient method or not
            
            if (conjugate_initial_counter == 0){
                dDirectionxs1 = new double[S1.size()];
                dDirectionxs2 = new double[S2.size()];
                for (int i = 0; i < dDirectionxs1.length; i++){
                    dDirectionxs1[i] = -gradientxs1[i];
                }
                for (int i = 0; i < dDirectionxs2.length; i++){
                    dDirectionxs2[i] = -gradientxs2[i];
                }
            } else {
                double akNum = 0;
                for (int i = 0; i < gradientxs1.length; i++){
                    akNum += gradientxs1[i]*(gradientxs1[i] - gradientxs1Prev[i]);
                }
                for (int i = 0; i < gradientxs2.length; i++){
                    akNum += gradientxs2[i]*(gradientxs2[i] - gradientxs2Prev[i]);
                }
                
                double ak = akNum/akDenom;
                
                for (int i = 0; i < gradientxs1.length; i++){
                    dDirectionxs1[i] = ak*dDirectionxs1[i] - gradientxs1[i];
                }
                for (int i = 0; i < gradientxs2.length; i++){
                    dDirectionxs2[i] = ak*dDirectionxs2[i] - gradientxs2[i];
                }
            }
            
            //Computing the complete change vector
            
            double[] dDirection1 = new double[n];
            
            for (int i = 0; i < S1.size(); i++){
                dDirection1[S1.get(i)] = dDirectionxs1[i];
                dDirection1[B1.get(OE1.getBackMap(S1.get(i)))] += -dDirectionxs1[i];
            }
            
            double[] dDirection2 = new double[n];
            
            for (int i = 0; i < S2.size(); i++){
                dDirection2[S2.get(i)] = dDirectionxs2[i];
                dDirection2[B2.get(OE2.getBackMap(S2.get(i)))] += -dDirectionxs2[i];
            }
            
            
            //Determining the closed set for tau, in order to mantain all edges with positive size. 
            double tau_max = 0;
            double tau_min = 0;
            boolean tauNeedsChange = true;
            
            Vector<PhyloTreeEdge> EdgesT1 = T1.getEdges();
            Vector<PhyloTreeEdge> EdgesT2 = T2.getEdges();
            
            Vector<Integer> potentialN1 = new Vector<Integer>();
            Vector<Integer> potentialN2 = new Vector<Integer>();
            
            for (int i = 0; i < EdgesT1.size(); i++){
                if (dDirection1[i] < 0){
                    if (tauNeedsChange || (-EdgesT1.get(i).getNorm()/dDirection1[i] < tau_max)){
                        tau_max = -EdgesT1.get(i).getNorm()/dDirection1[i];
                        
                        if (!N1.contains(i)){
                            potentialN1.clear();
                            potentialN1.add(i);
                        } else {
                            System.out.println("An element on N1 sneaked in (situation 1): "+ i);
                        }
                        tauNeedsChange = false;
                    } else if (-EdgesT1.get(i).getNorm()/dDirection1[i] == tau_max){
                        if (!N1.contains(i)){
                            potentialN1.add(i);
                        }else {
                            System.out.println("An element on N1 sneaked in (situation 2): "+ i);
                        }  
                    }
                }
            }
            
            for (int i = 0; i < EdgesT2.size(); i++){
                if (dDirection2[i] < 0){
                    if (tauNeedsChange || (-EdgesT2.get(i).getNorm()/dDirection2[i] < tau_max)){
                        tau_max = -EdgesT2.get(i).getNorm()/dDirection2[i];
                        
                        if (!N2.contains(i)){
                            potentialN1.clear();
                            potentialN2.clear();
                            potentialN2.add(i);
                        } else{
                            System.out.println("An element on N2 sneaked in (situation 1): "+ i);
                        }
                        tauNeedsChange = false;
                    } else if (-EdgesT2.get(i).getNorm()/dDirection2[i] == tau_max){
                        if (!N2.contains(i)){
                            potentialN2.add(i);
                        } else{
                            System.out.println("An element on N2 sneaked in (situation 2): "+ i);
                        }
                    }
                }
            }
            
            
            // We will look for the tau that minimizes f(x + tau* dDirection) between tau_min and tau_max.
            //We will first check if the minimum is the actual tau_max
            
            //Defining new values of the trees to compute geodesic and find the derivative; 
            Vector<PhyloTreeEdge> conjEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
            Vector<PhyloTreeEdge> conjEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
            
            
            //Computing the new values of the interior edges of the trees by moving in the direction of change
            for (int i = 0; i < EdgesT1.size(); i++){
                double[] tempVecEA = {EdgesT1.get(i).getNorm() + (tau_max-0.0000000000000001)*dDirection1[i]};
                EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                conjEdgesT1.get(i).setAttribute(tempEA);
            }
            
            for (int i = 0; i < EdgesT2.size(); i++){
                double[] tempVecEA = {EdgesT2.get(i).getNorm() + (tau_max-0.0000000000000001)*dDirection2[i]};
                EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                conjEdgesT2.get(i).setAttribute(tempEA);
            }
            
            PhyloTree conjT1 = new PhyloTree(conjEdgesT1, T1.getLeaf2NumMap(), T1.getLeafEdgeAttribs(), false);
            PhyloTree conjT2 = new PhyloTree(conjEdgesT2, T2.getLeaf2NumMap(), T2.getLeafEdgeAttribs(), false);
            
            //Computing geodesic in between these trees. 
            Geodesic conjGeode = parPolyMain.getGeodesic(conjT1, conjT2);
            
            RatioSequence conjRSeq = conjGeode.getRS();//The derivative will depend on the ratio sequence 
            
            //And on which common edges they have
            Vector<PhyloTreeEdge> conjCEs = conjGeode.getCommonEdges(); 
            
            
            double derivTau = 0;//Where the final derivative for tau will be saved
            
            //For each ratio we compute the contribution of the expression relating to the ratio in derivative with respect to tau_max
            Iterator<Ratio> conjRSIter = conjRSeq.iterator();
            while(conjRSIter.hasNext()){
                Ratio rat = (Ratio) conjRSIter.next();
                
                //Values that will contribute to the derivative of the ratio
                
                if (rat.getELength() > 0){
                    double ENum = 0;
                    for (PhyloTreeEdge e : rat.getEEdges()){
                        int eID = conjT1.getEdges().indexOf(e);
                        ENum += dDirection1[eID]*e.getNorm();
                    }
                    derivTau += ENum*(1 + (rat.getFLength()/rat.getELength()));
                }
                
                if (rat.getFLength() > 0){
                    double FNum = 0;
                    for (PhyloTreeEdge e : rat.getFEdges()){
                        int eID = conjT2.getEdges().indexOf(e);
                        FNum += dDirection2[eID]*e.getNorm();
                    }
                    derivTau += FNum*(1 + (rat.getELength()/rat.getFLength()));
                }
            }
            
            //For each common edge, we compute the contribution to the derivatives in the gradient. 
        
            for(PhyloTreeEdge e : conjCEs){
                int eID1 = conjT1.getSplits().indexOf(e.asSplit());
                int eID2 = conjT2.getSplits().indexOf(e.asSplit());
                
            
                derivTau += (dDirection1[eID1] - dDirection2[eID2])*(conjT1.getEdge(eID1).getAttribute().get(0) - conjT2.getEdge(eID2).getAttribute().get(0)); //The edge attribute in this case is the value in Tree 1 minus the value in Tree 2. 
            } 
            
            double tau = 0;
            
            if (derivTau <= 0){//In this case the minimum is reached right at the tau_max limit and the search is over.
                //System.out.println("   So it hitted a face");
                tau = tau_max;
                boolean ChangeInIndexMade = false;
                
                //We have hitted a boundary face, so we need to reclasify some variable to N1 or N2. 
                
                if (potentialN1.size() > 0){
                    //We want to give priority to indexes that have not been non-basic variables yet to try and avoid cycling. 
                    int[] IndexListOrdered = new int[potentialN1.size()];
                    int leftInd = 0;
                    int rightInd = potentialN1.size() - 1;
                    boolean AllN1Already = true;
                    for (int i = 0; i < potentialN1.size(); i++){
                        if (alreadyN1.contains(potentialN1.get(i))){
                            IndexListOrdered[rightInd] = i;
                            rightInd--;
                        } else {
                            IndexListOrdered[leftInd] = i;
                            leftInd++;
                            AllN1Already = false;
                        }
                    }
                    if (AllN1Already){
                        Collections.shuffle(potentialN1);
                    }
                    for (int i : IndexListOrdered){
                        if (S1.contains(potentialN1.get(i))){
                            N1.add(potentialN1.get(i));
                            alreadyN1.add(potentialN1.get(i));
                            S1.remove(Integer.valueOf(potentialN1.get(i)));
                            ChangeInIndexMade = true;
                            break;
                        } else if (B1.contains(potentialN1.get(i))){
                            int rowIndexTemp = B1.indexOf(potentialN1.get(i));
                            int newB1element = -1; 
                            for (int j : OE1.getMapList().get(rowIndexTemp)){
                                if (S1.contains(j)){
                                    newB1element = j;
                                    break;
                                }
                            }
                            if (newB1element != -1){
                                B1.set(rowIndexTemp, newB1element);
                                S1.remove(Integer.valueOf(newB1element));
                                N1.add(potentialN1.get(i));
                                alreadyN1.add(potentialN1.get(i));
                                ChangeInIndexMade = true;
                                break;
                            }
                        
                        }
                    }
                } else if (potentialN2.size() > 0){
                    //We want to give priority to indexes that have not been non-basic variables yet to try and avoid cycling. 
                    int[] IndexListOrdered = new int[potentialN2.size()];
                    int leftInd = 0;
                    int rightInd = potentialN2.size() - 1;
                    boolean AllN2Already = true;
                    for (int i = 0; i < potentialN2.size(); i++){
                        if (alreadyN2.contains(potentialN2.get(i))){
                            IndexListOrdered[rightInd] = i;
                            rightInd--;
                        } else {
                            IndexListOrdered[leftInd] = i;
                            leftInd++;
                            AllN2Already = false;
                        }
                    }
                    if (AllN2Already){
                        Collections.shuffle(potentialN2);
                    }
                    for (int i : IndexListOrdered){
                        if (S2.contains(potentialN2.get(i))){
                            N2.add(potentialN2.get(i));
                            alreadyN2.add(potentialN2.get(i));
                            S2.remove(Integer.valueOf(potentialN2.get(i)));
                            ChangeInIndexMade = true;
                            break;
                        } else if (B2.contains(potentialN2.get(i))){
                            int rowIndexTemp = B2.indexOf(potentialN2.get(i));
                            int newB2element = -1; 
                            for (int j : OE2.getMapList().get(rowIndexTemp)){
                                if (S2.contains(j)){
                                    newB2element = j;
                                    break;
                                }
                            }
                            if (newB2element != -1){
                                B2.set(rowIndexTemp, newB2element);
                                S2.remove(Integer.valueOf(newB2element));
                                N2.add(potentialN2.get(i));
                                alreadyN2.add(potentialN2.get(i));
                                ChangeInIndexMade = true;
                                break;
                            }
                        }
                    }
                }
                if (!ChangeInIndexMade){
                    System.out.println("ERROR: Although a variable should be reclassified as non-basic, it did not happen.");
                    break;
                }
                
                conjugate_initial_counter = 0; // We are re-initializing the conjugate gradient method in a new face;
                
                //Defining the new trees to go back to the main while loop: 
                
                Vector<PhyloTreeEdge> newEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
                Vector<PhyloTreeEdge> newEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
                double[] newEdgesValuesT1 = new double[n];
                double[] newEdgesValuesT2 = new double[n];
                
                for (int i = 0; i < B1.size(); i++){
                    newEdgesValuesT1[B1.get(i)] = OE1.getFixedLengths(i);
                }
                for (int i = 0; i < S1.size(); i++){
                    newEdgesValuesT1[S1.get(i)] = EdgesT1.get(S1.get(i)).getNorm() + tau*dDirection1[S1.get(i)];
                    newEdgesValuesT1[B1.get(OE1.getBackMap(S1.get(i)))] -= newEdgesValuesT1[S1.get(i)];
                }
                
                for (int i = 0; i < B2.size(); i++){
                    newEdgesValuesT2[B2.get(i)] = OE2.getFixedLengths(i);
                }
                for (int i = 0; i < S2.size(); i++){
                    newEdgesValuesT2[S2.get(i)] = EdgesT2.get(S2.get(i)).getNorm() + tau*dDirection2[S2.get(i)];
                    newEdgesValuesT2[B2.get(OE2.getBackMap(S2.get(i)))] -= newEdgesValuesT2[S2.get(i)];
                }
                
                
            
                //Computing the new values of the interior edges of the trees 
                for (int i = 0; i < newEdgesT1.size(); i++){
                    double[] tempVecEA = {newEdgesValuesT1[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newEdgesT1.get(i).setAttribute(tempEA);
                }
            
                for (int i = 0; i < newEdgesT2.size(); i++){
                    double[] tempVecEA = {newEdgesValuesT2[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newEdgesT2.get(i).setAttribute(tempEA);
                }
            
                T1 = new PhyloTree(newEdgesT1, T1.getLeaf2NumMap(), T1.getLeafEdgeAttribs(), false);
                T2 = new PhyloTree(newEdgesT2, T2.getLeaf2NumMap(), T2.getLeafEdgeAttribs(), false);
                
            
                tempGeode = parPolyMain.getGeodesic(T1, T2);
                
            } else {//We still need to find the optimum tau for this case. 
                //System.out.println("We are in the second option");
                int counterWhile = 0;
                tau = 0.1;
                if (tau > tau_max/2){
                    //System.out.println("tau initial changed");
                    tau = tau_max/2;
                }
                
                //System.out.println("tau max: "+ tau_max);
                //System.out.println("tau min: "+ tau_min);
                //System.out.println("tau value before while: "+ tau);
                while(((derivTau < -0.0000000000000001) || (derivTau > 0.0000000000000001)) && (((tau_max - tau_min) > 0.0000000001))){ //&&(counterWhile < 50)
                    counterWhile++;
                    //System.out.println("   Inside the tau while loop "+counterWhile);
                    tau = (tau_max + tau_min)/2;
                    
                    conjEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
                    conjEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
                    
                    //Computing the new values of the interior edges of the trees by moving in the direction of change
                    for (int i = 0; i < EdgesT1.size(); i++){
                        double[] tempVecEA = {EdgesT1.get(i).getNorm() + tau*dDirection1[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        conjEdgesT1.get(i).setAttribute(tempEA);
                    }
            
                    for (int i = 0; i < EdgesT2.size(); i++){
                        double[] tempVecEA = {EdgesT2.get(i).getNorm() + tau*dDirection2[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        conjEdgesT2.get(i).setAttribute(tempEA);
                    }
                    
                    conjT1 = new PhyloTree(conjEdgesT1, T1.getLeaf2NumMap(), T1.getLeafEdgeAttribs(), false);
                    conjT2 = new PhyloTree(conjEdgesT2, T2.getLeaf2NumMap(), T2.getLeafEdgeAttribs(), false);
                    
                    //Computing geodesic in between these trees. 
                    conjGeode = parPolyMain.getGeodesic(conjT1, conjT2);
            
                    conjRSeq = conjGeode.getRS();//The derivative will depend on the ratio sequence 
            
                    //And on which common edges they have
                    conjCEs = conjGeode.getCommonEdges(); 
            
                    derivTau = 0;//Where the final derivative for tau will be saved
                    
                    //For each ratio we compute the contribution of the expression relating to the ratio in derivative with respect to tau_max
                    conjRSIter = conjRSeq.iterator();
                    while(conjRSIter.hasNext()){
                        Ratio rat = (Ratio) conjRSIter.next();
                        
                        //Values that will contribute to the derivative of the ratio
                
                        if (rat.getELength() > 0){
                            double ENum = 0;
                            for (PhyloTreeEdge e : rat.getEEdges()){
                                int eID = conjT1.getEdges().indexOf(e);
                                ENum += dDirection1[eID]*e.getNorm();
                            }
                            derivTau += ENum*(1 + (rat.getFLength()/rat.getELength()));
                        }
                
                        if (rat.getFLength() > 0){
                            double FNum = 0;
                            for (PhyloTreeEdge e : rat.getFEdges()){
                                int eID = conjT2.getEdges().indexOf(e);
                                FNum += dDirection2[eID]*e.getNorm();
                            }
                            derivTau += FNum*(1 + (rat.getELength()/rat.getFLength()));
                        }
                    }
            
                    //For each common edge, we compute the contribution to the derivatives in the gradient. 
                    for(PhyloTreeEdge e : conjCEs){
                        int eID1 = conjT1.getSplits().indexOf(e.asSplit()); 
                        int eID2 = conjT2.getSplits().indexOf(e.asSplit());
            
                        derivTau += (dDirection1[eID1] - dDirection2[eID2])*(conjT1.getEdge(eID1).getAttribute().get(0) - conjT2.getEdge(eID2).getAttribute().get(0)); //The edge attribute in this case is the value in Tree 1 minus the value in Tree 2. 
                    }
                    
                    if (derivTau <= 0){// This would mean the minimum is between tau and tau_max
                        tau_min = tau;
                    } else {
                        tau_max = tau;
                    }
                }
                
                conjugate_initial_counter++; // Keeping count on how many loops we have done in this face. 
                
                //Defining the new trees to go back to the main while loop: 
                
                Vector<PhyloTreeEdge> newEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
                Vector<PhyloTreeEdge> newEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
                
                //System.out.println("   tau after while: " + tau);
                //System.out.println("   derivTau after while: " + derivTau);
                double[] newEdgesValuesT1 = new double[n];
                double[] newEdgesValuesT2 = new double[n];
                
                for (int i = 0; i < B1.size(); i++){
                    newEdgesValuesT1[B1.get(i)] = OE1.getFixedLengths(i);
                }
                for (int i = 0; i < S1.size(); i++){
                    newEdgesValuesT1[S1.get(i)] = EdgesT1.get(S1.get(i)).getNorm() + tau*dDirection1[S1.get(i)];
                    
                    newEdgesValuesT1[B1.get(OE1.getBackMap(S1.get(i)))] -= newEdgesValuesT1[S1.get(i)];
                }
                
                for (int i = 0; i < B2.size(); i++){
                    newEdgesValuesT2[B2.get(i)] = OE2.getFixedLengths(i);
                }
                for (int i = 0; i < S2.size(); i++){
                    newEdgesValuesT2[S2.get(i)] = EdgesT2.get(S2.get(i)).getNorm() + tau*dDirection2[S2.get(i)];
                    newEdgesValuesT2[B2.get(OE2.getBackMap(S2.get(i)))] -= newEdgesValuesT2[S2.get(i)];
                }
                
                
            
                //Computing the new values of the interior edges of the trees 
                for (int i = 0; i < newEdgesT1.size(); i++){
                    double[] tempVecEA = {newEdgesValuesT1[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newEdgesT1.get(i).setAttribute(tempEA);
                }
            
                for (int i = 0; i < newEdgesT2.size(); i++){
                    double[] tempVecEA = {newEdgesValuesT2[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newEdgesT2.get(i).setAttribute(tempEA);
                }
            
                T1 = new PhyloTree(newEdgesT1, T1.getLeaf2NumMap(), T1.getLeafEdgeAttribs(), false);
                T2 = new PhyloTree(newEdgesT2, T2.getLeaf2NumMap(), T2.getLeafEdgeAttribs(), false);
                
                tempGeode = parPolyMain.getGeodesic(T1, T2);
                
            }
            
        }
        
        //Getting the final values after the gradient descent has been performed. 
        Tree1 = T1;
        Tree2 = T2;
        Geodesic FinalGeode = tempGeode;
        Distance = FinalGeode.getDist();
        IterCount = iterCount;
        
    }// end of Constructor1
    
    
    public OrthExtDistance(OrthExt OE1, OrthExt OE2){
        O1ID = OE1.getOID();
        O2ID = OE2.getOID();
        Constructor1(OE1, OE2);
    }
    
    public OrthExtDistance(orthantExtPair orthEP){
        O1ID = orthEP.getOrthE1().getOID();
        O2ID = orthEP.getOrthE2().getOID();
        Constructor1(orthEP.getOrthE1(), orthEP.getOrthE2());
    }
    
    //The following function was used in Constructor2, but now it has been added directly inside the Constructor2 
    
    /*private PhyloTree[] NewMutualTrees(OrthExt OE1, OrthExt OE2){
        PhyloTree T1 = new PhyloTree(OE1.getStartTree());
        PhyloTree T2 = new PhyloTree(OE2.getStartTree());
        
        EdgeAttribute[] EAT1 = T1.getCopyLeafEdgeAttribs();
        EdgeAttribute[] EAT2 = T2.getCopyLeafEdgeAttribs();
        
        for (int i = 0; i < OE1.getCompleteLeafSet().size(); i++){
            if ((OE1.getCompLeaves2orgLeaves(i) == -1) && (OE2.getCompLeaves2orgLeaves(i) != -1)){
                EAT1[i] = EAT2[i].clone();
            }
            if ((OE2.getCompLeaves2orgLeaves(i) == -1) && (OE1.getCompLeaves2orgLeaves(i) != -1)){
                EAT2[i] = EAT1[i].clone();
            }
        }
        
        //T1.setLeafEdgeAttribs(EAT1);
        //T2.setLeafEdgeAttribs(EAT2);
        
        Vector<PhyloTreeEdge> T1Edges = polyAlg.Tools.myVectorClonePhyloTreeEdge(T1.getEdges());
        Vector<PhyloTreeEdge> T2Edges = polyAlg.Tools.myVectorClonePhyloTreeEdge(T2.getEdges());
        
        int orgNumLeaves1 = OE1.getOriginalLeaves().cardinality();
        int orgNumLeaves2 = OE2.getOriginalLeaves().cardinality();
        
        for (int i = orgNumLeaves1; i < OE1.getBackMap().length; i++){
            if(OE1.getBackMap(i) == -1){
                for (int j = 0; j < T2Edges.size(); j++){
                    PhyloTreeEdge e = T2Edges.get(j);
                    if (e.sameBipartition(OE1.getOrthantAxis(i-orgNumLeaves1))){
                        cur1Axis2Edges[i-orgNumLeaves1] = T1Edges.size();
                        cur1Edges2Axis.add(Integer.valueOf(i-orgNumLeaves1));
                        ET2toET1.put(Integer.valueOf(j), Integer.valueOf(T1Edges.size()));
                        PhyloTreeEdge eCl = e.clone();
                        eCl.setOriginalID(T1Edges.size());
                        T1Edges.add(eCl);
                    }
                }
            }
        }
        
        for (int i = orgNumLeaves2; i < OE2.getBackMap().length; i++){
            if(OE2.getBackMap(i) == -1){
                for (int j = 0; j < T1Edges.size(); j++){
                    PhyloTreeEdge e = T1Edges.get(j);
                    if (e.sameBipartition(OE2.getOrthantAxis(i-orgNumLeaves2))){
                        cur2Axis2Edges[i-orgNumLeaves2] = T2Edges.size();
                        cur2Edges2Axis.add(Integer.valueOf(i-orgNumLeaves2));
                        ET1toET2.put(Integer.valueOf(j), Integer.valueOf(T2Edges.size()));
                        PhyloTreeEdge eCl = e.clone();
                        eCl.setOriginalID(T2Edges.size());
                        T2Edges.add(eCl);
                    }
                }
            }
        }*/
        
        /*for(int i = 0; i < OE1.getOrthantAxis().size(); i++){
            if(cur1Axis2Edges[i] == -1){
                double[] tempVecEA = {0};
                EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                PhyloTreeEdge tempE = new PhyloTreeEdge(OE1.getOrthantAxis(i), tempEA, T1Edges.size());
                T1Edges.add(tempE);
            }
        }
        
        for(int i = 0; i < OE2.getOrthantAxis().size(); i++){
            if(cur2Axis2Edges[i] == -1){
                double[] tempVecEA = {0};
                EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                PhyloTreeEdge tempE = new PhyloTreeEdge(OE2.getOrthantAxis(i), tempEA, T2Edges.size());
                T2Edges.add(tempE);
            }
        }
        
        //T1.setEdges(T1Edges);
        //T2.setEdges(T2Edges);
        
        PhyloTree[] NewArray = new PhyloTree[2];
        NewArray[0] = new PhyloTree(T1Edges, T1.getLeaf2NumMap(), EAT1, false);
        NewArray[1] = new PhyloTree(T2Edges, T2.getLeaf2NumMap(), EAT2, false);
        
        return(NewArray);
    }*/
    
    //Constructor 2
    private void Constructor2(OrthExt OE1, OrthExt OE2){  
        //Scanner scan = new Scanner(System.in);
        //String pause = scan.next();
        double TolLimit = 0.00000001;
        PhyloNicePrinter treePrinter = new PhyloNicePrinter();
        //We start by the starting trees in each orthant extension.
        
        int[] cur1Axis2Edges = OE1.getCloneAxis2Edges();
        int[] cur2Axis2Edges = OE2.getCloneAxis2Edges();
        
        Vector<Integer> cur1Edges2Axis = OE1.getCloneEdges2Axis();
        Vector<Integer> cur2Edges2Axis = OE2.getCloneEdges2Axis();
        
        Map<Integer, Integer> ET1toET2 = new HashMap<>(); //These HashMaps serve to point consequential edges towards 
        Map<Integer, Integer> ET2toET1 = new HashMap<>(); // the common inconsequential in the other tree.
        
        //We will construct new starting trees that are mutually restricted
        EdgeAttribute[] EAT1 = OE1.getStartTree().getCopyLeafEdgeAttribs();
        EdgeAttribute[] EAT2 = OE2.getStartTree().getCopyLeafEdgeAttribs();
        
        for (int i = 0; i < OE1.getCompleteLeafSet().size(); i++){
            if ((OE1.getCompLeaves2orgLeaves(i) == -1) && (OE2.getCompLeaves2orgLeaves(i) != -1)){
                EAT1[i] = EAT2[i].clone();
            }
            if ((OE2.getCompLeaves2orgLeaves(i) == -1) && (OE1.getCompLeaves2orgLeaves(i) != -1)){
                EAT2[i] = EAT1[i].clone();
            }
        }
        
        Vector<PhyloTreeEdge> T1Edges = polyAlg.Tools.myVectorClonePhyloTreeEdge(OE1.getStartTree().getEdges());
        Vector<PhyloTreeEdge> T2Edges = polyAlg.Tools.myVectorClonePhyloTreeEdge(OE2.getStartTree().getEdges());
        
        int orgNumLeaves1 = OE1.getOriginalLeaves().cardinality();
        int orgNumLeaves2 = OE2.getOriginalLeaves().cardinality();
        
        for (int i = orgNumLeaves1; i < OE1.getBackMap().length; i++){
            if(OE1.getBackMap(i) == -1){
                for (int j = 0; j < T2Edges.size(); j++){
                    PhyloTreeEdge e = T2Edges.get(j);
                    if (e.sameBipartition(OE1.getOrthantAxis(i-orgNumLeaves1))){
                        cur1Axis2Edges[i-orgNumLeaves1] = T1Edges.size();
                        cur1Edges2Axis.add(Integer.valueOf(i-orgNumLeaves1));
                        ET2toET1.put(Integer.valueOf(j), Integer.valueOf(T1Edges.size()));
                        PhyloTreeEdge eCl = e.clone();
                        eCl.setOriginalID(T1Edges.size());
                        T1Edges.add(eCl);
                    }
                }
            }
        }
        
        for (int i = orgNumLeaves2; i < OE2.getBackMap().length; i++){
            if(OE2.getBackMap(i) == -1){
                for (int j = 0; j < T1Edges.size(); j++){
                    PhyloTreeEdge e = T1Edges.get(j);
                    if (e.sameBipartition(OE2.getOrthantAxis(i-orgNumLeaves2))){
                        cur2Axis2Edges[i-orgNumLeaves2] = T2Edges.size();
                        cur2Edges2Axis.add(Integer.valueOf(i-orgNumLeaves2));
                        ET1toET2.put(Integer.valueOf(j), Integer.valueOf(T2Edges.size()));
                        PhyloTreeEdge eCl = e.clone();
                        eCl.setOriginalID(T2Edges.size());
                        T2Edges.add(eCl);
                    }
                }
            }
        }
        
        PhyloTree T1 = new PhyloTree(T1Edges, OE1.getCompleteLeafSet(), EAT1, false);
        PhyloTree T2 = new PhyloTree(T2Edges, OE2.getCompleteLeafSet(), EAT2, false);
        
        
        /////////////////////////////////////////
        
        
        //Find the the geodesic in between these trees. 
        Geodesic tempGeode = parPolyMain.getGeodesic(T1, T2);
        
        //System.out.println("");
        //System.out.println("PROCESSING THE DISTANCE ALGORITHM... Starting at distance "+ tempGeode.getDist());
        //System.out.println("With T1: " + treePrinter.toString(T1) + "\n");
        //System.out.println("With T2: " + treePrinter.toString(T2) + "\n");
        
       // System.out.println("The number of ratios is: " + tempGeode.getRS().size());
        
        //Some useful counters
        int k1 = OE1.getDim();//Dimension of the first Orthant Extension Space
        int k2 = OE2.getDim();//Dimension of the second Orthant Extension Space
        int m1 = OE1.getFixedLengths().length;
        int m2 = OE2.getFixedLengths().length;
        
        int[] IVI1 = new int[OE1.getBackMap().length];//Help us to keep track of the indices of the actual variables.
        int[] IVI2 = new int[OE2.getBackMap().length];
        
        int idCount1 = 0;
        for (int i = 0; i < IVI1.length; i++){
            if(OE1.getBackMap(i) != -1){
                IVI1[i] = idCount1;
                idCount1++;
            } else {
                IVI1[i] = -1;
            }
        }
        
        int idCount2 = 0;
        for (int i = 0; i < IVI2.length; i++){
            if(OE2.getBackMap(i) != -1){
                IVI2[i] = idCount2;
                idCount2++;
            } else {
                IVI2[i] = -1;
            }
        }
        
        int ol1 = OE1.getOrgLeaves2compLeaves().length; //Variable where the number of original leaves for the first Extension is saved.
        int ol2 = OE2.getOrgLeaves2compLeaves().length; //number of original leaves for the second Extension
        
        int n = OE1.getOrthantAxis().size(); //The number of interior edges in binary trees with the complete leaf set and should coincide for both extension spaces. TO DO: verify it does coincide?  
        
        //TO DO: add code to verify both orthant extensions are in fact inside the same BHV tree space. For now, I just assume every user will be careful about this. 
        
        
        //The following while will perform reduced gradient method algorithm, with a conjugate gradient method in each classification of variables. In each iteration the gradient of the "active variables" (those clasified into S1 and S2) function from the current trees is computed, the optimal descent direction is selected following the conjugate gradient method, and the minimum in that direction is computed. If we hit a boundary, we reclasify variables in order to increment those forced to be zero. We continue until finding a semi-stationary point, and corroborate this is the optimum or add new non-basic variables otherwise.
        
        //Initializing the indexes sets B, S and N, with some extra structures to easy change. 
        //THIS COULD POTENTIALLY BE A PART OF OrthExt class TO AVOID IT BEING COMPUTED EVERY TIME A DISTANCE IS COMPUTED
        
        Vector<Integer> B1 = new Vector<Integer>();
        Vector<Integer> B2 = new Vector<Integer>();
        
        Vector<Integer> S1 = new Vector<Integer>();
        Vector<Integer> S2 = new Vector<Integer>();
        
        Vector<Integer> N1 = new Vector<Integer>();
        Vector<Integer> N2 = new Vector<Integer>();
        
        //We will keep a vector of indexes that have already been non-basic variables, to give priority to new potential non-basic variables with possible, trying to prevent cycling. 
        Vector<Integer> alreadyN1 = new Vector<Integer>();
        Vector<Integer> alreadyN2 = new Vector<Integer>();
        
        for (int i = 0; i < m1; i++){//For each row in the map matrix
            Vector<Integer> tempVect = OE1.getMapList().get(i); //Get the edges that merge into the final edge in the original tree
            B1.add(tempVect.get(0)); //Add the first entry of this list of edges into B1
            S1.addAll(tempVect.subList(1,tempVect.size())); //The rest is added to S1
        }
        
        for (int i = 0; i < m2; i++){//For each row in the map matrix for the second extension
            Vector<Integer> tempVect = OE2.getMapList().get(i); //Get the edges that merge into the final edge in the original tree
            B2.add(tempVect.get(0)); //Add the first entry of this list of edges into B1
            S2.addAll(tempVect.subList(1,tempVect.size())); //The rest is added to S1
        }
        
        //Some values before the iterations start
        
        int iterCount = 0; //Counter of the number of iterations performed.  
        
        boolean optimNotReached = true;//We will stop the loop when the gradient is small enough to guarantee we have reach the minimum. 
        
        int conjugate_initial_counter = 0; //Counter for re-initialization of the conjugate gradient method
        
        //We need to keep track on gradients and change directions
        
        double[] gradientxs1 = new double[S1.size()];
        double[] gradientxs2 = new double[S2.size()];
        
        double[] dDirectionxs1 = new double[S1.size()];
        double[] dDirectionxs2 = new double[S2.size()];
        
        /*System.out.println("ABOUT TO ENTER THE MAIN LOOP");
        System.out.println("IVI1: " + Arrays.toString(IVI1));
        System.out.println("IVI2: " + Arrays.toString(IVI2));
        System.out.println("");
        */
        
        while ((optimNotReached)){ // && (iterCount<2)
            iterCount++;
            
            if (iterCount == 10000){
                this.Message = "Warning: Number of iterations reached second threshold (10000) and the process was stopped. The last gradient for free variables was: " + Arrays.toString(gradientxs1) + " " + Arrays.toString(gradientxs2);
                break;
            }
            
            /*System.out.println(":::: ITERATION " + iterCount + "::::");
            System.out.println("   T1: \n" + treePrinter.toString(T1)+"\n \n");
            System.out.println("   T2: \n" + treePrinter.toString(T2)+"\n \n");
            System.out.println("   B1 = " + B1);
            System.out.println("   S1 = " + S1);
            System.out.println("   N1 = " + N1);
            System.out.println("   B2 = " + B2);
            System.out.println("   S2 = " + S2);
            System.out.println("   N2 = " + N2);
            System.out.println("");
            pause = scan.next();/*/
            
            if (conjugate_initial_counter > 15){
                conjugate_initial_counter = 0;
            }
            //System.out.println("Iteration number " + iterCount);
            double[] gradient1 = new double[m1 + k1];//Changing to the number of non-zero columns
            double[] gradient2 = new double[m2 + k2];
            
            //Supporting structures for the subgradient condition.
            Vector<Double> LooseGradient1 = new Vector<Double>();//Gradients associated with zero length supports
            Vector<Double> LooseGradient2 = new Vector<Double>();
            int[] map2LooseGradient1 = new int[m1 + k1];
            int[] map2LooseGradient2 = new int[m2 + k2];
            int numLooseGradient1 = 0;
            int numLooseGradient2 = 0;
            Map<Integer, Vector<Integer>> map2Gradients1 = new HashMap<Integer, Vector<Integer>>();//The vector that we build for each loose gradient, hopefully unit or lower. 
            Map<Integer, Vector<Integer>> map2Gradients2 = new HashMap<Integer, Vector<Integer>>();
            
   
            RatioSequence currentRSeq = tempGeode.getRS();//The derivaties will depend on the ratio sequence in the geodesic of the geodesic between current trees T1 and T2. 
            
            //System.out.println("The geodesic summary is: " + treePrinter.toString(tempGeode, OE1.getCompleteLeafSet()));
            
            /*System.out.println("   T1: \n" + treePrinter.toString(T1)+"\n \n");
            //System.out.println("   T2: \n" + treePrinter.toString(T2)+"\n \n");
            
            System.out.println("Reviewing which edge 'crosses'  the BG edge: ");
            PhyloTreeEdge edgeTemp = T2.getEdge(3);
            System.out.println(edgeTemp.toString());
            System.out.println(edgeTemp.isCompatibleWith(T1.getSplits(), OE1.getCompleteLeafSet().size()));
            System.out.println("!(disjointFrom(e) || this.contains(e) || e.contains(this))");
            System.out.println("");
            for (int i = 0; i<3; i++){
                System.out.println(edgeTemp.crosses(T1.getEdge(i), OE1.getCompleteLeafSet().size()));
                System.out.println(T1.getEdge(i).toString());
                System.out.println(edgeTemp.disjointFrom(T1.getEdge(i)));
                System.out.println(edgeTemp.contains(T1.getEdge(i)));
                System.out.println(T1.getEdge(i).contains(edgeTemp));
                System.out.println(" ");
            }*/
            
            //And it also depends on which common edges they have
            Vector<PhyloTreeEdge> currentECEs = tempGeode.geteCommonEdges(); 
            Vector<PhyloTreeEdge> currentFCEs = tempGeode.getfCommonEdges();
            
            //For each ratio we compute the contribution of the expression relating to the ratio in the final geodesic length in the derivative with respect to the edge. 
            Iterator<Ratio> rsIter = currentRSeq.iterator();
            while(rsIter.hasNext()){
                Ratio rat = (Ratio) rsIter.next();
                //System.out.println(" One ratio in beginning: " + rat.toStringVerbose(OE1.getCompleteLeafSet()));
                Vector<Integer> VforElength = new Vector<Integer>();
                Vector<Integer> VforFlength = new Vector<Integer>();
                
                if (rat.getELength() == 0){
                    LooseGradient1.add(Double.valueOf(rat.getFLength()));
                }
                if (rat.getFLength() == 0){
                    LooseGradient2.add(Double.valueOf(rat.getELength()));
                }
                for (PhyloTreeEdge e : rat.getEEdges()){
                    //int eID = e.getOriginalID();  
                    int eID = T1.getEdges().indexOf(e);
                    if (rat.getELength() == 0){
                        //Subgradient case, we replace gradient entry for zero instead.
                        gradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] += 0;//rat.getFLength()/Math.sqrt((double) rat.getEEdges().size()); //We add ol1 since in the unrestricted case, all internal edges are pushed to the end, because the original leave edges are all first. 
                        map2LooseGradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] = numLooseGradient1;
                        VforElength.add(Integer.valueOf(IVI1[cur1Edges2Axis.get(eID) + ol1]));
                    } else {
                        gradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] += e.getNorm()*(1 + (rat.getFLength()/rat.getELength()));
                        map2LooseGradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] = -1;
                    }       
                }
                for (PhyloTreeEdge e : rat.getFEdges()){
                    int eID = T2.getEdges().indexOf(e);
                    if (rat.getFLength() == 0){
                        //Subgradient case, we replace gradient entry for zero instead.
                        gradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] += 0;//rat.getELength()/Math.sqrt((double) rat.getFEdges().size());
                        map2LooseGradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] = numLooseGradient2;
                        VforFlength.add(Integer.valueOf(IVI2[cur2Edges2Axis.get(eID) + ol2]));
                    } else {
                        gradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] += e.getNorm()*(1 + (rat.getELength()/rat.getFLength()));
                        map2LooseGradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] = -1;
                    }   
                }
                if (rat.getELength() == 0){
                    map2Gradients1.put(Integer.valueOf(numLooseGradient1), VforElength);
                    numLooseGradient1++;
                }
                
                if (rat.getFLength() == 0){
                    map2Gradients2.put(Integer.valueOf(numLooseGradient2), VforFlength);
                    numLooseGradient2++;
                }
            }
            
            //For each common edge, we compute the contribution to the derivatives in the gradient.
            
            /*System.out.println("Here are the curent ECEs at beginning:");
            for(PhyloTreeEdge e : currentECEs){
                System.out.println("    " + treePrinter.toString(e, OE1.getCompleteLeafSet()));
            }
            
            System.out.println("Here are the curent FCEs at beginning:");
            for(PhyloTreeEdge e : currentFCEs){
                System.out.println("    " +  treePrinter.toString(e, OE1.getCompleteLeafSet()));
            }*/
            
            for(PhyloTreeEdge e : currentECEs){
                int eID = T1.getEdges().indexOf(e);
                if (eID == -1){
                    continue;
                }
                if (IVI1[cur1Edges2Axis.get(eID) + ol1] != -1){
                    EdgeAttribute T2EAtt = T2.getAttribOfSplit(e.asSplit());
                    if (T2EAtt == null){
                        Bipartition eClone = e.getOriginalEdge().clone();
                        eClone.complement(OE2.getCompleteLeafSet().size());
                        T2EAtt = T2.getAttribOfSplit(eClone);
                    }
                    if (T2EAtt == null){ //If T2EAtt is still null, then the "common" edge is actually not present.
                        gradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] += (e.getNorm());
                        map2LooseGradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] = -1;
                    } else {
                        gradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] += (e.getNorm() - T2EAtt.norm());
                        map2LooseGradient1[IVI1[cur1Edges2Axis.get(eID) + ol1]] = -1;
                    }
                }
                
            } 
            for(PhyloTreeEdge e : currentFCEs){
                int eID = T2.getEdges().indexOf(e);
                if (eID == -1){
                    continue;
                }
                if (IVI2[cur2Edges2Axis.get(eID) + ol2] != -1){
                    EdgeAttribute T1EAtt = T1.getAttribOfSplit(e.asSplit());
                    if (T1EAtt == null){
                        Bipartition eClone = e.getOriginalEdge().clone();
                        eClone.complement(OE1.getCompleteLeafSet().size());
                        T1EAtt = T1.getAttribOfSplit(eClone);
                    }
                    if (T1EAtt == null){//If T1EAtt is still null, then the "common" edge is actually not present.
                        gradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] += (e.getNorm());
                        map2LooseGradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] = -1;
                    } else {
                        gradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] += (e.getNorm() - T1EAtt.norm());
                        map2LooseGradient2[IVI2[cur2Edges2Axis.get(eID) + ol2]] = -1;
                    }
                    
                }
            } 
            
            
            //In the unrestricted case, lenghts of external edges to the original leaves are also potential variables, and are treated similarly to common edges.
            for (int i = 0; i < ol1; i++){
                gradient1[IVI1[i]] += (T1.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(i)].get(0) - T2.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(i)].get(0));
                 map2LooseGradient1[IVI1[i]] = -1;
            }
            for (int i = 0; i < ol2; i++){
                gradient2[IVI2[i]] += (T2.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(i)].get(0) - T1.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(i)].get(0));
                map2LooseGradient2[IVI2[i]] = -1;
            }
            
            
            //System.out.println("The gradient1: " + Arrays.toString(gradient1));
            //System.out.println("The gradient2: " + Arrays.toString(gradient2));
            //pause = scan.next();
            
            //System.out.println("Loose gradients1: " + LooseGradient1);
            //System.out.println("Loose gradients2: " + LooseGradient2);
            //System.out.println("map2LG1: " + Arrays.toString(map2LooseGradient1));
            //System.out.println("map2LG2: " + Arrays.toString(map2LooseGradient2));
            //pause = scan.next();
            
            //Using the gradients for each "variable" (the values of the edges for each current tree) we compute the gradients of the free variables in the reduced gradient method. But first, we need to save the previous values if we are not in the first iteration of a re-initialization of the conjugate gradient method. 
            
            double[] gradientxs1Prev = gradientxs1.clone();
            double[] gradientxs2Prev = gradientxs2.clone();
            
            double akDenom = 0;
            
            if (conjugate_initial_counter > 0){
                for (int i = 0; i < gradientxs1Prev.length; i++){
                    akDenom += gradientxs1Prev[i]*gradientxs1Prev[i];
                }
                for (int i = 0; i < gradientxs2Prev.length; i++){
                    akDenom += gradientxs2Prev[i]*gradientxs2Prev[i];
                }
            }
            
            boolean gradient_small = true; // as we compute the new gradient, we assess if the size is big enough to justify another loop or we have arrive to an stationary point. 
            
            gradientxs1 = new double[S1.size()];
            gradientxs2 = new double[S2.size()];
        
            //dDirectionxs1 = new double[S1.size()];
            //dDirectionxs2 = new double[S2.size()];
            
            if (iterCount == 2000){
                this.Message = "Warning: Number of iterations reached first threshold (2000) and the tolerance for the gradient was increased to 0.001";
                TolLimit = 0.001;
            }
            
            for (int i = 0; i < S1.size(); i++){
                gradientxs1[i] = gradient1[IVI1[S1.get(i)]] - gradient1[IVI1[B1.get(OE1.getBackMap(S1.get(i)))]];
                if((gradientxs1[i] < -TolLimit) || (gradientxs1[i] > TolLimit)){
                    gradient_small = false;
                }
            }
            
            for (int i = 0; i < S2.size(); i++){
                gradientxs2[i] = gradient2[IVI2[S2.get(i)]] - gradient2[IVI2[B2.get(OE2.getBackMap(S2.get(i)))]];
                if((gradientxs2[i] < -TolLimit) || (gradientxs2[i] > TolLimit)){
                    gradient_small = false;
                }
            }

            
            //System.out.println("The gradientxs1: " + Arrays.toString(gradientxs1));
            //System.out.println("The gradientxs2: " + Arrays.toString(gradientxs2));
            //pause = scan.next();
            
            //We use continue; in case we have arrived to an stationary point in the current face being considered.
            
            if(gradient_small){
                //If the gradient is small, we have arrived to an semi-stationary point. We will check if it holds the condition to be the optimum or we need to shuffle things around to find the potential one. 
                //System.out.println("   It entered the gradient small if...");
                Vector<Integer> promisingEN1 = new Vector<Integer>();
                Vector<Integer> promisingEN2 = new Vector<Integer>();
                
                optimNotReached = false; //Assume at first that the current semi-stationary point is in fact the optimum. 
                double[] gradientFactor1 = new double[m1 + k1];
                double[] gradientFactor2 = new double[m2 + k2];
                
                
                //System.out.println("It entered the gradient small if");
                //System.out.println("IVI1: " + Arrays.toString(IVI1));
                //System.out.println("IVI2: " + Arrays.toString(IVI2));
                
                
                for (int i = 0; i < N1.size(); i++){
                    if ((map2LooseGradient1[IVI1[N1.get(i)]] == -1) && (map2LooseGradient1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]] == -1)){
                        if ((gradient1[IVI1[N1.get(i)]] - gradient1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]]) < 0){
                            promisingEN1.add(N1.get(i));
                        }
                    } else if (map2LooseGradient1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]] == -1){
                        if (LooseGradient1.get(map2LooseGradient1[IVI1[N1.get(i)]]) != 0){
                            gradientFactor1[IVI1[N1.get(i)]] = Math.max(0, gradient1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]]/LooseGradient1.get(map2LooseGradient1[IVI1[N1.get(i)]]).doubleValue());
                        } else {
                            if (gradient1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]] > 0){
                                promisingEN1.add(N1.get(i));
                            }
                        }
                    } else if (map2LooseGradient1[IVI1[N1.get(i)]] == -1) {
                        if (LooseGradient1.get(map2LooseGradient1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]]) != 0){
                            gradientFactor1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]] = Math.min(gradientFactor1[IVI1[B1.get(OE1.getBackMap(N1.get(i)))]],Math.min(gradient1[IVI1[N1.get(i)]]/LooseGradient1.get(IVI1[B1.get(OE1.getBackMap(N1.get(i)))]).doubleValue(),0));
                        } else {
                            if (gradient1[IVI1[N1.get(i)]] < 0){
                                promisingEN1.add(N1.get(i));
                            }
                        }
                    }
                }
                
                for (int i = 0; i < N2.size(); i++){
                    if ((map2LooseGradient2[IVI2[N2.get(i)]] == -1) && (map2LooseGradient2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]] == -1)){
                        if ((gradient2[IVI2[N2.get(i)]] - gradient2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]]) < 0){
                            promisingEN2.add(N2.get(i));
                        }
                    } else if (map2LooseGradient2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]] == -1){
                        if (LooseGradient2.get(map2LooseGradient2[IVI2[N2.get(i)]]) != 0){
                            gradientFactor2[IVI2[N2.get(i)]] = Math.max(0, gradient2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]]/LooseGradient2.get(map2LooseGradient2[IVI2[N2.get(i)]]).doubleValue());
                        } else {
                            if (gradient2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]] > 0){
                                promisingEN2.add(N2.get(i));
                            }
                        }
                    } else if (map2LooseGradient2[IVI2[N2.get(i)]] == -1) {
                        if (LooseGradient2.get(map2LooseGradient2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]]) != 0){
                            gradientFactor2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]] = Math.min(gradientFactor2[IVI2[B2.get(OE2.getBackMap(N2.get(i)))]],Math.min(gradient2[IVI2[N2.get(i)]]/LooseGradient2.get(IVI2[B2.get(OE2.getBackMap(N2.get(i)))]).doubleValue(),0));
                        } else {
                            if (gradient2[IVI2[N2.get(i)]] < 0){
                                promisingEN2.add(N2.get(i));
                            }
                        }
                    }    
                }
                
                //System.out.println("Factors gradient1: " + Arrays.toString(gradientFactor1));
                //System.out.println("Factors gradient2: " + Arrays.toString(gradientFactor2));
                //pause = scan.next();
                
                
                for (Integer ks : map2Gradients1.keySet()){
                    double FactorLength = 0;
                    for (Integer t : map2Gradients1.get(ks)){
                        FactorLength += gradientFactor1[t.intValue()]*gradientFactor1[t.intValue()];
                    }
                    
                    if (FactorLength > 1){
                        for (int i = 0; i < N1.size(); i++){
                            if (map2LooseGradient1[IVI1[N1.get(i)]] == ks.intValue()){
                                promisingEN1.add(N1.get(i));
                            }
                        }
                    }
                }
                
                for (Integer ks : map2Gradients2.keySet()){
                    double FactorLength = 0;
                    for (Integer t : map2Gradients2.get(ks)){
                        FactorLength += gradientFactor2[t.intValue()]*gradientFactor2[t.intValue()];
                    }
                    
                    if (FactorLength > 1){
                        for (int i = 0; i < N2.size(); i++){
                            if (map2LooseGradient2[IVI2[N2.get(i)]] == ks.intValue()){
                                promisingEN2.add(N2.get(i));
                            }
                        }
                    }
                }
                
                if ((promisingEN1.size()>0) || (promisingEN2.size()>0)){
                    
                    //System.out.println("It not optimum because: ");
                    //System.out.println("Promising EN1: " + promisingEN1);
                    
                    //System.out.println("Promising EN2: " + promisingEN2);
                    //pause = scan.next();
                    N1.removeAll(promisingEN1);
                    S1.addAll(promisingEN1);
                    
                    N2.removeAll(promisingEN2);
                    S2.addAll(promisingEN2);
                    
                    conjugate_initial_counter = 0;
                    
                    
                    optimNotReached = true;
                    
                }
                
                continue;//We go back to the beginning of the loop. 
            }
            
            //We know need to determine the best direction of change depending on whether we are in the first iteration of a re=initialization of the conjutage gradient method or not
            
            if (conjugate_initial_counter == 0){
                //System.out.println("Should be equal");
                dDirectionxs1 = new double[S1.size()];
                dDirectionxs2 = new double[S2.size()];
                for (int i = 0; i < dDirectionxs1.length; i++){
                    dDirectionxs1[i] = -gradientxs1[i];
                }
                for (int i = 0; i < dDirectionxs2.length; i++){
                    dDirectionxs2[i] = -gradientxs2[i];
                }
            } else {
                double akNum = 0;
                for (int i = 0; i < gradientxs1.length; i++){
                    akNum += gradientxs1[i]*(gradientxs1[i] - gradientxs1Prev[i]);
                }
                for (int i = 0; i < gradientxs2.length; i++){
                    akNum += gradientxs2[i]*(gradientxs2[i] - gradientxs2Prev[i]);
                }
                
                double ak = akNum/akDenom;
                //System.out.println("ak is " + ak);
                
                for (int i = 0; i < gradientxs1.length; i++){
                    dDirectionxs1[i] = ak*dDirectionxs1[i] - gradientxs1[i];
                }
                for (int i = 0; i < gradientxs2.length; i++){
                    dDirectionxs2[i] = ak*dDirectionxs2[i] - gradientxs2[i];
                }
            }
            
            //Computing the complete change vector
            
            double[] dDirection1 = new double[m1 + k1];
            
            for (int i = 0; i < S1.size(); i++){
                dDirection1[IVI1[S1.get(i)]] = dDirectionxs1[i];
                dDirection1[IVI1[B1.get(OE1.getBackMap(S1.get(i)))]] += -dDirectionxs1[i];
            }
            
            double[] dDirection2 = new double[m2 + k2];
            
            for (int i = 0; i < S2.size(); i++){
                dDirection2[IVI2[S2.get(i)]] = dDirectionxs2[i];
                dDirection2[IVI2[B2.get(OE2.getBackMap(S2.get(i)))]] += -dDirectionxs2[i];
            }
            
            //System.out.println("The dDirectionxs1: " + Arrays.toString(dDirectionxs1));
            //System.out.println("The dDirectionxs2: " + Arrays.toString(dDirectionxs2));
            //System.out.println("The dDirection1: " + Arrays.toString(dDirection1));
            //System.out.println("The dDirection2: " + Arrays.toString(dDirection2));
            //pause = scan.next();
            
            //Determining the closed set for tau, in order to mantain all edges with positive size. 
            double tau_max = 0;
            double tau_min = 0;
            boolean tauNeedsChange = true;
            
            Vector<PhyloTreeEdge> EdgesT1 = T1.getEdges();
            Vector<PhyloTreeEdge> EdgesT2 = T2.getEdges();
            
            Vector<Integer> potentialN1 = new Vector<Integer>();
            Vector<Integer> potentialN2 = new Vector<Integer>();
            
            //We need to check the external edges in the original tree to find tau_max as well
           
            for (int i = 0; i < ol1; i++){
                if (dDirection1[IVI1[i]] < 0){
                    //System.out.println("For ol1 " + i + "it found dDirection1 < 0");
                    if(tauNeedsChange || (-T1.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(i)].get(0)/dDirection1[IVI1[i]] < tau_max)){
                        tau_max = -T1.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(i)].get(0)/dDirection1[IVI1[i]];
                        //System.out.println("Tau max changed in ol1 by " + i);
                        if (!N1.contains(i)){
                            potentialN1.clear();
                            potentialN1.add(i);
                        } else {
                            System.out.println("An element on N1 sneaked in (situation 1): "+ i);
                        }
                        tauNeedsChange = false;
                    } else if(-T1.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(i)].get(0)/dDirection1[IVI1[i]] == tau_max){
                        if (!N1.contains(i)){
                            potentialN1.add(i);
                        }else {
                            System.out.println("An element on N1 sneaked in (situation 2): "+ i);
                        } 
                    }
                }
            }
            for (int i = 0; i < ol2; i++){
                if (dDirection2[IVI2[i]] < 0){
                    //System.out.println("For ol2 " + i + "it found dDirection2 < 0");
                    if(tauNeedsChange || (-T2.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(i)].get(0)/dDirection2[IVI2[i]] < tau_max)){
                        tau_max = -T2.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(i)].get(0)/dDirection2[IVI2[i]];
                        //System.out.println("Tau max changed in ol2 by " + i);
                        if (!N2.contains(i)){
                            potentialN2.clear();
                            potentialN2.add(i);
                        } else {
                            System.out.println("An element on N2 sneaked in (situation 1): "+ i);
                        }
                        tauNeedsChange = false;
                    } else if(-T2.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(i)].get(0)/dDirection2[IVI2[i]] == tau_max){
                        if (!N2.contains(i)){
                            potentialN2.add(i);
                        }else {
                            System.out.println("An element on N2 sneaked in (situation 2): "+ i);
                        } 
                    }
                }
            }
            
            for (int i = 0; i < cur1Edges2Axis.size(); i++){
                if(IVI1[cur1Edges2Axis.get(i) + ol1] != -1){
                    if (dDirection1[IVI1[cur1Edges2Axis.get(i) + ol1]] < 0){
                        //System.out.println("For cur1Edges2Axis " + i + "it found dDirection1 < 0");
                        if(tauNeedsChange || (-EdgesT1.get(i).getNorm()/dDirection1[IVI1[cur1Edges2Axis.get(i) + ol1]] < tau_max)){
                            tau_max = -EdgesT1.get(i).getNorm()/dDirection1[IVI1[cur1Edges2Axis.get(i) + ol1]];
                            //System.out.println("Tau max changed in cur1Edges2axis by " + cur1Edges2Axis.get(i));
                            if (!N1.contains(cur1Edges2Axis.get(i)+ol1)){
                                potentialN1.clear();
                                potentialN1.add(cur1Edges2Axis.get(i)+ol1);
                            } else {
                                //System.out.println("An element on N1 sneaked in (situation 1): "+ (cur1Edges2Axis.get(i)+ol1));
                            }
                            tauNeedsChange = false;
                        } else if (-EdgesT1.get(i).getNorm()/dDirection1[IVI1[cur1Edges2Axis.get(i) + ol1]] == tau_max){
                            if (!N1.contains(cur1Edges2Axis.get(i)+ol1)){
                                potentialN1.add(cur1Edges2Axis.get(i)+ol1);
                            }else {
                                //System.out.println("An element on N1 sneaked in (situation 2): "+ (cur1Edges2Axis.get(i)+ol1));
                            }  
                        }
                    }
                }
            }
            
            for (int i = 0; i < cur2Edges2Axis.size(); i++){
                if(IVI2[cur2Edges2Axis.get(i) + ol2] != -1){
                    if (dDirection2[IVI2[cur2Edges2Axis.get(i) + ol2]] < 0){
                        //System.out.println("For cur2Edges2Axis " + i + "it found dDirection2 < 0");
                        if (tauNeedsChange || (-EdgesT2.get(i).getNorm()/dDirection2[IVI2[cur2Edges2Axis.get(i) + ol2]] < tau_max)){
                            tau_max = -EdgesT2.get(i).getNorm()/dDirection2[IVI2[cur2Edges2Axis.get(i) + ol2]];
                            //System.out.println("Tau max changed in cur2Edges2axis by " + cur1Edges2Axis.get(i));
                            if (!N2.contains(cur2Edges2Axis.get(i) + ol2)){
                                potentialN1.clear();
                                potentialN2.clear();
                                potentialN2.add(cur2Edges2Axis.get(i) + ol2);
                            } else{
                                //System.out.println("An element on N2 sneaked in (situation 1): "+ (cur2Edges2Axis.get(i) + ol2));
                            }
                            tauNeedsChange = false;
                        } else if (-EdgesT2.get(i).getNorm()/dDirection2[IVI2[cur2Edges2Axis.get(i) + ol2]] == tau_max){
                            if (!N2.contains(cur2Edges2Axis.get(i) + ol2)){
                                potentialN2.add(cur2Edges2Axis.get(i) + ol2);
                            } else{
                                //System.out.println("An element on N2 sneaked in (situation 2): "+ (cur2Edges2Axis.get(i) + ol2));
                            }
                        }
                    }
                }
            }
            
            // We will look for the tau that minimizes f(x + tau* dDirection) between tau_min and tau_max.
            //We will first check if the minimum is the actual tau_max
            
            //System.out.println("The tau_max is: " + tau_max);
            
            //Defining new values of the trees to compute geodesic and find the derivative; 
            Vector<PhyloTreeEdge> conjEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
            Vector<PhyloTreeEdge> conjEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
            
            
            //Computing the new values of the interior edges of the trees by moving in the direction of change
            for (int i = 0; i < cur1Edges2Axis.size(); i++){
                if(IVI1[cur1Edges2Axis.get(i)+ol1] != -1){
                    double[] tempVecEA = {EdgesT1.get(i).getNorm() + (tau_max-0.0000000000001)*dDirection1[IVI1[cur1Edges2Axis.get(i)+ol1]]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    conjEdgesT1.get(i).setAttribute(tempEA);
                    if(ET1toET2.containsKey(Integer.valueOf(i))){
                        conjEdgesT2.get(ET1toET2.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                    }
                }
            }
            
            for (int i = 0; i < cur2Edges2Axis.size(); i++){
                if(IVI2[cur2Edges2Axis.get(i)+ol2] != -1){
                    double[] tempVecEA = {EdgesT2.get(i).getNorm() + (tau_max-0.0000000000001)*dDirection2[IVI2[cur2Edges2Axis.get(i)+ol2]]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    conjEdgesT2.get(i).setAttribute(tempEA);
                    if(ET2toET1.containsKey(Integer.valueOf(i))){
                        conjEdgesT1.get(ET2toET1.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                    }
                }
            }
            
            //We also need to change the values in the Leaf Edge attribs in the unrestricted case. 
            
            EdgeAttribute[] T1LeafEdgeAtt = T1.getCopyLeafEdgeAttribs();
            EdgeAttribute[] T2LeafEdgeAtt = T2.getCopyLeafEdgeAttribs();
            
            
            for (int i = 0; i < ol1; i++){
                double[] tempVecEA = {T1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].get(0) + (tau_max-0.0000000000001)*dDirection1[IVI1[i]]};
                EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                T1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                if(OE2.getCompLeaves2orgLeaves(OE1.getOrgLeaves2compLeaves(i)) == -1){
                    T2LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                }
            }
            for (int i = 0; i < ol2; i++){
                double[] tempVecEA = {T2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].get(0) + (tau_max-0.0000000000001)*dDirection2[IVI2[i]]};
                EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                T2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                if(OE1.getCompLeaves2orgLeaves(OE2.getOrgLeaves2compLeaves(i)) == -1){
                    T1LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                }
            }
            
            PhyloTree conjT1 = new PhyloTree(conjEdgesT1, T1.getLeaf2NumMap(), T1LeafEdgeAtt, false);
            PhyloTree conjT2 = new PhyloTree(conjEdgesT2, T2.getLeaf2NumMap(), T2LeafEdgeAtt, false);
            
            //Computing geodesic in between these trees. 
            Geodesic conjGeode = parPolyMain.getGeodesic(conjT1, conjT2);
            
            //System.out.println("   conjT1: \n" + treePrinter.toString(conjT1)+"\n \n");
            //System.out.println("   conjT2: \n" + treePrinter.toString(conjT2)+"\n \n");
            //System.out.println("The geodesic summary is: " + treePrinter.toString(conjGeode, OE1.getCompleteLeafSet()));
            
            RatioSequence conjRSeq = conjGeode.getRS();//The derivative will depend on the ratio sequence 
            
            //And on which common edges they have
            Vector<PhyloTreeEdge> conjCEs = conjGeode.getCommonEdges(); 
            
            double derivTau = 0;//Where the final derivative for tau will be saved
            
            //For each ratio we compute the contribution of the expression relating to the ratio in derivative with respect to tau_max
            Iterator<Ratio> conjRSIter = conjRSeq.iterator();
            while(conjRSIter.hasNext()){
                Ratio rat = (Ratio) conjRSIter.next();
                //Values that will contribute to the derivative of the ratio
                
                if (rat.getELength() > 0){
                    double ENum = 0;
                    for (PhyloTreeEdge e : rat.getEEdges()){
                        int eID = conjT1.getEdges().indexOf(e);
                        ENum += dDirection1[IVI1[cur1Edges2Axis.get(eID) + ol1]]*e.getNorm();
                    }
                    derivTau += ENum*(1 + (rat.getFLength()/rat.getELength()));
                }
                
                if (rat.getFLength() > 0){
                    double FNum = 0;
                    for (PhyloTreeEdge e : rat.getFEdges()){
                        int eID = conjT2.getEdges().indexOf(e);
                        FNum += dDirection2[IVI2[cur2Edges2Axis.get(eID) + ol2]]*e.getNorm();
                    }
                    derivTau += FNum*(1 + (rat.getELength()/rat.getFLength()));
                }
            }
            
            for(PhyloTreeEdge e : conjCEs){
                int eID1 = conjT1.getSplits().indexOf(e.asSplit()); 
                int eID2 = conjT2.getSplits().indexOf(e.asSplit());
                
                if ((eID1 == -1) && (eID2 != -1)){
                    //System.out.println("Warning 1.1: " + treePrinter.toString(e, OE1.getCompleteLeafSet()));
                    derivTau += dDirection2[IVI2[cur2Edges2Axis.get(eID2) + ol2]]*(conjT2.getEdge(eID2).getAttribute().get(0)); 
                }else if ((eID1 != -1) && (eID2 == -1)){
                    //System.out.println("Warning 1.2: " + treePrinter.toString(e, OE2.getCompleteLeafSet()));
                    derivTau += dDirection1[IVI1[cur1Edges2Axis.get(eID1) + ol1]]*(conjT1.getEdge(eID1).getAttribute().get(0)); 
                } else if ((eID1 != -1) && (eID2 != -1)){
                    if ((IVI1[cur1Edges2Axis.get(eID1) + ol1] != -1) && (IVI2[cur2Edges2Axis.get(eID2) + ol2] != -1)){
                        derivTau += (dDirection1[IVI1[cur1Edges2Axis.get(eID1) + ol1]] - dDirection2[IVI2[cur2Edges2Axis.get(eID2) + ol2]])*(conjT1.getEdge(eID1).getAttribute().get(0) - conjT2.getEdge(eID2).getAttribute().get(0)); 
                    }
                }   
                
            } 
            
            //In the unrestricted case, we need to also consider the contribution of the external edge to the gradient. 
            
            for (int i = 0; i < ol1; i++){
                derivTau += dDirection1[IVI1[i]]*(T1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].get(0) - T2LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].get(0));
            }
            for (int i = 0; i < ol2; i++){
                derivTau += dDirection2[IVI2[i]]*(T2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].get(0) - T1LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].get(0));
            }

            double tau = 0;
            
            //System.out.println("   Tau derivative for tau_max ended being "+ derivTau);
            
            if (derivTau <= 0){//In this case the minimum is reached right at the tau_max limit and the search is over.
                //System.out.println("   So it hitted a face");
                //System.out.println("    Potential N1 is " + potentialN1);
                //System.out.println("    Potential N2 is " + potentialN2);
                
                //System.out.println("    Already N1 is " + alreadyN1);
                //System.out.println("    Already N2 is " + alreadyN2);
                tau = tau_max;
                boolean ChangeInIndexMade = false;
                
                //We have hitted a boundary face, so we need to reclasify some variable to N1 or N2. 
                
                if (potentialN1.size() > 0){
                    //We want to give priority to indexes that have not been non-basic variables yet to try and avoid cycling. 
                    int[] IndexListOrdered = new int[potentialN1.size()];
                    int leftInd = 0;
                    int rightInd = potentialN1.size() - 1;
                    boolean AllN1Already = true;
                    for (int i = 0; i < potentialN1.size(); i++){
                        if (alreadyN1.contains(potentialN1.get(i))){
                            IndexListOrdered[rightInd] = i;
                            rightInd--;
                        } else {
                            IndexListOrdered[leftInd] = i;
                            leftInd++;
                            AllN1Already = false;
                        }
                    }
                    if (AllN1Already){
                        Collections.shuffle(potentialN1);
                    }
                    for (int i : IndexListOrdered){
                        if (S1.contains(potentialN1.get(i))){
                            N1.add(potentialN1.get(i));
                            alreadyN1.add(potentialN1.get(i));
                            S1.remove(Integer.valueOf(potentialN1.get(i)));
                            ChangeInIndexMade = true;
                            break;
                        } else if (B1.contains(potentialN1.get(i))){
                            int rowIndexTemp = B1.indexOf(potentialN1.get(i));
                            int newB1element = -1; 
                            for (int j : OE1.getMapList().get(rowIndexTemp)){
                                if (S1.contains(j)){
                                    newB1element = j;
                                    break;
                                }
                            }
                            if (newB1element != -1){
                                B1.set(rowIndexTemp, newB1element);
                                S1.remove(Integer.valueOf(newB1element));
                                N1.add(potentialN1.get(i));
                                alreadyN1.add(potentialN1.get(i));
                                ChangeInIndexMade = true;
                                break;
                            }
                        }
                    }
                } else if (potentialN2.size() > 0){
                    //We want to give priority to indexes that have not been non-basic variables yet to try and avoid cycling. 
                    int[] IndexListOrdered = new int[potentialN2.size()];
                    int leftInd = 0;
                    int rightInd = potentialN2.size() - 1;
                    boolean AllN2Already = true;
                    for (int i = 0; i < potentialN2.size(); i++){
                        if (alreadyN2.contains(potentialN2.get(i))){
                            IndexListOrdered[rightInd] = i;
                            rightInd--;
                        } else {
                            IndexListOrdered[leftInd] = i;
                            leftInd++;
                            AllN2Already = false;
                        }
                    }
                    if (AllN2Already){
                        Collections.shuffle(potentialN2);
                    }
                    for (int i : IndexListOrdered){
                        if (S2.contains(potentialN2.get(i))){
                            N2.add(potentialN2.get(i));
                            alreadyN2.add(potentialN2.get(i));
                            S2.remove(Integer.valueOf(potentialN2.get(i)));
                            ChangeInIndexMade = true;
                            break;
                        } else if (B2.contains(potentialN2.get(i))){
                            int rowIndexTemp = B2.indexOf(potentialN2.get(i));
                            int newB2element = -1; 
                            for (int j : OE2.getMapList().get(rowIndexTemp)){
                                if (S2.contains(j)){
                                    newB2element = j;
                                    break;
                                }
                            }
                            if (newB2element != -1){
                                B2.set(rowIndexTemp, newB2element);
                                S2.remove(Integer.valueOf(newB2element));
                                N2.add(potentialN2.get(i));
                                alreadyN2.add(potentialN2.get(i));
                                ChangeInIndexMade = true;
                                break;
                            }
                        }
                    }
                }
                if (!ChangeInIndexMade){
                    System.out.println("ERROR: Although a variable should be reclassified as non-basic, it did not happen.");
                    break;
                }
                
                conjugate_initial_counter = 0; // We are re-initializing the conjugate gradient method in a new face;
                
                //Defining the new trees to go back to the main while loop: 
                
                /**System.out.println("    Just before new trees definition tau is " + tau);
                
                System.out.println("    The edges of T1 are " + EdgesT1);
                
                System.out.println("    The edges of T2 are " + EdgesT2);
                
                System.out.println("    S1 = " + S1);
                System.out.println("    B1 = " + B1);
                System.out.println("    S2 = " + S2);
                System.out.println("    B2 = " + B2);
                
                System.out.println("    And the directions are " + Arrays.toString(dDirection1) + " and " + Arrays.toString(dDirection2));*/
                
                
                Vector<PhyloTreeEdge> newEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
                Vector<PhyloTreeEdge> newEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
                double[] newEdgesValuesT1 = new double[EdgesT1.size()];
                double[] newEdgesValuesT2 = new double[EdgesT2.size()];
                double[] newLeafEdgesvaluesT1 = new double[ol1];
                double[] newLeafEdgesvaluesT2 = new double[ol2];
                
                for (int i = 0; i < B1.size(); i++){
                    if (B1.get(i) < ol1){
                        newLeafEdgesvaluesT1[B1.get(i)] = OE1.getFixedLengths(i);
                    } else {
                        newEdgesValuesT1[cur1Axis2Edges[B1.get(i) - ol1]] = OE1.getFixedLengths(i);
                    }
                }
                for (int i = 0; i < S1.size(); i++){
                    if (S1.get(i) < ol1){
                        newLeafEdgesvaluesT1[S1.get(i)] = T1.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(S1.get(i))].get(0) + tau*dDirection1[IVI1[S1.get(i)]];
                        if (B1.get(OE1.getBackMap(S1.get(i))) < ol1){
                            newLeafEdgesvaluesT1[B1.get(OE1.getBackMap(S1.get(i)))] -= newLeafEdgesvaluesT1[S1.get(i)];
                        } else {
                            newEdgesValuesT1[cur1Axis2Edges[B1.get(OE1.getBackMap(S1.get(i))) - ol1]] -= newLeafEdgesvaluesT1[S1.get(i)];
                        }
                    } else {
                        newEdgesValuesT1[cur1Axis2Edges[S1.get(i) - ol1]] = EdgesT1.get(cur1Axis2Edges[S1.get(i)-ol1]).getNorm() + tau*dDirection1[IVI1[S1.get(i)]];
                        if (B1.get(OE1.getBackMap(S1.get(i))) < ol1){
                            newLeafEdgesvaluesT1[B1.get(OE1.getBackMap(S1.get(i)))] -= newEdgesValuesT1[cur1Axis2Edges[S1.get(i) - ol1]];
                        } else {
                            newEdgesValuesT1[cur1Axis2Edges[B1.get(OE1.getBackMap(S1.get(i))) - ol1]] -= newEdgesValuesT1[cur1Axis2Edges[S1.get(i) - ol1]];
                        }
                    }
                }
                
                for (int i = 0; i < B2.size(); i++){
                    if (B2.get(i) < ol2){
                        newLeafEdgesvaluesT2[B2.get(i)] = OE2.getFixedLengths(i);
                    } else {
                        newEdgesValuesT2[cur2Axis2Edges[B2.get(i) - ol2]] = OE2.getFixedLengths(i);
                    }
                }
                for (int i = 0; i < S2.size(); i++){
                    if (S2.get(i) < ol2){
                        newLeafEdgesvaluesT2[S2.get(i)] = T2.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(S2.get(i))].get(0) + tau*dDirection2[IVI2[S2.get(i)]];
                        if (B2.get(OE2.getBackMap(S2.get(i))) < ol2){
                            newLeafEdgesvaluesT2[B2.get(OE2.getBackMap(S2.get(i)))] -= newLeafEdgesvaluesT2[S2.get(i)];
                        } else {
                            newEdgesValuesT2[cur2Axis2Edges[B2.get(OE2.getBackMap(S2.get(i))) - ol2]] -= newLeafEdgesvaluesT2[S2.get(i)];
                        }
                    } else {
                        newEdgesValuesT2[cur2Axis2Edges[S2.get(i) - ol2]] = EdgesT2.get(cur2Axis2Edges[S2.get(i)-ol2]).getNorm() + tau*dDirection2[IVI2[S2.get(i)]];
                        if (B2.get(OE2.getBackMap(S2.get(i))) < ol2){
                            newLeafEdgesvaluesT2[B2.get(OE2.getBackMap(S2.get(i)))] -= newEdgesValuesT2[cur2Axis2Edges[S2.get(i) - ol2]];
                        } else {
                            newEdgesValuesT2[cur2Axis2Edges[B2.get(OE2.getBackMap(S2.get(i))) - ol2]] -= newEdgesValuesT2[cur2Axis2Edges[S2.get(i) - ol2]];
                        }
                    }
                }
            
                //Computing the new values of the interior edges of the trees 
                for (int i = 0; i < cur1Edges2Axis.size(); i++){
                    if(IVI1[cur1Edges2Axis.get(i) + ol1] != -1){
                        double[] tempVecEA = {newEdgesValuesT1[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        newEdgesT1.get(i).setAttribute(tempEA);
                        if(ET1toET2.containsKey(Integer.valueOf(i))){
                            newEdgesT2.get(ET1toET2.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                        }
                    }
                }
            
                for (int i = 0; i < cur2Edges2Axis.size(); i++){
                    if(IVI2[cur2Edges2Axis.get(i) + ol2] != -1){
                        double[] tempVecEA = {newEdgesValuesT2[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        newEdgesT2.get(i).setAttribute(tempEA);
                        if(ET2toET1.containsKey(Integer.valueOf(i))){
                            newEdgesT1.get(ET2toET1.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                        }
                    }
                }
                
                //We also need to change the values in the Leaf Edge attribs in the unrestricted case. 
                
                
                EdgeAttribute[] newT1LeafEdgeAtt = T1.getCopyLeafEdgeAttribs();
                EdgeAttribute[] newT2LeafEdgeAtt = T2.getCopyLeafEdgeAttribs();
            
                for (int i = 0; i < ol1; i++){
                    double[] tempVecEA = {newLeafEdgesvaluesT1[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newT1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    if(OE2.getCompLeaves2orgLeaves(OE1.getOrgLeaves2compLeaves(i)) == -1){
                        newT2LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    }
                }
                for (int i = 0; i < ol2; i++){
                    double[] tempVecEA = {newLeafEdgesvaluesT2[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newT2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    if(OE1.getCompLeaves2orgLeaves(OE2.getOrgLeaves2compLeaves(i)) == -1){
                        newT1LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    }
                }
            
            
                T1 = new PhyloTree(newEdgesT1, T1.getLeaf2NumMap(), newT1LeafEdgeAtt, false);
                T2 = new PhyloTree(newEdgesT2, T2.getLeaf2NumMap(), newT2LeafEdgeAtt, false);
            
                tempGeode = parPolyMain.getGeodesic(T1, T2);
                
            } else {//We still need to find the optimum tau for this case. 
                //System.out.println("   So we are still in the same face");
                int counterWhile = 0;
                tau = 0.1;
                if (tau > tau_max/2){
                    tau = tau_max/2;
                }
                //System.out.println("    Prev tau = " + tau);
                while(((derivTau < -0.0000000000000001) || (derivTau > 0.0000000000000001)) && (((tau_max - tau_min) > 0.0000000001))){ //&&(counterWhile < 50)
                    counterWhile++;
                    //System.out.println("   INSIDE the tau while loop " + counterWhile);
                    tau = (tau_max + tau_min)/2;
                    //System.out.println("   with Taus: [" + tau_min + " < " + tau + " < " + tau_max + "]");
                    
                    conjEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
                    conjEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
                    
                    //Computing the new values of the interior edges of the trees by moving in the direction of change
                    for (int i = 0; i < cur1Edges2Axis.size(); i++){
                        if(IVI1[cur1Edges2Axis.get(i) + ol1] != -1){
                            double[] tempVecEA = {EdgesT1.get(i).getNorm() + tau*dDirection1[IVI1[cur1Edges2Axis.get(i)+ol1]]};
                            EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                            conjEdgesT1.get(i).setAttribute(tempEA);
                            if(ET1toET2.containsKey(Integer.valueOf(i))){
                                conjEdgesT2.get(ET1toET2.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                            }
                        }
                            
                    }
            
                    for (int i = 0; i < cur2Edges2Axis.size(); i++){
                        if(IVI2[cur2Edges2Axis.get(i) + ol2] != -1){
                            double[] tempVecEA = {EdgesT2.get(i).getNorm() + tau*dDirection2[IVI2[cur2Edges2Axis.get(i)+ol2]]};
                            EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                            conjEdgesT2.get(i).setAttribute(tempEA);
                            if(ET2toET1.containsKey(Integer.valueOf(i))){
                                conjEdgesT1.get(ET2toET1.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                            }
                        }
                    }
            
                    //We also need to change the values in the Leaf Edge attribs in the unrestricted case. 
            
                    T1LeafEdgeAtt = T1.getCopyLeafEdgeAttribs();
                    T2LeafEdgeAtt = T2.getCopyLeafEdgeAttribs();
            
                    for (int i = 0; i < ol1; i++){
                        double[] tempVecEA = {T1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].get(0) + tau*dDirection1[IVI1[i]]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        T1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                        if(OE2.getCompLeaves2orgLeaves(OE1.getOrgLeaves2compLeaves(i)) == -1){
                            T2LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                        }
                    }
                    
                    for (int i = 0; i < ol2; i++){
                        double[] tempVecEA = {T2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].get(0) + tau*dDirection2[IVI2[i]]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        T2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                        if(OE1.getCompLeaves2orgLeaves(OE2.getOrgLeaves2compLeaves(i)) == -1){
                            T1LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                        }
                    }
                    
                    conjT1 = new PhyloTree(conjEdgesT1, T1.getLeaf2NumMap(), T1LeafEdgeAtt, false);
                    conjT2 = new PhyloTree(conjEdgesT2, T2.getLeaf2NumMap(), T2LeafEdgeAtt, false);
                    
                    /*System.out.println("   With tau =" + tau);
                    System.out.println("    The original trees are: ");
                    System.out.println("      T1: \n" + treePrinter.toString(T1)+"\n");
                    System.out.println("      T2: \n" + treePrinter.toString(T2)+"\n");
                    System.out.println("    The direction of change is: ");
                    System.out.println("      Dd1: \n" + Arrays.toString(dDirection1)+"\n");
                    System.out.println("      Dd2: \n" + Arrays.toString(dDirection2)+"\n");
                    System.out.println("    So the 'new' trees are: ");
                    System.out.println("      T1_new: \n" + treePrinter.toString(conjT1)+"\n");
                    System.out.println("      T2_new: \n" + treePrinter.toString(conjT2)+"\n");*/
                    
                    //Computing geodesic in between these trees. 
                    conjGeode = parPolyMain.getGeodesic(conjT1, conjT2);
            
                    conjRSeq = conjGeode.getRS();//The derivative will depend on the ratio sequence 
            
                    //And on which common edges they have
                    conjCEs = conjGeode.getCommonEdges(); 
                    
                    derivTau = 0;//Where the final derivative for tau will be saved
                    
                    //For each ratio we compute the contribution of the expression relating to the ratio in derivative with respect to tau_max
                    conjRSIter = conjRSeq.iterator();
                    while(conjRSIter.hasNext()){
                        Ratio rat = (Ratio) conjRSIter.next();
                        //Values that will contribute to the derivative of the ratio
                
                        if (rat.getELength() > 0){
                            double ENum = 0;
                            for (PhyloTreeEdge e : rat.getEEdges()){
                                int eID = conjT1.getEdges().indexOf(e);
                                ENum += dDirection1[IVI1[cur1Edges2Axis.get(eID) + ol1]]*e.getNorm();
                            }
                            derivTau += ENum*(1 + (rat.getFLength()/rat.getELength()));
                        }
                
                        if (rat.getFLength() > 0){
                            double FNum = 0;
                            for (PhyloTreeEdge e : rat.getFEdges()){
                                int eID = conjT2.getEdges().indexOf(e);
                                FNum += dDirection2[IVI2[cur2Edges2Axis.get(eID) + ol2]]*e.getNorm();
                            }
                            derivTau += FNum*(1 + (rat.getELength()/rat.getFLength()));
                        }
                    }
                    
                    //For each common edge, we compute the contribution to the derivatives in the gradient. 
                    for(PhyloTreeEdge e : conjCEs){
                        int eID1 = conjT1.getSplits().indexOf(e.asSplit());
                        int eID2 = conjT2.getSplits().indexOf(e.asSplit());
                        
                        if ((eID1 == -1) && (eID2 != -1)){
                            //System.out.println("Warning 2.1: " + treePrinter.toString(e, OE1.getCompleteLeafSet()));
                            derivTau += dDirection2[IVI2[cur2Edges2Axis.get(eID2) + ol2]]*(conjT2.getEdge(eID2).getAttribute().get(0)); 
                        }else if ((eID1 != -1) && (eID2 == -1)){
                            //System.out.println("Warning 2.2: " + treePrinter.toString(e, OE2.getCompleteLeafSet()));
                            derivTau += dDirection1[IVI1[cur1Edges2Axis.get(eID1) + ol1]]*(conjT1.getEdge(eID1).getAttribute().get(0)); 
                        } else if ((eID1 != -1) && (eID2 != -1)){
                            if ((IVI1[cur1Edges2Axis.get(eID1) + ol1] != -1) && (IVI2[cur2Edges2Axis.get(eID2) + ol2] != -1)){
                                derivTau += (dDirection1[IVI1[cur1Edges2Axis.get(eID1) + ol1]] - dDirection2[IVI2[cur2Edges2Axis.get(eID2) + ol2]])*(conjT1.getEdge(eID1).getAttribute().get(0) - conjT2.getEdge(eID2).getAttribute().get(0)); 
                            }
                        } 
                    }
                    
                    //In the unrestricted case, we need to also consider the contribution of the external edge to the gradient. 
                    
                    for (int i = 0; i < ol1; i++){
                        derivTau += dDirection1[IVI1[i]]*(T1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].get(0) - T2LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].get(0));
                    }
                    for (int i = 0; i < ol2; i++){
                        derivTau += dDirection2[IVI2[i]]*(T2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].get(0) - T1LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].get(0));
                    }
                    
                    //System.out.println("        DerivTau = " + derivTau);
                    if (derivTau <= 0){// This would mean the minimum is between tau and tau_max
                        tau_min = tau;
                    } else {
                        tau_max = tau;
                    }
                }
                
                conjugate_initial_counter++; // Keeping count on how many loops we have done in this face. 
                
                //Defining the new trees to go back to the main while loop: 
                
                Vector<PhyloTreeEdge> newEdgesT1 = Tools.myVectorClonePhyloTreeEdge(EdgesT1);
                Vector<PhyloTreeEdge> newEdgesT2 = Tools.myVectorClonePhyloTreeEdge(EdgesT2);
                double[] newEdgesValuesT1 = new double[EdgesT1.size()];
                double[] newEdgesValuesT2 = new double[EdgesT2.size()];
                double[] newLeafEdgesvaluesT1 = new double[ol1];
                double[] newLeafEdgesvaluesT2 = new double[ol2];
                
                for (int i = 0; i < B1.size(); i++){
                    if (B1.get(i) < ol1){
                        newLeafEdgesvaluesT1[B1.get(i)] = OE1.getFixedLengths(i);
                    } else {
                        newEdgesValuesT1[cur1Axis2Edges[B1.get(i) - ol1]] = OE1.getFixedLengths(i);
                    }
                }
                for (int i = 0; i < S1.size(); i++){
                    if (S1.get(i) < ol1){
                        newLeafEdgesvaluesT1[S1.get(i)] = T1.getLeafEdgeAttribs()[OE1.getOrgLeaves2compLeaves(S1.get(i))].get(0) + tau*dDirection1[IVI1[S1.get(i)]];
                        if (B1.get(OE1.getBackMap(S1.get(i))) < ol1){
                            newLeafEdgesvaluesT1[B1.get(OE1.getBackMap(S1.get(i)))] -= newLeafEdgesvaluesT1[S1.get(i)];
                        } else {
                            newEdgesValuesT1[cur1Axis2Edges[B1.get(OE1.getBackMap(S1.get(i))) - ol1]] -= newLeafEdgesvaluesT1[S1.get(i)];
                        }
                    } else {
                        newEdgesValuesT1[cur1Axis2Edges[S1.get(i) - ol1]] = EdgesT1.get(cur1Axis2Edges[S1.get(i)-ol1]).getNorm() + tau*dDirection1[IVI1[S1.get(i)]];
                        if (B1.get(OE1.getBackMap(S1.get(i))) < ol1){
                            newLeafEdgesvaluesT1[B1.get(OE1.getBackMap(S1.get(i)))] -= newEdgesValuesT1[cur1Axis2Edges[S1.get(i) - ol1]];
                        } else {
                            newEdgesValuesT1[cur1Axis2Edges[B1.get(OE1.getBackMap(S1.get(i))) - ol1]] -= newEdgesValuesT1[cur1Axis2Edges[S1.get(i) - ol1]];
                        }
                    }
                }
                
                for (int i = 0; i < B2.size(); i++){
                    if (B2.get(i) < ol2){
                        newLeafEdgesvaluesT2[B2.get(i)] = OE2.getFixedLengths(i);
                    } else {
                        newEdgesValuesT2[cur2Axis2Edges[B2.get(i) - ol2]] = OE2.getFixedLengths(i);
                    }
                }
                for (int i = 0; i < S2.size(); i++){
                    if (S2.get(i) < ol2){
                        newLeafEdgesvaluesT2[S2.get(i)] = T2.getLeafEdgeAttribs()[OE2.getOrgLeaves2compLeaves(S2.get(i))].get(0) + tau*dDirection2[IVI2[S2.get(i)]];
                        if (B2.get(OE2.getBackMap(S2.get(i))) < ol2){
                            newLeafEdgesvaluesT2[B2.get(OE2.getBackMap(S2.get(i)))] -= newLeafEdgesvaluesT2[S2.get(i)];
                        } else {
                            newEdgesValuesT2[cur2Axis2Edges[B2.get(OE2.getBackMap(S2.get(i))) - ol2]] -= newLeafEdgesvaluesT2[S2.get(i)];
                        }
                    } else {
                        newEdgesValuesT2[cur2Axis2Edges[S2.get(i) - ol2]] = EdgesT2.get(cur2Axis2Edges[S2.get(i)-ol2]).getNorm() + tau*dDirection2[IVI2[S2.get(i)]];
                        if (B2.get(OE2.getBackMap(S2.get(i))) < ol2){
                            newLeafEdgesvaluesT2[B2.get(OE2.getBackMap(S2.get(i)))] -= newEdgesValuesT2[cur2Axis2Edges[S2.get(i) - ol2]];
                        } else {
                            newEdgesValuesT2[cur2Axis2Edges[B2.get(OE2.getBackMap(S2.get(i))) - ol2]] -= newEdgesValuesT2[cur2Axis2Edges[S2.get(i) - ol2]];
                        }
                    }
                }
                
                //Computing the new values of the interior edges of the trees 
                for (int i = 0; i < cur1Edges2Axis.size(); i++){
                    if(IVI1[cur1Edges2Axis.get(i) + ol1] != -1){
                        double[] tempVecEA = {newEdgesValuesT1[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        newEdgesT1.get(i).setAttribute(tempEA);
                        if(ET1toET2.containsKey(Integer.valueOf(i))){
                            newEdgesT2.get(ET1toET2.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                        }
                    }
                }
            
                for (int i = 0; i < cur2Edges2Axis.size(); i++){
                    if(IVI2[cur2Edges2Axis.get(i) + ol2] != -1){
                        double[] tempVecEA = {newEdgesValuesT2[i]};
                        EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                        newEdgesT2.get(i).setAttribute(tempEA);
                        if(ET2toET1.containsKey(Integer.valueOf(i))){
                            newEdgesT1.get(ET2toET1.get(Integer.valueOf(i)).intValue()).setAttribute(tempEA);
                        }
                    }
                }
                
                
                //We also need to change the values in the Leaf Edge attribs in the unrestricted case. 
            
                EdgeAttribute[] newT1LeafEdgeAtt = T1.getCopyLeafEdgeAttribs();
                EdgeAttribute[] newT2LeafEdgeAtt = T2.getCopyLeafEdgeAttribs();
            
                for (int i = 0; i < ol1; i++){
                    double[] tempVecEA = {newLeafEdgesvaluesT1[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newT1LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    if(OE2.getCompLeaves2orgLeaves(OE1.getOrgLeaves2compLeaves(i)) == -1){
                        newT2LeafEdgeAtt[OE1.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    }
                }
                for (int i = 0; i < ol2; i++){
                    double[] tempVecEA = {newLeafEdgesvaluesT2[i]};
                    EdgeAttribute tempEA = new EdgeAttribute(tempVecEA);
                    newT2LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    if(OE1.getCompLeaves2orgLeaves(OE2.getOrgLeaves2compLeaves(i)) == -1){
                        newT1LeafEdgeAtt[OE2.getOrgLeaves2compLeaves(i)].setEdgeAttribute(tempEA);
                    }
                }
            
            
                T1 = new PhyloTree(newEdgesT1, T1.getLeaf2NumMap(), newT1LeafEdgeAtt, false);
                T2 = new PhyloTree(newEdgesT2, T2.getLeaf2NumMap(), newT2LeafEdgeAtt, false);
            
                tempGeode = parPolyMain.getGeodesic(T1, T2);
                //System.out.println("The final tau was " + tau);
                
            }
            
        }
        
        //Getting the final values after the gradient descent has been performed. 
        Tree1 = T1;
        Tree2 = T2;
        Geodesic FinalGeode = tempGeode;
        Distance = FinalGeode.getDist();
        IterCount = iterCount;
        
        /*System.out.println("   Tree 1: \n" + treePrinter.toString(Tree1)+"\n \n");
        System.out.println("   Tree 2: \n" + treePrinter.toString(Tree2)+"\n \n");
        System.out.println(" With distance " + Distance);*/
        
    }// end of Constructor 2
    
    public OrthExtDistance(OrthExt OE1, OrthExt OE2, boolean restricted){
        O1ID = OE1.getOID();
        O2ID = OE2.getOID();
        if (restricted){
            Constructor1(OE1, OE2);
        } else {
            Constructor2(OE1, OE2);
        }
    }
    
    public OrthExtDistance(orthantExtPair orthEP, boolean restricted){
        O1ID = orthEP.getOrthE1().getOID();
        O2ID = orthEP.getOrthE2().getOID();
        if (restricted){
            Constructor1(orthEP.getOrthE1(), orthEP.getOrthE2());
        } else {
            Constructor2(orthEP.getOrthE1(), orthEP.getOrthE2());
        }
    }
    
    
    //Getters & Printers
    public PhyloTree getFirstTree(){
        return Tree1;
    }
    
    public PhyloTree getSecondTree(){
        return Tree2;
    }
    
    public double getDistance(){
        return Distance;
    }
    
    public int getIterCount(){
        return IterCount;
    }
    
    public int getO1ID(){
        return this.O1ID;
    }
    
    public int getO2ID(){
        return this.O2ID;
    }
    
    public String getWarning(){
        return this.Message;
    }
    
    public void PrintSummary(){
        System.out.println("The distance between the orthant extension spaces is " + Distance);
        PhyloNicePrinter treePrinter = new PhyloNicePrinter();
        System.out.println("Best Tree 1: ");
        System.out.println(treePrinter.toString(Tree1));
        System.out.println("");
        System.out.println("Best Tree 2: ");
        System.out.println(treePrinter.toString(Tree2));
        System.out.println("");
        System.out.println("Number of iterations for Computation: " + this.IterCount);
    }
}

