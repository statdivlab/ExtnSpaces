/** Part of the package BHVExtMinDistance

This class describes matrices with integer entries, which will represent linear constraints on the Extension spaces
restricted to an orthant; with the intention to generate useful computations for those matrices. 
*/

package BHVExtMinDistance;

import java.text.DecimalFormat;
import java.util.*;
import distanceAlg1.*;

public class PhyloNicePrinter{
    
    //Null constructor
    public PhyloNicePrinter(){
        
    }
    
    //The only functions we care about
    public String toString(PhyloTree T){
        Vector<String> Leaf2Num =  T.getLeaf2NumMap();
        EdgeAttribute[] leafEA =  T.getLeafEdgeAttribs();
        
        Vector<PhyloTreeEdge> tEdges = T.getEdges();
        
        Vector<String> LeafString = new Vector<String>();
        for(int i = 0; i < Leaf2Num.size(); i++){
            LeafString.add("{"  + Leaf2Num.get(i) + "} " + leafEA[i].toString());
        }
        
        Vector<String> tEdgesString = new Vector<String>();
        for(PhyloTreeEdge e : tEdges){
            Bipartition eClone = e.getOriginalEdge().clone();
            eClone.complement(Leaf2Num.size());
            tEdgesString.add("{" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "|" + eClone.toStringVerbose(Leaf2Num) + "} "+e.getAttribute().toString());
            //if(e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
            //    Bipartition eClone = e.getOriginalEdge().clone();
            //    eClone.complement(Leaf2Num.size());
            //    tEdgesString.add("{" + eClone.toStringVerbose(Leaf2Num) + "} "+e.getAttribute().toString());
            //} else {
            //    tEdgesString.add("{" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "} "+e.getAttribute().toString());
            //}
            
        }
        
        return("  Leaves: " + LeafString + ";\n  Edges: " + tEdgesString + ";");
        
    }
    
    public String toString(Bipartition bip, Vector<String> Leaf2Num){
        Bipartition comp = bip.clone();
        comp.complement(Leaf2Num.size());
        return("{"+bip.toStringVerbose(Leaf2Num)+"|"+comp.toStringVerbose(Leaf2Num)+"}");
    }
    
    public String toString(PhyloTreeEdge e, Vector<String> Leaf2Num){
        Bipartition eClone = e.getOriginalEdge().clone();
        eClone.complement(Leaf2Num.size());
        return("{" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "|" + eClone.toStringVerbose(Leaf2Num) + "} "+e.getAttribute().toString());
    }
    
    public String toString(Vector<PhyloTreeEdge> commonEdges, Vector<String> Leaf2Num){
        DecimalFormat d8o = new DecimalFormat("#0.########");
        String Print = "The Common Edges are: \n";
        
        for (PhyloTreeEdge e : commonEdges) {
            Bipartition eClone = e.getOriginalEdge().clone();
            eClone.complement(Leaf2Num.size());
            Print = Print + "   {" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "|" + eClone.toStringVerbose(Leaf2Num) + "} " + e.getAttribute().toString() + "\n";
            /**if (e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
                Bipartition eClone = e.getOriginalEdge().clone();
                eClone.complement(Leaf2Num.size());
                Print = Print + "   {" + eClone.toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString() + "\n";
            } else {
                Print = Print + "   {" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "} " + e.getAttribute().toString() + "\n";
            }*/
		}
        
        return Print;
    }
    
    public String toString(RatioSequence RS, Vector<String> Leaf2Num){
        DecimalFormat d8o = new DecimalFormat("#0.########");
        String Print = "The ratio sequence is: \n";
        for (int i = 0; i < RS.size(); i++){
            Ratio rat = RS.getRatio(i);
            Print = Print + "   [";
            Vector<PhyloTreeEdge> rEEdges = rat.getEEdges();
            for (int j = 0; j < rEEdges.size(); j++){
                PhyloTreeEdge e = rEEdges.get(j);
                Bipartition eClone = e.getOriginalEdge().clone();
                eClone.complement(Leaf2Num.size());
                Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+ "|" + eClone.toStringVerbose(Leaf2Num) +"} " + e.getAttribute().toString();
                /**if (e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
                    Bipartition eClone = e.getOriginalEdge().clone();
                    eClone.complement(Leaf2Num.size());
                    Print = Print + "{" + eClone.toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                } else {
                    Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                }*/
                if (j < (rEEdges.size()-1)){
                    Print = Print + ", ";
                } else {
                    Print = Print + "] ";
                }
            }
            
             Print = Print + d8o.format(rat.getELength()) + "/" + d8o.format(rat.getFLength()) + " [";
            
            double ratVal = rat.getELength()/rat.getFLength();
            
            Vector<PhyloTreeEdge> rFEdges = rat.getFEdges();
            for (int j = 0; j < rFEdges.size(); j++){
                PhyloTreeEdge e = rFEdges.get(j);
                Bipartition eClone = e.getOriginalEdge().clone();
                eClone.complement(Leaf2Num.size());
                Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+ "|" + eClone.toStringVerbose(Leaf2Num) +"} " + e.getAttribute().toString();
                /**if (e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
                    Bipartition eClone = e.getOriginalEdge().clone();
                    eClone.complement(Leaf2Num.size());
                    Print = Print + "{" + eClone.toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                } else {
                    Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                }*/
                if (j < (rFEdges.size()-1)){
                    Print = Print + ", ";
                } else {
                    Print = Print + "] --> t = " + d8o.format(ratVal/(1+ratVal)) +" \n";
                }
            }
            
        }
        
        return Print;
    }
    
    public String toString(Geodesic Geo, Vector<String> Leaf2Num){
        DecimalFormat d8o = new DecimalFormat("#0.########");
        Vector<PhyloTreeEdge> commonEdges = Geo.getCommonEdges();
        RatioSequence RSeq = Geo.getRS();
        
        
        String Print = "Geodesic of size " + d8o.format(Geo.getDist()) + "\n";
        
        Print = Print + "\n The Common Edges are: \n";
        
        for (PhyloTreeEdge e : commonEdges) {
            Bipartition eClone = e.getOriginalEdge().clone();
            eClone.complement(Leaf2Num.size());
            Print = Print + "   {" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "|" + eClone.toStringVerbose(Leaf2Num) + "} " + e.getAttribute().toString() + "\n";
            /**if (e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
                Bipartition eClone = e.getOriginalEdge().clone();
                eClone.complement(Leaf2Num.size());
                Print = Print + "   {" + eClone.toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString() + "\n";
            } else {
                Print = Print + "   {" + e.getOriginalEdge().toStringVerbose(Leaf2Num) + "} " + e.getAttribute().toString() + "\n";
            }*/
		}
        
        Print = Print + "\n And the Ratio Sequence is: \n";
        
        for (int i = 0; i < RSeq.size(); i++){
            Ratio rat = RSeq.getRatio(i);
            Print = Print + "   [";
            Vector<PhyloTreeEdge> rEEdges = rat.getEEdges();
            for (int j = 0; j < rEEdges.size(); j++){
                PhyloTreeEdge e = rEEdges.get(j);
                Bipartition eClone = e.getOriginalEdge().clone();
                eClone.complement(Leaf2Num.size());
                Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+ "|" + eClone.toStringVerbose(Leaf2Num) +"} " + e.getAttribute().toString();
                /**if (e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
                    Bipartition eClone = e.getOriginalEdge().clone();
                    eClone.complement(Leaf2Num.size());
                    Print = Print + "{" + eClone.toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                } else {
                    Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                }*/
                if (j < (rEEdges.size()-1)){
                    Print = Print + ", ";
                } else {
                    Print = Print + "] ";
                }
            }
            
             Print = Print + d8o.format(rat.getELength()) + "/" + d8o.format(rat.getFLength()) + " [";
            
            double ratVal = rat.getELength()/rat.getFLength();
            
            Vector<PhyloTreeEdge> rFEdges = rat.getFEdges();
            for (int j = 0; j < rFEdges.size(); j++){
                PhyloTreeEdge e = rFEdges.get(j);
                Bipartition eClone = e.getOriginalEdge().clone();
                eClone.complement(Leaf2Num.size());
                Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+ "|" + eClone.toStringVerbose(Leaf2Num) +"} " + e.getAttribute().toString();
                /**if (e.getOriginalEdge().getPartition().cardinality() > (Leaf2Num.size()/2)){
                    Bipartition eClone = e.getOriginalEdge().clone();
                    eClone.complement(Leaf2Num.size());
                    Print = Print + "{" + eClone.toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                } else {
                    Print = Print + "{" + e.getOriginalEdge().toStringVerbose(Leaf2Num)+"} " + e.getAttribute().toString();
                }*/
                if (j < (rFEdges.size()-1)){
                    Print = Print + ", ";
                } else {
                    Print = Print + "] --> t = " + d8o.format(ratVal/(1+ratVal)) +" \n";
                }
            }
            
        }
        
        return Print;
    }
    

}