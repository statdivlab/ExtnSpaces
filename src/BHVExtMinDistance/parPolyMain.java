/** This file is part of GTP, a program for computing the geodesic distance between phylogenetic trees,
 * and sturmMean, a program for computing the Frechet mean between phylogenetic trees.
    Copyright (C) 2008-2012  Megan Owen, Scott Provan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>. */

package BHVExtMinDistance;

import java.util.Vector; 
import java.util.Stack;

import distanceAlg1.EdgeAttribute;
import distanceAlg1.Geodesic;
import distanceAlg1.PhyloTree;
import distanceAlg1.PhyloTreeEdge;
import distanceAlg1.RatioSequence;
import distanceAlg1.Ratio;

import polyAlg.Tools;
import polyAlg.BipartiteGraph;

//import org.biojavax.bio.phylo.io.nexus.*;
//import org.biojava.bio.seq.io.ParseException;

public class parPolyMain {
    
/** Returns the geodesic between t1 and t2, which are assumed to have no common edges.
 *  Does not assume t1 and t2 have the same number of edges.
 *  Does not take into account the leaf edges.
 *  Uses polynomial algorithm.
 *  XXX: how to deal with multifurcating trees
 * 
 *  Returns:  a Geodesic with just the ratio sequence set 
 */
public static Geodesic getGeodesicNoCommonEdges(PhyloTree t1, PhyloTree t2 ) {
	int numEdges1 = t1.getEdges().size(); // number of edges in tree 1
	int numEdges2 = t2.getEdges().size(); // number of edges in tree 2
	RatioSequence rs = new RatioSequence();
	int [] aVertices, bVertices;
	Vector<Ratio> queue = new Vector<Ratio>();
	Ratio ratio;
	int[][] cover;
	
	if (numEdges1 == 0 && numEdges2 == 0) {
		return new Geodesic(new RatioSequence());
	}
		
	// double-check no common edges
	Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1, t2);
	if (commonEdges.size() != 0) {
		System.out.println("Exiting: tried to compute geodesic between subtrees that should not have common edges, but do!  t1 = " + t1 + " and t2 = " + t2);
		System.exit(1);
	}
	
	// double-check that both trees have splits.  Otherwise didn't remove a common edge.
	if (numEdges1 ==0 || numEdges2 == 0) {
		System.out.println("Exiting: tried to compute geodesic between subtrees that should not have common/compatible edges, but do!  t1 = " + t1 + " and t2 = " + t2);
		System.exit(1);
	}
	
	// if we can't split the ratio because it has too few edges in either the numerator or denominator
	if ((numEdges1 ==1) || (numEdges2 ==1)) {
		rs.add( new Ratio(t1.getEdges(), t2.getEdges()) );
		return new Geodesic(rs);
	}
	
	// initialize BipartiteGraph
	boolean[][] incidenceMatrix = Tools.getIncidenceMatrix(t1.getEdges(), t2.getEdges());

	
	BipartiteGraph bg = new BipartiteGraph(incidenceMatrix, t1.getIntEdgeAttribNorms(), t2.getIntEdgeAttribNorms());
	
	queue.add(new Ratio(t1.getEdges(), t2.getEdges() ));
	
	while(queue.size() >0) {
		ratio = queue.remove(0);
		
		aVertices = new int[ratio.getEEdges().size()];
		bVertices = new int[ratio.getFEdges().size()];
		
		// convert the ratio to what we pass to vertex cover
		for (int i = 0; i < ratio.getEEdges().size(); i++) {
			aVertices[i] = t1.getEdges().indexOf(ratio.getEEdges().get(i));
		}
		for (int i = 0; i < ratio.getFEdges().size(); i++) {
			bVertices[i] = t2.getEdges().indexOf(ratio.getFEdges().get(i));
		}
		
		// get the cover
		cover = bg.vertex_cover(aVertices, bVertices);
		
		// check if cover is trivial
		if ( (cover[0][0] == 0) || (cover[0][0] == aVertices.length) ){
			// add ratio to geodesic
			rs.add(ratio);

			
		}
		else {  // cover not trivial
			// make two new ratios
			Ratio r1 = new Ratio();
			Ratio r2 = new Ratio();
			
			int j = 0;  // for index in cover array
			
			// split the ratio based on the cover
			for (int i = 0; i < aVertices.length; i++) {
				if ( (j < cover[2].length) && (aVertices[i] == cover[2][j]) ) {
					r1.addEEdge( t1.getEdge(aVertices[i]) );
					j++;
				}
				else { // the split is not in the cover, and hence dropped first
					r2.addEEdge( t1.getEdge(aVertices[i]) );
				}
			}
			
			j = 0;   // reset index
			// split the ratio based on the cover
			for (int i = 0; i < bVertices.length; i++) {
				if ( (j < cover[3].length) && (bVertices[i] == cover[3][j]) ) {	
					r2.addFEdge( t2.getEdge(bVertices[i]) );
					j++;
				}
				else { // the split is not in the cover, and hence dropped first
					r1.addFEdge( t2.getEdge(bVertices[i]) );
				}
			}
			
			// add ratios to the queue
			queue.add(0, r2);
			queue.add(0, r1);
		}
	}
	
	return new Geodesic(rs);
}

	
/** Returns the distance between t1 and t2, accounting for any common edges and leaf edges.
 *  Calls recursive getGeodesic
 *  Does not assume t1 and t2 have the same number of edges.
 *  Pass in null for geoFile to not write to a file.
 * 
 */   
public static Geodesic getGeodesic(PhyloTree t1, PhyloTree t2) {
	double leafContributionSquared = 0;
	EdgeAttribute [] t1LeafEdgeAttribs = t1.getLeafEdgeAttribs();
	EdgeAttribute [] t2LeafEdgeAttribs = t2.getLeafEdgeAttribs();
	Geodesic geo = new Geodesic(new RatioSequence(),t1.getLeafEdgeAttribs(),t2.getLeafEdgeAttribs());
	
	// get the leaf contributions
	for(int i = 0; i < t1.getLeaf2NumMap().size(); i++ ) {
		if ( !(t1.getLeaf2NumMap().get(i).equals(t2.getLeaf2NumMap().get(i)) ) ) {
			System.out.println("Error getting geodesic: trees do not have the same sets of leaves");
			System.out.println("Starting tree leaves: " + t1.getLeaf2NumMap());
			System.out.println("Target tree leaves: " + t2.getLeaf2NumMap());
			
			System.out.println("Starting tree: " + t1.getNewick(true));
			System.out.println("Target tree: " + t2.getNewick(true));
			
			System.exit(1);
		}
//		System.out.println("leaf: " + t1.getLeaf2NumMap().get(i) + " | " + t1LeafEdgeLengths[i] + " - " + t2LeafEdgeLengths[i] + "| = " + (t1LeafEdgeLengths[i] - t2LeafEdgeLengths[i]) );

//		leafContributionSquared = leafContributionSquared+ Math.pow(Math.abs(t1LeafEdgeLengths[i] - t2LeafEdgeLengths[i]), 2);
		leafContributionSquared = leafContributionSquared+ Math.pow(EdgeAttribute.difference(t1LeafEdgeAttribs[i],t2LeafEdgeAttribs[i]).norm(), 2);
	}
	geo.setLeafContributionSquared(leafContributionSquared);
	
	Vector<PhyloTree> aTreesNoCommonEdges = new Vector<PhyloTree>();
	Vector<PhyloTree> bTreesNoCommonEdges = new Vector<PhyloTree>();
	
	// get the pairs of trees with no common edges put into aTreesNoCommonEdges and bTreesNoCommonEdges
	//  aTreesNoCommonEdges.get(i) goes with bTreesNoCommonEdges.get(i)
	//splitOnCommonEdge(t1,t2);
	
    //WARNING: This is intended as a replacement to the splitCommonEdge function, but not recursive anymore.
    {
        Stack<PhyloTree> stack1 = new Stack<PhyloTree>();
        Stack<PhyloTree> stack2 = new Stack<PhyloTree>();
        
        stack1.push(t1);
        stack2.push(t2);
        
        while(!stack1.empty()){
            PhyloTree currentTree1 = stack1.pop();
            PhyloTree currentTree2 = stack2.pop();
            
            Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(currentTree1, currentTree2);
            
            if(commonEdges.size() == 0){
                aTreesNoCommonEdges.add(currentTree1);
                bTreesNoCommonEdges.add(currentTree2);
            }else{
                PhyloTreeEdge commonEdge = commonEdges.get(0);
                
                Vector<PhyloTreeEdge> edgesA1 = new Vector<PhyloTreeEdge>();
                Vector<PhyloTreeEdge> edgesA2 = new Vector<PhyloTreeEdge>();
                Vector<PhyloTreeEdge> edgesB1 = new Vector<PhyloTreeEdge>();
                Vector<PhyloTreeEdge> edgesB2 = new Vector<PhyloTreeEdge>();
                
                for (PhyloTreeEdge e:  currentTree1.getEdges()) {
                    // tree A is the tree under the common edge (i.e. the common edge is the root)
                    // tree B is the tree above the common edge (i.e. think of all leaves below the common edge as one big leaf)
                    if (commonEdge.properlyContains(e)) {
                        edgesA1.add(e.clone());
                    } else if (!e.sameBipartition(commonEdge)) {
                        edgesB1.add(e.clone());
                    }
                }
                
                for (PhyloTreeEdge e : currentTree2.getEdges()) {
                    // tree A is the tree under the common edge (i.e. the common edge is the root)
                    // tree B is the tree above the common edge (i.e. think of all leaves below the common edge as one big leaf)
                    if (commonEdge.properlyContains(e)) {
                        edgesA2.add(e.clone());
                    }
                    else if (!e.sameBipartition(commonEdge)) {
                        edgesB2.add(e.clone());
                    }
                }
                
                // make the 4 trees
                PhyloTree tA1 = new PhyloTree(edgesA1, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
                PhyloTree tB1 = new PhyloTree(edgesB1, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
                PhyloTree tA2 = new PhyloTree(edgesA2, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
                PhyloTree tB2 = new PhyloTree(edgesB2, Tools.myVectorCloneString(t1.getLeaf2NumMap()), t1.isRooted());
                
                stack1.push(tA1);
                stack2.push(tA2);
                stack1.push(tB1);
                stack2.push(tB2);
            }
        }
    }
    
    
	//set the common edges
	Vector<PhyloTreeEdge> commonEdges = PhyloTree.getCommonEdges(t1,t2);
	geo.setCommonEdges(commonEdges);
	
	// set the t1 and t2 attributes of the common edges
	Vector<PhyloTreeEdge> eCommonEdges = Tools.myVectorClonePhyloTreeEdge(commonEdges);
	Vector<PhyloTreeEdge> fCommonEdges = Tools.myVectorClonePhyloTreeEdge(commonEdges);
	for (int i = 0; i < eCommonEdges.size(); i++) {
		eCommonEdges.get(i).setAttribute(t1.getAttribOfSplit(eCommonEdges.get(i)));
		fCommonEdges.get(i).setAttribute(t2.getAttribOfSplit(fCommonEdges.get(i)));
	}
	geo.seteCommonEdges(eCommonEdges);
	geo.setfCommonEdges(fCommonEdges);
	
		
	// find the geodesic between each pair of subtrees found by removing the common edges
	for(int i = 0; i < aTreesNoCommonEdges.size(); i++) {
		PhyloTree subTreeA = aTreesNoCommonEdges.get(i);
		PhyloTree subTreeB = bTreesNoCommonEdges.get(i);
		
		Geodesic newGeo = getGeodesicNoCommonEdges(subTreeA, subTreeB);
		
		geo.setRS(RatioSequence.interleave(geo.getRS(), newGeo.getRS()));
	}

	return geo;
}


public static double calcGeoDist(PhyloTree t1, PhyloTree t2) {
	return getGeodesic(t1, t2).getDist();
}

}

