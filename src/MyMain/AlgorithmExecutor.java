package MyMain;

import java.util.*;
import distanceAlg1.*;
import BHVExtMinDistance.*;
import polyAlg.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.FileNotFoundException;  // Import this class to handle errors
import java.util.Scanner; // Import the Scanner class to read text files
import java.util.Vector;
import java.util.Collections;
import java.io.PrintWriter;


import static polyAlg.PolyMain.getGeodesic;


public class AlgorithmExecutor{
    public static String formatDuration(long durationMillis) {
        long seconds = durationMillis / 1000;
        long minutes = seconds / 60;
        long hours = minutes / 60;
        long days = hours / 24;

        String formattedDuration = String.format("%d:%02d:%02d:%02d:%03d",
                days, hours % 24, minutes % 60, seconds % 60, durationMillis % 1000);

        return formattedDuration;
    }
    
    public static String escapeSpecialCharacters(String data){
        if (data == null) {
            throw new IllegalArgumentException("Input data cannot be null");
        }
        String escapedData = data.replaceAll("\\R", " ");
        if (data.contains(",") || data.contains("\"") || data.contains("'")) {
            data = data.replace("\"", "\"\"");
            escapedData = "\"" + data + "\"";
        }
        return escapedData;
    }
    
    private static int nextIndex(String t, int i, String s) {
		int minIndex = t.length();
		int tempIndex = -1;
		
		for (int j = 0; j < s.length(); j++) {
			tempIndex = t.substring(i+1).indexOf(s.charAt(j));
			if ((tempIndex != -1)  && (tempIndex + i + 1 < minIndex)) {
				minIndex = tempIndex + i +1;
			}
		}
		return minIndex;
	}
    
     public static boolean unrootedTree(String tNewick){
         String t = tNewick;
			
         //pull off ';' if at end
         if (t.charAt(t.length()-1) == ';') {
             t = t.substring(0, t.length() -1);
         }
		
         // pull off the first and last brackets (and a root length, between the last bracket and ;, if there is one.
         t = t.substring(t.indexOf('(') + 1);
         t = t.substring(0, t.lastIndexOf(')') );
         
         int i = 0;
         int OuterGroups = 0;
         boolean ActiveGroup = false;
         int innerParent = 0;
         
         Stack<Integer> stackBifurcations = new Stack<Integer>();
         int numBifurcations = 0;
         
         while (i < t.length() && i > -1){
             switch(t.charAt(i)) {
                 case '(':
                     if(!ActiveGroup){
                         ActiveGroup = true;
                     } else {
                         innerParent++;
                     }
                     stackBifurcations.push(0);
                     i++;
                     break;
                     
                 case ')':
                     if (ActiveGroup){
                         if (innerParent > 0){
                             innerParent--;
                         } else {
                             ActiveGroup = false;
                             OuterGroups++;
                         }
                     } else {
                         System.err.println("Reading error: invalid Newick string in tree " + t );
                         System.exit(1);
                     }
                     numBifurcations = stackBifurcations.pop();
                     if (numBifurcations == 0){
                         System.err.println("Reading error: node of degree two in Newick format for tree " + t );
                         System.exit(1);
                     }
                     i = nextIndex(t, i, ",)"); 
                     break;
                     
                 case ',': 
                     if(ActiveGroup){
                         numBifurcations = stackBifurcations.pop();
                         numBifurcations++;
                         stackBifurcations.push(numBifurcations);   
                     }
                     i++;
                     break;
                
                 default:
                     if (!ActiveGroup) {
                         OuterGroups++;
                     }
                     i = nextIndex(t, i, ",)");
                     break;
                     
             }
         }
         
         if (OuterGroups <= 2){
             return(false);
         } else {
             return(true);
         } 
	}
    
    public static void ExtensionDistance(PhyloTree FirstTree, PhyloTree SecondTree, String fileName, int Tnum, boolean TableProd, int nLinesTable, double Tole1, double Tole2){
        
        PhyloNicePrinter nicePrint = new PhyloNicePrinter();
        
        //Initializing Files to Write results from the code
        
        try {
            File outReport = new File(fileName+"_Report.txt");
            outReport.createNewFile();
            
            if(TableProd){
                File outTable = new File(fileName+"_Table.csv");
                outTable.createNewFile();
            }
            
        } catch (IOException e) {
            System.out.println("An error occurred opening the output documents.");
            e.printStackTrace();
        }
        
        Set<String> tempLeafSet = new HashSet<>();
        tempLeafSet.addAll(FirstTree.getLeaf2NumMap());
        tempLeafSet.addAll(SecondTree.getLeaf2NumMap());
        Vector<String> completeLeafSet = new Vector<String>(tempLeafSet);
        Collections.sort(completeLeafSet);
        
        ExtensionSpace firstES = new ExtensionSpace(FirstTree, completeLeafSet, false);
        ExtensionSpace secondES = new ExtensionSpace(SecondTree, completeLeafSet, false);
        
        long Start = 0;
        long End = 0;
        
        ExtensionSpaceDistance ESdistance;
        
        if ((Tnum == 0) || (Tnum == 1)){
            Start = System.currentTimeMillis();
            ESdistance = new ExtensionSpaceDistance(firstES, secondES, false, Tole1, Tole2);//, 6);
            End = System.currentTimeMillis();
        } else {
            Start = System.currentTimeMillis();
            ESdistance = new ExtensionSpaceDistance(firstES, secondES, false, Tnum, Tole1, Tole2);
            End = System.currentTimeMillis();
        }
        
        long TimeMilliSeconds = (End - Start);
        
        try {
            if (TableProd){
                FileWriter myWriter = new FileWriter(fileName+"_Report.txt");
                PrintWriter WriterTable = new PrintWriter(new FileWriter(fileName+"_Table.csv"));

                String header = "Orthant 1, Orthant 2, Distance, Iterations, Tree 1, Tree 2, Message"; 
                WriterTable.println(header); 

                //Computing some extra information
                int OrthPairsNum = ESdistance.getOOED().size();
                double MeanIts = 0;
                int NumberOptimals = 1;
                int nLT = 0;
                if (nLinesTable == 0){
                    nLT = OrthPairsNum;
                } else {
                    nLT = Math.min(nLinesTable, OrthPairsNum);
                }
                while((NumberOptimals < OrthPairsNum) && (ESdistance.getDistance() == ESdistance.getOOED().get(NumberOptimals).getDistance())){
                    NumberOptimals++;
                }
                String ListWarnings = "";
                for (int j = 0; j < ESdistance.getOOED().size(); j++){
                    OrthExtDistance OED = ESdistance.getOOED().get(j);
                    MeanIts = MeanIts*((double)j/(j+1)) + (double)OED.getIterCount()/(j+1);
                    if (OED.getWarning() != null){
                        ListWarnings = ListWarnings + "For orthant pair (" + OED.getO1ID() + ", " + OED.getO2ID() + ") " + OED.getWarning() + "\n";
                    }
                    if (j < nLT){
                        String Line = OED.getO1ID() + "," + OED.getO2ID() + "," + OED.getDistance() + "," + OED.getIterCount() + "," + escapeSpecialCharacters(OED.getFirstTree().getNewick(true)) +","+ escapeSpecialCharacters(OED.getSecondTree().getNewick(true)) +"," + OED.getWarning();

                        WriterTable.println(Line);
                    }
                }

                myWriter.write("The first original tree is: \n"+ nicePrint.toString(FirstTree)+"\n");
                myWriter.write("The second original tree is: \n"+ nicePrint.toString(SecondTree)+"\n");
                myWriter.write("The complete leaf set is: "+ completeLeafSet + "\n \n");

                myWriter.write("\n The Running Time is " + TimeMilliSeconds + " mili-seconds, which is: "+ formatDuration(TimeMilliSeconds) +" \n");
                myWriter.write("The number of orthant pairs is " + OrthPairsNum + "\n");
                myWriter.write("The mean number of iterations is " + MeanIts + "\n");
                myWriter.write("The distance found is " + ESdistance.getDistance() + "\n");
                myWriter.write("A total of " + NumberOptimals + " tree pairs achieved this distance. \n This are: \n \n");
                for (int j = 0; j < NumberOptimals; j++){
                    myWriter.write("\n ---------- \n");
                    myWriter.write("For orthant pair (" + ESdistance.getOOED().get(j).getO1ID() + ", " + ESdistance.getOOED().get(j).getO2ID() + ") the trees are \n");
                    
                    myWriter.write("\n T'_1:"+nicePrint.toString(ESdistance.getOOED().get(j).getFirstTree())+"\n");
                    
                    myWriter.write("\n T'_2:"+nicePrint.toString(ESdistance.getOOED().get(j).getSecondTree())+"\n \n");
                }
                
                myWriter.write("\n ---------- \n");
                myWriter.write("\n \n WARNINGS PRODUCED BY CODE: \n");
                myWriter.write(ListWarnings);

                myWriter.close();
                WriterTable.close();
            } else {
                FileWriter myWriter = new FileWriter(fileName+"_Report.txt");

                //Computing some extra information
                int OrthPairsNum = ESdistance.getOOED().size();
                double MeanIts = 0;
                int NumberOptimals = 1;
                while((NumberOptimals < OrthPairsNum) && (ESdistance.getDistance() == ESdistance.getOOED().get(NumberOptimals).getDistance())){
                    NumberOptimals++;
                }
                String ListWarnings = "";
                for (int j = 0; j < ESdistance.getOOED().size(); j++){
                    OrthExtDistance OED = ESdistance.getOOED().get(j);
                    MeanIts = MeanIts*((double)j/(j+1)) + (double)OED.getIterCount()/(j+1);
                    if (OED.getWarning() != null){
                        ListWarnings = ListWarnings + "For orthant pair (" + OED.getO1ID() + ", " + OED.getO2ID() + ") " + OED.getWarning() + "\n";
                    }
                }

                myWriter.write("The first original tree is: \n"+ nicePrint.toString(FirstTree)+"\n");
                myWriter.write("The second original tree is: \n"+ nicePrint.toString(SecondTree)+"\n");
                myWriter.write("The complete leaf set is: "+ completeLeafSet + "\n \n");

                myWriter.write("\n The Running Time is " + TimeMilliSeconds + " mili-seconds, which is: "+ formatDuration(TimeMilliSeconds) +" \n");
                myWriter.write("The number of orthant pairs is " + OrthPairsNum + "\n");
                myWriter.write("The mean number of iterations is " + MeanIts + "\n");
                myWriter.write("The distance found is " + ESdistance.getDistance() + "\n");
                myWriter.write("A total of " + NumberOptimals + " tree pairs achieved this distance. \n This are: \n \n");
                for (int j = 0; j < NumberOptimals; j++){
                    myWriter.write("\n ---------- \n");
                    myWriter.write("For orthant pair (" + ESdistance.getOOED().get(j).getO1ID() + ", " + ESdistance.getOOED().get(j).getO2ID() + ") the trees are \n");
                    
                    myWriter.write("\n T'_1:"+nicePrint.toString(ESdistance.getOOED().get(j).getFirstTree())+"\n");
                    
                    myWriter.write("\n T'_2:"+nicePrint.toString(ESdistance.getOOED().get(j).getSecondTree())+"\n \n ");
                   
                }
                
                myWriter.write("\n ---------- \n");
                myWriter.write("\n \n WARNINGS PRODUCED BY CODE: \n");
                myWriter.write(ListWarnings);

                myWriter.close();
            }
            System.out.println("Successfully wrote the report.");
        } catch (IOException e) {
            System.out.println("An error occurred writing on the documents.");
            e.printStackTrace();
        }
        
    }
    
    public static void main(String[] args) {
        try {
            String InputDocName = args[0];
            try{
                File InputFile = new File(InputDocName);
                Scanner myReader = new Scanner(InputFile);
                int nThreads = 0;
                String dirName = "ExtDistResults";
                String OutputDocName = "ExtensionDistance";
                int CountTreePairs = 0;
                double Tole1 = 0.00000001;
                double Tole2 = 0.0000000000000001;
                
                boolean ProduceTable = false;
                int NumberLinesTable = 0;
                
                for (int j = 1; j < args.length; j++){
                    if ("-d".equals(args[j])){
                        if (args.length > j+1){
                            dirName = args[j+1];
                        }
                    }
                    if ("-o".equals(args[j])){
                        if (args.length > j+1){
                            OutputDocName = args[j+1];
                        }
                    }
                    if ("-k".equals(args[j])){
                        if (args.length > j+1){
                            if (args[j+1].matches("\\d+")){
                                nThreads = Integer.parseInt(args[j+1]);
                            } else {
                                System.out.println("Warning: Invalid number of Threads. Fixed at 1 instead.");
                            }
                        }
                    }
                    if ("-Tol1".equals(args[j])){
                        if (args.length > j+1){
                            try{
                                Tole1 = Double.parseDouble(args[j+1]);
                            } catch (NumberFormatException e) {
                                System.out.println("Warning: Invalid number for Tol1. Fixed at 1E-8.");
                            }
                        }
                    }
                    if ("-Tol2".equals(args[j])){
                        if (args.length > j+1){
                            try{
                                Tole2 = Double.parseDouble(args[j+1]);
                            } catch (NumberFormatException e) {
                                System.out.println("Warning: Invalid number for Tol2. Fixed at 1E-16.");
                            }
                        }
                    }
                    if ("-T".equals(args[j])){
                        ProduceTable = true;
                    }
                    
                    if ("-l".equals(args[j])){
                        if (args.length > j+1){
                            if (args[j+1].matches("\\d+")){
                                NumberLinesTable = Integer.parseInt(args[j+1]);
                            } else {
                                System.out.println("Warning: Invalid number of rows for table. The report table will contain all orthant pairs instead.");
                            }
                        }
                    }
                }
                
                File theDir = new File(dirName);
                if (!theDir.exists()){
                    theDir.mkdirs();
                }
                
                
                while (myReader.hasNextLine()) {
                    CountTreePairs++;
                    boolean TreesAreOk = true;
                    String PairName = "Pair"+CountTreePairs;
                    String Tree1String = myReader.nextLine();
                    if (myReader.hasNextLine()){
                        try{
                            if(!unrootedTree(Tree1String)){
                                TreesAreOk = false;
                            }
                            PhyloTree FirstTree = new PhyloTree(Tree1String, false);
                            
                        } catch(Exception e){
                            if (!Tree1String.trim().isEmpty()){
                                PairName = new String(Tree1String.trim());
                            }
                            Tree1String = myReader.nextLine();
                        }
                        if (myReader.hasNextLine()){
                            String Tree2String = myReader.nextLine();
                            PhyloTree FirstTree = new PhyloTree(Tree1String, false);
                            if(!unrootedTree(Tree2String)){
                                TreesAreOk = false;
                            }
                            PhyloTree SecondTree = new PhyloTree(Tree2String, false);
                            if (!TreesAreOk){
                                System.out.println("Warning: The Newick format in the tree pair number " + CountTreePairs + " is not appropriate for unrooted trees (it might be describing a rooted tree instead). This pair will be skiped by the code.");
                                
                            } else {
                                ExtensionDistance(FirstTree, SecondTree, dirName+"/"+OutputDocName+ "_" +PairName, nThreads, ProduceTable, NumberLinesTable, Tole1, Tole2);
                                if (myReader.hasNextLine()){
                                    String EmptyLine = myReader.nextLine();
                                }
                            }
                            
                        } else {
                            System.out.println("Warning: The last tree does not have a second tree for comparison. It was ignored.");
                        }
                        
                    } else {
                        System.out.println("Warning: The last tree does not have a second tree for comparison. It was ignored.");
                    }
                }
                myReader.close();
            } catch (FileNotFoundException e) {
                System.out.println("Input File not found");
                e.printStackTrace();
            }
        }
        catch (ArrayIndexOutOfBoundsException e){
            System.out.println("Error: Argument for Input File missing");
        }
    }
}