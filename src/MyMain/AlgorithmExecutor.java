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
    
    public static void ExtensionDistance(PhyloTree FirstTree, PhyloTree SecondTree, String fileName, int Tnum){
        
        PhyloNicePrinter nicePrint = new PhyloNicePrinter();
        
        //Initializing Files to Write results from the code
        
        try {
            File outReport = new File(fileName+"_Report.txt");
            File outTable = new File(fileName+"_Table.csv");
            outReport.createNewFile();
            outTable.createNewFile();
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
            ESdistance = new ExtensionSpaceDistance(firstES, secondES, false);//, 6);
            End = System.currentTimeMillis();
        } else {
            Start = System.currentTimeMillis();
            ESdistance = new ExtensionSpaceDistance(firstES, secondES, false, Tnum);
            End = System.currentTimeMillis();
        }
        
        long TimeMilliSeconds = (End - Start);
        
        try {
            FileWriter myWriter = new FileWriter(fileName+"_Report.txt");
            PrintWriter WriterTable = new PrintWriter(new FileWriter(fileName+"_Table.csv"));
            
            String header = "Orthant 1, Orthant 2, Distance, Iterations, Tree 1, Tree 2, Message"; 
            WriterTable.println(header); 
            
            //Computing some extra information
            int OrthPairsNum = ESdistance.getOOED().size();
            double MeanIts = 0;
            String ListWarnings = "";
            for (int j = 0; j < ESdistance.getOOED().size(); j++){
                OrthExtDistance OED = ESdistance.getOOED().get(j);
                MeanIts = MeanIts*((double)j/(j+1)) + (double)OED.getIterCount()/(j+1);
                if (OED.getWarning() != null){
                    ListWarnings = ListWarnings + "For orthant pair (" + OED.getO1ID() + ", " + OED.getO2ID() + ") " + OED.getWarning() + "\n";
                }
                String Line = OED.getO1ID() + "," + OED.getO2ID() + "," + OED.getDistance() + "," + OED.getIterCount() + "," + escapeSpecialCharacters(OED.getFirstTree().getNewick(true)) +","+ escapeSpecialCharacters(OED.getSecondTree().getNewick(true)) +"," + OED.getWarning();
                WriterTable.println(Line);
            }
             
            myWriter.write("The first tree is: \n"+ nicePrint.toString(FirstTree)+"\n");
            myWriter.write("The second tree is: \n"+ nicePrint.toString(SecondTree)+"\n");
            myWriter.write("The complete leaf set is: "+ completeLeafSet + "\n \n");
             
            myWriter.write("\n The Running Time is " + TimeMilliSeconds + " mili-seconds \n");
            myWriter.write("Which is " + formatDuration(TimeMilliSeconds) + " \n");
            myWriter.write("Number of orthant pairs is " + OrthPairsNum + "\n");
            myWriter.write("The mean number of iterations is " + MeanIts + "\n");
            myWriter.write("The best orthant pair is (" + ESdistance.getOOED().get(0).getO1ID() + ", " + ESdistance.getOOED().get(0).getO2ID() + ") \n");
            myWriter.write("The distance found is " + ESdistance.getDistance() + "\n");
             
            myWriter.write(ListWarnings);
             
            myWriter.close();
            WriterTable.close();
            System.out.println("Successfully wrote the report.");
        } catch (IOException e) {
            System.out.println("An error occurred writing on the document.");
            e.printStackTrace();
        }
        
    }
    
    public static void main(String[] args) {
        try {
            String InputDocName = args[0];
            String OutputDocName = "ExtensionDistance";
            try{
                File InputFile = new File(InputDocName);
                Scanner myReader = new Scanner(InputFile);
                int nThreads = 0;
                String dirName = "ExtDistResults";
                int CountTreePairs = 0;
                
                if (args.length > 1){
                    if ((args.length == 2) && (args[1].matches("-?\\d+"))){
                        nThreads = Integer.parseInt(args[1]);
                    } else {
                        OutputDocName = args[1];
                    }
                }
                if (args.length == 3){
                     if (args[2].matches("-?\\d+")){
                         nThreads = Integer.parseInt(args[2]);
                     } else {
                         dirName = args[2];
                     }
                } else if (args.length > 3){
                    if (args[2].matches("-?\\d+")){
                         nThreads = Integer.parseInt(args[2]);
                         dirName = args[3];
                     } else {
                         dirName = args[2];
                         if(!args[3].matches("-?\\d+")){
                             System.out.println("Warning: Invalid number of Threads. Fixed at 1 instead.");
                             nThreads = 0;
                         } else {
                             nThreads = Integer.parseInt(args[3]);
                         }
                         
                     }
                }
                
                if (nThreads < 0){
                    System.out.println("Warning: Invalid number of Threads. Fixed at 1 instead.");
                    nThreads = 0;
                }
                
                
                File theDir = new File(dirName);
                if (!theDir.exists()){
                    theDir.mkdirs();
                }
                
                
                while (myReader.hasNextLine()) {
                    CountTreePairs++;
                    String Tree1String = myReader.nextLine();
                    if (myReader.hasNextLine()){
                        String Tree2String = myReader.nextLine();
                        PhyloTree FirstTree = new PhyloTree(Tree1String, false);
                        PhyloTree SecondTree = new PhyloTree(Tree2String, false);
                        ExtensionDistance(FirstTree, SecondTree, dirName+"/"+OutputDocName+"Pair"+CountTreePairs, nThreads);
                        if (myReader.hasNextLine()){
                            String EmptyLine = myReader.nextLine();
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