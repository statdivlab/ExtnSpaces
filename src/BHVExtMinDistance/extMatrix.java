/** Part of the package BHVExtMinDistance

This class describes matrices with integer entries, which will represent linear constraints on the Extension spaces
restricted to an orthant; with the intention to generate useful computations for those matrices. 
*/

package BHVExtMinDistance;

import java.util.*;

public class extMatrix{
    private int[][] mat; //Matrix entries
    private int nrow;//number of rows in the matrix
    private int ncol;// number of columns in the matrix
    
    
    //Constructor
    
    //We create a blank matrix with desired dimensions. All entries are set to zero.
    public extMatrix(int rows, int columns){
        mat = new int[rows][columns];

        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                mat[i][j] = 0;
            }
        }
        nrow = rows;
        ncol = columns;
    }
    
    //Constructor clone
    public extMatrix(extMatrix cMat){
        this.nrow = cMat.getNRow();
        this.ncol = cMat.getNCol();
        
        this.mat = new int[this.nrow][this.ncol];
        for(int i = 0; i < this.nrow; i++){
            for(int j = 0; j < this.ncol; j++){
                this.mat[i][j] = cMat.element(i,j); 
            }
        }
    }
    
    //Accessors.
    
    //get the number of rows
    public int getNRow(){return nrow;}
    
    //get the number of columns
    public int getNCol(){return ncol;}
    
    //get the entry of the matrix in certain position.
    public int element(int row, int column){return this.mat[row][column];}
    
    //this function returns a whole row from the Matrix.
    public int[] getRow(int i){return this.mat[i];}
    
    public void PrintMat(){
        for (int i = 0; i<nrow; i++){
            for (int j=0; j<ncol; j++){
                System.out.print(" " + this.mat[i][j]);
            }
            System.out.println("");
        }
    }
    
    //Setters 
    
    //Sets a particular entry equal to value in the matrix.
    public void setItem(int i, int j, int value){
        this.mat[i][j] = value;
    }
    
    //Adds 1 to the value in a particular entry.
    public void plus1(int i, int j){
         this.mat[i][j]++;
    }
    
    //Useful computations

}