#include "variables.h"     

void ReadInput()
{

 int I, J, K;
 FILE *read, *ic;
 int atoms=0;

 read = fopen("input","r");
 if(read==NULL){printf("error opening 'input' file to read\n"); exit(1);}
 ic = fopen("IC.xyz","w");
 if(ic==NULL){printf("error opening 'IC.xyz' file to write\n"); exit(1);}

 for(K=0;K<Nmax;K++){
     if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
          for(I=0;I<Nx;I++){
              for(J=0;J<Ny;J++){
                  fscanf(read,"%d ",&Box[I][J][K]); }
                  fscanf(read,"\n");   }
                  fscanf(read,"\n");   }}

 
// Calculating height of each site
 Height_Calculation();   
 
// Getting .xyz file for visualization purpose 
 for(K=0;K<Nmax;K++){
    if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
      for(I=0;I<Nx;I++){
         for(J=0;J<Ny;J++){
            if(Box[I][J][K]!=0) { atoms++; }}}}}
               fprintf(ic,"%d\n\n",atoms);
	       
 for(K=0; K<stp_bulk; K++){
     if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
          for(I=0; I<Nx; I++){
              for(J=0;J<Ny;J++){
                  if(Box[I][J][K]==1)      { fprintf(ic,"Ga %d %d %d\n",I,J,K); } 
                  else if(Box[I][J][K]==2) { fprintf(ic,"N %d %d %d\n",I,J,K); }
                  else if(Box[I][J][K]==3) { fprintf(ic,"AdGa %d %d %d\n",I,J,K); } 
                     }}}}
     
 for(K=stp_bulk;K<Nmax;K++){
     if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
          for(I=0; I<Nx; I++){
              for(J=0;J<Ny;J++){
                  if(Box[I][J][K]==1)      { fprintf(ic,"S %d %d %d\n",I,J,K); }
                  else if(Box[I][J][K]==2) { fprintf(ic,"O %d %d %d\n",I,J,K); } 
                  else if(Box[I][J][K]==3 ){ fprintf(ic,"C %d %d %d\n",I,J,K); }
                      }}}}
fclose(read);
fclose(ic);
}   // end of this program 
