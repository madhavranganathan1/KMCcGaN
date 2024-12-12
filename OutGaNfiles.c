#include "variables.h"
void OutGaNfiles()
{


 int I, J, K;
 FILE *output;
 char outfile[256];
 sprintf(outfile,"outputGaN%dx%d",outNfiles,rank);
// Write final output in some file//  // this format can be read by my ReadInput.c program
 output = fopen(outfile,"w");
 if(output==NULL){printf("error in opening outputGaN to write\n"); exit(1);}

 for(K=0; K<Nmax; K++){
     for(I=0; I<Nx; I++){
	 if(((K%8)==0) || (((K+5)%8)==0) || (((K+4)%8)==0) || (((K+1)%8)==0)){
              for(J=0; J<Ny; J++){
                 fprintf(output,"%d ",Box[I][J][K]); }
                 fprintf(output,"\n"); }
                 fprintf(output,"\n"); }}
fclose(output);
}   // end of this program
