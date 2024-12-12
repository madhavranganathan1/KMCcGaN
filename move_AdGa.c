// Jan 21, 2022 : Subroutine for diffusion of AdGa atom. AdGa atom can diffuse over N atom (having atleast one Ga atom as neighbor) and form epitaxil Ga at H= H+5.  But if Ga neighbor is not found for diffusing AdGa site, then AdGa will get deposit ove N at H= H-4( w.r.t N diffusing site).
// Jan 21, 2022 : AdGa atom can alo diffuse over Ga atoms (atleast having one N atom as a direct bond) with Ga/AdGa as NN leading to the formation of AdGa atoms. Probabilty is calcuated based on the number of direct bonds to which Ga atom is attached.


#include"variables.h"

int move_AdGa(int X, int Y)            
{ 
  int Xp1 = (X+1)%Nx,   Xm1 = (X-1+Nx)%Nx,  Yp1 = (Y+1)%Ny,  Ym1 = (Y-1+Ny)%Ny ;	
    
      cntGa = 0, Prob = 0.0, M_Z=0; nhxn = 0; nhxp = 0;
      
  /*  if((X+1)>=Nx){ nhxp= ch_stp;} else if((X-1)<0)  { nhxn= -ch_stp;}This if..else condition has to be uncommented for vicinal surface*/ 

      if(((H[X][Y]+1)%8) ==0){
                              if((Box[Xp1][Y]  [H[X][Y]-nhxp+1]) == 1)      { cntGa++; }
                              if((Box[Xm1][Ym1][H[X][Y]-nhxn+1]) == 1)      { cntGa++; }
                              if((Box[Xm1][Yp1][H[X][Y]-nhxn+1]) == 1)      { cntGa++; }}

      else if(((H[X][Y]+1)% 8)!=0){
                                   if((Box[Xm1][Y]  [H[X][Y]-nhxn+1]) ==1)     { cntGa++; }
                                   if((Box[Xp1][Ym1][H[X][Y]-nhxp+1]) ==1)     { cntGa++; }
                                   if((Box[Xp1][Yp1][H[X][Y]-nhxp+1]) ==1)     { cntGa++; }}

      if(cntGa==0) { /* Prob = 1.0 */ if(((H[X][Y]+1) % 8) ==0) { H7_AdGa( X, Y) ; }  else if (((H[X][Y]+1) % 8)!=0 ) { H3_AdGa( X, Y) ; } }     
	 			      	
      else {
            if     (cntGa==1) { Prob = Prob_AdGa1Ga; }	      
	    else if(cntGa==2) { Prob = Prob_AdGa2Ga; }
            else if(cntGa==3) { Prob = Prob_AdGa3Ga; } 
           
            double Racc = gsl_rng_uniform_pos(racc);    
	    if(Racc < Prob) {  
		             if(((H[X][Y]+1) % 8) ==0)        { H7_AdGa( X, Y) ; } 
 			     else if (((H[X][Y]+1) % 8) !=0 ) { H3_AdGa( X, Y) ; } }

	    else { /* failed diffusion */  NoAdGadiff++; }
           }
 return 0;
}                          
