// Jan 21, 2022 : Subroutine for diffusion of Ga atom. Ga atom can diffuse over N atom (having atleast one Ga atom as neighbor) and form epitaxil Ga at H= H+5. 
// Jan 21, 2022 : Ga atom can alo diffuse over Ga atoms (atleast having one N atom as a direct bond) with Ga/AdGa as NN leading to the formation of AdGa atoms. Probabilty is calcuated based on the number of direct bonds to which Ga atom is attached. 
// Jan 21, 2022 : Formation of AdGa will take place only for NNN and NNNN neighbor as dffusion to  NN neighbor will give rise N atom and that is not a stable configuration right now in code.
// Jan 21, 2022 : After Ga has diffused there will be a check for its NN neighbor (w.r.t N atom (X,Y)) for AdGa and isolated N atom. If AdGa found then shifting of AdGa at H= H+5 while if any isolated N ad atom having Ga at H= H+5 found, then Ga atom at H = H+5 will fall over H =H+5 w.r.t N(X,Y).
/// Jan 21, 2022: This Ga atom move can disturb the newly formed AdGa configuration at NNN and so again a check will be formed for Neighbor N atoms (w.r.t N(X,Y))  and if AdGa present as NN  neighbor then AdGa shifts over N atom.
// Jan 21, 2022 : Subsurfce diffusion of N will be attempted if Ga atom moves fails due to no site avalaible and Probabilty too high.

#include"variables.h"

int move_Ga (int X,int Y)           
{
 int Xp1 = (X+1)%Nx,  Xm1 = (X-1+Nx)%Nx,  Yp1 = (Y+1)%Ny,  Ym1 = (Y-1+Ny)%Ny; 
     
     cntN = 0;   cntAdGa = 0;   Prob = 0.0;  nhxn = 0; nhxp = 0;  
  
  /* if((X+1)>=Nx){ nhxp= ch_stp;} else if((X-1)<0)  { nhxn= -ch_stp;} //This if..else condition has to be uncommented for vicinal surface*/ 
      
     if((H[X][Y] % 8) == 0){                           
                            if((Box[Xm1][Y]  [H[X][Y]-nhxn-1]) ==2) { cntN++; }
                            if((Box[Xp1][Ym1][H[X][Y]-nhxp-1]) ==2) { cntN++; }
                            if((Box[Xp1][Yp1][H[X][Y]-nhxp-1]) ==2) { cntN++; }
                            if((Box[Xm1][Y]  [H[X][Y]-nhxn-1]) ==3) { cntAdGa++; }
                            if((Box[Xp1][Ym1][H[X][Y]-nhxp-1]) ==3) { cntAdGa++; }
                            if((Box[Xp1][Yp1][H[X][Y]-nhxp-1]) ==3) { cntAdGa++; }}

     else if((H[X][Y] % 8)!= 0){
                                if((Box[Xp1][Y]  [H[X][Y]-nhxp-1]) ==2) { cntN++; }
                                if((Box[Xm1][Ym1][H[X][Y]-nhxn-1]) ==2) { cntN++; }
                                if((Box[Xm1][Yp1][H[X][Y]-nhxn-1]) ==2) { cntN++; }
                                if((Box[Xp1][Y]  [H[X][Y]-nhxp-1]) ==3) { cntAdGa++; }
                                if((Box[Xm1][Ym1][H[X][Y]-nhxn-1]) ==3) { cntAdGa++; }
                                if((Box[Xm1][Yp1][H[X][Y]-nhxn-1]) ==3) { cntAdGa++; }} 

    
     if     ((cntN==0) && (cntAdGa==0))             { Prob = Prob_Gaad         ; }
     else if((cntN==0) && (cntAdGa==1))             { Prob = Prob_Ga1AdGa      ; }
     else if((cntN==0) && (cntAdGa==2))             { Prob = Prob_Ga2AdGa      ; }
     else if((cntN==0) && (cntAdGa==3))             { Prob = Prob_Ga3AdGa      ; }
     else if((cntN==1) && (cntAdGa==1))             { Prob = Prob_Ga1N_1AdGa   ; }
     else if((cntN==1) && (cntAdGa==2))             { Prob = Prob_Ga1N_2AdGa   ; }
     else if((cntN==2) && (cntAdGa==1))             { Prob = Prob_Ga2N_1AdGa   ; }
     else if((cntN==1) && (cntAdGa==0))             { Prob = Prob_Ga1N         ; }
     else if((cntN==2) && (cntAdGa==0))             { Prob = Prob_Ga2N         ; }
     else if((cntN==3) && (cntAdGa==0))             { NoGadiff++;      return 0; }
             
           double Racc = gsl_rng_uniform_pos(racc);                     
	  
	   if(Racc < Prob) {
		            if((H[X][Y] % 8) == 0)      { H8_Ga( X, Y) ; } 
			    else if((H[X][Y] % 8) !=0 ) { H4_Ga( X, Y) ; }  
	                   }
	   else { NoGadiff++; move_sbN( X, Y); }
	     	                   
 return 0;
}                  
