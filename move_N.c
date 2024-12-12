// Jan 21, 2022 : Subroutine for diffusion of N atom. N atom can diffuse over Ga atom (having atleast one N atom as neighbor) at H= H+5. If AdGa is found as NN of final diffused site, then shifting of AdGa at H= H+5.
// Jan 21, 2022 : N atom can also diffuse over AdGa atoms, leading to the formation of Ga atoms at H = H+3. Probabilty is calcuated based on the number of direct bonds to which Ga atom is attached.
// Jan 21, 2022 : After N has diffused there will be a check for its NN neighbor (w.r.t Ga atom (X,Y)) for isolated Ga atom. If any isolated Ga  ad atom having N at H= H+3 found, then N atom at H = H+3 will fall over H =H+3 w.r.t Ga(X,Y).


#include"variables.h"

int move_N(int X, int Y)            
{	
 int Xp1 = (X+1)%Nx,   Xm1 = (X-1+Nx)%Nx,  Yp1 = (Y+1)%Ny,  Ym1 = (Y-1+Ny)%Ny ;    	

     cntGa = 0, Prob = 0.0, N4N = 0 ; nhxn = 0; nhxp = 0;

   /* if((X+1)>=Nx){ nhxp= ch_stp;} else if((X-1)<0)  { nhxn= -ch_stp;} //This if..else condition has to be uncommented for vicinal surface*/

     if(((H[X][Y]+1) % 8) ==0){
                               if((Box[Xp1][Y]  [H[X][Y]-nhxp+1]) == 1)      { cntGa++; }
                               if((Box[Xm1][Ym1][H[X][Y]-nhxn+1]) == 1)      { cntGa++; }
                               if((Box[Xm1][Yp1][H[X][Y]-nhxn+1]) == 1)      { cntGa++; }}

     else if(((H[X][Y]+1) % 8) !=0){
                                    if((Box[Xm1][Y]  [H[X][Y]-nhxn+1]) ==1)     { cntGa++; }
                                    if((Box[Xp1][Ym1][H[X][Y]-nhxp+1]) ==1)     { cntGa++; }
                                    if((Box[Xp1][Yp1][H[X][Y]-nhxp+1]) ==1)     { cntGa++; }}                                       
     
     if( cntGa==0){
                   if(((H[X][Y]+1) % 8) ==0){
                                             if((Box[Xp1][Y]  [H[X][Y]-nhxp-4]) == 2)      { N4N++; }
                                             if((Box[Xm1][Ym1][H[X][Y]-nhxn-4]) == 2)      { N4N++; }
                                             if((Box[Xm1][Yp1][H[X][Y]-nhxn-4]) == 2)      { N4N++; }}
                                           
                   else if(((H[X][Y]+1) % 8) !=0){
                                                  if((Box[Xm1][Y]  [H[X][Y]-nhxn-4]) == 2)      { N4N++; }
                                                  if((Box[Xp1][Ym1][H[X][Y]-nhxp-4]) == 2)      { N4N++; }
                                                  if((Box[Xp1][Yp1][H[X][Y]-nhxp-4]) == 2)      { N4N++; }}                                  

                   if (N4N==1)     { Prob = Prob_N0Ga1N; }
                   else if(N4N==2) { Prob = Prob_N0Ga2N; }
                   else if(N4N==3) { Prob = Prob_N0Ga3N; }
  
                   double Racc = gsl_rng_uniform_pos(racc);
                   if(Racc < Prob) {
                                    if(((H[X][Y]+1) % 8) ==0)      { H7_N( X, Y) ; }
                                    else if (((H[X][Y]+1) % 8)!=0) { H3_N( X, Y) ; } }

                   else { /* failed diffusion */  NoNdiff++; }
                  }
 
        else if( cntGa==1) { /* Prob = Prob_nrot || Prob = Prob_nhop (in diff subroutine) ; */ if(((H[X][Y]+1) % 8) ==0) { H7_N_rot( X, Y) ; }  else if (((H[X][Y]+1) % 8)!=0) { H3_N_rot( X, Y) ; } }     
	  			      	
     else if((cntGa==2)||(cntGa==3)) { NoNdiff++;   return 0; }

 return 0;
}                              
