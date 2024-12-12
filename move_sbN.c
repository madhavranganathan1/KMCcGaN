// Jan 21, 2022 : Subroutine for diffusion of sub surface N atom. Sub surface atom can diffuse under AdGa atom and prefer to remain sub surface N atom only by pushing AdGa at H= H+5. If the sub surface N atom is liked to N atom (NNN) on same height that it is part of island, then this move is forbidden.
// Jan 21, 2022 : Ga atom(at H= H+5 w.r.t to sub surface N atom  should has Ga/AdGa only as NN, so AdGa forming after sub surface diffusion is preserved. Probabilty is calcuated based on the number of direct bonds to which Ga atom is attached.
// Jan 21, 2022 : After AdGa has diffused there will be a check for its NN neighbor (w.r.t AdGa atom (X,Y)) for isolated N atom. If any of N atom (NN) is found to be isolated, then this move is rejected. Therefore I m checking n starting that each of N atom (NN) at H= H+4 w.r.t to AdGa (X,Y) should atleast be havingmore then two Ga as a neighbor.


#include"variables.h"

int move_sbN (int X,int Y)           
{   
 int Xp1 = (X+1)%Nx,     Xm1 = (X-1+Nx)%Nx,    Yp1 = (Y+1)%Ny,     Ym1 = (Y-1+Ny)%Ny,    Xp2 = (X+2)%Nx,     Xm2 = (X-2+Nx)%Nx,
     Yp2 = (Y+2)%Ny,     Ym2 = (Y-2+Ny)%Ny,    Xp3 = (X+3)%Nx,     Xm3 = (X-3+Nx)%Nx,    Yp3 = (Y+3)%Ny,     Ym3 = (Y-3+Ny)%Ny;

     Prob = 0.0;   NN_sb = 0, M_Z =0; nhxn = 0; nhxp = 0;
     
   /*if((X+2)>=Nx){ nhxp= ch_stp;} else if((X-2)<0)  { nhxn= -ch_stp;} // This if..else condition has to be uncommented for vicinal surface*/  

     if((H[X][Y]-5) > 16){ 	     
                         stable_AdGa(X, Y, H[X][Y]);

			 if(Ga_Ga == 3){
                                        if((Box[Xm2][Ym1] [H[X][Y]-nhxn-5]) == 2)  { NN_sb++; }
                                        if((Box[Xm2][Yp1] [H[X][Y]-nhxn-5]) == 2)  { NN_sb++; }
                                        if((Box[X][Ym2]   [H[X][Y]-5])      == 2)  { NN_sb++; }
                                        if((Box[Xp2][Ym1] [H[X][Y]-nhxp-5]) == 2)  { NN_sb++; }
                                        if((Box[X][Yp2]   [H[X][Y]-5])      == 2)  { NN_sb++; }
                                        if((Box[Xp2][Yp1] [H[X][Y]-nhxp-5]) == 2)  { NN_sb++; }
       
                                        if(NN_sb == 0){				 
                                                       double Racc = gsl_rng_uniform_pos(racc);
                                                       Prob = Prob_NGa_Ga ;

                                                       if(Racc < Prob) {
                                                                        if((H[X][Y] % 8) == 0)     { H3_sbN( X, Y); }
                                                                        else if((H[X][Y] % 8)!= 0) { H7_sbN( X, Y); }}

                                                       else { NosbNdiff++; }}
                                        else {sbN_Island++; }}
			 else {Ga3_sbN++; }							
                        }                                                                       
return 0;
}                  
