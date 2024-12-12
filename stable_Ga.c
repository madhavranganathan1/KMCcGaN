// Jan 22, 2022: Code is written to check if N atom has zero direct neighbor. In that case Ga get deposited over N.


#include "variables.h"

int stable_Ga(int X, int Y, int Z)  
{
	
 int Xp1 = (X+1)%Nx,
     Xm1 = (X-1+Nx)%Nx,
     Yp1 = (Y+1)%Ny,
     Ym1 = (Y-1+Ny)%Ny;
     

     if((Z+1)%8==0){
                    if     ((Box[Xp1][Y]  [H[Xp1][Y]])   ==2)      { H[Xp1][Y]   = H[Xp1][Y]   + 5;  Box[Xp1][Y]  [H[Xp1][Y]]   = 1; Gaatoms++; dcnt++; }
		    else if((Box[Xm1][Ym1][H[Xm1][Ym1]]) ==2)      { H[Xm1][Ym1] = H[Xm1][Ym1] + 5;  Box[Xm1][Ym1][H[Xm1][Ym1]] = 1; Gaatoms++; dcnt++; }
		    else if((Box[Xm1][Yp1][H[Xm1][Yp1]]) ==2)      { H[Xm1][Yp1] = H[Xm1][Yp1] + 5;  Box[Xm1][Yp1][H[Xm1][Yp1]] = 1; Gaatoms++; dcnt++; }
                   }

     else if((Z+1)%8!=0){
	                 if     ((Box[Xm1][Y]  [H[Xm1][Y]])   ==2)  { H[Xm1][Y]   = H[Xm1][Y]   + 5;  Box[Xm1][Y]  [H[Xm1][Y]]   = 1; Gaatoms++; dcnt++; }
                         else if((Box[Xp1][Ym1][H[Xp1][Ym1]]) ==2)  { H[Xp1][Ym1] = H[Xp1][Ym1] + 5;  Box[Xp1][Ym1][H[Xp1][Ym1]] = 1; Gaatoms++; dcnt++; }
                         else if((Box[Xp1][Yp1][H[Xp1][Yp1]]) ==2)  { H[Xp1][Yp1] = H[Xp1][Yp1] + 5;  Box[Xp1][Yp1][H[Xp1][Yp1]] = 1; Gaatoms++; dcnt++; }
                        }
return 0;

}

        


