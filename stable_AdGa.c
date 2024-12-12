// Sept 28, 2020: Code is written to check Ga and AdGa neighbors to Ga site selected for diffusing an  Ga atom. If Ga_Ga = 3, then only Ga will be deposited.


#include "variables.h"

int stable_AdGa(int X, int Y, int Z)  
{
                               
   int Xp1 = (X+1)%Nx,
       Xm1 = (X-1+Nx)%Nx,
       Yp1 = (Y+1)%Ny,
       Ym1 = (Y-1+Ny)%Ny;

       Ga_Ga =0; 
      
       if(Z%8==0){
                  if(((Box[Xm1][Y]  [H[Xm1][Y]])   ==1)||((Box[Xm1][Y]  [H[Xm1][Y]])   ==3))      { Ga_Ga++; }
                  if(((Box[Xp1][Ym1][H[Xp1][Ym1]]) ==1)||((Box[Xp1][Ym1][H[Xp1][Ym1]]) ==3))      { Ga_Ga++; }
                  if(((Box[Xp1][Yp1][H[Xp1][Yp1]]) ==1)||((Box[Xp1][Yp1][H[Xp1][Yp1]]) ==3))      { Ga_Ga++; }		  
                 }

       else if(Z%8!=0){
                       if(((Box[Xp1][Y]  [H[Xp1][Y]])   ==1)||((Box[Xp1][Y]  [H[Xp1][Y]])   ==3)) { Ga_Ga++; }
                       if(((Box[Xm1][Ym1][H[Xm1][Ym1]]) ==1)||((Box[Xm1][Ym1][H[Xm1][Ym1]]) ==3)) { Ga_Ga++; }
                       if(((Box[Xm1][Yp1][H[Xm1][Yp1]]) ==1)||((Box[Xm1][Yp1][H[Xm1][Yp1]]) ==3)) { Ga_Ga++; }		       
                      }
return 0;
}

        


