//Jan 21, 2022 : Subroutine for deposition atom. 
// Jan 21, 2022 :Ga atom can deposit over N atom (having atleast one Ga atom as neighbor) and form epitaxil Ga at H= H+5. But if Ga neighbor is not found for diffusing N site, then Ga will get deposit ove N at H= H-4( w.r.t N diffusing site). Ga atom can alo diffuse over Ga atoms (atleast having one N atom as a direct bond) with Ga/AdGa as NN leading to the formation of AdGa atoms.
// Jan 21, 2022 : N atom can deposit over Ga atom (having atleast one N atom as neighbor) at H= H+3. But if N neighbor is not found for diffusing Ga site, then N will get deposit ove Ga at H= H-4( w.r.t N diffusing site). If AdGa is foundat deposition site, then shifting of AdGa at H= H+5 leading to the formation of Ga atom.


#include "variables.h"
void External_Flux()
{                          
  int Neighb, Xc, Yc, AdatomX, AdatomY, Rdepsite; 
      
      Rdepsite = (int) ((double)(nx*ny)*gsl_rng_uniform_pos(rdep));                  
      //printf("Rdepsite=%d",Rdepsite);
       
      Xc = Rdepsite / ny;             
      Yc = Rdepsite % ny;
   
      /* finding actual coordinates */
       AdatomY = Yc;
       if((Yc%2)==0){
         if((Xc%2)==0){AdatomX = (2*Xc);}
           else if((Xc%2)!=0){AdatomX = (2*Xc+1);}}
            
       else if((Yc%2)!=0){
           if((Xc%2)==0){AdatomX = (2*Xc+1);}
             else if((Xc%2)!=0){AdatomX = (2*Xc);}}

  int Xp1 = (AdatomX+1)%Nx,
      Xm1 = (AdatomX-1+Nx)%Nx,
      Yp1 = (AdatomY+1)%Ny,
      Ym1 = (AdatomY-1+Ny)%Ny;


  double Rtype;       
  Rtype = gsl_rng_uniform_pos(rflux);  

//*******************************************************  Ga DEPOSITION *****************************************************************
  if(Rtype < Flux_Ratio)
    {   
     Attempt_Ga++;	
 
     if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]) == 2){                         
         if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-3]) == 1){                 
             Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
             if(Neighb!=0){
                           H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 5;
                           Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 1;  
                           Gaatoms++;     dcnt++;   }
	     else if(Neighb==0){ stable_Ga(AdatomX,AdatomY,H[AdatomX][AdatomY]); }}}                                                                                           

     else if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]) == 1)                   
             {                                                   
              if(((H[AdatomX][AdatomY])%8) == 0)
                {               
                 if((Box[Xm1][AdatomY][H[Xm1][AdatomY]])==2){
                    if((Box[Xm1][AdatomY][H[Xm1][AdatomY]-3])==1){ 
                        Neighb = countNeighb(Xm1,AdatomY,H[Xm1][AdatomY]);
                        if(Neighb!=0){
                                      H[Xm1][AdatomY] = H[Xm1][AdatomY]  + 5;
                                      Box[Xm1][AdatomY][H[Xm1][AdatomY]] = 1;
                                      Gaatoms++;     dcnt++; }
			else if(Neighb==0){ stable_Ga(Xm1,AdatomY,H[Xm1][AdatomY]); }}}		 
                                                                                                                                                     
                 else if((Box[Xp1][Ym1][H[Xp1][Ym1]])==2){
                        if((Box[Xp1][Ym1][H[Xp1][Ym1]-3])==1){
                            Neighb = countNeighb(Xp1,Ym1,H[Xp1][Ym1]);
                            if(Neighb!=0){
                                          H[Xp1][Ym1] = H[Xp1][Ym1]  + 5;
                                          Box[Xp1][Ym1][H[Xp1][Ym1]] = 1;
                                          Gaatoms++;     dcnt++;   } 
                            else if(Neighb==0){ stable_Ga(Xp1,Ym1,H[Xp1][Ym1]); }}}			    
			   
                 else if((Box[Xp1][Yp1][H[Xp1][Yp1]])==2){
                        if((Box[Xp1][Yp1][H[Xp1][Yp1]-3])==1){
                            Neighb = countNeighb(Xp1,Yp1,H[Xp1][Yp1]);
                            if(Neighb!=0){
                                          H[Xp1][Yp1] = H[Xp1][Yp1]  + 5;
                                          Box[Xp1][Yp1][H[Xp1][Yp1]] = 1;
                                          Gaatoms++;     dcnt++;   }
			    else if(Neighb==0){ stable_Ga(Xp1,Yp1,H[Xp1][Yp1]); }}}
                                                  
                 else{
                      Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
                      if(Neighb!=0){
                                    stable_AdGa(AdatomX,AdatomY,H[AdatomX][AdatomY]) ;
                                    if(Ga_Ga == 3){
                                                   H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 3;
                                                   Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 3;
                                                   AdGaatoms++;     dcnt++;  }}}			      
		}                  
        
              else if(((H[AdatomX][AdatomY])%8) != 0)
                     {                           
                      if((Box[Xp1][AdatomY][H[Xp1][AdatomY]])==2){
                        if((Box[Xp1][AdatomY][H[Xp1][AdatomY]-3])==1){
                            Neighb = countNeighb(Xp1,AdatomY,H[Xp1][AdatomY]);
                            if(Neighb!=0){
                                          H[Xp1][AdatomY] = H[Xp1][AdatomY]  + 5;
                                          Box[Xp1][AdatomY][H[Xp1][AdatomY]] = 1;
                                          Gaatoms++;  dcnt++; }
			    else if(Neighb==0){ stable_Ga(Xp1,AdatomY,H[Xp1][AdatomY]); }}}
		      
                      else if((Box[Xm1][Ym1][H[Xm1][Ym1]])==2){
                             if((Box[Xm1][Ym1][H[Xm1][Ym1]-3])==1){
                                 Neighb = countNeighb(Xm1,Ym1,H[Xm1][Ym1]);
                                 if(Neighb!=0){
                                               H[Xm1][Ym1] = H[Xm1][Ym1]  + 5;
                                               Box[Xm1][Ym1][H[Xm1][Ym1]] = 1;
                                               Gaatoms++; dcnt++; }
	                         else if(Neighb==0){ stable_Ga(Xm1,Ym1,H[Xm1][Ym1]);}}}			 

                      else if((Box[Xm1][Yp1][H[Xm1][Yp1]])==2){
                             if((Box[Xm1][Yp1][H[Xm1][Yp1]-3])==1){
                                 Neighb = countNeighb(Xm1,Yp1,H[Xm1][Yp1]);
                                 if(Neighb!=0){
                                               H[Xm1][Yp1] = H[Xm1][Yp1]  + 5;
                                               Box[Xm1][Yp1][H[Xm1][Yp1]] = 1;
                                               Gaatoms++; dcnt++; }
		                 else if(Neighb==0){ stable_Ga(Xm1,Yp1,H[Xm1][Yp1]);}}}	     
	             
		      else{
                           Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
                           if(Neighb!=0){
                                         stable_AdGa(AdatomX,AdatomY,H[AdatomX][AdatomY]) ;
                                         if(Ga_Ga == 3){
                                                        H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 3;
                                                        Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 3;
                                                        AdGaatoms++;     dcnt++;  }}}	
		     }  
             }   
    }

// **********************************************************  N DEPOSITION  *****************************************************************

  else if((Flux_Ratio < Rtype) && (Rtype < 1.00))
          { 
           Attempt_N++;		  

           if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]) == 1){                     
             if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-5]) == 2){
                 Neighb = countNeighb(AdatomX,AdatomY,H[AdatomX][AdatomY]);
                 if(Neighb!=0){			       
                               H[AdatomX][AdatomY] = H[AdatomX][AdatomY]  + 3;
                               Box[AdatomX][AdatomY][H[AdatomX][AdatomY]] = 2;
			       Natoms++;    dcnt++;  stable_N(AdatomX,AdatomY); }
	         else if(Neighb==0){ stable_N_OGa(AdatomX,AdatomY,H[AdatomX][AdatomY]); }}}

           else if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]])==2)                        
                  {
                   if(((H[AdatomX][AdatomY]+1)%8)==0)  
                     {        	 	 
                      if((Box[Xp1][AdatomY][H[Xp1][AdatomY]])==1){                         
                        if((Box[Xp1][AdatomY][H[Xp1][AdatomY]-5])==2){
		            Neighb = countNeighb(Xp1,AdatomY,H[Xp1][AdatomY]);
                            if(Neighb!=0){ 
                                          H[Xp1][AdatomY] = H[Xp1][AdatomY]  + 3;
                                          Box[Xp1][AdatomY][H[Xp1][AdatomY]] = 2;
					  Natoms++;   dcnt++;   stable_N(Xp1,AdatomY); } 
			    else if(Neighb==0){ stable_N_OGa(Xp1,AdatomY,H[Xp1][AdatomY]); }}} 		      
		                           
                      else if((Box[Xm1][Ym1][H[Xm1][Ym1]])==1){
                             if((Box[Xm1][Ym1][H[Xm1][Ym1]-5])==2){
                                 Neighb = countNeighb(Xm1,Ym1,H[Xm1][Ym1]);
                                 if(Neighb!=0){
                                               H[Xm1][Ym1] = H[Xm1][Ym1]  + 3;
                                               Box[Xm1][Ym1][H[Xm1][Ym1]] = 2;
					       Natoms++;   dcnt++;   stable_N(Xm1,Ym1); } 
				 else if(Neighb==0){ stable_N_OGa(Xm1,Ym1,H[Xm1][Ym1]); }}}

                      else if((Box[Xm1][Yp1][H[Xm1][Yp1]])==1){
                             if((Box[Xm1][Yp1][H[Xm1][Yp1]-5])==2){
                                 Neighb = countNeighb(Xm1,Yp1,H[Xm1][Yp1]);
                                 if(Neighb!=0){
                                               H[Xm1][Yp1] = H[Xm1][Yp1]  + 3;
                                               Box[Xm1][Yp1][H[Xm1][Yp1]] = 2;
					       Natoms++;   dcnt++;   stable_N(Xm1,Yp1); } 
				 else if(Neighb==0){ stable_N_OGa(Xm1,Yp1,H[Xm1][Yp1]); }}}				 
		     } 

		   else if(((H[AdatomX][AdatomY]+1)%8)!=0)
                          { 
                           if((Box[Xm1][AdatomY][H[Xm1][AdatomY]])==1){                  
                             if((Box[Xm1][AdatomY][H[Xm1][AdatomY]-5])==2){  
                                 Neighb = countNeighb(Xm1,AdatomY,H[Xm1][AdatomY]);
                                 if(Neighb!=0){
                                               H[Xm1][AdatomY] = H[Xm1][AdatomY]  + 3;
                                               Box[Xm1][AdatomY][H[Xm1][AdatomY]] = 2;
					       Natoms++;   dcnt++;  stable_N(Xm1,AdatomY); } 
				 else if(Neighb==0){stable_N_OGa(Xm1,AdatomY,H[Xm1][AdatomY]); }}}
                    
                           else if((Box[Xp1][Ym1][H[Xp1][Ym1]])==1){
                                  if((Box[Xp1][Ym1][H[Xp1][Ym1]-5])==2){
                                      Neighb = countNeighb(Xp1,Ym1,H[Xp1][Ym1]);
                                      if(Neighb!=0){
                                                    H[Xp1][Ym1] = H[Xp1][Ym1]  + 3;
                                                    Box[Xp1][Ym1][H[Xp1][Ym1]] = 2;
                                                    Natoms++;   dcnt++;  stable_N(Xp1,Ym1); } 
				      else if(Neighb==0){ stable_N_OGa(Xp1,Ym1,H[Xp1][Ym1]); }}}   			   

                           else if((Box[Xp1][Yp1][H[Xp1][Yp1]])==1){
                                  if((Box[Xp1][Yp1][H[Xp1][Yp1]-5])==2){
                                      Neighb = countNeighb(Xp1,Yp1,H[Xp1][Yp1]);
                                      if(Neighb!=0){
                                                    H[Xp1][Yp1] = H[Xp1][Yp1]  + 3;
                                                    Box[Xp1][Yp1][H[Xp1][Yp1]] = 2;
						    Natoms++;   dcnt++;  stable_N(Xp1,Yp1); }
				      else if(Neighb==0){ stable_N_OGa(Xp1,Yp1,H[Xp1][Yp1]); }}}				      
                          } 
                  }		   

	   else if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]])==3) {
		       if((Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-3])==1){
                                       H[AdatomX][AdatomY] = H[AdatomX][AdatomY]    + 5 ;
				       Box[AdatomX][AdatomY][H[AdatomX][AdatomY]]   = 1 ;
				       Box[AdatomX][AdatomY][H[AdatomX][AdatomY]-5] = 2 ;
				       Natoms++; N_AdGa++;  dcnt++;   }}           	   
       } 

}    
