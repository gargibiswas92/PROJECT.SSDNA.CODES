//********************************************************************************************/
// electrostatics energy term
//
//********************************************************************************************/
#include <stdio.h>
#include <math.h>
void debyehuckel_(int IC[],
		  int JC[],
		  float Q1Q2[],
		  float X[],
		  float Y[],
		  float Z[],
		  float Fx[],
		  float Fy[],
		  float Fz[] ,
		  int *NE,
                  float *cutoffDist,
		  float *Eelec,
		  float *Eprot,
		  float *Edna,
		  float *Eprot_dna,
		  float *DebyeHuckelPotentials,
		  float *DebyeHuckelForces,
                  int *firstDNABead) 
{

  int i,intDistSquared;
  float helper;

  for(i =0;i < *NE;i++)
  {
     
      int C1 = IC[i]-1;
      int C2 = JC[i]-1;
      float dx = X[C1] - X[C2];
      float dy = Y[C1] - Y[C2];
      float dz = Z[C1] - Z[C2];

      float r2 = pow(dx,2) + pow(dy,2) + pow(dz,2);
      
      

      if (r2< *cutoffDist)
      {
      helper=100*r2;   //make sure to multiple at the inverse of interval in dhEnergyTable.c
      intDistSquared=(int) helper;
      *Eelec+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
       

       if ((C1 < *firstDNABead) && (C2 < *firstDNABead))
       {
            *Eprot+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
       }
 
       else if ((C1 >= *firstDNABead) && (C2 >= *firstDNABead))
       {
            *Edna+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
       }
	
      else
      {	*Eprot_dna+=DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
      }





       // force in the direction C1 to C2 ,devided by r (is used to compute forces in x,y,z directions)
       float F_over_r = DebyeHuckelForces[intDistSquared]*Q1Q2[i];
		
       Fx[C1]+=  F_over_r*dx;
       Fx[C2]-=  F_over_r*dx;
       Fy[C1]+=  F_over_r*dy;
       Fy[C2]-=  F_over_r*dy;
       Fz[C1]+=  F_over_r*dz;
       Fz[C2]-=  F_over_r*dz;
      }

  }

}

void debyehuckelfactor_(float *sigma,
			float *deConstant,
			float *screeningFactor,
			float *saltCoefficient,
			float *esEnergy)
{
  float K = 332.0; //?? not sure this is the correct number.
  *esEnergy = K*(*saltCoefficient)*exp(-(*screeningFactor)*(sqrt(*sigma)))/((*deConstant)*(sqrt(*sigma)));
}

