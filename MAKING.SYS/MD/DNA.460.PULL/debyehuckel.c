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
		  float Ene[],
		  float	EneElec_resi[],
                  int idna_elec[], 
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

  int i,intDistSquared,iprot;
  float helper;

  int Elec_resi[14]={1060,1283,1141,1093,1091,1119,1261,1032,784,800,882,540,552,638};
 
  float Emin_elec[14]={[0 ... 13]=100000};


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

       if ((C1 < *firstDNABead-1) && (C2 < *firstDNABead-1))
       {
            *Eprot+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
//	     printf("%4s,%4d,%4d,%10.3f\n","prot",C1, C2, *Eprot);
       }
 
       else if ((C1 >= *firstDNABead-1) && (C2 >= *firstDNABead-1))
       {
            *Edna+= DebyeHuckelPotentials[intDistSquared]*Q1Q2[i];
//	     printf("%3s,%4d,%4d,%10.3f\n","dna",C1, C2, *Edna);
       }
	
      else
      {	
	float En=DebyeHuckelPotentials[intDistSquared]*Q1Q2[i]; 	
	*Eprot_dna+=En;


//        printf("%4d,%4d,%10.3f\n",C1, C2, En);

	if (C1 > C2)
	{
//	printf("%4d,%15.6f,%15.6f\n",C1,En,Ene[C1]);
	Ene[C1]+=En;
//	printf("%4d,%15.6f\n",C1,Ene[C1]);
	}
	if (C1 < C2){
//	printf("%4d,%15.6f,%15.6f\n",C2,En,Ene[C2]);
	Ene[C2]+=En;	
//	printf("%4d,%15.6f\n",C2,Ene[C2]);
	}

//  total eneergy of 14 charge protein residues with each DNA phosphate //
//	for(iprot =0;iprot < 14;iprot++)
//  {	
//	if (C1 == Elec_resi[iprot] || C2 == Elec_resi[iprot])
//	{
//	printf("%4d,%4d,%4d,%4d,%15.6f,%15.6f\n",iprot,C1,C2,Elec_resi[iprot],En,EneElec_resi[iprot]);
//	EneElec_resi[iprot]+=En;
//	printf("%4d,%4d,%4d,%4d,%15.6f,%15.6f\n",iprot,C1,C2,Elec_resi[iprot],En,EneElec_resi[iprot]);
//	}
//  }



	for(iprot =0;iprot < 14;iprot++)
   {
       if (C1 == Elec_resi[iprot] || C2 == Elec_resi[iprot])
    {


    if (C1 > C2 && En < Emin_elec[iprot])
    {idna_elec[iprot]=C1+1;
    Emin_elec[iprot]=En;}
  

    if (C1 < C2 && En < Emin_elec[iprot])
    {idna_elec[iprot]=C2+1;
    Emin_elec[iprot]=En;}

//    printf("%4d,%4d,%4d,%4d,%4d\n", C1, C2, iprot,Elec_resi[iprot],idna_elec[iprot]);
    }
    

   }



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
//       printf("%4d,%4d,%15.6f,%15.6f,%15.6f\n",C1, C2, *Eprot,*Edna,*Eprot_dna);


  }

     	
//     for(iprot =0;iprot < 14;iprot++)
//     {printf("%4d,%4d,%4d\n",iprot,Elec_resi[iprot],idna_elec[iprot]);
//     }

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

