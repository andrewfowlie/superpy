/* Supersymmetric contribution to muon g-2
      
Author: A. Cottrant, june 2001             
Ref: A. Cottrant,  "Le moment anomal du muon en supersymetrie".
Rapport de DEA, Universite de Savoie.
These formulas  agree with the ones in S. Martin, J. Wells, hep-ph/0103067

Corrections:
28 march 2002: Correct factor of 1/2. in double result  line 87 
*/

#include"../../sources/micromegas.h"
#include "pmodel.h"
#include "matrix.h"

static double integr(int number ,double var) 

/*
Computes the integer being used in the calculation of g-2
with the notation of Martin & Wells :
integr(0, var) = -Fn1 (var)
integr(2,var) = Fn2 (var)/2
integr(1,var) = Fc1 (var)
integr(3,var) = Fc2 (var)
*/
	{double result;
	switch (number){
	case 0:
		{
		if ((var < 1.01)&&(var>0.99)) result = - 1; 
		else result = -2/(pow((1-var),4))*(1-6*var+3*pow(var,2)+2*pow(var,3)-6*pow(var,2)*log(var));	
		return(result);
		break;}
	case 2:
		{
		if ((var < 1.01)&&(var>0.99)) result =1/2.;
		else result = 3/(2*pow((1-var),3))*(1-pow(var,2)+2*var*log(var));
		return(result);
		break;}
	case 1:
		{if ((var < 1.01)&&(var>0.99)) result =1; 
		else result = 2/(pow((1-var),4))*(2+3*var-6*pow(var,2)+pow(var,3)+6*var*log(var));	
		return(result);
		break;}
	case 3:
		{if ((var < 1.01)&&(var>0.99)) result =1; 
		else result = -3/(2*pow((1-var),3))*(3-4*var+pow(var,2)+2*log(var));	
		return(result);
		break;}
			}
	return(0);
	}

static double calcgmuon (int csboson,double mmuon, double mbosino, double msfermion,double coeff1, double coeff2)
	/*compute the contribution for the g-2 from 1 graph
	The parameters are 
	csbosons : charge of the vector boson
	mmuon : muon's mass
	mbosino : sboson's mass
	msfermion : sfermion's mass
	coeff1 : 1st coefficient (for the term prop to mmuon )
	coeff2 : 2nd coefficient (for the term prop to mbosino )
	*/
	
	{
	double x,result;
	x = pow ((mbosino/msfermion),2);
	result = mmuon/(16*pow(M_PI,2))*(mmuon/(12*pow(msfermion,2))*coeff1*integr(csboson,x)+(2*mbosino)/(3*pow(msfermion,2))*coeff2*integr(csboson+2,x));
	return(result);
	}

double gmuon_(void)
{
/*compute the total g-2 from supersymmetry. Parameters are :
Mmu : muon mass
Mz  : Z-boson's mass
Mw  : W-boson's mass
s2thw : sin^2(thw)
e   : alpha
modM1 & argM1 : modulus and argument of M1		
modM2 & argM2 : modulus and argument of M2
modmu & argmu : modulus and argument of mu		
modAmu & argAmu : modulus and argument of Amu	
*/
/*declaration des variables*/
double thw, cb,cw,sw, ymu, g1,masssneutr, g2,gmuo = 0,beta,Mz,Mw,Mmu,e,tbeta;
matrix massmatrc,massmatrn,massmatrsmu;
matrix Ncomp, Nconj,N,Ucomp,Uconj,U,Vcomp,Vconj,V,X,Xconj,Xcomp;
vecteur massen,massec,massesmu;

double coeff1 ,coeff2;

matrix matrixtemp,matrixtemp2,matrixtemp3,matrixtemp4,matrixtemp5;
my_complex nc1,nc2;

int i,j;
/*Memory Allocation*/

massmatrn=initm(4);
massmatrc=initm(2);
massmatrsmu=initm(2);
Ncomp=initm(4);
Nconj=initm(4);
N=initm(4);

Ucomp=initm(2);
Vcomp=initm(2);
X=initm(2);
U=initm(2);
V=initm(2);
Uconj=initm(2);
Vconj=initm(2);
Xcomp=initm(2);
Xconj=initm(2);
massen=initv(4);
massesmu=initv(2);
massec=initv(2);

matrixtemp=initm(2);
matrixtemp2=initm(2);
matrixtemp3=initm(4);
matrixtemp4=initm(4);
matrixtemp5=initm(4);

/*Get the values for the masses and mixing matrices*/

e=sqrt(4*M_PI*0.00781653);
sw=0.48076;
Mz=91.1876;
Mmu=0.1057;

findVal("tB",&tbeta);


beta = atan(tbeta);
cb = cos(beta);
thw=asin(sw);
cw=cos(thw);
Mw=Mz*cw;
g1 = e/cw;
g2 = e/sw;
ymu = g2 * Mmu / (sqrt(2)*cb*Mw);

findVal("MNE1",&massen.e[0].r);
findVal("MNE2",&massen.e[1].r);
findVal("MNE3",&massen.e[2].r);
findVal("MNE4",&massen.e[3].r);

findVal("MC1",&massec.e[0].r);
findVal("MC2",&massec.e[1].r);

findVal("MSnm",&masssneutr);

findVal("Zn11",&Nconj.e[0][0].r); findVal("Zn12",&Nconj.e[1][0].r); findVal("Zn13",&Nconj.e[2][0].r); findVal("Zn14",&Nconj.e[3][0].r);
findVal("Zn21",&Nconj.e[0][1].r); findVal("Zn22",&Nconj.e[1][1].r); findVal("Zn23",&Nconj.e[2][1].r); findVal("Zn24",&Nconj.e[3][1].r);
findVal("Zn31",&Nconj.e[0][2].r); findVal("Zn32",&Nconj.e[1][2].r); findVal("Zn33",&Nconj.e[2][2].r); findVal("Zn34",&Nconj.e[3][2].r);
findVal("Zn41",&Nconj.e[0][3].r); findVal("Zn42",&Nconj.e[1][3].r); findVal("Zn43",&Nconj.e[2][3].r); findVal("Zn44",&Nconj.e[3][3].r);

findVal("Zu11",&Uconj.e[0][0].r); findVal("Zu12",&Uconj.e[1][0].r);
findVal("Zu21",&Uconj.e[0][1].r); findVal("Zu22",&Uconj.e[1][1].r);

findVal("Zv11",&Vcomp.e[0][0].r); findVal("Zv12",&Vcomp.e[0][1].r);
findVal("Zv21",&Vcomp.e[1][0].r); findVal("Zv22",&Vcomp.e[1][1].r);


{  double MSmuLL,MSmuRR,MSmuLR,Am,mu,MSmuth;

   MSmuLL=findValW("MSmL"); MSmuLL*=MSmuLL;
   MSmuRR=findValW("MSmR"); MSmuRR*=MSmuRR;
   Am=findValW("Am");
   mu=findValW("mu");
   MSmuLR=(Am-mu*tbeta)*Mmu;
 
   massesmu.e[0].r=sqrt((MSmuLL+MSmuRR-sqrt((MSmuLL-MSmuRR)*(MSmuLL-MSmuRR)+4*MSmuLR*MSmuLR))/2);
   massesmu.e[1].r=sqrt((MSmuLL+MSmuRR+sqrt((MSmuLL-MSmuRR)*(MSmuLL-MSmuRR)+4*MSmuLR*MSmuLR))/2);

   MSmuth=atan2(-2*MSmuLR, -MSmuLL+MSmuRR)/2;

   X.e[0][0].r=cos(MSmuth);
   X.e[0][1].r=sin(MSmuth);
   X.e[1][0].r=-X.e[0][1].r;
   X.e[1][1].r=X.e[0][0].r;
}




for(i=0;i<4;i++)
	if (massen.e[i].r<0)
		{for(j=0;j<4;j++)
		     {Nconj.e[j][i]=prod(icomp,Nconj.e[j][i]);
		     }
		 massen.e[i].r=-massen.e[i].r;
		}


for(i=0;i<2;i++)
	{if (massec.e[i].r<0)
		{for(j=0;j<2;j++)
		     { Uconj.e [j][i]=prod(icomp,Uconj.e[j][i]);
		     Vcomp.e [i][j]=prod(icomp,Vcomp.e[i][j]);
		     }
		 massec.e[i].r=-massec.e[i].r;
		}
	if (massesmu.e[i].r<0)
		{for(j=0;j<2;j++)
		     {X.e [i][j]=prod(icomp,X.e[i][j]);
		     
		     }
		massesmu.e[i].r=-massesmu.e[i].r;
		   
		}
	}

N=adj(Nconj);
U=adj(Uconj);

/*Compute the coefficients entering the formula for the neutral  and chargino
contribution to g-2_muon and calculates gmuon*/

/*neutralinos*/

for ( i=0;i<4;i=i+1)
	for ( j=0;j<2;j=j+1)
		{nc1=somme(prodscal(sqrt(2)*g1,prod(N.e[i][0],X.e[j][1])),prodscal(ymu,prod(N.e[i][2],X.e[j][0])));		
		nc2=diff(prodscal(1/sqrt(2),prod(somme(prodscal(g2,N.e[i][1]),prodscal(g1,N.e[i][0])),conjug(X.e[j][0]))),prodscal(ymu,prod(N.e[i][2],conjug(X.e[j][1]))));
		
		
		coeff1 =(somme(prod(nc1,conjug(nc1)),prod(nc2,conjug(nc2)))).r;
		coeff2 = prod(nc1,nc2).r;
			
		if (massesmu.e[j].r<1e-8) 
			{return 0;
			printf("erreur : Mass smuons nul\n");}		
		gmuo=gmuo+calcgmuon (0,Mmu, massen.e[i].r, massesmu.e[j].r,coeff1, coeff2);
	}

/*charginos*/

for (j=0;j<2;j=j+1)
	{
	nc1=prodscal(ymu,U.e[j][1]);
	nc2=prodscal(-g2,conjug(Vcomp.e[j][0]));
	coeff1 =(somme(prod(nc1,conjug(nc1)),prod(nc2,conjug(nc2)))).r;
	coeff2 = prod(nc1,nc2).r;
	if (masssneutr<1e-8) 
		{return 0;
		printf("erreur : Mass sneutrinos nul\n");}
	gmuo=gmuo+calcgmuon (1,Mmu, massec.e[j].r, masssneutr,coeff1, coeff2);	
	}

 return gmuo;

}

