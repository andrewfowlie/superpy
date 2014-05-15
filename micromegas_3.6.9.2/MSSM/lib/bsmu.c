/* b->s mu+ mu-  in MSSM with main tanbeta NLO contributions

Notation of Bobeth et al. hep-ph/0104284
*/

#include "../../sources/micromegas.h"
#include "pmodel.h"
#include "pmodel_aux.h"

static double D2(double x, double y) { 
if (fabs(y-1.)<.01)
{return
 1.0/(x-y)*(x*log(x)/(1.0-x)+1.);
}
if (fabs(x-y)<.0001)
{return
1.0/(1.-x)*(1.+log(x)/(1.0-x));
}
return 1.0/(x-y)*(x*log(x)/(1.0-x)-y*log(y)/(1.0-y));}
					  
static double D4(double x) {
if (fabs(x-1.)<.0001)
{return -.5;
}
return 1.0/(1.0-x)*(1.0+log(x)/(1.0-x));}

static double D3(double x) {return (x*log(x)/(1.0-x));}
				
static double D1(double x, double y, double z) {if (fabs(z-1)<.0001){
return
 (x*log(x)/(1.0-x)/(x-y)/(x-z)+y*log(y)/(1.0-y)/(y-x)/(y-z)
					   -1./(z-y)/(z-x));
}
 return 
					   (x*log(x)/(1.0-x)/(x-y)/(x-z)+y*log(y)/(1.0-y)/(y-x)/(y-z)
					   +z*log(z)/(1.0-z)/(z-y)/(z-x));}

#ifdef notUsed					  
static double H2(double x, double y) { return x*log(x)/(1.0-x)/(x-y)
				+y*log(y)/(1.0-y)/(y-x);}					 

static double B1(double x, double y) { return .5*(.5+1./(1.-y/x)+log(y/x)/(1-y/x)/(1-y/x)
				-log(y));}					   
					   				

static double alps6(double mu, double mz, double alpsmz)
{
	return alpsmz/(1.+7./2./M_PI*alpsmz*log(mu/mz));
}

	

static	double U[3][3], V[3][3], T[3][3], mc[3], mst[3], msq1, mglu, At, Ab;
static double B[3][3], msb[3], mss[3], mn[5], N[5][5], msnmu;
static	double m1, m2, mu, Att;
static	double MW, MZ, Mt, Mb, Mhc, Vts,Vtb,Vcb;

static	double alphem =  1.0/137.036;
static	double alphew = 0.007819;
static double GF=1.16639e-5; 
static double fB_s=0.230;
static double TauB_s=2.222e+12/*Lifetime of B_s=1.464e-12sec*/;
#endif

static double yb;
static	double mmu=.1059/* mass of muon*/;
static double MBs=5.369/* Mass B_s*/;
static double alpsmz=0.1185; /* alpha_s(Mz) */
static	double sq2, b, sb, cb;
static void Cpen_sp(struct read_param_tag* param, double *, double *,  double);
static void Cep_sp(struct read_param_tag* param, double *, double *,  double);
static void Cbox_sp(struct read_param_tag* param, double *, double *,  double);



static double alps(double mu, double mz, double alpsmz)
{ double v=1.-23.*alpsmz*log(mz/mu)/6./M_PI;
  return alpsmz*(1-116.*alpsmz*log(v)/v/4./M_PI/23.)/v;
}


double bsmumu_(void)
{
	double xtw,ytw;
	double beta0,beta1,g0m,g1m,mb;
	double mtp,mtmt,mtmw,ytmt,k,csk,cpk,bsmuk;
	double bsmu;
	double alpsmt,alpsmw,delmt,mtmw2;
	double C10sm,Yxt,rhct,CsHc;
	double Cbox_S,Cbox_P,Cpen_S,Cpen_P,Cep_S,Cep_P,Cs_susy,Cp_susy;
	double facbsmu,betamu;
	struct read_param_tag param;
	double eps_b=0, eps_bp=0, eps_tps=0;
	
	if(read_prm(&param))
	{
		puts("bsgammanew: can not read parameters.");
		return 0.0;
	}
	
/*	dump_prm(&param);
*/
        alpsmz=findValW("alfSMZ"); 
	sq2=sqrt(2.0);
	b=atan(param.tb);
	sb=sin(b);
	cb=cos(b);
	mb=MbPole;
	/*mc=1.392;*/
	mtp=findValW("Mtp");
	beta0=23./3.;
	beta1=116./3.;
	g0m=8.;
	g1m=404./3.-200./9.;
	/*mbs=50.;
	lam2=.12;
	del=.9;
	mub=4.8;
	mubar=4.8;
	ckmf=sqrt(.971);*/
	
	
	alpsmt=alps(mtp,param.MZ,alpsmz);
	alpsmw=alps(param.MW,param.MZ,alpsmz);
	delmt=alpsmt*g0m*(g1m/g0m-beta1/beta0)*(alpsmw/alpsmt-1.0)/8.0/beta0/M_PI;
		
	mtmw=mtp*(1-4*alpsmt/3./M_PI)*pow(alpsmw/alpsmt,12./23.)*(1.+delmt);
	mtmw2=mtp*mtp*(1-8.*alpsmt/3./M_PI)*pow(alpsmw/alpsmt,24./23.)
	*(1.+2.*delmt);
	xtw=mtmw2/param.MW/param.MW;
	ytw=mtmw2/param.Mhc/param.Mhc; 
	
	mtmt=mtp*(1.-4./3.*alpsmt/M_PI);
	ytmt=mtmt/sq2/param.MW*param.ee/param.sw/sb;
	yb=mb/sq2/param.MW*param.ee/param.sw/cb;
	
	

/* Standard model NLO */

	Yxt=.997*pow(mtmt/166.,1.56);
	C10sm=-Yxt/param.sw/param.sw;
	

/* Susy Contributions*/
	rhct=param.Mhc*param.Mhc/mtmt/mtmt;	
	CsHc=param.tb*param.tb*log(rhct)/(1.0-rhct);
	
/* Penguin counterterms leading*/

     Cep_sp(&param, &Cep_S, &Cep_P, mtmt);
/* Box contribution*/
    Cbox_sp(&param, &Cbox_S, &Cbox_P, mtmt);

/* Penguin direct*/

	Cpen_sp(&param, &Cpen_S, &Cpen_P, mtmt);
/*Large tan beta */
/* Calculates the eps_b eps_b' ..*/	
	calc_eps(&param, &eps_b, &eps_bp, &eps_tps);
		
       k=1.0/(1.0+eps_b*param.tb);	
	
	Cs_susy=Cbox_S+Cpen_S+Cep_S+CsHc;
	Cp_susy=Cbox_P+Cpen_P+Cep_P-CsHc;
	
	csk=Cs_susy*k*k;
	cpk=Cp_susy*k*k;
	
	betamu=pow((1.-4.*mmu*mmu/MBs/MBs),.5);
	facbsmu=7.827*pow(10.,-10)*betamu;
/*	printf("betamu %.3E %.3E  \n",betamu,facbsmu);
	printf("Box,pen,ep Hc SM %.3E %.3E %.3E %.5E %.3E %.5E \n",Cbox_S,Cpen_S,Cep_S,Cs_susy,CsHc,Yxt);
	printf("Box,pen,epP %.3E %.3E %.3E %.3E\n",Cbox_P,Cpen_P,Cep_P,Cp_susy);
	
	printf("Cs-susy Cp-susy %.3E %.3E \n",Cs_susy,Cp_susy);
	printf("k,eps_b %.3E %.3E %.3E\n",k,eps_b);
	*/
	/*Cs_susy=0.;
	Cp_susy=0.;
	*/
	
	bsmu=facbsmu*(betamu*betamu*(Cs_susy*Cs_susy)*
	MBs*MBs*MBs*MBs/16./param.MW/param.MW/param.MW/param.MW+
	(MBs*MBs/4./param.MW/param.MW*Cp_susy-2.*Yxt)*(MBs*MBs/4./param.MW/param.MW*Cp_susy-2.*Yxt));
 /* for test with keps_b factor*/	
	 bsmuk=facbsmu*(betamu*betamu*(csk*csk)*
	MBs*MBs*MBs*MBs/16./param.MW/param.MW/param.MW/param.MW+
	(MBs*MBs/4./param.MW/param.MW*cpk-2.*Yxt)*(MBs*MBs/4./param.MW/param.MW*cpk-2.*Yxt));
	
	return bsmu;
}



#define CT (param->T[1][1])
#define ST (param->T[1][2])
#define CB (param->B[1][1])
#define SB (param->B[1][2])



/* Leading order chargino contribution contains the leading higher order effect assumes that sparticles have mass  1 TeV with possibly lighter charginos and Higgs*/

static void Cep_sp(struct read_param_tag* param, double *Cep_S, double *Cep_P,  double mtmt)

{
		
	
	int a;
	double yst1,yst2;
	double t=0.;

	for(a=1;a<=2;a++) 
	{yst1=param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a];
	 yst2=param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a];
	 t+=param->mc[a]*ST*CT*mtmt*param->V[a][2]*param->U[a][2]*(D3(yst1)-D3(yst2));
	 t+=-param->mc[a]*sq2*param->MW*param->V[a][1]*param->U[a][2]*(CT*CT*D3(yst1)+ST*ST*D3(yst2)- 			  								
	 D3(param->msq1*param->msq1/param->mc[a]/param->mc[a]));
	 
	 
	}



	*Cep_S=t*param->tb*param->tb*param->tb/(param->Mhc*param->Mhc-param->MW*param->MW);	

	
	*Cep_P=-t*param->tb*param->tb*param->tb/(param->Mhc*param->Mhc-param->MW*param->MW);	
}

static void Cbox_sp(struct read_param_tag* param, double *Cbox_S, double *Cbox_P,  double mtmt)

{
		
	int i,j;
	double yst1,yst2,ysq,xli,f1,f2,fq,f1p,f2p,fqp,zji;
	
	double t=0.;
	double tp=0.;
	
	for(i=1;i<=2;i++) 
	{yst1=param->mst[1]*param->mst[1]/param->mc[i]/param->mc[i];
	 yst2=param->mst[2]*param->mst[2]/param->mc[i]/param->mc[i];
	 ysq=param->msq1*param->msq1/param->mc[i]/param->mc[i];
	 xli= param->msnmu*param->msnmu/param->mc[i]/param->mc[i];
	 
	 	for(j=1;j<=2;j++) 
	 	{zji=param->mc[j]*param->mc[j]/param->mc[i]/param->mc[i];
		f1=(yst1*param->U[j][2]*param->V[i][1]+param->mc[j]/param->mc[i]*param->U[i][2]*param->V[j][1])*D1(xli,yst1,zji);
		f2=(yst2*param->U[j][2]*param->V[i][1]+param->mc[j]/param->mc[i]*param->U[i][2]*param->V[j][1])*D1(xli,yst2,zji);
		fq=(ysq*param->U[j][2]*param->V[i][1]+param->mc[j]/param->mc[i]*param->U[i][2]*param->V[j][1])*D1(xli,ysq,zji);
		
		
		t+=param->U[j][2]/param->mc[i]/param->mc[i]*ST*CT*mtmt*param->V[i][2]*(f1-f2);
	 	t+=-param->U[j][2]/param->mc[i]/param->mc[i]*sq2*param->MW*param->V[i][1]*(CT*CT*f1+ST*ST*f2-fq);
		f1p=(yst1*param->U[j][2]*param->V[i][1]-param->mc[j]/param->mc[i]*param->U[i][2]*param->V[j][1])*D1(xli,yst1,zji);
		f2p=(yst2*param->U[j][2]*param->V[i][1]-param->mc[j]/param->mc[i]*param->U[i][2]*param->V[j][1])*D1(xli,yst2,zji);
		fqp=(ysq*param->U[j][2]*param->V[i][1]-param->mc[j]/param->mc[i]*param->U[i][2]*param->V[j][1])*D1(xli,ysq,zji);
		tp+=param->U[j][2]/param->mc[i]/param->mc[i]*ST*CT*mtmt*param->V[i][2]*(f1p-f2p);
	 	tp+=-param->U[j][2]/param->mc[i]/param->mc[i]*sq2*param->MW*param->V[i][1]*(CT*CT*f1p+ST*ST*f2p-fqp);
	 	}
		
	}
	/*printf("Box %.3E %.3E %.3E \n",zji,t,tp);
	printf("Box %.3E %.3E %.3E \n",f1,f2,fq);
	printf("Box %.3E %.3E %.3E \n",f1p,f2p,fqp);
	*/
	
	*Cbox_S=t*param->tb*param->tb*param->MW/sq2;	
	*Cbox_P=-tp*param->tb*param->tb*param->MW/sq2;	
	/*printf("CT ST CB SB %.3E %.3E %.3E %.3E\n",CT,ST,CB,SB);
	*/
}
 

static void Cpen_sp(struct read_param_tag* param, double *Cpen_S, double *Cpen_P,  double mtmt)
{
		
	
	int a,i,j;
	double yst1,yst2,zij,ysq;
	double g1,g2,gq,g1p,g2p,gqp;
	double t=0.;
	double tp=0.;
	double ts2=0.;
	double tp2=0.; 
	
	
	for(i=1;i<=2;i++) 
	{
	
	 for(j=1;j<=2;j++) 
	 {yst1=param->mst[1]*param->mst[1]/param->mc[j]/param->mc[j];
	  yst2=param->mst[2]*param->mst[2]/param->mc[j]/param->mc[j];
	  ysq=param->msq1*param->msq1/param->mc[j]/param->mc[j];
	  zij=param->mc[i]*param->mc[i]/param->mc[j]/param->mc[j];
         g1=(yst1*param->U[j][2]*param->V[i][1]+param->mc[i]/param->mc[j]*param->U[i][2]*param->V[j][1])
		 *D2(yst1,zij);
	  g2=(yst2*param->U[j][2]*param->V[i][1]+param->mc[i]/param->mc[j]*param->U[i][2]*param->V[j][1])
		*D2(yst2,zij);
	  gq=(ysq*param->U[j][2]*param->V[i][1]+param->mc[i]/param->mc[j]*param->U[i][2]*param->V[j][1])*D2(ysq,zij);
	  
	  	
	  t+=-param->U[j][2]*ST*CT*mtmt*param->V[i][2]*(g1-g2);
	  
	
	  t+=param->U[j][2]*sq2*param->MW*param->V[i][1]*(CT*CT*g1+ST*ST*g2-gq);
	  
	  
	  g1p=(yst1*param->U[j][2]*param->V[i][1]-param->mc[i]/param->mc[j]*param->U[i][2]*param->V[j][1])
		*D2(yst1,zij);
	  g2p=(yst2*param->U[j][2]*param->V[i][1]-param->mc[i]/param->mc[j]*param->U[i][2]*param->V[j][1])
		*D2(yst2,zij);
	  gqp=(ysq*param->U[j][2]*param->V[i][1]-param->mc[i]/param->mc[j]*param->U[i][2]*param->V[j][1])
	  	*D2(ysq,zij);
	  
	  tp+=-param->U[j][2]*ST*CT*mtmt*param->V[i][2]*(g1p-g2p);
	  tp+=param->U[j][2]*sq2*param->MW*param->V[i][1]*(CT*CT*g1p+ST*ST*g2p-gqp);
	 }
		
	}

	for(a=1;a<=2;a++) 
	{yst1=param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a];
	 yst2=param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a];
	 ts2+=param->U[a][2]/param->mc[a]*(-(sq2*param->MW*param->V[a][1]*CT*CT-mtmt*param->V[a][2]*ST*CT)*
	 D4(yst1));
	 ts2+=param->U[a][2]/param->mc[a]*(+(sq2*param->MW*param->V[a][1]*ST*ST+mtmt*param->V[a][2]*ST*CT)*
	 D4(yst2));
	 ts2+=param->U[a][2]/param->mc[a]*(sq2*param->MW*param->V[a][1]*(CT*CT-ST*ST)*D2(yst1,yst2));	 
	}

	for(a=1;a<=2;a++) 
	{yst1=param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a];
	 yst2=param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a];
	 tp2+=mtmt*mtmt*param->mu*param->U[a][2]*param->V[a][2]/param->mc[a]*D2(yst1,yst2);	 
	}

  

	*Cpen_S=param->tb*param->tb/(param->Mhc*param->Mhc-param->MW*param->MW)
		*(sq2*param->MW*t+2.*CT*ST*mtmt*param->mu*ts2);
	*Cpen_P=-param->tb*param->tb/(param->Mhc*param->Mhc-param->MW*param->MW)
		*(sq2*param->MW*tp+tp2);  
}
