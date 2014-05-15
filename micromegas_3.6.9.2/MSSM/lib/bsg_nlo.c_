/* b->s gamma in MSSM with main tanbeta NLO contributions

Complete desription of the formulae used in 
" Implementation of b-> sgamma in micrOMEGAs "


For the standard model the formulas are taken from 
SM1: A. Kagan and M. Neubert hep-ph/9805303.
with some formulas from
SM2: Chetyrkin, Misiak, Munz hep-ph/9612313
see also
P. Gambino, M. Misiak hep-ph/0104034

The BSM corrections are from the papers:
1- G.Degrassi, P.Gambino, G.F. Giudice, JHEP12(2000) 009
2- M.Ciuchini, G. Degrassi, P.Gambino, G.F.Giudice, hep-ph/9710335
*/


#include "pmodel.h"
#include "pmodel_aux.h"
#include "../../sources/micromegas.h"

int dMQcorrections=1;

/*
struct read_param_tag
{ double U[3][3], V[3][3], T[3][3], mc[3], mst[3], msq1, mglu, At, Ab;
 double B[3][3], msb[3], mss[3], mn[5], N[5][5],msnmu, tb;
  double m1, m2, mu,MW, MZ, Mt, Mb, Mhc, Vts,Vtb,Vcb,ee, sw, alphS_MZ;
};
*/ 

/* 
void calc_eps(struct read_param_tag* param, double *eps_b, double *eps_bp, double *eps_tps);
int read_prm(struct read_param_tag* param);
*/  

// Default values of B parameters
#define Br_BXenu 0.1064
#define Csl 0.546  /* See  Eq. 2,3,7 in Gambino and Giordano, 0805.0271 */
#define ckmf 0.9613 /* ckmf= abs(Vts*Vtb/Vcb)^2, ckmf=0.9613 corresponds in the Wolfenstein parametrisation to 
                       A=0.808, lamda=0.2253, rhobar=0.132, etabar=0.341 , PDG 2010 */
#define lam2 0.12 		       
#define mbs 50. /* mb/ms*/
#define alphem 0.00729735 /* 1.0/137.036*/
#define alphew 0.007819/*1.0/127.9*/



static double H2(double x, double y) { return x*log(x)/(1.0-x)/(x-y)
				+y*log(y)/(1.0-y)/(y-x);}
#ifdef DEBUG

static double fga1(double x) {double x1=x-1.0; return 
					   (7.0-5.0*x-8.0*x*x)/36.0/x1/x1/x1
					  +x*(3.0*x-2.0)/6.0/x1/x1/x1/x1*log(x);}
					  
static double fga2(double x) {double x1=x-1.0; return
	           (3.0-5.0*x)/6.0/x1/x1+(3.0*x-2.0)/3.0/x1/x1/x1*log(x);}


static double fg1(double x) {double x1=x-1.0; return
		       (2.0+5.0*x-x*x)/12.0/x1/x1/x1-x/2/x1/x1/x1/x1*log(x);}

static double fg2(double x) {double x1=x-1.0; return
		       (3.0-x)/2.0/x1/x1-log(x)/x1/x1/x1;}

static double fga3(double x) {return (1.0-x)*fga1(x)-x*fga2(x)/2.0-23.0/36.0;}

static double fg3(double x) {return (1.0-x)*fg1(x)-x*fg2(x)/2.0-1.0/3.0;}
				
/*  B0 is really B0(0, x,y) and x=m^2/q^2 with q some scale*/
static double B0(double x, double y) { return 1+.5*((x+y)/(x-y)*log(y/x)
				-log(x*y));}				
static double B1(double x, double y) { return .5*(.5+1./(1.-y/x)+log(y/x)/(1-y/x)/(1-y/x)
				-log(y));}				
#endif				
/* Dilogarithm with truncated serie order 12 is better than per mil for x<.6*/

static double Li2(double z) 
	{double z4=z*z*z*z;
 	double z8=z4*z4;
	double x=1/z;
	double x4=x*x*x*x;
	double x8=x4*x4;
	
		if(z> -1)
		{return z+z*z/4.+z*z*z/9.+z4/16.+z4*z/25.+z4*z*z/36.
	+z4*z*z*z/49.+z8/64.+z8*z/81.+z8*z*z/100.+z8*z*z*z/121.+z8*z4/144.
	+z8*z4*z/13./13.+z8*z4*z*z/14./14.+z8*z4*z*z*z/225.+z8*z8/16./16.;
		}
		else
		{return -.5*log(-x)*log(-x)-M_PI*M_PI/6.-(x+x*x/4.+x*x*x/9.+x4/16.+x4*x/25.+x4*x*x/36.
	+x4*x*x*x/49.+x8/64.+x8*x/81.+x8*x*x/100.+x8*x*x*x/121.+x8*x4/144.
	+x8*x4*x/13./13.+x8*x4*x*x/14./14.+x8*x4*x*x*x/225.+x8*x8/16./16.);
		}	
	}
/* Definitions of function, see  Degrassi et al., JHEP12(2000) 009 Eq. 2.4*/
				
static double F71(double x)
	{double x1=(x-1);
        return    x*(7.-5.*x-8.*x*x)/24./x1/x1/x1
	 +x*x*(3.*x-2.)/4./x1/x1/x1/x1*log(x);
	 }

static double F81(double x)
	{double x1=(x-1);
        return x*(2.+5.*x-x*x)/8./x1/x1/x1-3.*x*x/4./x1/x1/x1/x1*log(x);
	}		
static double F72(double x)
	{double x1=(x-1);
       return x*(3.-5.*x)/12./x1/x1+x*(3.*x-2.)/6./x1/x1/x1*log(x);
       }
static double F82(double x)
	{double x1=(x-1);
       return x*(3.-x)/4./x1/x1-x/2./x1/x1/x1*log(x);
	}     	

static double E7(double x) 
	{ double x1=(x-1);
	return x*(-18+11*x+x*x)/12./(x1*x1*x1)+
	x*x*(15.-16.*x+4.*x*x)*log(x)/6./(x1*x1*x1*x1)-2*log(x)/3.;}
	
static double G71(double y)
 	{double y1=(y-1);
	return 4*(-3+7*y-2*y*y)*Li2(1.-1./y)+(8.-14.*y-3.*y*y)*log(y)*log(y)/y1
	+2*(-3-y+12*y*y-2*y*y*y)*log(y)/y1+3.*(7.-13.*y+2*y*y);
	}
	
static double G72(double y)
 	{double y1=(y-1);	
	return y*(18.-37.*y+8.*y*y)*Li2(1.-1./y)+
	y*(-14.+23.*y+3.*y*y)*log(y)*log(y)/y1+
	(-50+251*y-174*y*y-192*y*y*y+21*y*y*y*y)*log(y)/9./(y-1.)
	+(797.-5436.*y+7569.*y*y-1202*y*y*y)/108.;
	}
	
static double G7H(double y, double tb)
	{double y1=(y-1);
	return -4.*y*G71(y)/9./(y1*y1*y1)+2*y*G72(y)/(y1*y1*y1*y1)/9./tb/tb;
	}

static double D71(double y)
 	{double y1=(y-1);	
	return y1*(21.-47.*y+8.*y*y)+2.*(-8.+14.*y+3.*y*y)*log(y);
	}

static double D72(double y)
 	{double y1=(y-1);	
	return	y1*(-31.-18.*y+135.*y*y-14.*y*y*y)/6.
	+y*(14.-23.*y-3.*y*y)*log(y);
	}
	
static double D7H(double y, double tb)
	{double y1=(y-1);
	return	2*y/9./(y1*y1*y1*y1)*(-D71(y)+D72(y)/y1/tb/tb);
	}
	
static double EH(double y, double tb)
	{double y1=(y-1);
	return	y*(y1*(16.-29.*y+7.*y*y)+
	6.*(3.*y-2.)*log(y))/36./(y1*y1*y1*y1)/tb/tb;
	}
	
static double G81(double y)
 	{double y1=(y-1);
	return .5*(-36+25*y-17*y*y)*Li2(1.-1./y)+
	(19+17*y)*log(y)*log(y)/y1+.25*(-3.-187.*y+12.*y*y-14.*y*y*y)*log(y)/y1
	+3.*(143.-44.*y+29.*y*y)/8.;
	}


static double G82(double y)
 	{double y1=(y-1);
	return y*(30.-17.*y+13.*y*y)*Li2(1.-1./y)
	-y*(31.+17.*y)*log(y)*log(y)/(y1)+
	(-226.+817.*y+1353.*y*y+318.*y*y*y+42.*y*y*y*y)*log(y)/36./y1
	+(1130.-18153.*y+7650.*y*y-4451.*y*y*y)/216.;
	}


static double G8H(double y, double tb)
	{double y1=(y-1);
	return -y*G81(y)/3./(y1*y1*y1)+ y*G82(y)/(y1*y1*y1*y1)/6./tb/tb;
	}

static double D81(double y)
 	{double y1=(y-1);
	return y1*(81.-16.*y+7.*y*y)-2.*(19.+17.*y)*log(y);
	}
	
static double D82(double y)
 	{double y1=(y-1);
	return y1*(-38.-261.*y+18.*y*y-7.*y*y*y)/6.+y*(31.+17.*y)*log(y);
	}
	
static double D8H(double y, double tb)
	{double y1=(y-1);
	return (y/6./(y1*y1*y1*y1))*(-D81(y)+D82(y)/y1/tb/tb);	
	}


static void c78_su(struct read_param_tag* param, double *, double *, double, double);
static	double  yt, yb;



#ifdef DEBUG
static double alpssu, alph3 = 0.12;
#endif

static	double sq2, b, sb, cb;

static double f(double z)
	{return 1.-8.*z+8.*z*z*z-z*z*z*z-12*z*z*log(z);}
	
static double f77(double z)
	{return (10.+z-2.*z*z/3.+(z-4.)*log(z))*z/3.;}
	
static double f78(double z)
	{return 8*(Li2(1.-z)-M_PI*M_PI/6.-z*log(z)+9.*z/4.-z*z/4.+z*z*z/12.)/9.;
	}

static double f88(double z, double b)
	{return (4.*Li2(1.-z)-2.*M_PI*M_PI/3.-z*(2.+z)*log(z)+8.*log(1.-z)
	+7.*z+3.*z*z-2.*z*z*z/3.-2.*(2.*z+z*z+4.*log(1-z))*log(b))/27.;
	}

/* For f22 and f27 approximation only valid for delta=.9*/	
static double f22(double z)
	{return (.107636-.208484*sqrt(z)-.156146*z);
	}
	
static double f27(double z)
	{return (-.190805+.948865*sqrt(z)-.787805*z);
	}	
	
		
static double alps(double mu, double mz, double alpsmz)
{	double v=1.-23.*alpsmz*log(mz/mu)/6./M_PI;
	return alpsmz*(1-116.*alpsmz*log(v)/v/4./M_PI/23.)/v;
}	

static double alps6(double mu, double mz, double alpsmz)
{
	return alpsmz/(1.+7./2./M_PI*alpsmz*log(mu/mz));
}


double deltaMb(void)
{
	double xtw,ytw;
	double beta0,beta1,g0m,g1m,mb;
		
	double alpshp;
	double alpssu,mtp,mtmt,mtmw,mtsu,ytmt;

		
	double alpsmt,alpsmw,alpsmz,delmt,mtmw2,eta_wmu;
	double c70smw,c70Hw,c80smw,c80Hw;
	double mub,mubar,del;

#ifdef DEBUG
    double eps_b,eps_bp,eps_tps;
    struct read_param_tag param;

	double al3Wb=0.548; /* alpha_s(MW)/alpha_s(Mb) */
    double x,y;
	double br,br_old;
	int err=0;
	double phsp,lamsl, r=0.316;     /* mc/mb */
	double z0=.0841; /* z0=.0841 (mc/mb)^2=.29^2 Ratio of pole masses*/
//	double z1=.22*.22;/*z1=.0484; *(mc/mb)^2=.22^2*  MSbar running mc(mb) */
    double z1=.29*.29; /* Change the value to reproduce NNLO result */
	double c7_epsb,c8_epsb,c7_epsbh,c8_epsbh;
	double c7_chi_mw, c7_su, c8_chi_mw, c8_su;
	double c2_sum,c7_hp_mw, c8_hp_mw, eta_h;
    double eta_s;
    double c70_w,c80_w,c20,c2_tum
	double c71Hw,c81Hw,c71_w_sm,c81_w_sm;
	double c80_mu,c70_mu,c71_w,c81_w,c71_sum,c81_tum,c71_wmu;
	double c7em_mu,cem8_mu,cem2_mu;
    double kap,rer2,rer8,r7,
	double K77,K27,K78,K22,K88,K28,kslem_mub;
	double kap1,ga77,g27,g87,S_mub,facbsg,bsg_lo,Knlo,bsg_nlo;
    double c70suw,c80suw,c78_sm_tbcorr,c78_hp_tbcorr;
    double c7_hp, c8_hp,c7_mw_old, c8_mw_old;
    double  c7_mw, c7_mb, c8_mw;
	double test1,test2,test3;
	double B22,B27,B77,B28,B78,B88,D27,D28,D77,D88,D78,D87;
#endif
	double eps_b, eps_bp, eps_tps;
	struct read_param_tag param;
	
if(!dMQcorrections) return 0;	
	if(read_prm(&param))
	{
		puts("bsgammanew: can not read parameters.");
		return 0.0;
	}

        alpsmz=param.alphS_MZ;

#ifdef DEBUG	
	dump_prm(&param);
#endif	
	sq2=sqrt(2.0);
	b=atan(param.tb);
	sb=sin(b);
	cb=cos(b);
//	mb=MbPole;
	mb=findValW("MbMb");
//	 printf("MbPole =%.3e\n",MbPole);
	mtp=findValW("Mtp");
	beta0=23./3.;
	beta1=116./3.;
	g0m=8.;
	g1m=404./3.-200./9.;
	del=.9;
/* Set the renormalisation scale mub, this value can be changed between Mbpole/2 to 2 MbPole*/
	mub=MbPole;
/* Set the renormalisation scale mubar, this value is the one relevant for semileptonic B decays, can be changed between Mbpole/2 to 2 MbPole*/
	mubar=MbPole;
	

/* Evolution of alpha_s with 6 flavors*/		
	alpssu=alps6(fabs(param.mglu),param.MZ,alpsmz);
	alpshp=alps6(param.Mhc,param.MZ,alpsmz);
		
/* Evolution of alpha_s with 5 flavors*/	
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
	mtsu=mtmw*pow(alpssu/alpsmt,4.0/7.0)*pow(alpsmt/alpsmz,12.0/23.0)
	      /sqrt(1.+9.*ytmt*ytmt/(8.*M_PI*alpsmt)*
	      (pow(alpssu/alpsmt,1.0/7.0)-1.));
	yt=mtsu/sq2/param.MW*param.ee/param.sw/sb;
	
	eta_wmu=alpsmw/alps(mub,param.MZ,alpsmz);
	

/* Standard model and Higgs at Leading order */
	
	c70smw=F71(xtw);
	c70Hw=F71(ytw)/3/param.tb/param.tb+F72(ytw);
	c80smw=F81(xtw);
	c80Hw=F81(ytw)/3./param.tb/param.tb+F82(ytw);
/*-----------------------------------------------------------*/

   calc_eps(&param, &eps_b, &eps_bp, &eps_tps);

   return  eps_b*param.tb;
}

double deltaMs(void)
{   
   struct read_param_tag param;
   double q,alpssu,eps_s,alpsmz;
   double mglu, mLMSG,mRMSG;
   double MSdL=findValW("MSsL");   
   double MSuL=findValW("MScL");
if(!dMQcorrections) return 0;   
   read_prm(&param);
    alpsmz=param.alphS_MZ;
   mglu=param.mglu;     
   q=fabs(mglu);
   alpssu=alps6(q,param.MZ,alpsmz);
   mLMSG=MSdL/mglu;
   mRMSG=findValW("MSsR")/mglu;
     eps_s=-2.0/3.0*alpssu/M_PI*(param.mu-  findValW("Ad")/param.tb)/mglu*
   H2(mLMSG*mLMSG, mRMSG*mRMSG);
   if(fabs(param.mu*param.mu/param.m2/param.m2-1)<.01)
   {eps_s+=alphew/param.sw/param.sw/4.0/M_PI*(-param.mu*param.m2)*
    (MSuL*MSuL/(param.m2*param.m2-MSuL*MSuL)/(param.mu*param.mu-MSuL*MSuL)*
    (log(MSuL*MSuL/param.m2/param.m2)+(param.m2*param.m2/MSuL/MSuL-1.0))
    +MSdL*MSdL/(param.m2*param.m2-MSdL*MSdL)/(param.mu*param.mu-MSdL*MSdL)*
    (log(MSdL*MSdL/param.m2/param.m2)+(param.m2*param.m2/MSdL/MSdL-1.0))/2.);
   }
   else
   { 
     eps_s+=alphew/param.sw/param.sw/4.0/M_PI*(param.mu*param.m2)*(
      H2(param.m2*param.m2/MSuL/MSuL,param.mu*param.mu/MSuL/MSuL)/MSuL/MSuL
     +H2(param.m2*param.m2/MSdL/MSdL,param.mu*param.mu/MSdL/MSdL)/MSdL/MSdL/2.);
   }
      
   return  eps_s*param.tb;    
}

double deltaMd(void)
{   
   struct read_param_tag param;
   double q,alpssu,eps_s,alpsmz;
   double mglu, mLMSG,mRMSG;
   double MSdL=findValW("MSdL");   
   double MSuL=findValW("MSuL");
   
if(!dMQcorrections) return 0;  
   read_prm(&param);
    alpsmz=param.alphS_MZ;
   mglu=param.mglu;     
   q=fabs(mglu);
   alpssu=alps6(q,param.MZ,alpsmz);
   mLMSG=findValW("MSdL")/mglu;
   mRMSG=findValW("MSdR")/mglu;
     eps_s=-2.0/3.0*alpssu/M_PI*(param.mu-  findValW("Ad")/param.tb)/mglu*
   H2(mLMSG*mLMSG, mRMSG*mRMSG);
 
    if(fabs(param.mu*param.mu/param.m2/param.m2-1)<.01)
   {eps_s+=alphew/param.sw/param.sw/4.0/M_PI*(-param.mu*param.m2)*
    (MSuL*MSuL/(param.m2*param.m2-MSuL*MSuL)/(param.mu*param.mu-MSuL*MSuL)*
    (log(MSuL*MSuL/param.m2/param.m2)+(param.m2*param.m2/MSuL/MSuL-1.0))
    +MSdL*MSdL/(param.m2*param.m2-MSdL*MSdL)/(param.mu*param.mu-MSdL*MSdL)*
    (log(MSdL*MSdL/param.m2/param.m2)+(param.m2*param.m2/MSdL/MSdL-1.0))/2.);
   }
   else
   { 
     eps_s+=alphew/param.sw/param.sw/4.0/M_PI*(param.mu*param.m2)*(
      H2(param.m2*param.m2/MSuL/MSuL,param.mu*param.mu/MSuL/MSuL)/MSuL/MSuL
     +H2(param.m2*param.m2/MSdL/MSdL,param.mu*param.mu/MSdL/MSdL)/MSdL/MSdL/2.);
   }
     
   return  eps_s*param.tb;    
}



    
double bsgnlo_(double * SMbsg)
{
	double x,xtw,ytw;
#ifdef DEBUG
	double al3Wb=0.548; /* alpha_s(MW)/alpha_s(Mb) */
    double y;
	double br,br_old;
	int err=0;
	double phsp,lamsl, r=0.316;     /* mc/mb */
#endif	
	double z0=.0841; /* z0=.0841 (mc/mb)^2=.29^2*/
//	double z1=.22*.22;/*z1=.0484; *(mc/mb)^2=.22^2*/
    double z1=.29*.29;	
	double beta0,beta1,g0m,g1m,mb,mcPole;
		
	double c7_epsb,c8_epsb,c7_epsbh,c8_epsbh;
	double c7_chi_mw, c7_su, c8_chi_mw, c8_su;
	double c7_hp_mw, c8_hp_mw, eta_h,alpshp,alpsmz;
	double c2_sum;
	double alpssu,mtp,mtmt,mtmw,mtsu,ytmt,eta_s;
    double bsg_nlo_sm;
		
	double alpsmt,alpsmw,delmt,mtmw2,eta_wmu;
	double c70smw,c70Hw,c70_w,c80smw,c80Hw,c80_w,c20,c2_tum;
	double c71Hw,c81Hw,c71_w_sm,c81_w_sm,c71sm_w,c81sm_w,c7emsm_mu;
	double c70sm_mu,c80sm_mu,c71sm_wmu,bsg_lo_sm;
	double c80_mu,c70_mu,c71_w,c81_w,c71_sum,c81_tum,c71_wmu;
	double c7em_mu,cem8_mu,cem2_mu;
	double kap,rer2,rer8,r7,mub,mubar,del;
	double K77,K27,K78,K22,K88,K28,kslem_mub;
	double kap1,ga77,g27,g87,S_mub,facbsg,bsg_lo,Knlo,bsg_nlo;
#ifdef DEBUG
    double c70suw,c80suw,c78_sm_tbcorr,c78_hp_tbcorr;
    double c7_hp, c8_hp,c7_mw_old, c8_mw_old;
    double  c7_mw, c7_mb, c8_mw;
	double test1,test2,test3;
	double B22,B27,B77,B28,B78,B88,D27,D28,D77,D88,D78,D87;
#endif
	double eps_b, eps_bp, eps_tps;
	struct read_param_tag param;
	
	if(read_prm(&param))
	{
		puts("bsgammanew: can not read parameters.");
		return 0.0;
	}
        alpsmz=param.alphS_MZ;
#ifdef DEBUG	
	dump_prm(&param);
#endif	
	sq2=sqrt(2.0);
	b=atan(param.tb);
	sb=sin(b);
	cb=cos(b);
//	mb=MbPole;
    mcPole=sqrt(z0)*MbPole;
	mb=findValW("MbMb");
	mtp=findValW("Mtp");
	beta0=23./3.;
	beta1=116./3.;
	g0m=8.;
	g1m=404./3.-200./9.;
	del=.9;
/* Set the renormalisation scale mub, this value can be changed between Mbpole/2 to 2 MbPole*/
	mub=MbPole;
/* Set the renormalisation scale mubar, this value is the one relevant for semileptonic B decays, can be changed between Mbpole/2 to 2 MbPole*/
	mubar=MbPole;
	

/* Evolution of alpha_s with 6 flavors*/		
	alpssu=alps6(fabs(param.mglu),param.MZ,alpsmz);
	alpshp=alps6(param.Mhc,param.MZ,alpsmz);
		
/* Evolution of alpha_s with 5 flavors*/	
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
	mtsu=mtmw*pow(alpssu/alpsmt,4.0/7.0)*pow(alpsmt/alpsmz,12.0/23.0)
	      /sqrt(1.+9.*ytmt*ytmt/(8.*M_PI*alpsmt)*
	      (pow(alpssu/alpsmt,1.0/7.0)-1.));
	yt=mtsu/sq2/param.MW*param.ee/param.sw/sb;
	
	eta_wmu=alpsmw/alps(mub,param.MZ,alpsmz);


/* Standard model and Higgs at Leading order */
	
	c70smw=F71(xtw);
	c70Hw=F71(ytw)/3/param.tb/param.tb+F72(ytw);
	c80smw=F81(xtw);
	c80Hw=F81(ytw)/3./param.tb/param.tb+F82(ytw);

/* Calculates the eps_b eps_b' ..*/	
	calc_eps(&param, &eps_b, &eps_bp, &eps_tps);

/* Large tan beta corrections to SM and Higgs*/	
	c7_epsb=(eps_b-eps_bp)/(1.+eps_b*param.tb)*param.tb*F72(xtw);
	c7_epsbh=-(eps_tps+eps_b)*param.tb/(1.+eps_b*param.tb)*F72(ytw);
	c8_epsb=(eps_b-eps_bp)/(1.+eps_b*param.tb)*param.tb*F82(xtw);
	c8_epsbh=-(eps_tps+eps_b)*param.tb/(1.+eps_b*param.tb)*F82(ytw);
	

	c70Hw+=c7_epsbh;
	c80Hw+=c8_epsbh;

/* Running from mh+ to Mw*/	
	eta_h=alpshp/alpsmz;
		
	c7_hp_mw=pow(eta_h,16.0/21.0)*c70Hw+8./3.*(
	pow(eta_h,2.0/3.0)-pow(eta_h,16.0/21.0))*c80Hw;
	c8_hp_mw=pow(eta_h,2.0/3.0)*c80Hw;
	
/* Leading order susy includes already large tan beta factors*/
	c7_su=0.0;
	c8_su=0.0;
	c78_su(&param, &c7_su, &c8_su, eps_b, mtsu);
	
	eta_s=alpssu/alpsmz;
	
/* Chargino contribution at mw scale , see  Degrassi et al., JHEP12(2000) 009  Eq.3.5*/
		
	c7_chi_mw=pow(eta_s,16.0/21.0)*c7_su+8./3.*(
	pow(eta_s,2.0/3.0)-pow(eta_s,16.0/21.0))*c8_su;
	c8_chi_mw=pow(eta_s,2.0/3.0)*c8_su;

	c70_w=c70smw+c7_epsb+c70Hw+c7_chi_mw;
	c80_w=c80smw+c8_epsb+c80Hw+c8_chi_mw;
		
	c20=(pow(eta_wmu,(-12./23.))+pow(eta_wmu,(6./23.)))/2.;
		
	c2_sum= 626126.0/272277.0*pow(eta_wmu,14.0/23.0)
			-56281.0/51730.0 *pow(eta_wmu,16.0/23.0)
			-3.0/7.0         *pow(eta_wmu,6.0/23.0)
			-1.0/14.0        *pow(eta_wmu,-12.0/23.0)
			-0.6494          *pow(eta_wmu,0.4086)
			-0.0380          *pow(eta_wmu,-0.423)
			-0.0186          *pow(eta_wmu,-0.8994)
			-0.0057          *pow(eta_wmu,0.1456);	
	
	c2_tum= -.9135*pow(eta_wmu,.4086)
			+.0873 *pow(eta_wmu,(-.4230))
			-.0571 *pow(eta_wmu,(-.8994))
			+.0209 *pow(eta_wmu,.1456);
			
/* Runnning from mw to mub  : SM+HIGGS+SUSY*/
	c80_mu=(c80_w+313063./363036.)*pow(eta_wmu,(14./23.))+c2_tum;
	c70_mu=c70_w*pow(eta_wmu,(16./23.))+8.*(pow(eta_wmu,(14./23.))								
	-pow(eta_wmu,(16./23.)))*c80_w/3.+c2_sum;
	
	

/* Runnning from mw to mub  : SM only */
	c80sm_mu=(c80smw+313063./363036.)*pow(eta_wmu,(14./23.))+c2_tum;
	c70sm_mu=c70smw*pow(eta_wmu,(16./23.))+8.*(pow(eta_wmu,(14./23.))								
	-pow(eta_wmu,(16./23.)))*c80smw/3.+c2_sum;
		
/* NLO Standard model*/	
	
	x=xtw;			
	c71_w=(-16.*x*x*x*x-122.*x*x*x+80.*x*x-8.*x)*Li2(1.-1./x)/9./
		pow((x-1.),4.)				
		+(6.*x*x*x*x+46.*x*x*x-28*x*x)*log(x)*log(x)/3./pow((x-1.),5.);
	c71_w+=(-102.*pow(x,5.)-588.*x*x*x*x-2262.*x*x*x+3244.*x*x-1364.*x+208.)
	*log(x)/81./pow((x-1),5.);
	c71_w+=(1646.*x*x*x*x+12205.*x*x*x-10740.*x*x+2509.*x-436.)/486./
	pow((x-1.),4.);

	c71_w_sm=c71_w;

	c81_w=(-4*x*x*x*x+40.*x*x*x+41.*x*x+x)*Li2(1.-1./x)/6./pow((x-1.),4.)+
	(-17.*x*x*x-31.*x*x)*log(x)*log(x)/2./pow((x-1.),5);
	c81_w+=(-210.*pow(x,5.)+1086.*x*x*x*x+4893.*x*x*x+2857.*x*x-1994.*x
	+280.)*log(x)/216./pow((x-1.),5.);
	c81_w+=(737.*x*x*x*x-14102.*x*x*x-28209.*x*x+610.*x-508.)/1296./
		pow((x-1.),4.);
	c81_w_sm=c81_w;

/*NLO Higgs**/

	c71Hw=G7H(ytw,param.tb)+D7H(ytw,param.tb)*log(param.MW*param.MW/param.Mhc/param.Mhc)-4*EH(ytw,param.tb)/9;
	c81Hw=G8H(ytw,param.tb)+D8H(ytw,param.tb)*log(param.MW*param.MW/param.Mhc/param.Mhc)-EH(ytw,param.tb)/6;

/* Sum of SM +Higgs*/
	c71_w+=c71Hw;
	c81_w+=c81Hw;
	
	c71_sum= (4661194/816831.*eta_wmu*E7(xtw)-17.3023+14.8088*eta_wmu)*pow(eta_wmu,14.0/23.0)
		+(-8516/2217.*eta_wmu*E7(xtw)+8.5027-10.8090*eta_wmu)*pow(eta_wmu,16.0/23.0)
+(4.5508-.8740*eta_wmu)*pow(eta_wmu,6.0/23.0)
+(0.7519+.4218*eta_wmu)*pow(eta_wmu,-12.0/23.0)
+(-1.9043*eta_wmu*E7(xtw)+2.0040-2.9347*eta_wmu)*pow(eta_wmu,0.4086)
+(-.1008*eta_wmu*E7(xtw)+.7476+.3971*eta_wmu)*pow(eta_wmu,-0.423)
+(.1216*eta_wmu*E7(xtw)-.5385+.1600*eta_wmu)*pow(eta_wmu,-0.8994)
+(.0183*eta_wmu*E7(xtw)+.0914+0.0225*eta_wmu)*pow(eta_wmu,0.1456);	

	c81_tum= 297664/14283.*pow(eta_wmu,16./23.)
		 -7164416/357075. *pow(eta_wmu,14./23.)
		 +256868/14283.*pow(eta_wmu,37./23.)
		 -6698884/357075. *pow(eta_wmu,39./23.);
		 
/* SM, running from Mw to mub*/

	c71sm_wmu=c71_w_sm*pow(eta_wmu,39/23.)+
		8.*(pow(eta_wmu,37/23.)-pow(eta_wmu,39/23.))*c81_w/3.
		+c71_sum+c81_tum*c80smw+
		37208.*(pow(eta_wmu,39/23.)-pow(eta_wmu,16/23.))*c70smw/4761.;
/* end SM*/

/* SM+Higgs, running from Mw to mub*/

	c71_wmu=c71_w*pow(eta_wmu,39/23.)+
		8.*(pow(eta_wmu,37/23.)-pow(eta_wmu,39/23.))*c81_w/3.
		+c71_sum+c81_tum*c80_w+
		37208.*(pow(eta_wmu,39/23.)-pow(eta_wmu,16/23.))*c70_w/4761.;
	
	
	
/*e-m corrections Formula 11 in Kagan-Neubert hep-ph/9805303, here c70_w and
c80_w contains SM+HIGGS+SUSY*/

	cem2_mu= -190/8073.*pow(eta_wmu,-35/23.)
			-359/3105. *pow(eta_wmu,-17./23.)
			+4276/121095. *pow(eta_wmu,-12./23.)
			+350531/1009125. *pow(eta_wmu,-9./23.)
			+2/4347. *pow(eta_wmu,-7./23.)
			-5956/15525.   *pow(eta_wmu,6./23.)
			+38380/169533. *pow(eta_wmu,14./23.)
			-748/8625. *pow(eta_wmu,16./23.);
	cem8_mu= -32/575.*pow(eta_wmu,-9/23.)
			+32/1449. *pow(eta_wmu,-7./23.)
			+640/1449. *pow(eta_wmu,14./23.)
			-704/1725. *pow(eta_wmu,16./23.);
	c7em_mu=(32/75.*pow(eta_wmu,-9/23.)-40/69.*pow(eta_wmu,-7/23.)
	+88/575.*pow(eta_wmu,16/23.))*c70_w+cem8_mu*c80_w+cem2_mu;

/* sm only*/	
	c7emsm_mu=(32/75.*pow(eta_wmu,-9/23.)-40/69.*pow(eta_wmu,-7/23.)
	+88/575.*pow(eta_wmu,16/23.))*c70smw+cem8_mu*c80smw+cem2_mu;
	
	

/**  The  K_ij functions */
	
	kap=3.382-4.14*(sqrt(z0)-.29);
	kap1=3.672-4.14*(sqrt(z1)-.22);
	rer2=-4.987+12.78*(sqrt(z1)-.22);
	rer8=44./9.-8.*M_PI*M_PI/27.;
	r7=-10./3.-8.*M_PI*M_PI/9.;
	

	ga77=32./3.;
	g27=416./81.;
	g87=-32./9.;
	S_mub=exp(-2.*alps(mub,param.MZ,alpsmz)*log(del)*(log(del)+7./2.)/3./M_PI);


	K77=S_mub*(1+alps(mub,param.MZ,alpsmz)*(r7+ga77*log(MbPole/mub)-16/3.)/2./M_PI
	 +6.*(pow((1-z0),4)/f(z0)-1)*lam2/MbPole/MbPole)+alps(mub,param.MZ,alpsmz)*
	 f77(del)/M_PI+S_mub*alps(mubar,param.MZ,alpsmz)*kap1/M_PI/2.;

/*	K27=S_mub*(alps(mub,MZ,alpsmz)*(rer2+g27*log(mb/mub))/2./M_PI
	 -lam2/9./z0/mb/mb)+  alps(mub,MZ,alpsmz)*f27(z0,del)/M_PI;
*/

/* In K27,K28,K88, the numerical values correspond to delta=.9 and z0=.29^2*/
	
	K27=S_mub*(alps(mub,param.MZ,alpsmz)*(rer2+g27*log(MbPole/mub))/2./M_PI
	 -lam2/9./mcPole/mcPole)+  alps(mub,param.MZ,alpsmz)*(f27(z1))/M_PI;

	K78=S_mub*alps(mub,param.MZ,alpsmz)*(rer8+g87*log(MbPole/mub))/2./M_PI
		+ alps(mub,param.MZ,alpsmz)*f78(del)/M_PI;
		
/*	K22=alps(mub,MZ,alpsmz)*f22(z0,del)/M_PI;*/
	K22=alps(mub,param.MZ,alpsmz)*f22(z1)/M_PI;
	K88=alps(mub,param.MZ,alpsmz)*f88(del,mbs)/M_PI;
	K28=alps(mub,param.MZ,alpsmz)*(-f27(z1))/3./M_PI;

	kslem_mub=2.*alps(mub,param.MZ,alpsmz)*log(param.MW/mub)/M_PI;
	
	facbsg=6.*Br_BXenu*alphem*ckmf/M_PI/Csl;
/* SM only*/
    bsg_lo_sm=facbsg*(c70sm_mu*c70sm_mu);
    Knlo=K22*c20*c20+ K77*c70sm_mu*c70sm_mu+K88*c80sm_mu*c80sm_mu+K27*c20*c70sm_mu 
	 +K28*c20*c80sm_mu+K78*c70sm_mu*c80sm_mu
	 +S_mub*alps(mub,param.MZ,alpsmz)*c71sm_wmu*c70sm_mu/2./M_PI
	 +S_mub*alphem*(2*c7em_mu*c70sm_mu-c70sm_mu*c70sm_mu*kslem_mub)/
	 alps(mub,param.MZ,alpsmz);
	 if(SMbsg) *SMbsg=facbsg*Knlo;	
	 
//	 printf("SM: bsg_lo=%.5e Knlo=%.5e bsg_nlo=%.4e\n",bsg_lo_sm,Knlo,bsg_nlo);
	 

/* SM+Higgs+SUSY */

	bsg_lo=facbsg*(c70_mu*c70_mu);
    Knlo=K22*c20*c20+ K77*c70_mu*c70_mu+K88*c80_mu*c80_mu+K27*c20*c70_mu 
	 +K28*c20*c80_mu+K78*c70_mu*c80_mu
	 +S_mub*alps(mub,param.MZ,alpsmz)*c71_wmu*c70_mu/2./M_PI
	 +S_mub*alphem*(2*c7em_mu*c70_mu-c70_mu*c70_mu*kslem_mub)/
	 alps(mub,param.MZ,alpsmz);
	 
	 
//	printf("bsg_lo=%.5e Knlo=%.5e\n",bsg_lo,Knlo);
//	printf("c70=%.5e c80=%.5e\n",c70_w,c80_w);
//	printf("c71=%.5e c81_sm=%.5e c81=%.5e\n",c71_w,c81_w_sm,c81_w);
	 
	bsg_nlo=facbsg*Knlo;	
	
	return bsg_nlo; 
}





static double F73(double x)
	{double x1=(x-1);
     return (5.-7.*x)/6./x1/x1+x*(3.*x-2.)/3./x1/x1/x1*log(x);}
static double F83(double x)
	{double x1=(x-1);
     return (1.+x)/2./x1/x1-x/x1/x1/x1*log(x);}
     
#define CT (param->T[1][1])
#define ST (param->T[1][2])
#define CB (param->B[1][1])
#define SB (param->B[1][2])
	

void calc_eps(struct read_param_tag* param, double *eps_b, double *eps_bp, double *eps_tps)
{
	int i;
#ifdef DEBUG 
        double test1,test2;
#endif        	
	double att,q,alpssu,alpsmz;
	att=0.;
	q=fabs(param->mglu);
	alpsmz=param->alphS_MZ;
	alpssu=alps6(q,param->MZ,alpsmz);
	
	
/* Add the full corrections to eps_b, does not include the terms subdominant in tanbeta*/			
		
/*	eps_b=-2.0/3.0*alpssu/M_PI*(mu-Ab/tb)/mglu*
			H2(msb[1]*msb[1]/mglu/mglu,msb[2]*msb[2]/mglu/mglu);
*/
	*eps_b=-2.0/3.0*alpssu/M_PI*(param->mu-param->Ab/param->tb)/param->mglu*
			H2(param->msb[1]*param->msb[1]/param->mglu/param->mglu,param->msb[2]*param->msb[2]/param->mglu/param->mglu);
/*	eps_b+=-yt*yt/16.0/M_PI/M_PI*(At/mu-1./tb)*
			H2(mst[1]*mst[1]/mu/mu,mst[2]*mst[2]/mu/mu);
*/	

*eps_b+=-yt*yt/16.0/M_PI/M_PI*(param->At-param->mu/param->tb)*
(param->U[1][2]*param->V[1][2]/param->mc[1]
*H2(param->mst[1]*param->mst[1]/param->mc[1]/param->mc[1],param->mst[2]*param->mst[2]/param->mc[1]/param->mc[1])
+param->U[2][2]*param->V[2][2]/param->mc[2]
*H2(param->mst[1]*param->mst[1]/param->mc[2]/param->mc[2],param->mst[2]*param->mst[2]/param->mc[2]/param->mc[2]));

	
/*	*eps_b-=1.0/3.0*alpssu/M_PI/param->tb*(B1(param->mglu*param->mglu/q/q,param->msb[1]*param->msb[1]/q/q)+	
		B1(param->mglu*param->mglu/q/q,param->msb[2]*param->msb[2]/q/q));
*/
		
	if(fabs(param->mu*param->mu/param->m2/param->m2-1)<.01)
	{*eps_b+=alphew/param->sw/param->sw/4.0/M_PI*(-param->mu*param->m2)*
	(CT*CT*param->mst[1]*param->mst[1]/(param->m2*param->m2-param->mst[1]*param->mst[1])/(param->mu*param->mu-param->mst[1]*param->mst[1])*
	(log(param->mst[1]*param->mst[1]/param->m2/param->m2)+(param->m2*param->m2/param->mst[1]/param->mst[1]-1.0))
	+ST*ST*param->mst[2]*param->mst[2]/(param->m2*param->m2-param->mst[2]*param->mst[2])/(param->mu*param->mu-param->mst[2]*param->mst[2])*
	(log(param->mst[2]*param->mst[2]/param->m2/param->m2)+(param->m2*param->m2/param->mst[2]/param->mst[2]-1.0))
	+CB*CB*param->msb[1]*param->msb[1]/(param->m2*param->m2-param->msb[1]*param->msb[1])/(param->mu*param->mu-param->msb[1]*param->msb[1])*
	(log(param->msb[1]*param->msb[1]/param->m2/param->m2)+(param->m2*param->m2/param->msb[1]/param->msb[1]-1.0))/2.
	+SB*SB*param->msb[2]*param->msb[2]/(param->m2*param->m2-param->msb[2]*param->msb[2])/(param->mu*param->mu-param->msb[2]*param->msb[2])*
	(log(param->msb[2]*param->msb[2]/param->m2/param->m2)+(param->m2*param->m2/param->msb[2]/param->msb[2]-1.0))/2.);	
	}
	else
	{
	*eps_b+=alphew/param->sw/param->sw/4.0/M_PI*(param->mu*param->m2)*(
	(CT*CT*H2(param->m2*param->m2/param->mst[1]/param->mst[1],param->mu*param->mu/param->mst[1]/param->mst[1])/param->mst[1]/param->mst[1]
	+ST*ST*H2(param->m2*param->m2/param->mst[2]/param->mst[2],param->mu*param->mu/param->mst[2]/param->mst[2])/param->mst[2]/param->mst[2])
	+CB*CB*H2(param->m2*param->m2/param->msb[1]/param->msb[1],param->mu*param->mu/param->msb[1]/param->msb[1])/param->msb[1]/param->msb[1]/2.
	+SB*SB*H2(param->m2*param->m2/param->msb[2]/param->msb[2],param->mu*param->mu/param->msb[2]/param->msb[2])/param->msb[2]/param->msb[2]/2.);	
	}
	

/*	yb=yb/(1.0+eps_b*tb);*/
	*eps_bp=-2.0/3.0*alpssu/M_PI*(param->mu-param->Ab/param->tb)/param->mglu*(
		CT*CT*CB*CB*H2(param->mst[1]*param->mst[1]/param->mglu/param->mglu,param->msb[2]*param->msb[2]/param->mglu/param->mglu)
	   +CT*CT*SB*SB*H2(param->mst[1]*param->mst[1]/param->mglu/param->mglu,param->msb[1]*param->msb[1]/param->mglu/param->mglu)
	   +ST*ST*CB*CB*H2(param->mst[2]*param->mst[2]/param->mglu/param->mglu,param->msb[2]*param->msb[2]/param->mglu/param->mglu)
	   +ST*ST*SB*SB*H2(param->mst[2]*param->mst[2]/param->mglu/param->mglu,param->msb[1]*param->msb[1]/param->mglu/param->mglu)
			);
	
		
	for(i=1;i<=4;i++)
	*eps_bp+=yt*yt/16.0/M_PI/M_PI*param->N[i][4]*param->N[i][3]*(param->At-param->mu/param->tb)/param->mn[i]*(
	    CT*CT*CB*CB*H2(param->mst[2]*param->mst[2]/param->mn[i]/param->mn[i],param->msb[1]*param->msb[1]/param->mn[i]/param->mn[i])
	   +CT*CT*SB*SB*H2(param->mst[2]*param->mst[2]/param->mn[i]/param->mn[i],param->msb[2]*param->msb[2]/param->mn[i]/param->mn[i])
	   +ST*ST*CB*CB*H2(param->mst[1]*param->mst[1]/param->mn[i]/param->mn[i],param->msb[1]*param->msb[1]/param->mn[i]/param->mn[i])
	   +ST*ST*SB*SB*H2(param->mst[1]*param->mst[1]/param->mn[i]/param->mn[i],param->msb[2]*param->msb[2]/param->mn[i]/param->mn[i])
			);
					
		

	if(fabs(param->mu*param->mu/param->m2/param->m2-1)<.01)
	{*eps_bp+=alphew/param->sw/param->sw/4.0/M_PI*(-param->mu*param->m2)*
	(CT*CT*param->mst[1]*param->mst[1]/(param->m2*param->m2-param->mst[1]*param->mst[1])/(param->mu*param->mu-param->mst[1]*param->mst[1])*
	(log(param->mst[1]*param->mst[1]/param->m2/param->m2)+(param->m2*param->m2/param->mst[1]/param->mst[1]-1.0))/2.
	+ST*ST*param->mst[2]*param->mst[2]/(param->m2*param->m2-param->mst[2]*param->mst[2])/(param->mu*param->mu-param->mst[2]*param->mst[2])*
	(log(param->mst[2]*param->mst[2]/param->m2/param->m2)+(param->m2*param->m2/param->mst[2]/param->mst[2]-1.0))/2.
	+CB*CB*param->msb[1]*param->msb[1]/(param->m2*param->m2-param->msb[1]*param->msb[1])/(param->mu*param->mu-param->msb[1]*param->msb[1])*
	(log(param->msb[1]*param->msb[1]/param->m2/param->m2)+(param->m2*param->m2/param->msb[1]/param->msb[1]-1.0))
	+SB*SB*param->msb[2]*param->msb[2]/(param->m2*param->m2-param->msb[2]*param->msb[2])/(param->mu*param->mu-param->msb[2]*param->msb[2])*
	(log(param->msb[2]*param->msb[2]/param->m2/param->m2)+(param->m2*param->m2/param->msb[2]/param->msb[2]-1.0))
	);	
	}
	else
	{
	*eps_bp+=alphew/param->sw/param->sw/4.0/M_PI*(param->mu*param->m2)*(
	(CT*CT*H2(param->m2*param->m2/param->mst[1]/param->mst[1],param->mu*param->mu/param->mst[1]/param->mst[1])/param->mst[1]/param->mst[1]
	+ST*ST*H2(param->m2*param->m2/param->mst[2]/param->mst[2],param->mu*param->mu/param->mst[2]/param->mst[2])/param->mst[2]/param->mst[2])/2.
	+(CB*CB*H2(param->m2*param->m2/param->msb[1]/param->msb[1],param->mu*param->mu/param->msb[1]/param->msb[1])/param->msb[1]/param->msb[1]
	+SB*SB*H2(param->m2*param->m2/param->msb[2]/param->msb[2],param->mu*param->mu/param->msb[2]/param->msb[2])/param->msb[2]/param->msb[2])
	);	
	}
	
	
/*	*eps_tpb=-2.0/3.0*alpssu/M_PI*(mu+At/tb)/mglu*(
		CT*CT*CB*CB*H2(mst[2]*mst[2]/mglu/mglu,msb[1]*msb[1]/mglu/mglu)
	   +CT*CT*SB*SB*H2(mst[2]*mst[2]/mglu/mglu,msb[2]*msb[2]/mglu/mglu)
	   +ST*ST*CB*CB*H2(mst[1]*mst[1]/mglu/mglu,msb[1]*msb[1]/mglu/mglu)
	   +ST*ST*SB*SB*H2(mst[1]*mst[1]/mglu/mglu,msb[2]*msb[2]/mglu/mglu)
			);
			
	for(i=1;i<=4;i++)
		*eps_tpb+=yb*yb/16.0/M_PI/M_PI*N[i][4]*N[i][3]*(Ab+mu/tb)/mn[i]*(
	    CT*CT*CB*CB*H2(mst[1]*mst[1]/mn[i]/mn[i],msb[2]*msb[2]/mn[i]/mn[i])
	   +CT*CT*SB*SB*H2(mst[1]*mst[1]/mn[i]/mn[i],msb[1]*msb[1]/mn[i]/mn[i])
	   +ST*ST*CB*CB*H2(mst[2]*mst[2]/mn[i]/mn[i],msb[2]*msb[2]/mn[i]/mn[i])
	   +ST*ST*SB*SB*H2(mst[2]*mst[2]/mn[i]/mn[i],msb[1]*msb[1]/mn[i]/mn[i])
			);
*/
			

	*eps_tps=-2.0/3.0*alpssu/M_PI*(param->mu+param->At/param->tb)/param->mglu*(
		CT*CT*H2(param->mst[2]*param->mst[2]/param->mglu/param->mglu,param->mss[1]*param->mss[1]/param->mglu/param->mglu)
	   +ST*ST*H2(param->mst[1]*param->mst[1]/param->mglu/param->mglu,param->mss[1]*param->mss[1]/param->mglu/param->mglu)
			);
						
			

	for(i=1;i<=4;i++)
		*eps_tps+=yb*yb/16.0/M_PI/M_PI*param->N[i][4]*param->N[i][3]*(param->mu/param->tb)/param->mn[i]*(
	    CT*CT*CB*CB*H2(param->mst[1]*param->mst[1]/param->mn[i]/param->mn[i],param->msb[2]*param->msb[2]/param->mn[i]/param->mn[i])
	   +CT*CT*SB*SB*H2(param->mst[1]*param->mst[1]/param->mn[i]/param->mn[i],param->msb[1]*param->msb[1]/param->mn[i]/param->mn[i])
	   +ST*ST*CB*CB*H2(param->mst[2]*param->mst[2]/param->mn[i]/param->mn[i],param->msb[2]*param->msb[2]/param->mn[i]/param->mn[i])
	   +ST*ST*SB*SB*H2(param->mst[2]*param->mst[2]/param->mn[i]/param->mn[i],param->msb[1]*param->msb[1]/param->mn[i]/param->mn[i])
			);
					
		
}


/* Leading order chargino contribution contains the leading higher order effect assumes that sparticles have mass  1 TeV with possibly lighter charginos and Higgs*/

static void c78_su(struct read_param_tag* param, double *c7_su, double *c8_su, double eps_b, double mtsu)
{
#ifdef DEBUG		
	double t1,t2;
        double 	eta,bet0,k1;
#endif        
	int a;	
	double k,t,t8;
	k=1.0/(1.0+eps_b*param->tb);
	
	t=0.0;

	for(a=1;a<=2;a++) /* c7; JHEP12(2000)009, Eq. 4.1 */
	{
	t+=2.0/3.0*param->MW*param->MW/param->msq1/param->msq1*param->V[a][1]*param->V[a][1]*
	   F71(param->msq1*param->msq1/param->mc[a]/param->mc[a]);
	  			  
	t-=2.0/3.0*(CT*param->V[a][1]-ST*param->V[a][2]*mtsu/sq2/sb/param->MW)*
	   (CT*param->V[a][1]-ST*param->V[a][2]*mtsu/sq2/sb/param->MW)*param->MW*param->MW/param->mst[1]/param->mst[1]				   
	   *F71(param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a]);
				
	t-=2.0/3.0*(ST*param->V[a][1]+CT*param->V[a][2]*mtsu/sq2/sb/param->MW)*
	   (ST*param->V[a][1]+CT*param->V[a][2]*mtsu/sq2/sb/param->MW)*param->MW*param->MW/param->mst[2]/param->mst[2]  		          
	   *F71(param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a]);
	   
	 				
	t+=k/cb*param->U[a][2]*param->V[a][1]*param->MW/sq2/param->mc[a]*(F73(param->msq1*param->msq1/param->mc[a]/param->mc[a])
				-CT*CT*F73(param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a])
				-ST*ST*F73(param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a]));
	t+=k/cb*CT*ST*param->U[a][2]*param->V[a][2]*mtsu/2.0/sb/param->mc[a]
	   *(F73(param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a])-F73(param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a]));
	  			
	}



	*c7_su=t;	
	t8=0.0;
	for(a=1;a<=2;a++) /* c8;  JHEP12(2000)009, Eq.4.1 */
	{
	t8+=2.0/3.0*param->MW*param->MW/param->msq1/param->msq1*param->V[a][1]*param->V[a][1]*F81(param->msq1*param->msq1/param->mc[a]/param->mc[a]);	 
	t8-=2.0/3.0*(CT*param->V[a][1]-ST*param->V[a][2]*mtsu/sq2/sb/param->MW)*
	    (CT*param->V[a][1]-ST*param->V[a][2]*mtsu/sq2/sb/param->MW)*param->MW*param->MW/param->mst[1]/param->mst[1]      
				*F81(param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a]);
				
	t8-=2.0/3.0*(ST*param->V[a][1]+CT*param->V[a][2]*mtsu/sq2/sb/param->MW)
	    *(ST*param->V[a][1]+CT*param->V[a][2]*mtsu/sq2/sb/param->MW)*param->MW*param->MW/param->mst[2]/param->mst[2]
				*F81(param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a]);				
	t8+=k/cb*param->U[a][2]*param->V[a][1]*param->MW/sq2/param->mc[a]*(F83(param->msq1*param->msq1/param->mc[a]/param->mc[a])
				-CT*CT*F83(param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a])
				-ST*ST*F83(param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a]));
	
	t8+=k/cb*ST*CT*param->U[a][2]*param->V[a][2]*mtsu/2.0/sb/param->mc[a]
				*(F83(param->mst[1]*param->mst[1]/param->mc[a]/param->mc[a])-
				F83(param->mst[2]*param->mst[2]/param->mc[a]/param->mc[a]));
	}
	*c8_su=t8;
}

int read_prm(struct read_param_tag* param)
{
	int err=0;
	
	/*param->MW=80.423;
	param->MZ=91.1876;
        param->ee=sqrt(4*M_PI*0.00781653);
	param->sw=0.48076;
   */ 
    err+=findVal("MZ",&param->MZ);
    { double alfEMZ;   
      err+=findVal("alfEMZ",&alfEMZ); 
      param->ee=sqrt(alfEMZ*4*M_PI);
    }  
    err+=findVal("SW",&param->sw);
    param->MW=findValW("MZ")*sqrt(1-param->sw*param->sw);
	err+=findVal("Mtp",&param->Mt);
	err+=findVal("MHc",&param->Mhc);
	err+=findVal("tb",&param->tb);
/*	err+=findVal("Vts",&param->Vts);
	err+=findVal("Vtb",&param->Vtb);
	err+=findVal("Vcb",&param->Vcb);

	param->Vts=-.0405;
	param->Vtb=.99915;
	param->Vcb=.0411;
*/
	err+=findVal("Zv11",&(param->V[1][1]));
	err+=findVal("Zv12",&(param->V[1][2]));
	err+=findVal("Zv21",&(param->V[2][1]));
	err+=findVal("Zv22",&(param->V[2][2]));
	err+=findVal("Zu11",&(param->U[1][1]));
	err+=findVal("Zu12",&(param->U[1][2]));
	err+=findVal("Zu21",&(param->U[2][1]));
	err+=findVal("Zu22",&(param->U[2][2]));

	err+=findVal("Zt11",&(param->T[1][1]));
	err+=findVal("Zt12",&(param->T[1][2]));
	err+=findVal("Zt21",&(param->T[2][1]));
	err+=findVal("Zt22",&(param->T[2][2]));

	err+=findVal("Zb11",&(param->B[1][1]));
	err+=findVal("Zb12",&(param->B[1][2]));
	err+=findVal("Zb21",&(param->B[2][1]));
	err+=findVal("Zb22",&(param->B[2][2]));

	err+=findVal("Zn11",&(param->N[1][1]));
	err+=findVal("Zn21",&(param->N[2][1]));
	err+=findVal("Zn31",&(param->N[3][1]));
	err+=findVal("Zn41",&(param->N[4][1]));
	err+=findVal("Zn12",&(param->N[1][2]));
	err+=findVal("Zn22",&(param->N[2][2]));
	err+=findVal("Zn32",&(param->N[3][2]));
	err+=findVal("Zn42",&(param->N[4][2]));
	err+=findVal("Zn13",&(param->N[1][3]));
	err+=findVal("Zn23",&(param->N[2][3]));
	err+=findVal("Zn33",&(param->N[3][3]));
	err+=findVal("Zn43",&(param->N[4][3]));
	err+=findVal("Zn14",&(param->N[1][4]));
	err+=findVal("Zn24",&(param->N[2][4]));
	err+=findVal("Zn34",&(param->N[3][4]));
	err+=findVal("Zn44",&(param->N[4][4]));
	
	err+=findVal("MC1",&param->mc[1]);
	err+=findVal("MC2",&param->mc[2]);
	err+=findVal("MNE1",&param->mn[1]);
	err+=findVal("MNE2",&param->mn[2]);
	err+=findVal("MNE3",&param->mn[3]);
	err+=findVal("MNE4",&param->mn[4]);
	err+=findVal("MSt1",&param->mst[1]);
	err+=findVal("MSt2",&param->mst[2]);
	err+=findVal("MSb1",&param->msb[1]);
	err+=findVal("MSb2",&param->msb[2]);
	err+=findVal("MSsL",&param->mss[1]);
	err+=findVal("MSsR",&param->mss[2]);
	err+=findVal("MSnm",&param->msnmu);
	
	err+=findVal("MG1",&param->m1);
	err+=findVal("MG2",&param->m2);
	err+=findVal("MSG",&param->mglu);
	err+=findVal("At",&param->At);
	err+=findVal("Ab",&param->Ab);
	err+=findVal("mu",&param->mu);
	
	err+=findVal("MbMb",&param->Mb);
	err+=findVal("alfSMZ",&param->alphS_MZ);

//  printf("alphasMZ =%.3e \n",param->alphS_MZ);
        { double MSuL;
	  double MZ=param->MZ;
          double SW=param->sw;
	  double tb=param->tb;
          double sb=tb/sqrt(1+tb*tb);
	   err+=findVal("MSuL",&MSuL);
           param->msq1=sqrt(MSuL*MSuL-(MZ*MZ*(0.5-2/3.*SW*SW)*(1-2*sb*sb)));
        }      

	return err;
}

#ifdef DEBUG
int dump_prm(struct read_param_tag* param)
{
	int err=0;
	printf("\n dump_parm \n");
	printf("MW= %.3E \n",param->MW);
	printf("MZ= %.3E \n",param->MZ);
        printf("ee= %.3E \n",param->ee);
	printf("sw= %.3E \n",param->sw);
	printf("Mt= %.3E \n",param->Mt);
	printf("Mhc= %.3E \n",param->Mhc);
	printf("tb= %.3E \n",param->tb);
	printf("Vts= %.3E \n",param->Vts);
	printf("Vtb= %.3E \n",param->Vtb);
	printf("Vcb= %.3E \n",param->Vcb);
	
	
	printf("Zv11= %.3E \n",(param->V[1][1]));
	printf("Zv12= %.3E \n",(param->V[1][2]));
	printf("Zv21= %.3E \n",(param->V[2][1]));
	printf("Zv22= %.3E \n",(param->V[2][2]));
	
	printf("Zu11= %.3E \n",(param->U[1][1]));
	printf("Zu12= %.3E \n",(param->U[1][2]));
	printf("Zu21= %.3E \n",(param->U[2][1]));
	printf("Zu22= %.3E \n",(param->U[2][2]));

	printf("Zt11= %.3E \n",(param->T[1][1]));
	printf("Zt12= %.3E \n",(param->T[1][2]));
	printf("Zt21= %.3E \n",(param->T[2][1]));
	printf("Zt22= %.3E \n",(param->T[2][2]));
	
	printf("Zb11= %.3E \n",(param->B[1][1]));
	printf("Zb12= %.3E \n",(param->B[1][2]));
	printf("Zb21= %.3E \n",(param->B[2][1]));
	printf("Zb22= %.3E \n",(param->B[2][2]));

	printf("Zn11= %.3E \n",(param->N[1][1]));
	printf("Zn21= %.3E \n",(param->N[2][1]));
	printf("Zn31= %.3E  \n",(param->N[3][1]));
	printf("Zn41= %.3E \n",(param->N[4][1]));
	printf("Zn12= %.3E \n",(param->N[1][2]));
	printf("Zn22= %.3E \n",(param->N[2][2]));
	printf("Zn32= %.3E \n",(param->N[3][2]));
	printf("Zn42= %.3E \n",(param->N[4][2]));
	printf("Zn13= %.3E \n",(param->N[1][3]));
	printf("Zn23= %.3E \n",(param->N[2][3]));
	printf("Zn33= %.3E \n",(param->N[3][3]));
	printf("Zn43= %.3E \n",(param->N[4][3]));
	printf("Zn14= %.3E \n",(param->N[1][4]));
	printf("Zn24= %.3E \n",(param->N[2][4]));
	printf("Zn34= %.3E \n",(param->N[3][4]));
	printf("Zn44= %.3E \n",(param->N[4][4]));
	
	printf("MC1= %.3E \n",param->mc[1]);
	printf("MC2= %.3E \n",param->mc[2]);
	printf("MNE1= %.3E \n",param->mn[1]);
	printf("MNE2= %.3E \n",param->mn[2]);
	printf("MNE3= %.3E \n",param->mn[3]);
	printf("MNE4= %.3E \n",param->mn[4]);
	printf("MStop1= %.3E \n",param->mst[1]);
	printf("MStop2= %.3E \n",param->mst[2]);
	printf("MSbot1= %.3E \n",param->msb[1]);
	printf("MSbot2= %.3E \n",param->msb[2]);
	printf("MSsL= %.3E \n",param->mss[1]);
	printf("MSsR= %.3E \n",param->mss[2]);
	printf("Mq1= %.3E \n",param->msq1);
	printf("MSnmu= %.3E \n",param->msnmu);
	
	printf("MG1= %.3E \n",param->m1);
	printf("MG2= %.3E \n",param->m2);
	printf("MG3= %.3E \n",param->mglu);
	printf("At= %.3E \n",param->At);
	printf("Ab= %.3E \n",param->Ab);
	printf("mu= %.3E \n",param->mu);
	
	
	return err;
}

#endif


