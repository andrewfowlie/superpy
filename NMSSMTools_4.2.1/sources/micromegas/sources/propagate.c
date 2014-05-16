#include<stdlib.h>
#include<math.h>
#include<stdio.h>

#include "micromegas.h"
#include "micromegas_aux.h"
extern double Zi(int i);

/*===================  micrOMEGAs Halo profile =====================*/

static double   Ha     =   20.;
static double   Halpha =   1;
static double   Hgam   =   1;
static double   Hbeta  =   3;

double hProfileZhao(double r)
{
   double rhomax,rho,a1,a2,a3; 
   a1=pow((Rsun/r),Hgam);
   a2=1+pow(Rsun/Ha,Halpha);
   a3=1+pow(r/Ha,Halpha);
   rho = a1 * pow(a2/a3,(Hbeta-Hgam)/Halpha );
   
   rhomax = 1e10;
   if(rho >= rhomax) return rhomax; else return rho;
}

void setProfileZhao(double Halpha_n, double Hbeta_n, double Hgam_n, double Ha_n)
{
   Ha        =  Ha_n;
   Halpha    =  Halpha_n;
   Hgam      =  Hgam_n;
   Hbeta     =  Hbeta_n;
   setHaloProfile(hProfileZhao);
}

static double alphaEin=0.17;

double hProfileEinasto(double r)
{ 
   return  exp(-2*(pow(r/Rsun,alphaEin)-1)/alphaEin);
} 

void setProfileEinasto(double alpha)
{  alphaEin=alpha;
   setHaloProfile(hProfileEinasto);
}
  
static double f_clump=0, rho_clump=0;

extern void setClumpConst(double f,double rho);
extern double rhoClumpsConst(double r);
extern void setRhoClumps(double (*cProfile)(double));

void setClumpConst(double f,double rho)
{  f_clump=f; rho_clump=rho; setRhoClumps(rhoClumpsConst); }

static double (*hProfile_)(double)=hProfileZhao;

double rhoClumpsConst(double r)
{ double v=(rho_clump - rhoDM*(*hProfile_)(r)*f_clump);
  if(v<=0) return 0;
  return f_clump*v;
}

static double (*rhoClumpEff_)(double)=rhoClumpsConst;

void setHaloProfile(double (*hProfile)(double)){ hProfile_=hProfile;}

void setRhoClumps(double (*cProfile)(double)) {   rhoClumpEff_=cProfile; }


#define CELERITY_LIGHT 2.99792458e10 /* cm.s{-1} */
#define DIGIT 15.0
#define CM_PER_KPC 3.0856775807e+21

#define Eps  (1.E-3)

/* ------------------ */
/* Conversion Factors */
/* ------------------ */
/* mass, cross section */
#define mb_to_cm2 1.e-27           /* mb to cm2 (i.e. 1mb = mb_to_cm2 * cm2) */
#define amu2GeV 931.494013e-3      /* amu (atomic mass unit) to GeV          */
#define GeV2kg 1.78264e-27         /* etc. */
/* length */
#define cm_to_kpc 3.2407793e-22    /* From PDG    */
#define km_to_kpc 3.2407793e-17
/* time */
#define s_to_Myr ( 1./(1.e6*31557600.) )  /* Sideral year (from PDG 2005) */
#define h_to_Myr ( s_to_Myr/3600. )
#define d_to_Myr ( h_to_Myr/24. )
#define yr_to_Myr 1.e-6
#define kyr_to_Myr 1.e-3
/* speed */
#define kmpers_to_kpcperMyr 1.022712e-3 
#define cmpers_to_kpcperMyr 1.022712e-8
/* Diffusion */
#define kpc2perMyr_to_cm2pers 3.01715e29

/* ----------- */
/* Constant... */
/* ----------- */
#define c_light (29979245800.*cmpers_to_kpcperMyr)  /* Light speed kpc/Myr */
#define m_p 0.938271998        /* Proton mass (GeV)   */
#define m_n 0.939565561        /* Proton mass (GeV)   */
#define m_D 1.875              /* Deuteron mass (GeV) */


/*--- Galaxy-related (r and z are cylindrical coordinates) */
#define h .1          /* Disk half-thickness (kpc): this parameter should not be modified */
#define DensityH  .9  /* All hydrogen (MilkyWay disc) in cm-3  */
#define DensityHe .1  /* All Helium (Milky Way disk) in cm-3   */

#define rHaloMin 0.000001

/*=================================  Photon propagation ====================================*/


static double fi_, dfi_,sn_;

static double xIntegrand(double x) 
{  double r=Rsun*sqrt(x*x+sn_*sn_);
   double pf;
   if(r<rHaloMin) r=rHaloMin;
   pf=hProfile_(r);
   return  pf*(pf +rhoClumpEff_(r)/rhoDM );
}

static double yIntegrand(double y)
{ double x,r,res,pf;
  if(y==0) return 0;
  x=1/y-1;
  r=Rsun*sqrt(x*x+sn_*sn_);
  if(r<rHaloMin) r=rHaloMin;
  pf=hProfile_(r);
  res= pf*(pf+rhoClumpEff_(r)/rhoDM)/y/y;

  return res;  
}

static double psiIntegrand(double psi)
{ double Int,a,b;
  double cs_= cos(fi_-dfi_) - sin(fi_)*sin(dfi_)*(1-sin(psi));
 
  a=sin(fi_-dfi_);
  b=sin(fi_)*sin(dfi_)*(1-sin(psi)); 
  sn_= sqrt( a*a+b*(2*cos(fi_-dfi_)-b)); 

  
  if(cs_>0)
  { double xmin;
    if(rHaloMin/Rsun>sn_) xmin=sqrt((rHaloMin/Rsun-sn_)*(rHaloMin/Rsun+sn_)); else xmin=0;
    Int=2*simpson(xIntegrand,xmin,cs_,1.E-4);
  } else Int=0;
  
  Int+= simpson(yIntegrand,0,1/(fabs(cs_)+1),1.E-4);
  return Int;
}

static double fiIntegrand(double fi)
{ 
  dfi_=fi;
  return sin(fi)*simpson(psiIntegrand,-M_PI/2,M_PI/2,1.E-4)/M_PI;
}  

double HaloFactor(double fi,double dfi)
{  double res;
   double sm_in_kpc=3.0856775807E21;
   double Norm;
   
   fi_=fi;
   if(dfi<=0) {dfi_=0;res=psiIntegrand(0);} 
   else res=simpson(fiIntegrand,0,dfi,1.E-4)/(1-cos(dfi));
   Norm=rhoDM/hProfile_(Rsun)/Mcdm;
   return res*Norm*Norm/(8*M_PI)*Rsun*sm_in_kpc;
}

void gammaFluxTab(double fi,double dfi, double sigmaV, double *Sp, double *Sobs)
{
  int i; 
  double hf=HaloFactor(fi, dfi)*sigmaV;
  if(dfi>0) hf*=(1-cos(dfi))*2*M_PI; 
  for(i=0;i<NZ;i++) Sobs[i]=hf*Sp[i]; 
}

double gammaFlux(double fi, double dfi,  double dSigmadE )
{
  return  dSigmadE*HaloFactor(fi, dfi)*(1-cos(dfi))*2*M_PI;
}


/*********************************************************************************/
/*********************************************************************************/
/*  Positrons propagation related routines */
/*********************************************************************************/
/*********************************************************************************/


static double propagator_quantic(double K0_tau, double z)
{
  int n;
  double sum=0,ds1=0,ds2=0, L=L_dif;

  for(n=0; ds1>=sum*0.03*Eps || ds2>=sum*0.03*Eps  ;n++)
  { double kn=0.5*M_PI*(1+2*n)/L;
    double en=cos(kn*z)*exp(-kn*kn*K0_tau);
    ds2=ds1;
    ds1=fabs(en);
    sum+=en;
  }
  return 2*sum/L;
}

static double propagator_image(double Kt, double zS)
{
  double sum;
  double x1,x2,sgn=1;
       
  if(zS<0) zS*=-1;
  if(zS>=L_dif) return 0;
  x1=zS, x2=2*L_dif-zS;
  for(sum=0;;x1+=2*L_dif,x2+=2*L_dif,sgn*=-1)
  {
    double dSum= exp(-x1*x1/(4*Kt)) -  exp(-x2*x2/(4*Kt)) ;
    sum+=sgn*dSum;
    if(dSum<0.03*Eps*sum) break;
  }      
  return sum/sqrt(M_PI*Kt); 
}


static double V(double E) { return pow(E,Delta_dif-1)/(1-Delta_dif);}

static double green_v(double Kt, double zS)
{
  if(zS>=L_dif) return 0;
  if(L_dif*L_dif/(4*Kt) >= 1) return propagator_image(Kt,zS);    /* [kpc^{-1}] */
      else                    return propagator_quantic(Kt,zS);  /* [kpc^{-1}] */
}

static double azimuthInt(double k)
{
/* return 1/pi*int_0^pi exp(- k sin^2(fi/2)) dfi == exp(-k/2)*bessi0(k/2)  */    

double X[18]={0.000000E+00,3.125000E-02,6.250000E-02,1.250000E-01,1.875000E-01,
 2.500000E-01,2.812500E-01,3.125000E-01,3.750000E-01,4.062500E-01,4.375000E-01,
 5.000000E-01,5.625000E-01,6.250000E-01,7.500000E-01,8.750000E-01,9.375000E-01,
 1.000000E+00};    
   
double Y[18]={3.544908E+00,3.660214E+00,3.785458E+00,4.072489E+00,4.424998E+00,
 4.880735E+00,5.174304E+00,5.534618E+00,6.449089E+00,6.942868E+00,7.404071E+00,
 8.105752E+00,8.444493E+00,8.465206E+00,7.930969E+00,7.108068E+00,6.687196E+00,
 6.283185E+00};
      
double z=1/(sqrt(k)+1);
       
return polint4(z,18,X,Y)*z/2/M_PI;
}

static double Kt_,dR_,z_;
#ifdef NEW
static double  rIntegrand(double r)
{  double I_angular_rS, prQ,rr;
   if(r==0) return 0;
   I_angular_rS=azimuthInt(Rsun*r /Kt_);
   rr=sqrt(r*r+z_*z_);
   prQ=hProfile_(rr);
   prQ*=(prQ+rhoClumpEff_(rr)/rhoDM);
   return r*I_angular_rS*exp(-(Rsun-r)*(Rsun-r)/(4*Kt_))*prQ;
}

static double zIntegrand(double z)
{ 
  double rMin,rMax;
  rMin   = Rsun - dR_;
  rMax   = Rsun + dR_;
  if (rMin <= 0.0  ) rMin = 0.0;
  if (rMin*rMin < rHaloMin*rHaloMin-z*z) rMin=sqrt(rHaloMin*rHaloMin-z*z);  
  if (rMax > Rdisk) rMax = Rdisk;

  z_=z;  
  return   green_v(Kt_,z)*simpson(rIntegrand,rMin,rMax, 1.E-4);
}        


static double integral_cal_In(double K0_tau)
{
   if(K0_tau<0.001) { double prQ=hProfile_(Rsun); prQ*=prQ+rhoClumpEff_(Rsun)/rhoDM;   return prQ/2;}
   Kt_=K0_tau;
   dR_ = sqrt(4.0 * Kt_ * log(10.0) * DIGIT);  

   return  simpson(zIntegrand,0,dR_, 1.E-3)/(4*Kt_);
}

#else

static double r_;

static double zIntegrand(double z)
{ double rr=sqrt(r_*r_+z*z), prQ;
  prQ=hProfile_(rr); prQ*=prQ+rhoClumpEff_(rr)/rhoDM;
  return   green_v(Kt_,z)*prQ; 
}        

static double  rIntegrand(double r)
{  double I_angular_rS,I_z,zMax,zMin;
   if(r==0) return 0;
   r_=r;
   if(dR_>L_dif) zMax=L_dif; else zMax=dR_;
   if(r>=rHaloMin)zMin=0; else zMin=sqrt(rHaloMin*rHaloMin -r*r );
   I_z=simpson(zIntegrand,zMin,zMax,0.1*Eps);
   I_angular_rS=azimuthInt(Rsun*r /Kt_);
   return  r*I_angular_rS*exp(-(Rsun-r)*(Rsun-r)/(4*Kt_))*I_z;
}

static double integral_cal_In(double K0_tau)
{
  double rMin,rMax;
  if(K0_tau<0.001) {double prQ; prQ=hProfile_(Rsun); prQ*=prQ+rhoClumpEff_(Rsun)/rhoDM; return prQ/2;}
  Kt_=K0_tau;
  
  dR_  = sqrt(4.0 * Kt_ * log(100./Eps));
  rMin   = Rsun - dR_;
  rMax   = Rsun + dR_;
  if (rMin <= 0.0  ) rMin = 0.0;
  if (rMax > Rdisk) rMax = Rdisk;

  return simpson(rIntegrand,rMin,rMax, 0.3*Eps)/(4*Kt_);
}


#endif
static double Eobs;
static double*tab_;


static double dKt(double Es,double Eo)
{
/*  double   epsilon_t = (0.1 * (3600. * 24. * 365.25));  */	 /* [sec] */
  double epsilon_t= 1.E-7; /* Myr */

  return   K_dif*( (Tau_dif*s_to_Myr)*(V(Eo)-V(Es)) +  epsilon_t);
}

static double PosifluxIntegrand(double Erun)
{  return  integral_cal_In(dKt(Erun,Eobs))*SpectdNdE(Erun,tab_);  }


double posiFlux(double E, double sigmav, double *tab)
{
  double flu;
  double rho0;
  if(E>=Mcdm)return 0;
  rho0=rhoDM/hProfile_(Rsun)/Mcdm;
  tab_=tab; 
  flu =Tau_dif/(E*E)*sigmav/(4*M_PI)*CELERITY_LIGHT;
  Eobs=E;
  return  flu*rho0*rho0*simpson(PosifluxIntegrand,E,Mcdm,Eps); 
}

static double* xa_,*ya_;
static int N_;

static double SpectIntegrand(double Erun)
{ return  polint4(dKt(Erun,Eobs),N_,xa_,ya_)*SpectdNdE(Erun,tab_); }

void posiFluxTab(double Emin, double sigmav, double *tab, double *tabOut)
{
  double flu,rho0;
  int i;
  double buff[NZ];
  tab_   =tab;
  rho0=rhoDM/hProfile_(Rsun)/Mcdm;
  flu = Tau_dif*sigmav/(4*M_PI)*CELERITY_LIGHT;

  buildInterpolation(integral_cal_In,0., dKt(Mcdm,Emin),-Eps,&N_,&xa_,&ya_);

//printf("N_=%d\n",N_);
  
  for(i=0;i<NZ;i++)
  {
     Eobs=Mcdm*exp(Zi(i));
     if(Eobs<Emin*0.9) buff[i]=0; 
     else buff[i]= (flu/Eobs)* simpson(SpectIntegrand, Eobs, Mcdm, Eps);
  }   

  for(i=0;i<NZ;i++) tabOut[i]=rho0*rho0*buff[i]; 
  free(xa_); free(ya_);    
}

/*********************************************************************************/
/*  End of positrons propagation related routines */
/*********************************************************************************/

/*********************************************************************************/
/*  Antiprotons propagation related routines                                     */
/**** Source code adapted from Maurin, Taillet, Combet, astro-ph/0609522 *********/
/************* PROPAGATOR FORMULA (constant wind model) **************************/
/* ------ */
/* MACROS */
/* ------ */
/* beta=v/c and rigidity=p/|Z| from Ek (for pbar) */
#define SQR(a)          ( (a)*(a) )
#define pbar_beta(Ekin)   ( sqrt(SQR(Ekin)+2.*m_p*Ekin)/(Ekin+m_p) )
#define pbar_rig(Ekin)    ( sqrt(SQR(Ekin)+Ekin*2.*m_p) )
#define Dbar_beta(Ekin)   ( sqrt(SQR(Ekin)+2.*m_D*Ekin)/(Ekin+m_D) )
#define Dbar_rig(Ekin)    ( sqrt(SQR(Ekin)+Ekin*2.*m_D) )

static double Leff=0;
static double sigma_inelastic_pbarH_TAN_and_NG(double EK_pbar)
{
  double result;
  
  if (EK_pbar<=0.0) return 0;
  result  = 1. + 0.584*pow(EK_pbar,-0.115) + 0.856*pow(EK_pbar,-0.566);
  result *= 24.7;
  return result;  /*[mb]*/
}

/********************************************************************************************/
static double sigma_inelastic_pH_TAN_and_NG(double EK_proton)
{
  double U,Cp,E_proton;
  double result;
  
  E_proton = EK_proton+m_p;
  result = 0.0;
  if (EK_proton<=0.0)
  {
    return result;
  }
  else if (EK_proton<0.3)
  {
    return result;
  }
  else if (EK_proton<3.0)
  {
    U = log(E_proton/200.0);
    result  = 1. + 0.0273*U;
    result *= 32.2e-27;
    Cp = 17.9 + 13.8*log(EK_proton) + 4.41*log(EK_proton)*log(EK_proton);
    result /= 1. + 2.62e-3*pow(EK_proton,-Cp);
    return result/mb_to_cm2; /*[mb]*/
  }
  else if (E_proton<200.0)
  {
    U = log(E_proton/200.0);
    result = 1. + 0.0273*U;
    result *= 32.2e-27;
    return result/mb_to_cm2; /*[mb]*/
  }
  else
  {
    U = log(E_proton/200.0);
    result = 1. + 0.0273*U + 0.01*U*U;
    result *= 32.2e-27;
    return result/mb_to_cm2; /*[mb]*/
  }
}

static double sigma_inelastic_NOANN_pbarH_TAN_and_NG(double EK_pbar)
{
  if (EK_pbar<=0.0)  return 0;
  else if (EK_pbar < 13.3)
  {
      double result = 661.*(1. + 0.0115 * pow(EK_pbar, -0.774) - 0.948 * pow(EK_pbar, 0.0151));
      return (sigma_inelastic_pbarH_TAN_and_NG(EK_pbar) - result);
  }
  else return sigma_inelastic_pH_TAN_and_NG(EK_pbar);
}

/********************************************************************************************/


static double sigma_inelastic_ANN_pbarH_TAN_and_NG(double T_pbar)
{
  double result;

  if (T_pbar<=0.0)  return 0; else
  result= (sigma_inelastic_pbarH_TAN_and_NG(T_pbar) - sigma_inelastic_NOANN_pbarH_TAN_and_NG(T_pbar));
  if (result<0)  result=0.;
  return result;
}


/**************************************************************************/
/************    HELIUM        ********************************************/
/**************************************************************************/
static double sigma_inelastic_pbarHe_TAN_and_NG(double T_pbar)
{
double result;
result =sigma_inelastic_pbarH_TAN_and_NG(T_pbar)*2.519842 ; /*[mb]*/
return result;
}


static double sigma_inelastic_NOANN_pbarHe_TAN_and_NG(double T_pbar)
{
double result;
result = sigma_inelastic_NOANN_pbarH_TAN_and_NG(T_pbar)*2.519842; /*[mb]*/
return result;
}

static double sigma_inelastic_ANN_pbarHe_TAN_and_NG(double T_pbar)
{
double result;
result=sigma_inelastic_pbarHe_TAN_and_NG(T_pbar)-
	 sigma_inelastic_NOANN_pbarHe_TAN_and_NG(T_pbar);
return result;
}

/*========================================*/

/************* Propagator for pbar (Constant wind model) *************/
static double Kdif, Gtot,Vdif;
static double pbar_PropagatorInfiniteR(double rSrc, double zSrc, double ek)
{
   double distance;
   double propagator,kspal,kv,kd,cn,knL;
   int n,i;

   if (rSrc==0.) return 0;

   distance= sqrt(SQR(rSrc) + SQR(zSrc));

//   n_terms_infinite_sum=1000000;

//   if (distance<L_dif/100.) return 1./(4.*M_PI*Kdif*distance);

   propagator=0.;
   kspal= 2.*h*Gtot/Kdif;
   kv=Vdif/(2.*Kdif);
   kd= (kspal + 2.*kv);
   for (n=0;  ; n++)
   {  double dPropagator;
//      if(n<=lastK) knL=kn[n]; else
      { 
          knL= (n+0.5)*M_PI; cn=1;
          if(kd)
          { int nn;
            if(kd>0) nn=n+1; else nn=n; 		
            for(i=0; i<10; ++i) knL = nn*M_PI - atan(2.*knL/(kd*L_dif));
          }
//          kn[n]=knL;
//          lastK=n; 
      }         
      cn= 1. - sin(2.*knL)/(2.*knL);
      dPropagator=bessk0(sqrt(SQR(knL/L_dif)+SQR(kv))*rSrc);
      propagator+=dPropagator*sin(knL*(L_dif-zSrc)/L_dif)*sin(knL)/cn;
      if(fabs(dPropagator) <= 0.1*Eps*fabs(propagator)) break; 
   }
   return propagator*exp(-kv*zSrc)/(2.*M_PI*Kdif*L_dif);
}

static double ek_,z_,r_;

#define rHalo1  0.1

static double rhoQ_2(double r){ double prQ=hProfile_(r); prQ*=prQ+rhoClumpEff_(r)/rhoDM;  return r*r*prQ;  }
static double rhoQ_3(double r){ double prQ=hProfile_(r); prQ*=prQ+rhoClumpEff_(r)/rhoDM;  return r*r*r*prQ;}

static double  thetaIntegrandP(double uTheta)
{
  double theta=uTheta*uTheta;
  double sinTh=sin(theta/2);
  double prQ;
  double d=sqrt(z_*z_ + (Rsun-r_)*(Rsun-r_) +4*r_*Rsun*sinTh*sinTh);
  if(d>=Rdisk) return 0;
  prQ=hProfile_(d); prQ*=prQ+rhoClumpEff_(d)/rhoDM;
  return  uTheta*prQ*(1-exp(-2*(Rdisk-d)/Leff));
}

static double rIntegrandP(double r)
{
  double fluxDM_r_z;
  double thetaMin,thetaMax,sinMax;
  double drQ  = rHalo1*rHalo1-z_*z_ -(Rsun-r)*(Rsun-r);
  double RQ = Rdisk*Rdisk  -z_*z_ -(Rsun-r)*(Rsun-r);
  double prQ;
  if(r==0)  return 0;  
  if(RQ<=0) return 0;
  if(drQ>0) thetaMin= 2*asin(sqrt(drQ/(4*Rsun*r))); else thetaMin=0;
  sinMax=sqrt(RQ/(4*Rsun*r));
  if(sinMax>=1) thetaMax=M_PI; else thetaMax=2*asin(sinMax);

  r_=r;
  prQ=hProfile_(rHalo1); prQ*=prQ+rhoClumpEff_(rHalo1)/rhoDM;
  
  fluxDM_r_z=2*2*simpson(thetaIntegrandP,sqrt(thetaMin),sqrt(thetaMax),
  0.1*Eps)+2*thetaMin*prQ;   	

  return r*fluxDM_r_z*pbar_PropagatorInfiniteR(r,z_,ek_);
}

static double zIntegrandP(double u)
{ 
  double z=u*u;
  double Ismall=0,a,b,r1,r2,p1,p2;
  if(u==0) return 0;
  if(z>L_dif-0.00001) return 0;
z_=z; 
  r1=Leff/25;  p1=rIntegrandP(r1)*sqrt(r1*r1+z*z)/r1; 
  r2=Leff/20;  p2=rIntegrandP(r2)*sqrt(r2*r2+z*z)/r2;
  
   a=(p1*r2*r2-p2*r1*r1)/(r2*r2-r1*r1);
   b=(p1-p2)/(r1*r1-r2*r2);

  Ismall=(a-b*z*z)*(sqrt(r1*r1+z*z)-z) +b/3*(pow(r1*r1+z*z,3./2.)-z*z*z);

//  return u*(Ismall+simpson(rIntegrandP,r1,10*Leff, Eps/3));
  return  u*(Ismall+simpson(rIntegrandP,r1,Rsun, Eps/3)
                     +simpson(rIntegrandP,Rsun,Rdisk+Rsun,Eps/3));
                       
  
  return  u*(Ismall+simpson(rIntegrandP,r1,Rsun, Eps/3)
                   +simpson(rIntegrandP,Rsun,Rdisk+Rsun,Eps/3));
}

int  Gtot_style=2;

static double pbarPropRate(double ek)
{
  double fluxDM, dFluxDM;
  double i2,i3;
  ek_=ek; 

  Kdif=pbar_beta(ek)*K_dif*pow(pbar_rig(ek),Delta_dif);
//printf("Kdif=%E\n",Kdif);
  switch(Gtot_style)
  { case 0: Gtot=0; break;
    case 1: Gtot= ( DensityH*sigma_inelastic_ANN_pbarH_TAN_and_NG(ek) 
                  + DensityHe*sigma_inelastic_ANN_pbarHe_TAN_and_NG(ek) );
            break;
            
   default: Gtot= 0.8*sigma_inelastic_pbarH_TAN_and_NG(ek)*(DensityH + DensityHe*2.519842);

  }   
  
  Gtot*= pbar_beta(ek)*c_light/cmpers_to_kpcperMyr/s_to_Myr * mb_to_cm2;

  Vdif=Vc_dif*kmpers_to_kpcperMyr;

  { double kv,kd,knL;
    int i;
    kv=Vdif/(2.*Kdif);
    kd= (2.*h*Gtot/Kdif    + 2.*kv);
    if(kd)
    { int nn;
      knL= (0.5)*M_PI;
      if(kd>0) nn=1; else nn=0; 		
      for(i=0; i<10; ++i) knL = nn*M_PI - atan(2.*knL/(kd*L_dif));
    } else knL= (0.5)*M_PI;
//    lastK=0;
//    kn[0]=knL;
    Leff=1/sqrt(SQR(knL/L_dif)+SQR(kv)); 
  }  
          
  fluxDM=4*simpson(zIntegrandP,0.,sqrt(L_dif),Eps);

  i2=simpson(rhoQ_2,rHaloMin, rHalo1, 1.E-2) -rHalo1*rhoQ_2(rHalo1)/3;
  i3=simpson(rhoQ_3,rHaloMin, rHalo1, 1.E-2) -rHalo1*rhoQ_3(rHalo1)/4; 
    
  dFluxDM=4*M_PI*i2*pbar_PropagatorInfiniteR(Rsun,0.5*i3/i2,ek);
  fluxDM+=dFluxDM;   
  return fluxDM*pbar_beta(ek)*c_light/cm_to_kpc/(4.*M_PI);
}

double pbarFlux(double E, double dSigmadE)
{ double rho0;
  rho0=rhoDM/hProfile_(Rsun)/Mcdm ;
  return 0.5*rho0*rho0*pbarPropRate(E)*dSigmadE;
}


static double logPbarRate(double x)
{ return log(0.5*pbarPropRate( Mcdm*exp(x)));}

void pbarFluxTab(double Emin, double sigmav, double *tab, double *tabOut)
{
  int i;
  int N;
  double * Egrid,*Fgrid;
  double tab2[NZ];
  double rho0;

  buildInterpolation(logPbarRate,log( Emin/Mcdm),log(0.9),0.01,&N,&Egrid,&Fgrid);
/*printf("Npbar=%d\n",N);*/
  rho0=rhoDM/hProfile_(Rsun)/Mcdm;
  for(i=0;i<NZ;i++)
  {  double z=Zi(i);
     double E=Mcdm*exp(z);
     if(E<Emin*0.9) tab2[i]=0; else tab2[i]= 
     sigmav*exp(polint4(z,N, Egrid, Fgrid) )*zInterp(z,tab);
  }   
  for(i=0;i<NZ;i++) tabOut[i]=rho0*rho0*tab2[i];
 
  free(Egrid); free(Fgrid);   
}

/*********************************************************************************/
/*  End of antiprotons propagation related routines */
/*********************************************************************************/

void solarModulation(double PHI, double mass, double * inTab, double * outTab)
{ double buff[NZ];
  int i;
  for(i=0;i<NZ;i++)
  { 
    double X,E,P,Xs,Es,Ps;
    X=exp(Zi(i));
    E=Mcdm*X;
    P=sqrt( E*(E+2*mass));
    Es=E+PHI/1000.;
    if(Es>=Mcdm) buff[i]=0; else
    {  Ps=sqrt( Es*(Es+2*mass));
       Xs=Es/Mcdm;
       buff[i]= (X/Xs)*zInterp(log(Xs),inTab)*pow(P/Ps,2);    
    }
  }
  for(i=0;i<NZ;i++) outTab[i]=buff[i];
}


double pBarBackgroundFlux(double E)  /* arXiv:astro-ph/0609522v3 */
{ const   double C[5]={-3.211,0.12145,-0.2728,-0.075265,-0.007162};
  const   double D[2]={-2.02735,1.16463};
  double s, x;
  if(E<0) return 0; 

  x=log(E);
  if(E<11) 
  { int i;
    for( i=3,s=C[4];i>=0;i--)  s=s*x+C[i];
  } else s=D[0]*pow(x,D[1]);
  if(E<11) return 1.E-4*exp(s);  /* m^2 => cm^2 */
  else return 8.621348E-01*1.E-4*exp(s);
}   

void pBarBackgroundTab(double *pBarTab)
{
int i;
for(i=0;i<NZ;i++) {double E=Mcdm*exp(Zi(i));pBarTab[i]= E*pBarBackgroundFlux(E);}

}