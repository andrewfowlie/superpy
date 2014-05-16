#include"sf_isr.h"
#include<math.h>
#include"tools.h"

#include"crt_util.h"
#include"interface.h"
#include"simpson.h"
#include"runVegas.h"
#include"subproc.h"
#include"4_vector.h"
#include"plot.h" 

#define NPOINTS 100
#define   EM     5.1099906e-4
#define   EGAM   0.5772156649          /* Euler constant */
#define  ALPHA   0.0072973530796448189
#define  mmToGeV 5.067E12

#define EPS 1.e-6
#define B2 (1./3.)
#define R3 (1./3.)


static double scale=1., xy_nm=560, z_mm=0.4, qTot=2.E10;
static int bOn=0;

static double beta, coeff, b_ncl,b_ng,b_ips,b_cips;
static double xi[NPOINTS+4],yi[NPOINTS+4];

/* E.A.Kuraev,V.S.Fadin:Sov.J.Nucl.Phys.41(1985)466 */
/* S.Jadach,B.F.L.Ward:Comp.Phys.Commun.56(1990)351 */




int mc_isr__(int i) 
{ int N;
  pinf_int(Nsub,i,NULL,&N);
  return N;
}

static double b_h(double etax)
{
  int  n;
  double s0=0., ds, g[3];

  g[0] = 0.37328217390739632;
  if(etax<=0.) return g[0]* gammai_(2, b_ng);
  g[1] = pow(etax, R3)/1.354117939426404/2.;
  g[2] = pow(etax, 2.*R3) / 6.;
    
  for(n=1;n<4;n++) s0+=g[n-1]* gammai_(n+1, b_ng);
    
  do
  {  int n3=(n-1)%3;
     g[n3]*= 3 *etax /(n*(n-1)*(n-2)*(n-3));
     ds = g[n3]*gammai_(++n, b_ng);
     s0 += ds;
  } while (ds > s0 * EPS/100); 

   return s0;
}


static double cfbeam(double x)
{
   if (x <= 0.) return 0.; else
   {  double k=2./(b_ips * 3.);
      double etax=k * (1/x - 1);
      if (etax > 50.)  return 0;   
      return  (1+x*(b_cips-1 ))*pow(k/x, R3)*exp(-etax)*b_h(etax*pow(1+x*(b_cips-1 ),3));
   }
} 


static double cfbeamLog(double y){return cfbeam(exp(-y))*pow(divy_(y),-2.*R3);} 

static double cfisr(double x)
{ return coeff*(x *x+1-beta*(log(x)*(3.*x*x + 1.)/ 2.+(1.-x) *(1.-x))/2.)/2;}

static double cfisrLog(double y){return cfisr(exp(-y))*pow(divy_(y),beta-1);} 

int p_isr__(int* pNum) { if(abs(pNum[0])==11 && pNum[1]==0) return 1; else return 0; }

void n_isr__(int i, char *name)
{char *c;
  if(scale<=1.)  sprintf(name,"Qisr=%.2E*sqrtS,Beamstr:",scale);
     else        sprintf(name,"Qisr=%.2E[GeV],Beamstr:",scale);

   if(bOn) sprintf(name+strlen(name),"%.1f,%.3f,%.1E",xy_nm,z_mm,qTot);
   else    strcat(name," OFF");
   c=strstr(name,"E+0"); if(c) for(c+=3; ; c++) { c[-2]=c[0]; if(c[0]==0) break;} 
   c=strstr(name,"E+0"); if(c) for(c+=3; ; c++) { c[-2]=c[0]; if(c[0]==0) break;}
}

static void calc_params(void)
{
  double   sqrt_S,Q;
/*  vinf_(0,NULL,&sqrt_S); */
  REAL Pcm;;
  incomkin(0.,0.,inP1,inP2,NULL,&Pcm,NULL);
  sqrt_S=2*Pcm;
   
  if(scale>1.) Q=scale; else Q=scale*sqrt_S;

  beta = ALPHA*(2*log(Q/EM)-1)/M_PI;
  coeff = exp(beta * (.75 - EGAM)- lgamma(1+beta)); 
  if(bOn)
  {
    b_ncl=  25*ALPHA*ALPHA*qTot/(12*EM*(xy_nm*1.E-6)*mmToGeV);
    b_ips=5*ALPHA*qTot*sqrt_S/(12*EM*EM*EM*z_mm*(xy_nm*1.E-6)*mmToGeV*mmToGeV);
    if(b_ips>=1) b_cips=sqrt(1+pow(b_ips,0.666666666666666666)); else b_cips=1;
    b_ng=b_ncl/b_cips;
  }
}


static double norm_test(double y) { return 3*cfbeam(1-y*y*y);}
 

static double cisrf_(double x) { return cfisr(1-pow(x,1/beta));}

static double norm_tot_func(double z)
{ double x=1-pow(z,1/beta);
  double y=pow(1-x,R3);
  return dinter_(y, NPOINTS+4, xi, yi);
}


int i_isr__(int ii,double * be,double * mass)
{
  int i;
  static int bOn_old=-1;
  static double beta_old=0, coeff_old=0, b_ncl_old=0, b_ips_old=0;
  double iRest;
  
  *mass=5.11E-4;
  calc_params();

  if(beta==beta_old && coeff==coeff_old &&  bOn==bOn_old)
  {  if(!bOn)  {*be= beta; return 1;}
     if(b_ncl== b_ncl_old && b_ips==b_ips_old) {*be=beta; return 1;}
  }

  for (i = 0; i < NPOINTS; ++i)  xi[i] = (double)(i)/NPOINTS;
  for(i= NPOINTS;i<NPOINTS+4;i++)  xi[i]= 1- (1-xi[i-1])*0.5;
  if(bOn) iRest= b_ng-simpson(norm_test,0.,1.,1.E-8);
  for (i = 0; i < NPOINTS+4; ++i)
  {  double x=xi[i];
     x=1-x*x*x;
     yi[i] = cfisr(x);
     if (bOn)
     {   double lx=-log(x);
         yi[i] = (yi[i] * iRest
  	         + pow(1-x, B2)* pow(divy_(lx), 1-beta-B2)
  	        * convol_(cfisrLog, cfbeamLog, beta, B2, lx, EPS)) / b_ng;
     }
  }
//#define CHECK_NORN
#ifdef CHECK_NORN
if(bOn)
 printf("Check of Beamstralung normalization  delta(1-x)*%E  + %E(Beamstralung) =1?\n",(1-exp(-b_ng))/b_ng, 
 simpson(norm_test,0.,1.,1.E-8)/b_ng);
 
 printf("ISR [+ beamstrahlung] normalization: %E\n", simpson(norm_tot_func,0.00001,0.99999,1.E-8)); 
  
#endif
    
  beta_old=beta;
  coeff_old=coeff;
  bOn_old=bOn;
  b_ncl_old=b_ncl;
  b_ips_old=b_ips;

  *be=beta;
  return 1;
}


int r_isr__(int i, char *name)
{  
   double   z0,z1,z2,z3;
   char buff[100];
   if (sscanf(name,"Qisr=%lf%*[^:]:%lf,%lf,%lf",&z0,&z1,&z2,&z3)==4)
   {  bOn=1;
      scale=z0;
      xy_nm=z1;
      z_mm =z2;
      qTot =z3;
   }else if(sscanf(name,"Qisr=%lf%*[^:]: %s",&z0,buff)==2)
   {   bOn=0;
       scale=z0;
   } else  return 0;
   
   return 1;
}
              
double c_isr__(int i, double x, double q)
{x=pow(1-x,R3); return dinter_(x, NPOINTS+4, xi, yi);}


  

int  m_isr__(int i, int*pString)
{
  void * pscr=NULL;
  int mode =1;

  for(;;)
  {  char strmen[]="\40" 
     " ISR scale      = XXX           "
     " Beamstralung          ON       "
     " Bunch x+y sizes (nm)= YYY      "
     " Bunch lenght (mm)   = ZZZ      "
     " Number of particles = NNN      "
     "          *  N_gamma = NG       "
     "          *  Upsilon = UPS      "
     " Beamstrahlung F(x) plot        "
     " Beamstrahlung F(x)*(1-x)^(2/3) ";

     if(scale<=1) improveStr(strmen,"XXX","%.2E*sqrtS",scale);
     else improveStr(strmen,"XXX","%.0fGeV",scale);

     calc_params(); 
     if(bOn) 
     {    
        if(b_cips<1) improveStr(strmen,"_gamma","_cl");
        improveStr(strmen,"YYY","%.1f",xy_nm);
        improveStr(strmen,"ZZZ","%.3f",z_mm);
        improveStr(strmen,"NNN","%.1e",qTot);
        improveStr(strmen,"NG","%.2f", b_ng);
        improveStr(strmen,"UPS","%.2f",b_ips);
     } else
     {
        improveStr(strmen,"ON","%3.3s","OFF");
        strmen[2*strmen[0]+1]=0;  
     }
     menu1(46,10,"",strmen,"n_sf_isr",&pscr,&mode);
    
     switch(mode)
     { 
       case 0: return 1;

       case 1: correctDouble(25,16,"Enter value:( <=1 in part of sqrt(S); >1 in GeV) ",&scale,1);break;
       case 2: bOn=!bOn; break;
       case 3: correctDouble(52,16,"Enter new value ",&xy_nm,1);break;
       case 4: correctDouble(52,16,"Enter new value ",&z_mm,1);break; 
       case 5: correctDouble(52,16,"Enter new value ",&qTot,1);break;    
       case 6:
       case 7: messanykey(10,10, "This parameter is a function of\n"
                                 "above ones and Sqrt(S)");    break;
       case 8: case 9:
       { double xMin=0, xMax=1;
         int nPoints=100;
         double f[100]; 
         void * screen;
         int i;
         get_text(1,1,maxCol(),maxRow(),&screen);
         for(i=0;i<nPoints;i++)
         {  double x=xMin+(i+0.5)*(xMax-xMin)/nPoints;
              if(mode==8)f[i]=cfbeam(x)/pow(1- x,2*R3)/b_ng;
              else       f[i]=cfbeam(x)/b_ng;
         }
         {
              char mess[100];
              double dx=(xMax-xMin)/(2*nPoints);
              if(mode==8) sprintf(mess,"F(x)");
              else sprintf(mess,"F(x)*(1-x)^(2/3)");
              plot_1(xMin+dx,xMax-dx,nPoints,f,NULL,"Beamstrahlung","x",mess);  
         }
         put_text(&screen);
//         printf("cfbeam(1)/b_ng=%E\n",cfbeam(1)/b_ng);
        
       }                    
     }
  }
  return 1;
}
