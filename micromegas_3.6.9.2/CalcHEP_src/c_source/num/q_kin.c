/*
 Copyright (C) 1997, 1999  Alexander Pukhov 
*/

#include"../../include/nType.h"

#include"q_kin.h"
#include"decay.h"
#include"4_vector.h"

#include<math.h>
#include "interface.h"
#include"cut.h"
#include"kinaux.h"
#include"strfun.h"
#include"regfunal.h"
#include"subproc.h"
#include"q_kin.h"
#include"tools.h"
#include"drandXX.h"


static void LorRot(REAL rapidity, int ntot,REAL*pvect)
{
    static REAL rapid___ = 0;
    static REAL sh = 0;
    static REAL ch = 1;
    REAL ee, pp;
    int nv;

    if (rapidity != rapid___) 
    {
	rapid___ = rapidity;
	sh = sinh(rapidity);
	ch = sqrt(sh * sh + 1);
    }
    if (rapidity) for(nv=4*ntot-4; nv>=0; nv -=4)
    {  ee = pvect[nv];
       pp = pvect[nv+3];
       pvect[nv] = ee * ch + pp * sh;
       pvect[nv+3] = ee * sh + pp * ch;
    }
} 

static void AzimuthRot(REAL fi, int ntot,REAL*pvect)
{
  REAL sn = sin(fi);
  REAL cn = cos(fi);
  int i;
    
  for(i=0;i<ntot;i++)
  {  int nx=4*i+1,ny=4*i+2;
     REAL X,Y;
     X=pvect[nx];
     Y=pvect[ny];
     pvect[nx]= X*cn+Y*sn;
     pvect[ny]=-X*sn+Y*cn;
  }
} 


static void Rot3D(double  cn0 , double fi, double psi, int ntot, REAL *pvect)
{
   REAL cn=cn0, sn=sqrt(fabs(1-cn*cn));
   REAL fi_1=fi;
   REAL E[3]={ cn, sn*sin(fi_1), sn*cos(fi_1)};
   REAL EPS[3][3]= { {  0,  -E[2], E[1]},
                     { E[2],   0 ,-E[0]},
                     {-E[1], E[0],  0  }  }; 

   REAL P[3][3]=  {  {1-E[0]*E[0],  -E[0]*E[1], -E[0]*E[2] },
                     { -E[1]*E[0], 1-E[1]*E[1], -E[1]*E[2] },
                     { -E[2]*E[0],  -E[2]*E[1],1-E[2]*E[2] } };
    
  REAL ROT[3][3],cn2,sn2;
  double x1,x2,x3,f3; 
  int i,j,k;


  x1=psi;
  x2=M_PI;
  for(;;)
  { x3=0.5*(x1+x2);
    f3=x3-psi-sin(x3);
    if(f3<0) x1=x3;else x2=x3;
    if(fabs(f3)<1.E-4) break;
  }
                             
  cn2=cos(x3),sn2=sin(x3);
  
  for(i=0;i<3;i++) for(j=0;j<3;j++) ROT[i][j]= (cn2-1)*P[i][j] + sn2*EPS[i][j]; 
      
  for(k=0;k<ntot;k++) 
  {  REAL * X= pvect+4*k+1;
     REAL Y[3]={X[0],X[1],X[2]}; 
     for(i=0;i<3;i++)  for(j=0;j<3;j++) Y[i]+=ROT[i][j]*X[j];  
     for(i=0;i<3;i++) X[i]=Y[i];
  }                       

}


#define DEPTH 10

typedef  struct
{  int ncsreg[2];
   int ncscut[2];
   int tcscut[2];
   char  lvpole[PLISTLEN];
   int  itypep;
   REAL sph_we;
}  sphereStr;



static REAL sqrt_S, rapidity, stop, sbot, pcm;
static double ssmin, ssmax;

static REAL tfact0;
static int nout1;

static REAL pm[PLISTLEN];

static int nvpos0, nvposx, nvposy;

static int nvout[DEPTH][2], nvin[DEPTH];
static int  lnkbab[DEPTH], lnkleg[DEPTH];
static REAL  summas[DEPTH][2];
static int nmsreg[DEPTH][2], nmscut[DEPTH][2], nss;	
static int nsph[DEPTH];

static REAL beta[2];
static  sphereStr  sph_inf[DEPTH][10]; 



static REAL  cmfun( REAL E,REAL PM1,REAL PM2)
{
  REAL e2,s,d;
  e2=E*E;
  s=PM1+PM2;
  d=PM1-PM2;  
  return sqrt((e2-s*s)*(e2-d*d))/(2*E);     
}

int mkmom(double *x, double *tfact, REAL* pvect)
{
  int i,k,l;
  int nx=0;
  double  fct0, fct1, fct2;
  REAL  pIn[2][4]={{0,0,0,0},{0,0,0,0}};
  REAL  pXY[2][4]={{0,1,0,0},{0,0,1,0}};
    
  REAL  xcos, xfi, parfi;
  double cosmin, cosmax, parcos;
  REAL ytilda=0; 
  
  int  i__2;
  REAL d__1;

  REAL ff, al;
  REAL  xx, bes; 
  double fct;
 
  double stilda;
  REAL rstilda,  pcmtilda, xtilda; 
  int  ns;

  REAL psy1, psy2, x1,x2;

  double  smin, smax;
  REAL  amass[DEPTH][2];

  int nvpole;
  int nvpos;

  REAL hsum[2], hdif;

  REAL fct_1__;
  int nsing;
  
  iDecay memDecay;
    
  sing_struct singar[200];

  *tfact = tfact0;
/* **   MOMENTS */
    if(nout_int==1)
    { 
       REAL s=pm[0]+pm[1], d=pm[0]-pm[1], p=sqrt((pm[2]-s)*(pm[2]+ s)*(pm[2]-d)*(pm[2]+d))/(2*pm[2]);
       REAL e1=sqrt(p*p+pm[0]*pm[0]), e2=sqrt(p*p+pm[1]*pm[1]);
       REAL y1= log( (pcm +sqrt(pm[0]*pm[0]+pcm*pcm))/(e1+p) );
       REAL y2=-log( (pcm +sqrt(pm[1]*pm[1]+pcm*pcm))/(e2+p) );
       if(pm[0]>0) { REAL y=  log( pm[0]/(e1+p)); if(y2<y) y2=y; }
       if(pm[1]>0) { REAL y= -log( pm[1]/(e2+p)); if(y1>y) y1=y; }
       for(i=0;i<12;i++)pvect[i]=0;
       pvect[0]=e1; pvect[3]=p;
       pvect[4]=e2; pvect[7]=-p;
       pvect[8]=pm[2];
       *tfact=389379660.0*M_PI/(2*p*pm[2]*sqrt_S*sqrt_S);
       if(sf_num[0] && sf_num[1]) { ytilda=x[0]*y1+ (1-x[0])*y2;    *tfact*=(y1-y2);}   
        else if(sf_num[0]) ytilda=y2; else if(sf_num[1]) ytilda=y1; else  return 1;
       LorRot(ytilda,3,pvect);
       x[0]= sf_num[0]?  pvect[3]/pcm : 1;
       x[1]= sf_num[1]? -pvect[7]/pcm : 1;
       
       LorRot(rapidity,3,pvect);
       return 0;
    }
    if (nin_int == 2) 
    {   REAL y1,y2;
	if (sf_num[0] || sf_num[1]) 
	{
	    nsing = 0;
	    getreg_(&nsing, singar, 0., 1., nss);
	    bes =  beta[0]+ beta[1];
	    if (bes >= 1) al = 1; else if (nsing)  al = 0.5; else al = 0;

	    xx = x[nx++];
	    if (xx < al) 
	    {
		xx /= al;
		regfun_(2,nsing,singar,ssmin,ssmax,xx,&stilda,&fct1);
	    } else 
	    {   
		if(xx==al) xx=1;else xx = (1 - xx) / (1 - al);
		stilda=stop-pow(xx*pow(stop-ssmin,bes)+(1-xx)*pow(stop-ssmax,bes),1/bes);
		regfct_(2,nsing,singar,ssmin,ssmax,stilda,&fct1);
	    }
	    fct1 /= stop-sbot;
	    if (bes < 1) 
	    {
		ff = pow(1 - ssmin/stop, bes) - pow(1 - ssmax/stop, bes);
		fct2 = ff / bes * pow( 1 - stilda/stop, 1 - bes);
		*tfact *= ff/ (al * fct2 / fct1 + (1 - al));
	    } else  *tfact *= fct1;
	    
            xtilda=(stilda-sbot)/(stop-sbot);

	    if (sf_num[0] && sf_num[1]) 
	    { REAL yy = log(1/xtilda); 
		xx = x[nx++];
		if (beta[0] < 1 && beta[1] < 1) 
		{
		    al = beta[1]/(beta[0]+beta[1]);
		    if (xx < al) 
		    {
			xx /= al;
			psy1 = pow(xx, 1 / beta[0]);
			psy2 = 1 - psy1;
		    } else 
		    {
			if(xx==al) xx=1; else xx = (1 - xx) / (1 - al);
			psy2 = pow(xx, 1 / beta[1]);
			psy1 = 1 - psy2;
		    }
		    y1 = yy * psy1;
		    y2 = yy * psy2;

		    *tfact *=pow(divy_(y1),beta[0]-1)*pow(divy_(y2),beta[1]-1)/
		    	     (pow(psy1,1-beta[0])+pow(psy2,1-beta[1]));

		    if (bes < 1)  *tfact *= pow( divy_(yy),1-bes);
		     else          *tfact *= bes*pow(yy,bes-1);
		} else if (beta[0] < 1) 
		{
		    psy1 = pow(xx,  1 / beta[0]);
		    y1 = yy * psy1;
		    y2 = yy * (1 - psy1);
		    *tfact *=  pow(yy, *beta)*pow(divy_(y1),beta[0]-1);
		} else if (beta[1] < 1) 
		{
		    psy2 = pow(xx,  1 / beta[1]);
		    y2 = yy * psy2;
		    y1 = yy * (1 - psy2);
		    *tfact *= pow(yy,beta[1])*pow(divy_(y2),beta[1]-1);
		} else 
		{
		    y1 = yy * xx;
		    y2 = yy - y1;
		    *tfact *= yy;
		}
		x1 = exp(-y1);
		x2 = exp(-y2);
                ytilda=(y2-y1)/2;
	    } else if (sf_num[0]) 
	    {   REAL e2=sqrt(pm[1]*pm[1]+pcm*pcm);
		x1 = xtilda;
		x2 = 1;
                ytilda=0.5*log((e2+(2*x1-1)*pcm)/(e2+pcm));
	    } else 
	    {   REAL e1=sqrt(pm[0]*pm[0]+pcm*pcm);
		x1 = 1;
		x2 = xtilda;
                ytilda=-0.5*log((e1+(2*x2-1)*pcm)/(e1+pcm));
	    }
	} else 
	{
	    x1 = 1;
	    x2 = 1;
	    stilda = sqrt_S*sqrt_S;
            ytilda=0;
	}
        rstilda = sqrt(stilda);
        pcmtilda=cmfun(rstilda,pm[0],pm[1]);
	pIn[0][3]=pcmtilda;
	pIn[1][3]=-pcmtilda;

     
/*  *sf_fact= *tfact/tfact0; */
    }
/* *  FILLING ZERO COMPONENTS FOR in-PARTICLES */

    for (k = 0; k < nin_int; ++k) pvFill(pm[k],pIn[k],k+1,pvect);
	 
    if (nin_int == 2)   *tfact /= 4*pcmtilda *rstilda; else 
    { rstilda = pm[0]; *tfact /= rstilda * 2; ytilda=0;} 
     
/* *    X & Y AXISES */

    pvFill(0,pXY[0],nvposx,pvect);
    pvFill(0,pXY[1],nvposy,pvect);  

    nvpos = nvpos0;
    nvpole =nvpos++;
    
    
    
/* *    MASS INTEGRATION */

    for (i = 0; i < nout1; ++i) 
    {   REAL sval= i? amass[lnkbab[i]][lnkleg[i]]: rstilda;   

       	for (k = 0; k < 2; ++k) 
	{
	    if (kinmtc_1[i].lvout[k][1]) 
	    { double sqmass;
	      REAL  xx=x[nx++];
	        d__1=  k?  sval - amass[i][0] : sval - summas[i][1];
		         
	        smax = d__1 * d__1;
		
		d__1 = summas[i][k];
		smin = d__1 * d__1;

		if (nmscut[i][k]) rancor_(&smin, &smax,0., 1., nmscut[i][k]);
		
		if (smin >= smax)  {*tfact = 0; return 0;}
		
		if (nmsreg[i][k] )
		{
		    nsing = 0;
		    getreg_(&nsing, singar, 0., 1.,nmsreg[i][k]);
		    regfun_(2,nsing,singar,smin,smax,xx,&sqmass,&fct);
		} else
		{
		    sqmass = xx * smax + (1 - xx) * smin;
		    fct = smax - smin;
		} 
		amass[i][k] = sqrt(sqmass);
		*tfact *= fct;
	    }
	    else amass[i][k]= summas[i][k];
	}
    }

    lvtonv(kinmtc_1[0].lvin, 0 , nvin[0],pvect); /*very stupid*/ 

        
    for (i = 0; i < nout1; ++i)  /*  MAIN CYCLE */
    {   int ns___=nsph[i]-1;
        double Emax[2];
	if (i == 0 && nin_int == 1)  xcos = 0.1 /* was fixed  0.1 */; 
                               else  xcos = x[nx++];
	al = 0;
	l = 0;
	if (i == 0 || (i == 1 && nin_int == 1)) 
	{
	    xfi = 0.1;  /* was fixed  0.1;    */
            for(;l<=ns___;l++)
	    {  al +=  sph_inf[i][l].sph_we;
	       if (xcos <= al) ns___ = l;
	    }
	    xcos = (al - xcos) / sph_inf[i][ns___].sph_we;
	} else 
	{
	    xfi = x[nx++];
            for(;l<=ns___;l++)
            {
	       al += sph_inf[i][l].sph_we;
	       if (xfi <= al) ns___ = l;
	    }
	    xfi = (al - xfi) /sph_inf[i][ns___].sph_we;
	}
	lvtonv( sph_inf[i][ns___].lvpole,nin_int, nvpole,pvect);

	decay_0(nvin[i], amass[i][0], amass[i][1], &fct0, Emax,&memDecay,pvect);
	if(fct0==0) {*tfact=0; return 0;}
	if(!isfinite(fct0))
	{ printf("mkmom infinite factor\n");
	  *tfact=0;
	  return 0;
	}  
	decay_1(nvpole, hsum, &hdif,&memDecay,pvect);

	cosmin = -1;
	cosmax = 1;
	nsing = 0;
	for (k = 0; k < 2; ++k) 
	{   int ncM=sph_inf[i][ns___].ncscut[k];
	    int ncT=sph_inf[i][ns___].tcscut[k];
	    d__1 = ((k << 1) - 1) / hdif;
            getreg_(&nsing,singar,hsum[k],d__1,sph_inf[i][ns___].ncsreg[k]);
            if(ncM) rancor_(&cosmin,&cosmax,hsum[k],d__1,ncM);
            if(ncT) rancor_t(&cosmax,hsum[k],d__1,Emax[k], pm[sph_inf[i][ns___].lvpole[0]-1],pcmtilda, 
                            amass[i][k], invcut_1[ncT-1].cvmin );
	    if (cosmin >= cosmax) {*tfact = 0; return 0;}
	}
	regfun_(sph_inf[i][ns___].itypep,nsing,singar,cosmin,cosmax,xcos,&parcos,&fct);
	fct_1__ = sph_inf[i][ns___].sph_we / fct;
	parfi = (xfi * 2 - 1) * M_PI;
	decay_3(nvposy, parcos, parfi, nvout[i][0], nvout[i][1],&memDecay,pvect);

	
	i__2 = nsph[i];
	for (ns = 0; ns < i__2; ++ns)  if (ns != ns___)
	{
	    lvtonv(sph_inf[i][ns].lvpole, nin_int, nvpole,pvect);
	    decay_1(nvpole, hsum, &hdif,&memDecay,pvect);
	    decay_2(nvout[i][1], &parcos,&memDecay,pvect);
	    cosmin = -1;
	    cosmax = 1;
	    nsing = 0;
	    for (k = 0; k < 2; ++k) 
	    {   int ncM=sph_inf[i][ns].ncscut[k];
                int ncT=sph_inf[i][ns].tcscut[k];
		d__1 = ((k << 1) - 1) / hdif;
		getreg_(&nsing,singar,hsum[k],d__1,sph_inf[i][ns].ncsreg[k]);
		if(ncM) rancor_(&cosmin,&cosmax,hsum[k],d__1,ncM);
                if(ncT) rancor_t(&cosmax,hsum[k],d__1,Emax[k], pm[sph_inf[i][ns___].lvpole[0]-1],pcmtilda, 
                            amass[i][k], invcut_1[ncT-1].cvmin );

		if (cosmin>=parcos || parcos>=cosmax){*tfact=0; return 0;}
	    }
	    regfct_(sph_inf[i][ns].itypep,nsing,singar,cosmin,cosmax,parcos, &fct);
	    fct_1__ += sph_inf[i][ns].sph_we/ fct;
	}
	*tfact = *tfact * fct0 / fct_1__;
    }
    if(nin_int==2)
    { if(sf_num[0]) x[0]= x1;
      if(sf_num[1]) x[1]= x2;
      LorRot(rapidity+ytilda,nin_int+nout_int,pvect);
      AzimuthRot(drandXX()*2*M_PI,nout_int, pvect+8);        
    } else 
    { 
        Rot3D(2*(drandXX()-0.5),drandXX()*2*M_PI,drandXX()*M_PI,nout_int,pvect+4); 
        LorRot(rapidity,nin_int+nout_int,pvect);
    }
    
    if(!isfinite(*tfact))
    {  fprintf(stderr,"mkmom: infinite factor\n");
       *tfact=0;
      return 0;
    }

    for(i=0;i<(nin_int+nout_int)*4;i++) if(!isfinite(pvect[i])) {*tfact=0; return 0;}

    return 0;
}    


int imkmom(double P1, double P2)
{
    int i, j, k, l,ns;
    char lvbuf[PLISTLEN];
    int ndim;
    physValRec * pList;

    beta[0]=sf_be[0];
    beta[1]=sf_be[1];
    
    if(nin_int==2)  
    {  if(nout_int==1) ndim=1; else
       { 
         ndim = nout_int * 3 - 5;
         if (sf_num[0]) ndim++;
         if (sf_num[1]) ndim++;
         tfact0 = 2*M_PI*389379660.0;
       }  
    }else { tfact0 = 2*M_PI; if(nout_int==2) ndim=1; else ndim = nout_int * 3 - 7;}
    
    for (i=0; i <  nin_int + nout_int; i++) pinf_int(Nsub,i+1,pm+i,NULL);

    nout1 = nout_int - 1;
    if(nout1>DEPTH) return 0;
    nvposx = nin_int + nout_int + 1;
    nvposy = nvposx + 1;
    nvpos0 =  nvposy + 1;

/* *  NVOUT( , ) FILLING */
    for (i = 0; i < nout1; ++i)  for (k = 0; k < 2; ++k) 
    {
       if (kinmtc_1[i].lvout[k][1]) nvout[i][k] = nvpos0++;
       else                         nvout[i][k] = kinmtc_1[i].lvout[k][0];   
    }
    
    nvin[0] = nvpos0++;
    for (i = 1; i < nout1; ++i) 
    {   nvin[i]=0;
	for (j = 0; j < i; ++j) for (k = 0; k < 2; ++k) 
	{
	   if (eqvect_(kinmtc_1[i].lvin, kinmtc_1[j].lvout[k])) 
	   {
               nvin[i] = nvout[j][k];
	       lnkbab[i] = j;
	       lnkleg[i] = k;
	   }
	}
	if(!nvin[i]) { fprintf(stderr,"Error in kinematics \n"); sortie(52); }
    }


    for (i = 0; i < nout1; ++i) for (k = 0; k < 2; ++k) 
    {   
      REAL ss = 0; 
      int pn;
      
      for(j=0; pn=kinmtc_1[i].lvout[k][j];j++) ss += pm[pn - 1];
      summas[i][k] = ss;
    }

    if (nin_int == 2) 
    {  REAL m1=sf_num[0]?sf_mass[0]:pm[0];
       REAL m2=sf_num[1]?sf_mass[1]:pm[1];
       
       incomkin(m1, m2, P1, P2,  &sqrt_S, &pcm, &rapidity);

       ssmin=pm[0]+pm[1]; 
       if(ssmin<summas[0][0]+summas[0][1]) ssmin=summas[0][0]+summas[0][1];
       ssmin*=ssmin;      
 
       if( sf_num[0] && sf_num[1] ) {sbot=0;stop=4*pcm*pcm;}
       else if(sf_num[0]){sbot=m2*m2; stop=sbot+2*pcm*(pcm+sqrt(pcm*pcm+sbot));}
       else if(sf_num[1]){sbot=m1*m1; stop=sbot+2*pcm*(pcm+sqrt(pcm*pcm+sbot));} 
       else  stop = sqrt_S*sqrt_S;
       ssmax=stop;
    } else 
    {   REAL m1=pm[0];
        rapidity=log((P1+sqrt(P1*P1+m1*m1))/m1 ) ;
    }
    for(i = 0; i < nout1; ++i) 
    {
       nsph[i] = 0;
       for(k=0;k<2;k++)
       {
          for(ns=0; ns<10;ns++)
          {
             sph_inf[i][ns].ncsreg[k] = 0;
             sph_inf[i][ns].ncscut[k] = 0;
             sph_inf[i][ns].tcscut[k]=0;
          }
	  nmsreg[i][k] = 0;
	  nmscut[i][k] = 0;
       }
    }
    
    nss=0;
    for(l=0; invreg_1[l].lvinvr[0]; l++) 
    {   int orig=1, ll=0;
        for(;ll<l;ll++) {if( invreg_1[ll].nextrg == l+1) { orig=0; break;}}	
	if(orig) 
	{
	   sngpos_(invreg_1[l].lvinvr, &i, &k, lvbuf);
	   if (i==0)  nss = l+1; else
	   {  i--; k--;
	      if (lvbuf[0] == 0)  nmsreg[i][k] = l+1;
	      else 
	      {
		 for (ns = 0; ns <nsph[i]; ++ns) 
                 if (eqvect_(lvbuf,sph_inf[i][ns].lvpole)) 
		 {
		    sph_inf[i][ns].ncsreg[k] = l+1;
		    break;
		 }
		 if(ns==nsph[i] && ns<10)
		 {
		    nsph[i]++;
		    strcpy(sph_inf[i][ns].lvpole,lvbuf);
		    if (spole_(invreg_1[l].lvinvr))sph_inf[i][ns].itypep = -2; 
		    else sph_inf[i][ns].itypep = -1;
		         
	            sph_inf[i][ns].ncsreg[k] = l+1;
	         }
	      }
	   }
	}
    }

    for(l=0;l<nCuts;l++) if( invcut_1[l].key[0] == 'M')
    for(pList=invcut_1[l].pLists;pList;pList=pList->next)
    {   char buff[20];
        strcpy(buff,pList->pstr);
        coninv_(buff);
	sngpos_(buff, &i, &k, lvbuf);
	if (i == 0) rancor_(&ssmin, &ssmax, 0., 1., l+1);
	
	else if (lvbuf[0] == 0)  nmscut[i-1][k-1] = l+1;
	else 
	{   i--; k--;
	    for (ns = 0; ns <  nsph[i ]; ++ns) 
	    {
		if (eqvect_(lvbuf, sph_inf[i][ns].lvpole)) 
		{
		    sph_inf[i][ns].ncscut[k]  = l+1;
		    break;
		}
	    }
	    if(ns==nsph[i] && ns <10 )
	    {
	       nsph[i]++;
	       strcpy(sph_inf[i][ns].lvpole,lvbuf);
	       sph_inf[i][ns].itypep = 2;
               sph_inf[i][ns].ncscut[k] = l+1;
	    }
	}
    }
    if(nin_int==2 && ssmin>=ssmax) return 0;

    if(nin_int==2) for(l=0;l<nCuts;l++)
    {  int m;
       invcut_ tc=invcut_1[l];
/* 
       if( tc.key[0]=='T' && tc.key[1]!='^' && tc.minon 
                          && tc.pLists      && tc.pLists->pstr[1]==0 )  
       for(pList=invcut_1[l].pLists;pList;pList=pList->next)for(m=1;m<=2;m++)
       {  char str[4];
          strcpy(str+1,pList->pstr);
          str[0]=m;
          sngpos_(str, &i, &k, lvbuf);
	  {  i--; k--;
	     for (ns = 0; ns <  nsph[i ]; ++ns) 
	     {
		if(eqvect_(lvbuf, sph_inf[i][ns].lvpole)
                   &&strcmp(kinmtc_1[i].lvout[k],str+1)==0)
		{
		    sph_inf[i][ns].tcscut[k]  = l+1;
		    break;
		}
	     }
	  }
       }
*/ 
    }
        
    for(i = 0; i < nout1; ++i) 
    {
	if (nsph[i] == 0) 
	{
	   nsph[i] = 1;
           sph_inf[i][0].lvpole[0] = (i == 0 && nin_int == 1 )?  nvposx:1;
	   sph_inf[i][0].lvpole[1] = 0;
	   sph_inf[i][0].itypep = 1;
	} else 
	{
	   ns = nsph[i]-1;
	   if(sph_inf[i][ns].ncsreg[0] || sph_inf[i][ns].ncsreg[1])  
	   sph_inf[i][0].itypep *=-1;
	}
    }

    for(i = 0; i<nout1; ++i) 
    {   REAL wesum = 0;
	for (ns = 0; ns<nsph[i]; ++ns) 
	{
	    int   nwe = 0;
	    for (k = 0; k < 2; ++k) 
	    for (l = sph_inf[i][ns].ncsreg[k];l;l=invreg_1[l-1].nextrg) ++nwe;
		    
            if (sph_inf[i][ns].itypep >= 0) ++nwe;
	    sph_inf[i][ns].sph_we = nwe;
	    wesum += nwe;
	}
	for (ns = 0; ns <nsph[i]; ++ns)   sph_inf[i][ns].sph_we /= wesum;
    }

#ifdef DEBUG   
    for (i = 0; i < nout1; ++i) 
    {   
        printf("Decay number %d     nmscut= (%d,%d) nmsreg = (%d,%d)\n",
        i, nmscut[i][0], nmscut[i][1],  nmsreg[i][0], nmsreg[i][1]);
        
        {int  l,c;
          printf("kinematics= (");
          for (l=0;c= kinmtc_1[i].lvin[l];l++) printf("%d",c);
          printf(")->(");
          for (l=0;c=kinmtc_1[i].lvout[0][l];l++)  printf("%d",c);
          printf(")+(");
          for (l=0;c=kinmtc_1[i].lvout[1][l];l++)  printf("%d",c);
          printf(")\n");
        }        
        printf(" summas=(%f,%f)\n",summas[i][0],summas[i][1]);

        for (ns = 0; ns < nsph[i]; ++ns) 
	{   int c;
	    printf("   Sphere number = %d  weight=%f type=%d \n", 
	     ns, sph_inf[i][ns].sph_we, sph_inf[i][ns].itypep);
            printf("     pole vector(");
	    for(k=0; c=sph_inf[i][ns].lvpole[k]; k++) printf("%d",c);
            printf(")\n");
		
	    printf("    ncsreg=(%d,%d) ncscut=(%d,%d) tcscut=(%d,%d) \n",
	    sph_inf[i][ns].ncsreg[0], sph_inf[i][ns].ncsreg[1], 
            sph_inf[i][ns].ncscut[0], sph_inf[i][ns].ncscut[1],
            sph_inf[i][ns].tcscut[0], sph_inf[i][ns].tcscut[1]); 
	}
    }    
#endif 
    return ndim;
} /* mkmom_ */
