#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "vegas.h"
#include "crt_util.h"
#include "drandXX.h"
#include "syst.h"
#define MAX_DIM   15
#define MAX_NDMX  50

long EventGrid=10000;

static void drand_arr(int dim, double * x)
{ int i;
  unsigned int umax= UINT_MAX;
  int rest=dim;
  unsigned int randpos=umax*drandXX();
  
  for(i=0;i<dim;i++) x[i]=-1; 
  
  for(i=0;i<dim;i++)
  {
     int pos=randpos%rest;
     int j;
     for(j=0;;j++) if(x[j]<0)
     {if(pos)pos--;else {x[j]=drandXX();break;}}
     randpos/=rest;
     umax/=rest;
     rest--;
     if(umax<rest){ umax= UINT_MAX;randpos=umax*drandXX();}
  }
}


static double amotry(double *p, double *y, int ndim,
	     double (*f)(double *,void*), int ilo, double fac,void*Aux)
{
   int i,j;
   double  ytry, fac1=(1.0-fac)/ndim;
   double * p_buff=p+(ndim+1)*ndim;
   double * p_ilo =p+ilo*ndim;
   
   for(j=0;j<ndim;j++) p_buff[j]=p_ilo[j]*fac;
   for(i=0;i<=ndim;i++)  if(i!=ilo) 
     {double *p_i=p+i*ndim;  for(j=0;j<ndim;j++) p_buff[j] +=p_i[j]*fac1;} 
   ytry=(*f)(p_buff,Aux);
   
   if (ytry > y[ilo]) 
   {  for(j=0;j<ndim;j++) p_ilo[j]=p_buff[j];
      y[ilo]=ytry;
   }
/*printf("amotry returns fac=%f %E\n",fac,  ytry);   */
   return ytry;
}

static double amoeba(double*p,double*y,int ndim,double (*f)(double *,void*), double eps, int *nCalls,void*Aux)
{
   int i,ilo,ihi,inlo,j;
   double ysave,ytry;

   for (;;) 
   {
      ihi=0;									     
      ilo = y[0]<y[1] ? (inlo=1,0) : (inlo=0,1);				     
      for (i=0;i<=ndim;i++)							     
      {										     
     	 if (y[i] >= y[ihi]) ihi=i;
     	 if (y[i] < y[ilo]) { inlo=ilo; ilo=i; } 
     	 else if (y[i] < y[inlo] && i != ilo) inlo=i;
      }										     

      if((*nCalls)<=0) break;
      { double d=0;
        for(i=0;i<ndim;i++) d+=(p[i+ihi*ndim]-p[i+ilo*ndim])*(p[i+ihi*ndim]-p[i+ilo*ndim]);
        d=sqrt(d);
        if(y[ihi]-y[ilo] < eps*d) break;
      }
       										     
      
//      if((*nCalls)<=0||2*(y[ihi]-y[ilo])/(fabs(y[ilo])+fabs(y[ihi]))<eps)break;
     										     
      ytry=amotry(p,y,ndim,f,ilo,-1.0,Aux); (*nCalls)--;				     
      if (ytry >= y[ihi]) {ytry=amotry(p,y,ndim,f,ilo,2.,Aux); (*nCalls)--;}	     
      else if (ytry <= y[inlo])							     
      {										     
         ysave=y[ilo];								     
     	 ytry=amotry(p,y,ndim,f,ilo,0.5,Aux);  (*nCalls)--;
     	 if (ytry <= ysave)
     	 {  
     	    for (i=0;i<=ndim;i++)
     	    {  double * p_ihi=p+ihi*ndim;
               if (i != ihi)
     	       {  double * p_i=p+i*ndim;
     		  for(j=0;j<ndim;j++) p_i[j]=0.5*(p_i[j]+p_ihi[j]);
     		  y[i]=(*f)(p_i,Aux);
     	       }
            }
/*printf("srink\n");            */
     	    (*nCalls) -= ndim;
         }									     
      }										     
   }
   return y[ihi];
}
/*========================== end of amoeba ================*/


static int Ng[MAX_DIM];
static int Ndim;
static double (*f_)(double*, double,double*);
static int Ndmx;
static double *Xgrid;
static double *Cgrid;

#define XG(j,i) Xgrid[(i)+(j)*(Ndmx)]


static int nroot(long N, int n)
{  int i,r;
   long N_;
   if(n==1) return N;
   r=pow(N,1./n);
   
   for(i=1,N_=r; i<n; i++) N_*=r; 
   
/*   printf("N_0=%d\n",N_);*/

   for(; N_<N; ) { r++; for(i=1,N_=r; i<n; i++) N_*=r;}
   for(; N_>N; ) { r--; for(i=1,N_=r; i<n; i++) N_*=r;} 
   return r;   
}
static void  generateVegasCubes(vegasGrid * vegPtr,long * nCubes) 
{  int i;
   double nCubes_=*nCubes;

   Ndim=vegPtr->ndim;
   Ndmx = vegPtr->ndmx+1;
   Xgrid= vegPtr->x_grid;
   Cgrid= vegPtr->c_grid;
   *nCubes=1;
   for(i=0;i<Ndim;i++)
   {
      Ng[i]   =nroot(nCubes_,Ndim -i);
      nCubes_ /= Ng[i];
      *nCubes *= Ng[i];
   }
}



static void Local2Global(int*Kg,double*XLOC,double*XGLOB,double*JACOB,int*GRID_LOC)
{ int j,n;
  double xlj,xn,xn_;                               
  *JACOB=1;          
  for (j = 0; j < Ndim; ++j)                       
  {  xlj = (Kg[j ] + XLOC[j])/Ng[j];               
     n=(int)(xlj*(Ndmx-1));                              
     if (n) xn_= XG(j,n);else  xn_=0;               
     xn = XG(j,n+1);
     XGLOB[j] = xn_ +(xn-xn_)*(xlj*(Ndmx-1)-n);          
     *JACOB *= (xn-xn_)*(Ndmx-1);                         
     if(GRID_LOC) GRID_LOC[j] = n;                 
  }
}


vegasGrid *  vegas_init(int dim,int nd)
{ 
   vegasGrid * vegPtr;

   if((dim>MAX_DIM)||(nd>MAX_NDMX)) return NULL;
   vegPtr=(vegasGrid * )malloc(sizeof(vegasGrid));
   if(vegPtr)
   { 
      vegPtr->ndmx=nd;
      vegPtr->ndim = dim;
      vegPtr->x_grid=malloc(dim*(nd+1)*sizeof(double));
      vegPtr->c_grid=malloc(dim*(nd+1)*sizeof(double));
      
      Xgrid=vegPtr->x_grid;
      Ndmx = vegPtr->ndmx+1;
      if(vegPtr->x_grid && vegPtr->c_grid)
      { int i, j;
        double * x_=vegPtr->x_grid;
        double * c_=vegPtr->c_grid;
        for(j=0;j<dim;j++) for(i=0;i<=nd;i++,x_++,c_++) 
        {  *x_=i/(double)nd;
           *c_=1./nd;
        }
      }else 
      { 
        if(vegPtr->x_grid) free(vegPtr->x_grid);
        if(vegPtr->c_grid) free(vegPtr->c_grid); 
        free(vegPtr); 
        return NULL;
      }
      vegPtr->nCubes=0;
      vegPtr->fMax=NULL;
   }
   return vegPtr;
}


void vegas_finish(vegasGrid * vegPtr) 
{   if(vegPtr)
    {  free(vegPtr->x_grid); 
       free(vegPtr->c_grid); 
       if(vegPtr->fMax) free(vegPtr->fMax);
       free(vegPtr); 
       vegPtr=NULL;
    }
}

/*     			*  VEGAS  *
      SUBROUTINE PERFORMS NDIM-DIMENSIONAL MONTE CARLO INTEG'N 
      - BY G.P. LEPAGE    SEPT 1976/(REV)AUG 1979 
      - ALGORITHM DESCRIBED IN J COMP PHYS 27,192(1978) 
*/

int vegas_int(vegasGrid * vegPtr, long ncall0, double alph, 
 double (*fxn)( double *,double), double *ti, double *tsi)
{
   int dim= vegPtr->ndim;

   double *d=malloc((vegPtr->ndmx+1)*dim*sizeof(double));
#define DD(j,i) d[(i)+(j)*Ndmx]

   double x[MAX_DIM];
   double xlocal[MAX_DIM];
   int    ia[MAX_DIM];
   int Kg[MAX_DIM],kg[MAX_DIM],ng[MAX_DIM];
   int i,j,first;
   double  f2,fb,f2b;
   int  npg=2;
   long nCubes=ncall0/npg;
   long cCube,l;
   int ret_code=0;
   float *pmax=NULL;
      
   if(alph==0)
   {
     double nCubes_=EventGrid;
     l=1;
     for(i=0;i<dim;i++)
     {
       ng[i]   =nroot(nCubes_,dim -i);
       nCubes_ /= ng[i];
       l *= ng[i];
     }   
     if(vegPtr->fMax && l!=vegPtr->nCubes) 
     { printf("remove old grid (it is a mistake!)\n");
       free(vegPtr->fMax); vegPtr->fMax=NULL;vegPtr->nCubes=0;
     }  
     if(vegPtr->fMax) 
     {
        pmax=vegPtr->fMax; 
        l=vegPtr->nCubes;
     } else 
     { 
       vegPtr->fMax=pmax=malloc(l*sizeof(float));
       vegPtr->nCubes=l;
       for(cCube=0;cCube<l;cCube++) pmax[cCube]=0;
     } 
   } else  if(vegPtr->fMax){free(vegPtr->fMax);vegPtr->fMax=NULL;vegPtr->nCubes=0;}
   
   generateVegasCubes(vegPtr,&nCubes);

   npg=ncall0/nCubes;
    
   *ti  = 0;
   *tsi = 0;
   for (j = 0; j < dim; ++j) { for (i = 0; i < Ndmx-1; ++i)  DD(j,i) = 0;}

/*    - MAIN INTEGRATION LOOP */
   for(i=0;i<dim;i++) Kg[i]=0;
   for(cCube=0; cCube<nCubes; cCube++) 
   {  
      if(informline(cCube,nCubes))  { ret_code=1; goto exi;}
      fb  = 0;
      f2b = 0;
      for(i = 0; i<npg; i++) 
      {   double f;
          drand_arr(dim,xlocal);
          Local2Global(Kg,xlocal,x, &f, ia);
          f *= (*fxn)(x,f);
          fb += f;
          f2= f*f;
          f2b += f2;
          for (j = 0;j<dim;++j) DD(j,ia[j]) += f2;

          if(pmax)
          { 
            for(j=0;j<dim;j++) { double xj= (Kg[j]+xlocal[j])/(double)Ng[j];    kg[j]=xj*ng[j];}
            for(j=1, l=kg[0];j<dim;j++) l=l*ng[j]+kg[j];
            f=fabs(f);
            if( pmax[l]<f ) pmax[l]=f;
          }  
      }
       
      f2b = sqrt(f2b/npg);
      fb /=npg;
       
      f2b = (f2b - fb) * (f2b + fb)/(npg-1);             
       
      *ti  += fb/nCubes;
      *tsi += f2b/((double)nCubes * (double)nCubes);   
    
      for(i=dim-1; i>=0; i--){if(++Kg[i]<Ng[i]) break; else Kg[i]=0;}
    }
    
    *tsi = sqrt(fabs(*tsi));

    if(*tsi < 1.E-6*fabs(*ti)) *tsi=1.E-6*fabs(*ti);

    if(alph>0)  /* REFINE GRID */
    { 
        double r[MAX_NDMX]; 					       
        double dt[MAX_DIM];
        double xin[MAX_NDMX];
        
        double  xn, xo,dr;
        int k,ndm= Ndmx - 2;
        
        for (j = 0; j < dim ; ++j)
        {
            xo = DD(j,0);
            xn = DD(j,1);
            DD(j,0) = (xo + xn) / 2;
            dt[j] = DD(j,0);
            for (i = 1; i < ndm; ++i)
            {
                DD(j,i) = xo + xn;
                xo = xn;
                xn = DD(j,i+1);
                DD(j,i) = (DD(j,i) + xn) / 3;
                dt[j] += DD(j,i);
            }
            DD(j,ndm) = (xn + xo) / 2;
            dt[j] += DD(j,ndm);
        }
        for (j = 0; j < dim; ++j)
        {  double rc = 0;
    	   for (i = 0; i <= ndm; ++i)
    	   {
    	      r[i] = 0;
    	      if (DD(j,i) > 0)
    	      {  double xoln = log(dt[j]/DD(j,i));
    	         if (xoln <= 70.f)  r[i] = pow( (1 - exp(-xoln))/xoln, alph);
    	         else               r[i] = pow(  1/xoln,               alph);
    	      }
    	      rc += r[i];  
    	   }

    	   rc /= (Ndmx-1);
           if(rc)
	   { 
             for(i=0,k=0,xn=0,dr=0;i<ndm;) 
       	     {
                do
                {  dr += r[k];
    	           xo = xn;
    	           xn = XG(j,k+1);
    	           k++;
    	         } while (rc > dr);
    	         do
    	         {  dr -= rc;
    	            xin[i] = xn-(xn-xo)*dr/r[k-1];
                    i++;
                 } while (rc<=dr);
    	      }
    	      for (i=0;i<ndm;++i)  XG(j,i+1) = xin[i];
    	      XG(j,ndm+1) = 1;
    	      XG(j,0)=0;
           }
        }
     }
exi:  free(d); return ret_code;
#undef DD
} 
/*
#define INCUB(x)   ((x)>0   ? (0.5+(x))/((x)+1.) : 0.5/(1.-(x)))
#define OUTCUB(x)  ((x)>0.5 ? (0.5-(x))/((x)-1.) : 1-0.5/x)
*/

#define OUTCUB(x) (x)

static double INCUB(double x) { while(x> 2)x-=2; if(x>1) x=2-x; 
                         while(x<-1)x+=2; if(x<0) x=-x;
                         return x;}


static double f_max(double* x,void *Kg)
{
  int i;
  double f;
  double xg[MAX_DIM];
  double xl[MAX_DIM];
  double tmp;

  for(i=0;i<Ndim;i++) xl[i]=INCUB(x[i]);
  Local2Global((int*)Kg,xl,xg, &f,NULL);

  tmp=f_(xg,f,NULL);
  f*=tmp;
  
  if(f<0) return -f; else return f;
}


static double run_amoeba(int ndim, int*Kg,double *xx, double *y, double step, double eps, int nCalls)
{
   int i,j;
   for (i=1; i<=ndim;++i)
   { double * x_i=xx+i*ndim; 
     for(j=0;j<ndim;++j) x_i[j]=xx[j];
     if(x_i[i-1] >0.5)x_i[i-1] -= step; else x_i[i-1] += step; 
     for(j=0;j<ndim;j++) x_i[j]=OUTCUB(x_i[j]);
     y[i]=f_max(x_i,Kg);
   }

   for(j=0;j<ndim;j++) xx[j]=OUTCUB(xx[j]);

/*printf(" test y[0]: %E %E \n",y[0], f_max(xx));*/

   {
      double r=amoeba(xx,y,ndim,f_max,eps,&nCalls,Kg);
      for(j=0;j<ndim*(ndim+1);j++) xx[j]=INCUB(xx[j]);
      return r; 
   }
}


int vegas_max(vegasGrid * vegPtr, long  nCubes, long nRandom, long nSimplex, double milk,
 double (*fxn)(double *,double,double*), double * eff)
{
   int dim= vegPtr->ndim;
   
   double x[MAX_DIM], xlocal[MAX_DIM], xx[(MAX_DIM+2)*MAX_DIM],y[(MAX_DIM+2)];
   double average =0,sum=0;
   float *fmax;
   int Kg[MAX_DIM];
   long  cCube;
   int i,ret_code=0;

   if(nRandom<=0) nRandom=2; 
   generateVegasCubes(vegPtr,&nCubes);
   if(vegPtr->fMax && vegPtr->nCubes !=nCubes)
   { free(vegPtr->fMax); vegPtr->fMax=NULL;}
   if(vegPtr->fMax) fmax=vegPtr->fMax; else 
   { vegPtr->nCubes=nCubes;    
     vegPtr->fMax=malloc(nCubes*sizeof(float));
     fmax=vegPtr->fMax;
     for(i=0;i<nCubes;i++) fmax[i]=0;
   } 
   f_=fxn;   // for  amoeba 

//   for(i=0;i<dim;i++) Kg[i]=0;
   for(cCube=0; cCube<nCubes; cCube++) 
   {  double cMax=fmax[cCube]; 
      long k,L0; 
      double f;

      xx[0]=-1;
      if(informline(cCube,nCubes)) 
      { free(vegPtr->fMax);
        vegPtr->fMax=NULL;  
        return 1;
      }
      L0=cCube; for (i =dim-1;i>=0; i--) {Kg[i]=L0%Ng[i]; L0=L0/Ng[i];}      

      for(k = 0; k<nRandom ; k++) 
      {  
         drand_arr(dim,xlocal);
         Local2Global(Kg,xlocal,x, &f, NULL);f *=(*fxn)(x,f,NULL); if(f<0) f=-f;
         average+=f;
         if(f>cMax){cMax=f;  for(i=0;i<dim;i++)xx[i]=xlocal[i]; y[0]=cMax; }
      }
      if(cMax && nSimplex && xx[0]>=0) cMax=run_amoeba(dim, Kg,xx, y, 0.1 , 0.5, nSimplex); 
      fmax[cCube]=cMax;
//      for(i=dim-1; i>=0; i--){if(++Kg[i]<Ng[i]) break; else Kg[i]=0;}
   }
   
   for(sum=0,cCube=0; cCube<nCubes; cCube++) sum+=fmax[cCube];
   if(sum==0)
   { free(vegPtr->fMax);
     vegPtr->fMax=NULL; 
     *eff=0;
   } else *eff=(average/nRandom)/sum;

/*
   if(milk)
   {   double minFmax=sum*milk/nCubes;
       for(cCube=0; cCube<nCubes; cCube++)
       if(fmax[cCube]<minFmax) fmax[cCube]=minFmax;
       *eff/=1+milk;  
   }
*/

//   printf("average=%E smax=%E eff=%e milk=%e \n",average,sum,*eff,milk);
exi:    
   return ret_code;
}


long vegas_events(vegasGrid * vegPtr,  long  nEvents, double gmax,
double (*fxn)( double *,double,double*), 
   void (*out)(long,int,char*,double*),int recalc,
   double * eff, int * nmax, int * mult, int * neg)
{
   int dim= vegPtr->ndim;
   double x[MAX_DIM], xlocal[MAX_DIM],xx[MAX_DIM*(MAX_DIM+2)],y[MAX_DIM+2];
   long  cEvent,cCube,L0,L1,Ntry=0, nCubes=vegPtr->nCubes;
   double rc,sum=0;
   float * smax;
   double pvect[40]; 

   struct  { double pvect[40]; double f; double sgn; }  *events=NULL;
   int i,nRecT=0;


   if(!(vegPtr->fMax)) return -1;
   smax=vegPtr->fMax;

   f_=fxn;
   
   generateVegasCubes(vegPtr,&nCubes);
   
   *eff=*nmax=*mult=*neg=0;
   
   for(cCube=1; cCube<nCubes; cCube++) smax[cCube]+=smax[cCube-1];
   sum=smax[nCubes-1];
   
   if(sum==0) 
   { messanykey(10,10,"Integrand is zero for all scanned points"); 
     if(blind) sortie(140); else return -1;
   }
   
   Ntry=0;
   for(cEvent=0; cEvent<nEvents; ) 
   {  long L;
      double f,ds;
      int n,sgn=0;
      char *drandXXstate ;      
      double oldMax,newMax;
      int Kg[MAX_DIM];

      Ntry++;

      rc=drandXX()*sum;
      L0=0;
      L1=nCubes-1;
      while(L0+1 < L1)     
      {  L=(L0+L1)/2;
         if(smax[L]<=rc) L0=L;  else L1=L;  
      }  
      if(smax[L0]>rc) L=L0; else L=L1;
      
 
      if(informline(cEvent,nEvents))  break;

      L0=L; for (i =dim-1;i>=0; i--) {Kg[i]=L0%Ng[i]; L0=L0/Ng[i];}
/*
 printf("L=%d dim=%d Ng={%d %d %d}, Kg={%d %d %d}\n",
        L,dim, Ng[0],Ng[1],Ng[2],Kg[0],Kg[1],Kg[2]);      
*/
      drandXXstate=seedXX(NULL);                                               
      drand_arr(dim,xlocal); 
      Local2Global(Kg,xlocal,x,&f,NULL);f*=(*fxn)(x,f,pvect);if(f<0){f=-f;sgn=-1;} else sgn=1;      
      (*eff)++;
      
      oldMax= L? smax[L]-smax[L-1] :smax[0];
      newMax=f;
      if(f<=oldMax*gmax && f>oldMax*gmax*drandXX())
      { cEvent++; 
        (*out)(L,sgn,drandXXstate,pvect);
        if(sgn<0) *neg++;  
      }  else  if(f>oldMax*gmax) 
      {   
           long dN=0;
           int i,k;
           int nRec=0;   
           char state[100]; 
           strcpy(state,drandXXstate); 
           
           if(!nRecT){ events=malloc(sizeof(*events));nRecT=1;}
           events[0].f=f;
           events[0].sgn=sgn;
           for(i=0;i<40;i++) events[0].pvect[i]=pvect[i];
           nRec=1;
           
           if(recalc)for(dN=1;dN < Ntry*newMax/sum/gmax;dN++)
           {
              drandXXstate=seedXX(NULL);                                               
              drand_arr(dim,xlocal);
              Local2Global(Kg,xlocal,x,&f,NULL);f*=(*fxn)(x,f,pvect);if(f<0){f=-f;sgn=-1;} else sgn=1;
              (*eff)++;      
           
              if(f>oldMax*gmax)  
              { 
                 if(nRecT<=nRec){ nRecT++;  events=realloc(events, sizeof(*events)*(nRecT));} 
                 events[nRec].f=f;
                 events[nRec].sgn=sgn;
                 for(i=0;i<40;i++) events[nRec].pvect[i]=pvect[i];
                 nRec++;
                 if(f>newMax) newMax=f;
              }                 
           }
           
           for(n=0,k=0;k<nRec&& cEvent < nEvents ;k++) if(events[k].f > newMax*drandXX())
           {  n++;
              (*out)(L,events[k].sgn,drandXXstate,events[k].pvect);
              if(events[k].sgn<0) *neg++;
              cEvent++;              
           }
           if(*nmax<n) *nmax=n;
//           printf("nRec=%d cor=%8d  old=%E new=%E   L=%8d  %s N=%8d dN=%8d  \n", nRec, n, oldMax,newMax,L,state, (int)(Ntry*oldMax/sum) ,dN  );
           if(n>1) *mult+=n-1;           
           Ntry+=dN*(newMax-oldMax)/sum;
           { ds=newMax-oldMax;
             for(cCube=L;cCube<nCubes;cCube++) smax[cCube]+=ds;
             sum+=ds;
           }
      }                                                                   
   }    
   *eff=nEvents/(*eff);
   { double mem=smax[0],mem_;
     for(cCube=1;cCube<nCubes;cCube++) { mem_=smax[cCube]; smax[cCube]-=mem; mem=mem_;} 
   } 
   free(events);
   return cEvent;
} /* vegas_ */
