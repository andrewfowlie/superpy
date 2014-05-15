#include"micromegas.h"
#include"micromegas_aux.h"

/* Numerical recipes codes */

static double*dym,*dyt,*yt,*dysav,*ysav,*ytemp;
static int RKQCprnFlag=1;

static void rk4(double*y, double*dydx, int n, double x,double h,double * yout,
    void (*derivs)(double,double*,double*))
{
   int i;
   double hh=h/2, h6=h/6, xh=x+hh;

   for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
   (*derivs)(xh,yt,dyt);
   for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
   (*derivs)(xh,yt,dym);
   for (i=0;i<n;i++) { yt[i]=y[i]+h*dym[i]; dym[i] += dyt[i]; }
   (*derivs)(x+h,yt,dyt);
   for (i=0;i<n;i++) yout[i]=(y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]));
}


#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666		/* 1.0/15.0 */
#define SAFETY 0.9
#define ERRCON 6.0e-4


static int rkqc(double * y, double * dydx, int n, double * x, double htry, 
     double eps,  double * yscal,  double* hdid, double* hnext,
     void (*derivs)(double,double *,double *))
{
   int i;
   double xsav=(*x),h=htry;

   for(i=0;i<n;i++) {ysav[i]=y[i]; dysav[i]=dydx[i];}

   for (;;) 
   {  double hh= 0.5*h, errmax=0;
      rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
      *x=xsav+hh; 
      (*derivs)(*x,ytemp,dydx);
      rk4(ytemp,dydx,n,*x,hh,y,derivs);
      *x=xsav+h;
      if (*x == xsav && RKQCprnFlag) 
      { printf("Step size too small in routine RKQC\n"); RKQCprnFlag=0;
        return 1;
      }
      rk4(ysav,dysav,n,xsav,h,ytemp,derivs);
      for (i=0;i<n;i++) 
      {  double temp;
         ytemp[i]=y[i]-ytemp[i];
         if(!isfinite( ytemp[i])) { errmax=ytemp[i]; break;} 
         temp= fabs(ytemp[i]/yscal[i]);
	 if (errmax < temp) errmax=temp;
      }
      if(!isfinite(errmax)){ h=h/10; continue;}
      errmax /= eps;

      if (errmax <= 1.0) 
      {
	 *hdid=h;
	 *hnext=((errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4*h));
	  break;
      }
      {  double h_=(SAFETY*h*exp(PSHRNK*log(errmax)));
//              if(h_/h>10 ) h=10*h;
//         else 
         if(h_/h<0.1) h=0.1*h; else  h=h_;
      }  
   }
   for (i=0;i<n;i++) y[i] += (double) (ytemp[i]*FCOR);
   return 0;
}

#define MAXSTP 100000
#define TINY 1.0e-30


int  odeint(double * ystart, int nvar, double x1, double x2, double eps, 
         double h1, void (*derivs)(double,double *,double *))
{
   int nstp,i;
   double x,hnext,hdid,h;

   double *yscal,*y,*dydx;


   double ** allAlloc[9]={NULL,NULL,NULL,&dym,&dyt,&yt,&dysav,&ysav,&ytemp};
   allAlloc[0]=&yscal;allAlloc[1]=&y; allAlloc[2]=&dydx;
   for(i=0;i<9;i++) *allAlloc[i]=(double*)malloc(nvar*sizeof(double)); 

   RKQCprnFlag=1; 
   x=x1;
   h=((x2 > x1) ? fabs(h1) : -fabs(h1));
   for (i=0;i<nvar;i++) y[i]=ystart[i];
   for (nstp=1;nstp<=MAXSTP;nstp++) 
   {
      (*derivs)(x,y,dydx);
      for (i=0;i<nvar;i++) yscal[i]=(fabs(y[i])+fabs(dydx[i]*h)+TINY);

      if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
      if(rkqc(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs))break;
      
      if ((x-x2)*(x2-x1) >= 0.) 
      {
         for (i=0;i<nvar;i++) ystart[i]=y[i];
         for(i=0;i<9;i++) free(*allAlloc[i]);
         return 0;
      }
      h=hnext;
   }
   for(i=0;i<9;i++) free(*allAlloc[i]);
   return 1;
}


static double bessi0(double  x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans= (1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))));
	} else {
		y=3.75/ax;
		ans= ((exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2)))))))));
	}
	return ans;
}

static double bessi1(double x)
{
	double ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=(ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))));
	} else {
		y=3.75/ax;
		ans=(0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
			-y*0.420059e-2)));
		ans= (0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*ans)))));
		ans *= ((exp(ax)/sqrt(ax)));
	}
	return x < 0.0 ? -ans : ans;
}

double bessk0(double x)
{
/*
   M.Abramowitz and I.A.Stegun, Handbook of Mathematical Functions,
   Applied Mathematics Series vol. 55 (1964), Washington.
*/
         
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return  ans;
}


double bessk1(double x)
{
	double y,ans;

	if (x <= 2.0) {
		y=x*x/4.0;
		ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y=2.0/x;
		ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return  ans;
}

double bessk2(double x)
{
	double bk,bkm,bkp,tox;

	tox= 2.0/x;
	bkm=bessk0(x);
	bk=bessk1(x);
	bkp=bkm+tox*bk;
	bkm=bk;
	bk=bkp;
	return bk;
}

double K2pol(double x)
{
   if(x<0.1) return 1+ 1.875*x*(1+0.4375*x*(1-0.375*x));
   else      return bessk2(1/x)*exp(1/x)*sqrt(2/M_PI/x);
}

double K1pol(double x)
{
  if(x<0.1) return 1+ 0.375*x*(1-0.3125*x*(1+0.875*x));
  else      return bessk1(1/x)*exp(1/x)*sqrt(2/M_PI/x); 
}

static double  polintN(double x, int n,  double *xa, double *ya)
{  double z[20];
   int i,m;
   for(i=0;i<n;i++) z[i]=ya[i];
   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}

static int  leftXN(int n,int dim,  double * xa, double x)
{  int k1,k2,k3;

   k1=n/2;                         
   k2=dim-(n+1)/2-1;

   if(xa[0]< xa[dim-1])
   {              
     if(x<=xa[k1]) return 0;
     if(x>=xa[k2]) return dim-n;                   
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]>x)k2=k3; else k1=k3;
     }
   } else 
   {  
     if(x>=xa[k1]) return 0;
     if(x<=xa[k2]) return dim-n;
     while(k2-k1>1)                
     { k3=(k1+k2)/2;               
       if(xa[k3]<x)k2=k3; else k1=k3;
     }
   }
   return k1+1-n/2;
}

static int  leftX(int dim, double * xa, double x)
{  int k1,k2,k3;
                                                                                
   if(x<=xa[0]) return 0;
   if(x>=xa[dim-3]) return dim-3;
                                                                                
   k1=0;
   k2=dim-3;
                                                                                
   while(k2-k1>1)
   { k3=(k1+k2)/2;
     if(xa[k3]>x)k2=k3; else k1=k3;
   }
   return k1;
}


double polint2(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(2,n, xa, x);
   return polintN(x,2,xa+shift, ya+shift);
}


double polint3(double x, int n,  double *xa, double *ya)
{ int shift=leftX(n, xa, x);
  double ar;
  ar=polintN(x,3,xa+shift, ya+shift); 
  if(shift==0) return ar;
  if(shift==n-3 &&(  (x> xa[n-2] && xa[0]<xa[n-1])|| (x<xa[n-2] && xa[0]>xa[n-1]) )) return ar;
  shift--;
  return 0.5*( ar+ polintN(x,3,xa+shift, ya+shift));
}


double polint4(double x, int n,  double *xa, double *ya)
{ int shift=leftXN(4,n, xa, x);
   return polintN(x,4,xa+shift, ya+shift);
}


static void del(int k, int * N, double *xa, double *ya)
{  int i;
    (*N)--;
    for(i=k;i<*N;i++) { xa[i]=xa[i+1]; ya[i]=ya[i+1];}    
} 

static void ins(int k,  double x, double y, int*N,double *xa,double *ya)
{ 
  int i;
  for(i=*N;i>k;i--) { xa[i]=xa[i-1]; ya[i]=ya[i-1];}
   xa[k]=x;ya[k]=y; (*N)++;        
}

void printInterpolation(int N,double *xa, double * ya)
{ int i; for(i=0;i<N;i++) printf("{ %E %E}\n",xa[i],ya[i]);
  printf("\n");
}

int buildInterpolation( double (*Fun)(double), double x1,double x2, double eps,
int * N_, double ** xa_, double **ya_)
{  int i,cnt,N,k;
   double *xa,*ya;
   double dx0;
   
   dx0=fabs(x2-x1)*0.01;      
   N=5;
   xa=malloc(N*sizeof(double));
   ya=malloc(N*sizeof(double));
   
   for(i=0;i<5;i++) {xa[i]=x1+ (x2-x1)/4*i; ya[i]=Fun(xa[i]);}  

   for(cnt=1;cnt;)
   { cnt=0; 
     for(i=0; i<N; i++)
     {  double x=xa[i], y=ya[i], yy;
        if(i<N-1 && fabs(xa[i+1]-xa[i]) < dx0) continue; else
        if(i>0   && fabs(xa[i]-xa[i-1]) < dx0) continue;
                           
        del(i,&N,xa,ya);
        yy=polint4(x, N, xa, ya);
        ins(i, x, y, &N,xa, ya);
        if( (eps>0 && fabs(yy-y) > eps) || (eps<0 && fabs(yy-y)> -eps*(y)))  
        {
           cnt=1;
           xa=realloc(xa,sizeof(double)*(N+1));
           ya=realloc(ya,sizeof(double)*(N+1));
           
           if(i==0)   k=1;  
           else if(i==N-1) k=N-1;
           else if(fabs(xa[i-1]-xa[i])< fabs(xa[i]-xa[i+1])) k=i+1;
           else  k=i;
                                                            
           x=(xa[k-1]+xa[k])/2;                                                                        
           y=Fun(x); 
           ins(k,x,y,&N,xa,ya);
           i++;     
        }
         
     }
   }   
   *N_=N;
   *xa_=xa;
   *ya_=ya;
   return 0;  
}


#define SQ(x)  ((x)*(x))

double   LintIk(int II,double MSQ,double MQ,double MNE)
{  
  double LAM,SPPM,SPMM,R1,R2,R3,del,CMD;
  double msq2=MSQ*MSQ, mq2=MQ*MQ, mne2=MNE*MNE;

  SPPM=  msq2+mq2-mne2, SPMM= msq2-mq2-mne2;
  R1  =(msq2-mq2)/mne2;
  R2  =(mq2-mne2)/msq2;
  R3  =(msq2-mne2)/mq2;

  del =2.*mne2*(mq2+msq2)-mne2*mne2-SQ(msq2-mq2);
  
  if(del>0) LAM=2.*atan(sqrt(del)/SPPM)/sqrt(del);
  if(del<0) LAM=log((SPPM+sqrt(-del))/(SPPM-sqrt(-del)))/sqrt(-del);

  switch(II)
  { 
    case 1:  
      CMD=1./del*(R2/3-2/3.*R3-5/3.+ (2*msq2-2/3.*mne2)*LAM);	
    break; 
    case 2: 
      CMD=(log(msq2/mq2)-SPMM*LAM)/2./mne2/mne2+
         ( ((mq2*mq2-mq2*msq2)/mne2-7/3.*mq2+2/3.*(mne2-msq2))*LAM+R2/3+R1+2/3.
         )/del;
    break;
    case 3:
      CMD=-3/SQ(del)*SPPM+LAM/del*(-1+6*mq2*msq2/del);
    break;	
    case 4:
      CMD=((log(msq2/mq2) - SPMM*LAM)/2/mne2-1/msq2 -mq2*SPMM/del*LAM)/mne2/mne2
      
         +( mq2/mne2/mne2-SQ(1-mq2/mne2)/msq2+0.5/mne2
            +3*mq2/del*(1 +  R1 + (-R1*mq2-2*mq2-msq2+mne2)*LAM)
          )/del;
    break;
    case 5:
     CMD=(log(msq2/mq2)-SPMM*LAM)/(2*mne2*mne2)-(LAM*(2*(msq2-mne2)+3*mq2+R1*mq2)-3-R1)/del;
     break;
    default: CMD=0.; 
  }
  return CMD;
}

double MaxGapLim(double x, double mu) 
/* S.Yellin, Phys.Rev. D66,032005(2002)

   There is a theoretical model which predicts homogenious event distribution 
   with everage number of events mu. Let experiment gets a gap bitween points 
   where according to theory x point are expected. Then the theoretical model 
   is non-confirmed with probability MaxGap   
*/  
{
  int k;
  double C0,kf;
  if(x>mu-1.E-5) return 1-exp(-mu);
  for(k=0,C0=0,kf=1;k<=mu/x; k++,kf*=k) {C0+= pow(k*x-mu,k)*exp(-k*x)*(1+k/(mu-k*x))/kf;}
  return C0;   
}

/*
int main(void)
{ double x;
  for(x=0.00001; x<1; x*=1.5)
  printf("x=%e bessk2=%e\n",x,  bessk2(x)*x*x); 

}
*/

#define BUFFSIZE 500

int readTable(char * fileName, int *Ncolumn, double **tab)
{  FILE *f;
   char buff[BUFFSIZE];   
   int nRec=0,nCol=0,nCom=0;

   f=fopen(fileName,"r");
   
   
   if(!f) return 0;

   while(fgets(buff,BUFFSIZE,f))
   { int i;
     char*ch;
     for(i=0; buff[i] && buff[i]==' ';i++);
     if(buff[i]==0 || buff[i]=='#') {nCom++; continue;}
     ch=strtok(buff," \n");
     if(ch[0]=='#' || ch[0]==0) continue;
     for(i=0;ch;i++,ch=strtok(NULL," \n"))
     { 
       if(nRec==0) {tab[i]=malloc(sizeof(double)); nCol++;} 
       else
       { if(i==nCol){fclose(f); for(i=0;i<nCol;i++) free(tab[i]);  return -(nRec+1+nCom);}
         tab[i]=realloc(tab[i],(nRec+1)*sizeof(double));
       }
       if(1!=sscanf(ch,"%lf",tab[i]+nRec)) break;
     } 
     nRec++;
   }
   fclose(f);
   if(Ncolumn) *Ncolumn=nCol;
   return nRec;
}

static double amotry(double *p, double *y, int ndim,
	     double (*f)(double *), int ilo, double fac)
{
   int i,j;
   double  ytry, fac1=(1.0-fac)/ndim;
   double * p_buff=p+(ndim+1)*ndim;
   double * p_ilo =p+ilo*ndim;
   
   for(j=0;j<ndim;j++) p_buff[j]=p_ilo[j]*fac;
   for(i=0;i<=ndim;i++)  if(i!=ilo) 
     {double *p_i=p+i*ndim;  for(j=0;j<ndim;j++) p_buff[j] +=p_i[j]*fac1;} 
   ytry=(*f)(p_buff);
   
   if (ytry > y[ilo]) 
   {  for(j=0;j<ndim;j++) p_ilo[j]=p_buff[j];
      y[ilo]=ytry;
   }
//printf("amotry returns fac=%f %E\n",fac,  ytry);   
   return ytry;
}

double amoeba(double *p, double * y, int ndim, double (*f)(double *), 
                                                    double eps, int *nCalls)
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

//printf("nCall=%d  ndim=%E\n",*nCalls,y[ilo]);   
     										     
      if((*nCalls)<=0||2*(y[ihi]-y[ilo])/(fabs(y[ilo])+fabs(y[ihi]))<eps)break;
     										     
      ytry=amotry(p,y,ndim,f,ilo,-1.0); (*nCalls)--;				     
      if (ytry >= y[ihi]) {ytry=amotry(p,y,ndim,f,ilo,2.); (*nCalls)--;}	     
      else if (ytry <= y[inlo])							     
      {										     
         ysave=y[ilo];								     
     	 ytry=amotry(p,y,ndim,f,ilo,0.5);  (*nCalls)--;
     	 if (ytry <= ysave)
     	 {  
     	    for (i=0;i<=ndim;i++)
     	    {  double * p_ihi=p+ihi*ndim;
               if (i != ihi)
     	       {  double * p_i=p+i*ndim;
     		  for(j=0;j<ndim;j++) p_i[j]=0.5*(p_i[j]+p_ihi[j]);
     		  y[i]=(*f)(p_i);
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
