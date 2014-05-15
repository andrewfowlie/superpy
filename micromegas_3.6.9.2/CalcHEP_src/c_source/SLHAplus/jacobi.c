#include"SLHAplus.h"

/* sMT is a macro  which  presents numeration of independent elements of symmetry matrix */
#define sMT(i,j) ((i)*dim-((i)*((i)+1))/2 +(j))

/* NMAX defines  maximum number of cycles used for matrix diagonalization  */
#define NMAX 50

/* leftMULTr, rightMULTr, leftMULTc, rightMULTc - macros for matrix multiplication 
   on 2-dimension rotation.  left/right    distinguish left/right multiplication,
   The 'r'/'c' final characters distinguish  real and complex matrices   
*/

#define leftMULTr(A,i1,i2,c1,s) { REAL a1,a2; int k; for(k=0;k<dim;k++)\
  { a1=A[i1*dim+k]; a2=A[i2*dim+k]; A[i1*dim+k]+=c1*a1+s*a2; A[i2*dim+k]+=c1*a2-s*a1;}}

#define rightMULTr(A,i1,i2,c1,s) { REAL a1,a2; int k; for(k=0;k<dim;k++)\
  { a1=A[k*dim+i1]; a2=A[k*dim+i2]; A[k*dim+i1]+=c1*a1-s*a2 ; A[k*dim+i2]+=c1*a2+s*a1;}}

#define leftMULTc(A,i1,i2,c1,s) { COMPLEX a1,a2,sc=conj(s); int k; for(k=0;k<dim;k++)\
  { a1=A[i1*dim+k]; a2=A[i2*dim+k]; A[i1*dim+k]+=c1*a1+s*a2 ; A[i2*dim+k]+=c1*a2-sc*a1;}}

#define rightMULTc(A,i1,i2,c1,s) { COMPLEX a1,a2,sc=conj(s); int k; for(k=0;k<dim;k++)\
  { a1=A[k*dim+i1]; a2=A[k*dim+i2]; A[k*dim+i1]+=c1*a1-sc*a2 ; A[k*dim+i2]+=c1*a2+s*a1;}}
  
/*
   Ar/Ac are used to create a working copies of incoming matrices.
   It allows to keep incoming matrices unchanged
*/


static REAL  *Ar=NULL;  static int dimr=0;
static COMPLEX *Ac=NULL;  static int dimc=0;


/* rJacobi, cJacobiH, rJacobiA, cJacobiA, cJacobiS  routines for matrix diagonalizing 
   based of Jacobi algorithm. 
   See specification of parameters in the manuscript (arXiv:1008.0181).
   See details of algorithm  in "Numerical Recipes in C" and arXiv:physics/0607103.   
*/

int rJacobi(REAL* as, int dim, REAL *d, REAL* v)
{
   int N,i1,i2,i,j;
   REAL tresh,theta,t,s,h,g,c,dc;
   
   if(dimr<dim){ Ar=realloc(Ar,dim*dim*sizeof(REAL)); dimr=dim; }  
   for(i=0;i<dim;i++) for(j=i;j<dim;j++) Ar[i*dim+j]=Ar[j*dim+i]=as[sMT(i,j)];
   for(i=0;i<dim*dim;i++) v[i]=0; 
   for(i=0;i<dim;i++) {v[i+dim*i]=1; d[i]=Ar[i+dim*i]; }

   for(N=0;N<NMAX;N++) 
   {  REAL sm=0;
      for(i1=0;i1<dim-1;i1++) for(i2=i1+1;i2<dim;i2++) sm += fabs(Ar[i1+dim*i2]); 
      if(sm == 0) break;
      if(i < 4) tresh= 0.2*sm/(dim*dim); else tresh=0;
      for(i1=0;i1< dim-1;i1++) for(i2=i1+1;i2<dim;i2++) 
      {
         g= fabs(Ar[i1+dim*i2]);
         if(i>4 &&  d[i1]+100*g== d[i1] && d[i2]+100*g==d[i2])  Ar[i1+dim*i2]=0;
         else if(g > tresh) 
         {
            h=d[i2]-d[i1];
            if((fabs(h)+g) == fabs(h)) t=(Ar[i1+dim*i2])/h; else 
            {
               theta= (0.5*h/(Ar[i1+dim*i2]));
               t= 1/(fabs(theta)+sqrt(1.0+theta*theta));
               if (theta < 0) t = -t;
            }
            c= 1/sqrt(1+t*t);
            dc=-t*t*c*c/(1+c);
            s=t*c;
            h=t*Ar[i1+dim*i2];
            d[i1] -= h;
            d[i2] += h;
            
            rightMULTr(Ar,i1,i2,dc,s)
            leftMULTr(Ar,i1,i2,dc,(-s)) 
            rightMULTr(v,i1,i2,dc,s)

            Ar[i1+dim*i2]=Ar[i1*dim+i2]=0;
            Ar[i1+dim*i1]=d[i1];
            Ar[i2*dim+i2]=d[i2];              
         }
      }
   }
  
   for(i=0;i<dim-1; )if(fabs(d[i])>fabs(d[i+1]))
   {  int i1=i+1;
      REAL mr=d[i]; d[i]=d[i1]; d[i1]=mr;
      for(j=0;j<dim;j++)
      {  
         mr=v[j*dim+i]; v[j*dim+i]=v[j*dim+i1]; v[j*dim+i1]=mr;
      }
      if(i)i--; else  i=1;
   } else i++;
   
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { REAL mem=v[j*dim+i]; v[j*dim+i]=v[i*dim+j];v[i*dim+j]=mem;}
   
   
         
   if(N<NMAX) return 0; else return 1;
}


int cJacobiH(COMPLEX* ah, int dim, REAL *d, COMPLEX* v)
{
   int N,i1,i2,i,j;
   REAL tresh,theta,sm,g,h,c,dc,t;
   COMPLEX s;

   if(dimc<dim) { Ac=realloc(Ac,dim*dim*sizeof(COMPLEX)); dimc=dim; }     
   for(i=0;i<dim;i++) for(j=i;j<dim;j++) Ac[i*dim+j]=ah[sMT(i,j)];
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) Ac[j*dim+i]=conj(ah[sMT(i,j)]);
   
   for(i=0;i<dim*dim;i++) v[i]=0; 
   for(i=0;i<dim;i++) {v[i+dim*i]=1; d[i]=Ac[i+dim*i];}
   for(N=0;N<NMAX;N++) 
   {
      sm=0;
      for(i1=0;i1<dim-1;i1++) for(i2=i1+1;i2<dim;i2++) sm += cabs(Ac[i1+dim*i2]); 
      if(sm == 0)  break; 
      if(i < 4) tresh= 0.2*sm/(dim*dim); else tresh=0;
      for(i1=0;i1< dim-1;i1++) for(i2=i1+1;i2<dim;i2++) 
      {  
        g= cabs(Ac[i1+dim*i2]);         
        if(i>4 && d[i1]+100*g==d[i1] && d[i2]+100*g==d[i2]) Ac[i1+dim*i2]=0;
        else if(g > tresh) 
        {
            h=d[i2]-d[i1];
            if(h+g == h) t=g/h; else 
            {
               theta= (0.5*h/g);
               t= 1/(fabs(theta)+sqrt(1.0+theta*theta));
               if (theta < 0) t = -t;
            }
            c= 1/sqrt(1+t*t);
            s=t*c*cexp(-I*carg(Ac[i1+dim*i2]));
            dc=-t*t/(1+t*t)/(1+c);
            h=t*g;
            d[i1] -= h;
            d[i2] += h;

            rightMULTc(Ac,i1,i2,dc,s) 
            leftMULTc(Ac,i1,i2, dc,(-s)) 
            rightMULTc(v,i1,i2,dc,s)  

            Ac[i1+dim*i2]=Ac[i1*dim+i2]=0;
            Ac[i1+dim*i1]=d[i1];
            Ac[i2+dim*i2]=d[i2];                    
        }
      }
   }

   for(i=0;i<dim-1;) if(fabs(d[i])>fabs(d[i+1]))
   {  int i1=i+1;
      REAL  mr=d[i];  d[i]=d[i1]; d[i1]=mr;
       for(j=0;j<dim;j++)
       {  
         COMPLEX mc=v[j*dim+i]; v[j*dim+i]=v[j*dim+i1]; v[j*dim+i1]=mc;
       }
       if(i) i--; else  i=1;
   } else i++;
   
   for(i=0;i<dim*dim;i++)   v[i]=conj(v[i]);
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { COMPLEX mem=v[j*dim+i]; v[j*dim+i]= v[i*dim+j];v[i*dim+j]=mem;}
     
   if(N<NMAX) return 0; else return 1;
}

int rJacobiA(REAL*aa, int dim, REAL*d, REAL*u, REAL*v)
{
  int N,i1,i2,i,j,l;
  REAL tresh,sm,g;

  if(dimr<dim) { Ar=realloc(Ar,dim*dim*sizeof(REAL)); dimr=dim; }
  for(i=0;i<dim*dim;i++) Ar[i]=aa[i];
  for(i1=0;i1<dim;i1++) for(i2=0;i2<dim;i2++) {l=i1+dim*i2; u[l]=v[l]=0;} 
  for(i1=0;i1<dim;i1++) {l=i1+dim*i1;u[l]=v[l]=1;}

  for(N=0;N<NMAX;N++) 
  { 
     sm=0;
     for(i1=0;i1<dim-1;i1++)for(i2=i1+1;i2<dim;i2++)sm+=fabs(Ar[i1+dim*i2])+fabs(Ar[i2+dim*i1]); 
     if(sm == 0.0) break;
     if(i < 4) tresh= 0.1*sm/(dim*dim); else tresh=0;
     for(i1=0;i1< dim-1;i1++) for(i2=i1+1;i2<dim;i2++) 
     {  REAL a11,a12,a21,a22;

        a11=Ar[i1*dim+i1];
        a12=Ar[i1*dim+i2];
        a21=Ar[i2*dim+i1];
        a22=Ar[i2*dim+i2];

        g= fabs(a12)+fabs(a21);
        if(i>4 && a11+100*g== a11 && a22+100*g==a22) { Ar[i1+dim*i2]=Ar[i2+dim*i1]=0; }
        else if( g  > 2*tresh) 
        {
           REAL m11,m12,m22,m21,tu,cu,tv,cv,su,sv,dcu,dcv;
           m11=a11*a11+a12*a12;    /* m=a*at */
           m22=a22*a22+a21*a21;
           m12=a11*a21+a22*a12;
           if(2*fabs(m12) < fabs(m11-m22))
           {  g=-2*m12/(m11-m22);
              tu=g/(1+sqrt(1+g*g)); 
           }            
           else  if(2*fabs(m12)> fabs(m11-m22))
           { g=-(m11-m22)/m12/2;
             tu=1/(fabs(g)+sqrt(1+g*g));
             if(g<0) tu=-tu;      
           } else
           { tu=1/(1+sqrt(2));
             if(m12*(m11-m22)>0) tu=-tu;
           }

           cu=1/sqrt(1+tu*tu);
           su=cu*tu;
           m11=a11-tu*a21;    /* m = ut * a */
           m12=a12-tu*a22;
           m21=a21+tu*a11;
           m22=a22+tu*a12;
           if(m11 != 0) {tv=-m12/m11; cv=1/sqrt(1+tv*tv); sv=tv*cv;}
           else if(m22 != 0) { tv=m21/m22; cv=1/sqrt(1+tv*tv);sv=tv*cv;}
           else {cv=0;sv=1;}

           dcu=-tu*tu*cu*cu/(1+cu);
           dcv=-tv*tv*cv*cv/(1+cv);

           rightMULTr(v,i1,i2,dcv,sv)     /* V' =  V v   */
           rightMULTr(u,i1,i2,dcu,su)     /* U' =  U u   */
 
           leftMULTr(Ar,i1,i2,dcu,(-su))   /* A' = ut A   */
           rightMULTr(Ar,i1,i2,dcv,sv)     /* A' =  A v   */
           
           Ar[i2*dim+i1]=0; Ar[i1*dim+i2]=0;
        }
     }
   }

   for(i=0; i<dim; i++) d[i]=Ar[i+dim*i];
   for(i=0; i<dim; i++) if(d[i]<0) {for(j=0;j<dim;j++) v[j*dim+i]*=-1; d[i]*=-1;}
  
   for(i=0;i<dim-1; ) if(d[i]>d[i+1])
   {   int i1=i+1;
       REAL mr=d[i]; d[i]=d[i1];d[i1]=mr;
       for(j=0;j<dim;j++)
       {  
         mr=v[j*dim+i]; v[j*dim+i]=v[j*dim+i1]; v[j*dim+i1]=mr;
         mr=u[j*dim+i]; u[j*dim+i]=u[j*dim+i1]; u[j*dim+i1]=mr;
       }
       if(i) i--; else i=1;
   } else i++;
   
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { REAL mem=v[j*dim+i]; v[j*dim+i]=v[i*dim+j];v[i*dim+j]=mem;}
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { REAL mem=u[j*dim+i]; u[j*dim+i]=u[i*dim+j];u[i*dim+j]=mem;} 

   if(N<NMAX) return 0; else return 1;
}

int cJacobiA(COMPLEX*aa, int dim, REAL*d, COMPLEX*u, COMPLEX*v)
{
  int N,i1,i2,i,j,l;
  REAL  tresh,sm;
  COMPLEX g;
  
  if(dimc<dim) { Ac=realloc(Ac,dim*dim*sizeof(COMPLEX)); dimc=dim; }  
  for(i=0;i<dim*dim;i++) Ac[i]=aa[i];   
  for(i1=0;i1<dim;i1++) for(i2=0;i2<dim;i2++) {l=i1+dim*i2; u[l]=v[l]=0;} 
  for(i1=0;i1<dim;i1++) {l=i1+dim*i1;u[l]=v[l]=1;}

  for(N=0;N<NMAX;N++) 
  { 
     sm=0;
     for(i1=0;i1<dim-1;i1++)for(i2=i1+1;i2<dim;i2++)sm+=cabs(Ac[i1+dim*i2])+cabs(Ac[i2+dim*i1]); 
     if(sm == 0) break;
     if(i < 4) tresh= 0.1*sm/(dim*dim); else tresh=0;
     for(i1=0;i1< dim-1;i1++) for(i2=i1+1;i2<dim;i2++) 
     {  COMPLEX a11,a12,a21,a22;

        a11=Ac[i1*dim+i1];
        a12=Ac[i1*dim+i2];
        a21=Ac[i2*dim+i1];
        a22=Ac[i2*dim+i2];

        g= cabs(a12)+cabs(a21);
        if(i>4 && a11+100*g==a11 && a22+100*g==a22) { Ac[i1+dim*i2]=Ac[i2+dim*i1]=0; }
        else if( cabs(g)  > tresh) 
        {
           COMPLEX m11,m12,m22,m21,cu,tu,tv,cv,su,sv,dcu,dcv;
                   
           m11=a11*conj(a11)+a12*conj(a12);    /* m=a*at */
           m22=a22*conj(a22)+a21*conj(a21);
           m12=a11*conj(a21)+a12*conj(a22);
           m21=conj(m12);
           if(4*cabs(m12) < cabs(m11-m22))
           {  g= 2*m12/(m11-m22);
              tu=-g/(1+csqrt(1+g*g*m21/m12)); 
           }                      
           else if(cabs(m12)!=0)
           {  COMPLEX sqr;
              g=(m11-m22)/(2*m12);
              sqr=csqrt(g*g+m21/m12);
              tu= (cabs(g+sqr) > cabs(g-sqr))? -1/(g+sqr):-1/(g-sqr);
           } else tu=0;           

           cu=1/sqrt(1+tu*conj(tu));
           dcu=-tu*conj(tu)*cu*cu/(1+cu);
           
           su=cu*tu;
           m11=a11-tu*a21;    /* m = ut * a */
           m12=a12-tu*a22;
           m21=a21+conj(tu)*a11;
           m22=a22+conj(tu)*a12;
           if(m11 != 0) {tv=-m12/m11; cv=1/sqrt(1+tv*conj(tv)); sv=tv*cv;}
           else if(m22 != 0) { tv=m21/m22; cv=1/csqrt(1+tv*conj(tv));sv=tv*cv;}
           else {cv=0;sv=1;}

           
           dcv=-tv*conj(tv)*cv*cv/(1+cv);

           rightMULTc(v,i1,i2,dcv,sv)
           rightMULTc(u,i1,i2,dcu,su)
           leftMULTc(Ac,i1,i2, dcu,(-su))
           rightMULTc(Ac,i1,i2,dcv,sv)
           Ac[i2*dim+i1]=0; Ac[i1*dim+i2]=0;
        }
     }
   }
   for(i=0;i<dim;i++)
   {   COMPLEX e=cexp(-I*carg(Ac[i*dim+i]));   
       for(j=0;j<dim;j++) v[j*dim+i]*=e;
       d[i]=Ac[i+dim*i]*e;
   }

     for(i=0;i<dim-1;) if(d[i]>d[i+1])
     { int i1=i+1;
       REAL mr=d[i]; d[i]=d[i1]; d[i1]=mr;
       for(j=0;j<dim;j++)
       {  
        COMPLEX mc=v[j*dim+i]; v[j*dim+i]=v[j*dim+i1]; v[j*dim+i1]=mc;
                mc=u[j*dim+i]; u[j*dim+i]=u[j*dim+i1]; u[j*dim+i1]=mc;
       }
       if(i) i--; else i=1;
     } else i++;
      

   for(i=0;i<dim*dim;i++)   v[i]=conj(v[i]);
   for(i=0;i<dim*dim;i++)   u[i]=conj(u[i]);
      
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { COMPLEX mem=v[j*dim+i]; v[j*dim+i]=v[i*dim+j];v[i*dim+j]=mem;}
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { COMPLEX mem=u[j*dim+i]; u[j*dim+i]=u[i*dim+j];u[i*dim+j]=mem;} 
            
      
   if(N<NMAX) return 0; else return 1;
}


int cJacobiS(COMPLEX*as, int dim, REAL*d, COMPLEX*v)
{
  int N,i1,i2,i,j;
  REAL tresh;

  if(dimc<dim) { Ac = realloc(Ac,dim*dim*sizeof(COMPLEX)); dimc=dim; }  
  for(i=0;i<dim;i++)   for(j=i;j<dim;j++)    Ac[i*dim+j]=as[sMT(i,j)];
  for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++)  Ac[j*dim+i]=Ac[i*dim+j];
  for(i=0;i<dim*dim;i++) v[i]=0; 
  for(i=0;i<dim;i++)  v[i*dim+i]=1;
  for(N=0;N<NMAX;N++) 
  {  REAL sm=0;
     for(i=0;i<dim-1;i++)for(j=i+1;j<dim;j++) sm+=cabs(Ac[i+dim*j]);
      
     if(sm == 0) break;
     if(i < 4) tresh= 0.2*sm/(dim*dim); else tresh=0;
     for(i1=0;i1< dim-1;i1++) for(i2=i1+1;i2<dim;i2++) 
     {  COMPLEX a11,a12,a22;
        a11=Ac[i1*dim+i1];  a12=(Ac[i1*dim+i2]+Ac[i2*dim+i1])/2;  
                            a22=Ac[i2*dim+i2]; 
        if(i>4 && a11+100*a12==a11 && a22+100*a12==a22){Ac[i1*dim+i2]=Ac[i2*dim+i1]=0;}
        else if(cabs(a12)  >  tresh) 
        {  REAL cv,dcv,tau,nn,fi,fa;
           COMPLEX ff,sv;
           nn=a11*conj(a11)-a22*conj(a22);
           ff=conj(a11)*a12+a22*conj(a12);
           fi=carg(ff);
           fa=cabs(ff);
             
           if( fabs(nn) < fa)
           { REAL R= 0.5*nn/fa;
             if(R>0) tau=-1/(R+sqrt(R*R+1));else tau=-1/(R-sqrt(R*R+1));
           }            
           else 
           { REAL r=2*fa/nn;
             tau=-r/(1+sqrt(1+r*r));
           }
           
           cv=1/sqrt(1+tau*tau);
           dcv=-tau*tau*cv*cv/(1+cv);
           sv=cv*cexp(I*fi)*tau;
                    
           rightMULTc(v,i1,i2,dcv,sv)              /* V' =  V v   */           
           leftMULTc(Ac,i1,i2,dcv,(-conj(sv)))     /* A' = ut A   */
           rightMULTc(Ac,i1,i2,dcv,sv)             /* A' =  A v   */
           Ac[i2*dim+i1]=0; Ac[i1*dim+i2]=0;
        }
     }
   }
   
   for(i=0;i<dim;i++)
   {   REAL fi=carg(Ac[i*dim+i]);
       COMPLEX e=cexp(-I*fi/2);   
       for(j=0;j<dim;j++) v[j*dim+i]*=e;
       d[i]=Ac[i+dim*i]*cexp(-I*fi);
   }

   for(i=0;i<dim-1; ) if(d[i]>d[i+1])
   { int i1=i+1;
     REAL mr=d[i]; d[i]=d[i1]; d[i1]=mr;
     for(j=0;j<dim;j++)
     {  
        COMPLEX mc=v[j*dim+i]; v[j*dim+i]=v[j*dim+i1]; v[j*dim+i1]=mc;
     }
     if(i) i--; else i=1;
   } else i++;

   for(i=0;i<dim*dim;i++)   v[i]=conj(v[i]);
   for(i=0;i<dim-1;i++) for(j=i+1;j<dim;j++) { COMPLEX mem=v[j*dim+i]; v[j*dim+i]=v[i*dim+j];v[i*dim+j]=mem;}

      
   if(N<NMAX) return 0; else return 1;
}
