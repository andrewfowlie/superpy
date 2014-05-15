#include<math.h>
#include"pmodel.h"

static void r_simpson( double(*func)(double),double * f,double a,double b, 
double eps, double * aEps, double * ans, double * aAns, int *deepness)
{
  double f1[5];
  int i;
  int d1=*deepness+1,d2=*deepness+1;
  double s1,s2,s3,e_err;

/*printf("a=%E b=%E d=%d\n",a,b,*deepness);*/

  s1=(f[0]+4*f[4]+f[8])/6;
  s2=(f[0]+4*f[2]+2*f[4]+4*f[6]+f[8])/12;
  s3=(f[0]+4*f[1]+2*f[2]+4*f[3]+2*f[4]+4*f[5]+2*f[6]+4*f[7]+f[8])/24;


  e_err=eps*fabs(s3);
  i=0;
  if( ( fabs(s3-s2) < e_err && fabs(s3-s1) < 16*e_err)) i=1; else
  if( fabs(s3-s2)*(b-a) < 0.1*(*aEps) && fabs(s3-s1)*(b-a) < 1.6*(*aEps)) 
  { i=1;  *aEps -= fabs((s3-s2)*(b-a));}
  
  if(i || *deepness>20)
  { *ans+=s3*(b-a);
    *aAns+=(fabs(f[0])+4*fabs(f[2])+2*fabs(f[4])+4*fabs(f[6])+fabs(f[8]))
          *fabs(b-a)/12;
    return ;
  }
  
  for(i=0;i<5;i++) f1[i]=f[4+i];
  for(i=8;i>0;i-=2)f[i]=f[i/2];

  for(i=1;i<8;i+=2) f[i]=(*func)(a+i*(b-a)/16);

  r_simpson(func,f,a,(a+b)/2,eps,aEps,ans,aAns,&d1);
  for(i=0;i<5;i++) f[2*i]=f1[i];
  for(i=1;i<8;i+=2) f[i]=(*func)((a+b)/2+i*(b-a)/16);
  r_simpson(func, f,(a+b)/2,b,eps,aEps,ans,aAns,&d2);
  if(d1>d2) *deepness=d1; else *deepness=d2;   
}

double simps( double (*func)(double),double a,double b, double  eps)
{
  double f[9];
  double aEps; /* absolute error  */
  int i;	
  
  aEps=0;
  if(a==b) return 0;
  for(i=0;i<9;i++) 
  { f[i]=(*func)(a+i*(b-a)/8); aEps +=fabs(f[i]);}
  if(aEps==0.)  return 0;
  eps=eps/2;
  aEps = eps*aEps*fabs(b-a)/9;

  for(;;)
  {  double ans=0., aAns=0.; 
     int deepness=1;
     r_simpson(func,f,a,b,eps,&aEps,&ans,&aAns,&deepness);
     if(5*aAns*eps > aEps) return ans;
     for(i=0;i<9;i++)  f[i]=(*func)(a+i*(b-a)/8);
     aEps=aAns*eps;
  }  
}
