#include"num_in.h"

REAL Helicity[2]={0.,0.};
REAL HelicityN[2]={0.,0.};

double Fmax=0; 
int  CalcConst=0;


#ifdef forMICROMEGAS
double (*loopFF__)(double,double, double,double)=NULL;
int WIDTH_FOR_OMEGA=0;
static int stype(int nin, char * momStr)
{ 
  int c,i;
  if(nin==1) return 1;
  for(i=0,c=0; momStr[i];i++) if (momStr[i]<=2) c++;
  if(c&1) return 0; else return 1;
}
     
#endif

int indx_(int k,int l)
{
  int i,j;
  if(k<l) {i=k;j=l;} else {i=l;j=k;}
  return i+(j*(j-1))/2;
}

void sprod_(int ntot, REAL * momenta,REAL*DP)
{
   int k, i,j;
   for (i = 0; i < ntot-1; ++i)
   {
      for (j = i + 1; j < ntot; ++j)
      {  REAL *sum=DP+indx_(i,j);
         REAL *v1=momenta+4*(i);
         REAL *v2=momenta+4*(j);
         *sum=*v1**v2;
         for(k=1;k<=3;k++) (*sum)-=v1[k]*v2[k];
      }
   }
}


static REAL sqrMom(int nin, char * momnum, REAL * lv)
{  char * ii;
   REAL s[4]={0,0,0,0};
   char nin_char;
   nin_char=nin;

   ii=momnum;
   while(*ii)
   {  int k;
      if(*ii>nin_char) for(k=0;k<4;k++) s[k]-=lv[k+4*(*ii-1)];
      else             for(k=0;k<4;k++) s[k]+=lv[k+4*(*ii-1)];
      ii++;
   }
   return (s[0]-s[1])*(s[0]+s[1])-s[2]*s[2]-s[3]*s[3];
}


int prepDen(int nden, int nin,double BWrange2, REAL*dmass,REAL*dwidth, char** Qtxt,REAL*momenta,
     REAL*Q0,COMPLEX*Q1,REAL*Q2)
{ static int ndenmem=0;
  static REAL computer_eps=1;
  double s0max; 
  int i;
  int err=0;

  if(computer_eps>0.5)
  {  REAL one=1, one_plus_eps;
     do{ computer_eps=computer_eps/2; one_plus_eps=one+computer_eps;}
     while( one_plus_eps !=one); computer_eps*=2;
  }

  for(i=0,s0max=0;i<nin;i++) s0max+=momenta[4*i];
  s0max=computer_eps*s0max*s0max;
  
  for(i=1;i<= nden;i++)
  { 
#ifdef forMICROMEGAS
   int sgn=0;
   if(loopFF__)
   {      if(strcmp(Qtxt[i],"\1\2")==0 ) sgn= 1;
     else if(strcmp(Qtxt[i],"\1\4")==0 ) sgn=-1;
   }
/* printf("mom=%d %d %d\n",Qtxt[i][0],Qtxt[i][1],Qtxt[i][2]); */
   if(sgn)
   { double mq2=momenta[0]*momenta[0],mne2=momenta[4]*momenta[4];
     int k;
     for(k=1;k<4;k++){ mq2 -=momenta[k]*momenta[k]; 
                       mne2-=momenta[k+4]*momenta[k+4];}
     Q1[i]=(*loopFF__)(sgn,sqrt(fabs(mq2)),dmass[i],sqrt(fabs(mne2)));
     Q2[i]=Q1[i]*Q1[i];
     Q0[i]=1;
     continue;
   } 
   if(WIDTH_FOR_OMEGA && !stype(nin,Qtxt[i])) dwidth[i]=0.01*fabs(dmass[i]);
   
#endif

  Q1[i]=dmass[i]*dmass[i]-sqrMom(nin,Qtxt[i],momenta);   
  if(dwidth[i])
  {  REAL w,w2, q2=Q1[i]*Q1[i];
     w=dmass[i]*dwidth[i];
     w2=w*w;
     if(q2>BWrange2*w2) {if(q2<(BWrange2+1)*w2) q2=(BWrange2+1) *w2; w2=0; }

     Q2[i]=1/(q2+w2);
     Q0[i]=Q2[i]*Q1[i]*Q1[i];
     Q1[i]=1/(Q1[i] + I*sqrt(w2));
  } else
  {  if(cabs(Q1[i]) < 10*s0max) err=2;
     if(!Q1[i]) Q1[i]=s0max;
     Q1[i]=1/Q1[i];
     Q2[i]=Q1[i]*Q1[i];
     Q0[i]=1;
  }
  }
  return err;
}
