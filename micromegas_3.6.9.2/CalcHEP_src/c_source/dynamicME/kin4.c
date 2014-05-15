#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include "../../include/VandP.h"
#include"SLHAplus.h"
#include"../ntools/include/vegas.h"
#include"../ntools/include/simpson.h"
#include"../num/include/alphas2.h"

#include"dynamic_cs.h"
#include "vp.h"

#define VVmassGap  70

int ForceUG=0;
decayTableStr* decayTable=NULL;


/*=============   decayPcm   and decayPcmW ================*/      

double  decayPcm(double am0,  double  am1,  double  am2)
{
  double  summ, diffm;
  summ = am1 + am2;
  diffm = am1 - am2;
  if(am0<summ) return 0;
  return sqrt((am0-summ)*(am0+ summ)*(am0-diffm)*(am0+diffm))/(am0*2);
}

static double w1_,w2_,m0_,m1_,m2_,m1_p,m2_p;
static int nGauss=0;

static double ME(double m0,double m1,double m2)
{ //return 1;
  //return  pow( m1*m2/(m1+m2)/(m1+m2),2);
  return  m2*m2*m1*m1 +(m0*m0 -m2*m2 -m1*m1)*(m0*m0 -m2*m2 -m1*m1)/8;
}

static double y11,y12,y21,y22;

static double intDecay2(double y2)
{ double m2;
  m2= m2_p*m2_p + w2_*m2_p*tan(y2);
  if(m2<=0) return 0;
  m2=sqrt(m2);
  return  ME(m0_,m1_,m2)*decayPcm(m0_,m1_,m2);
}   

static double intDecay2_(double x)
{ if(x<=0 || x>=1) return 0;
  return intDecay2(y21+x*x*x*(4-3*x)*(y22-y21))*12*x*x*(1-x)*(y22-y21);
}

static double  intDecay1(double y1)
{ 
   m1_=m1_p*m1_p + w1_*m1_p*tan(y1);
   if(m1_<0) return 0;
   m1_=sqrt(m1_);
   
   y21=atan(-m2_p/w2_);
   y22=atan( ((m0_-m1_)*(m0_-m1_)-m2_p*m2_p)/(m2_p*w2_));
   if(nGauss) return  gauss(intDecay2_,0,1,nGauss);
   else       return  simpson(intDecay2_,0,1,1.E-3);
//   return  simpson(intDecay2, atan(-m2_p/w2_), atan( ((m0_-m1_)*(m0_-m1_)-m2_p*m2_p)/(m2_p*w2_)), 1.E-3);
}

static double intDecay1_(double x)
{ if(x<=0 || x>=1) return 0;
   return intDecay1(y11+x*x*x*(4-3*x)*(y12-y11))*12*x*x*(1-x)*(y12-y11);
}
  

double decayPcmW(double m0,double m1,double m2,double w1,double w2, int N)
{  
  nGauss=N;
  m0_=m0;
  if(w1==0 && w2==0) return decayPcm(m0,m1,m2);
  else if(w1==0)
  { if(m1>m0) return 0;
    m1_=m1;
    m2_p=m2;
    w2_=w2;
    y21=atan(-m2/w2);
    y22=atan( ((m0-m1)*(m0-m1)  -m2*m2)/(m2*w2));
    if(nGauss) return gauss( intDecay2_,0,1,nGauss)/M_PI/ME(m0,m1,m2);
     else      return simpson( intDecay2_,0,1,1.E-3)/M_PI/ME(m0,m1,m2); 
  }
  else if(w2==0)
  { if(m2>m0) return 0;
    m1_=m2;
    m2_p=m1;
    w2_=w1;
    y21=atan(-m1/w1);
    y22=atan( ((m0-m2)*(m0-m2)-m1*m1)/(m1*w1));
    if(nGauss) return gauss( intDecay2_,0,1,nGauss)/M_PI/ME(m0,m1,m2);  
    else       return simpson(intDecay2_,0,1,1.E-3)/M_PI/ME(m0,m1,m2);
  }
  else 
  { w1_=w1;
    w2_=w2; 
    m1_p=m1;
    m2_p=m2; 
   y11=atan(-m1/w1); y12=atan( (m0*m0-m1*m1)/(m1*w1));

   if(nGauss) return gauss(intDecay1_,0,1,nGauss)/(M_PI*M_PI)/ME(m0,m1,m2);
   else       return simpson(intDecay1_,0,1,1E-3)/(M_PI*M_PI)/ME(m0,m1,m2);

  }
}



/*================  kinematic 1->3 and  1->4  =======================*/

static double kinematic_1_3(REAL *pmass, int i3, double m12, double xcos, REAL * P)
{ 
  double factor;
  REAL pout,mQ,chY,shY,xsin, E1,P12,P13,E2,P22,P23, m0,m1,m2,m3;
  int i,i1,i2;
  
  for(i=1;i<4;i++)if(i3!=i) {i1=i; break;}
  for(i++;i<4;i++)if(i3!=i) {i2=i; break;}
  
  m0=pmass[0];
  m1=pmass[i1];
  m2=pmass[i2];
  m3=pmass[i3];

  if(m12<=m1+m2) return 0;
  for(i=0;i<16;i++) P[i]=0;

  P[0]=m0; 
  factor=1/(64*M_PI*M_PI*M_PI*m0*m0);
  
  pout=decayPcm(m0,m12,m3);
  if(!pout) return 0;  
  P[i3*4]=sqrt(pout*pout+m3*m3); P[i3*4+3]=-pout; 

  factor*=pout;  
  
  shY=pout/m12;
  chY=sqrt(1+shY*shY);  
  pout=decayPcm(m12,m1,m2);
  if(!pout) return 0;
  factor*=pout;
  xsin=sqrt(1-xcos*xcos);
  E1=sqrt(m1*m1+pout*pout);    E2=sqrt(m2*m2+pout*pout);
  P13=xcos*pout;               P23=-P13;
  P12=xsin*pout;               P22=-P12;
  
  P[4*i1]  =chY*E1 + shY*P13;  P[4*i2]  =chY*E2 + shY*P23;
  P[4*i1+3]=shY*E1 + chY*P13;  P[4*i2+3]=shY*E2 + chY*P23;
  P[4*i1+2]=P12;               P[4*i2+2]=P22;
  
  return factor;
}

static double kinematic_1_4(REAL *pmass, double xm1, double xm2, double xcos1, double xcos2,double fi2,  REAL * P)
{ 
  double factor,M1,M2,Pcm,p1cm,p2cm,chY,shY,xsin;
  double p0,p3;
  int i,j;

  factor= 1./(pow(2*M_PI,8)*2*pmass[0]*16);
  
  M1= pmass[1]+pmass[2]+ xm1*(pmass[0]-pmass[1]-pmass[2]-pmass[3]-pmass[4]);
  M2= pmass[3]+pmass[4]+ xm2*(pmass[0]-M1-pmass[3]-pmass[4]); 

  factor*=(pmass[0]-pmass[1]-pmass[2]-pmass[3]-pmass[4])*(pmass[0]-M1-pmass[3]-pmass[4]);
  Pcm=decayPcm(pmass[0],M1,M2);            factor*=4*M_PI*Pcm/pmass[0];
  p1cm=decayPcm(M1,pmass[1],pmass[2]);     
  p2cm=decayPcm(M2,pmass[3],pmass[4]);     factor*=2*M_PI*p1cm*p2cm;
  
  
  P[0]=pmass[0]; P[1]=P[2]=P[3]=0;
  
  P[4+0]=sqrt(pmass[1]*pmass[1]+p1cm*p1cm);       P[8+0]=sqrt(pmass[2]*pmass[2]+p1cm*p1cm);
  P[4+1]=0;   
  P[4+2]=p1cm*sqrt(1-xcos1*xcos1);
  P[4+3]=p1cm*xcos1;
  
          
  P[12+0]=sqrt(pmass[3]*pmass[3]+p2cm*p2cm);       P[16+0]=sqrt(pmass[4]*pmass[4]+p2cm*p2cm);
  P[12+1]=sin(fi2)*p2cm*sqrt(1-xcos2*xcos2);
  P[12+2]=cos(fi2)*p2cm*sqrt(1-xcos2*xcos2); 
  P[12+3]=p2cm*xcos2;  
  
  for(i=1;i<=2;i++) for(j=1;j<=3;j++) P[i*8+j]=-P[i*8-4+j];
  shY=Pcm/M1;
  chY=sqrt(1+shY*shY);
  for(i=1;i<3;i++)
  { double  p0=P[4*i];
    double  p3=P[4*i+3];
    P[4*i]=chY*p0+shY*p3;
    P[4*i+3]=shY*p0+chY*p3;
  }

  shY=-Pcm/M2;            
  chY=sqrt(1+shY*shY);  
  for(i=3;i<5;i++)
  { double  p0=P[4*i];
    double  p3=P[4*i+3];
    P[4*i]=chY*p0+shY*p3;
    P[4*i+3]=shY*p0+chY*p3;
  }

//printf("Energy conservation\n");
/*
for(i=0;i<4;i++)
{ double sum=P[0+i]-P[4+i]-P[8+i]-P[12+i]-P[16+i];
  if(fabs(sum/P[0]) > 1.E-4)
  { printf("No Energy conservation %E i=%d  \n",sum/P[0],i);
    exit(22);
  }  
}    
*/
  for(i=0;i<5;i++)
  { double m;
    m=sqrt(fabs(P[4*i]*P[4*i]-P[4*i+1]*P[4*i+1]-P[4*i+2]*P[4*i+2]-P[4*i+3]*P[4*i+3]));
    if(fabs(m-pmass[i])>pmass[0]*1.E-5) { printf("wrong mass %d (%E != %E) \n",i,m,pmass[i]); exit(33);}
  }
  return factor;
}


/* ===========  Intergration ================ */

double (*sqme)(int nsub,double GG, REAL *pvect, int * err_code)=NULL;
static int  nsub_stat;
static REAL*Q=NULL;
static REAL Pmass[5];
static int i3_;
static double M_;
static double GG=1.2;

static double dWidthdCos(double xcos)
{
  double factor;
  REAL P[16];
  int err_code=0;


  factor=kinematic_1_3(Pmass,i3_,M_,xcos, P);
  
  if(factor==0) return 0;
  
//printf("xcos=%e factor=%E sqme=%e\n", xcos,factor,(*sqme)(nsub_stat,GG,P,&err_code));  
  return  factor*(*sqme)(nsub_stat,GG,P,&err_code);

}

static double dWidthdM(double M)
{ double r;

  M_=M; r= simpson(dWidthdCos,-1.,1.,1.E-4);
  
//  printf("M=%e dWidthdM= %E\n", M_,r);
  return r;
}

static double width13(numout * cc, int nsub, int * err) 
{
  int i;
  if(passParameters(cc)){ *err=4; return 0;}
  for(i=0;i<4;i++) cc->interface->pinf(nsub,1+i,Pmass+i,NULL);
  if(cc->SC ) GG=*(cc->SC); else  GG=sqrt(4*M_PI*alpha_2(Pmass[0]));
  *err=0;  
 i3_=1; 
  sqme=cc->interface->sqme;
  nsub_stat=nsub; 
  return simpson(dWidthdM,Pmass[2]+Pmass[3], Pmass[0]-Pmass[1],1.E-2);
}


static double wInt14(double *x, double w)
{
   REAL pvect[20];
   int err_code=0;
   double res;
   res=kinematic_1_4(Pmass,x[0], x[1], 2*(x[2]-0.5),2*(x[3]-0.5),2*M_PI*x[4],pvect); 
   if(res==0) return 0;  
   res*= (*sqme)(1,GG, pvect, &err_code);
   if(err_code) return 0;
   return res*8*M_PI;
}

static double width14(numout * cc, int * err) 
{
  int i;
  vegasGrid * vegPtr;
  long ncall0=10000;                     /* number of integrand calls */
  double alph=1.5;                       /* rate of grid improvement  */
  double ti;                             /* integral estimation */
  double tsi;                            /* standard deviation */
                          
  *err=0;

  if(passParameters(cc)){ *err=4; return 0;}
    
  for(i=0;i<5;i++) cc->interface->pinf(1,1+i,Pmass+i,NULL);  
  if(cc->SC ) GG=*(cc->SC); else  GG=sqrt(4*M_PI*alpha_2(Pmass[0]));
  *err=0;

  sqme=cc->interface->sqme;
  
  vegPtr=vegas_init(5,50);
  
  vegas_int(vegPtr,ncall0,alph,wInt14,&ti,&tsi); 
  vegas_int(vegPtr,ncall0,alph,wInt14,&ti,&tsi);
  vegas_finish(vegPtr);
   
  return ti;
}


                       
int pname2lib(char*pname, char * libname)
{
  int n,p;
  char buff[30];
  strcpy(buff,pname);
  p=strlen(buff)-1;
  if(buff[p]=='%') {buff[p]=0;p=1;} else p=0;
  n=pTabPos(buff);
  if(!n) {printf("Wrong particle name '''%s'''\n",pname); libname[0]=0; return 1;}
  if(p) { if(n>0) sprintf(libname,"pp%d",n); else sprintf(libname,"ap%d",-n);}
  else  { if(n>0) sprintf(libname,"p%d",n); else sprintf(libname,"a%d",-n);}   
  return 0;
}


static int decodeProcess(char *txt,int*inList,int*outList)
{ char name[20];
   char *ch_,*ch;
   int i,p;   
   ch_=strstr(txt,"->");
   if(!ch_) { inList[0]=0; ch=txt;}
   else
   { 
     for(p=0,ch=txt;; )
     { sscanf(ch," %[^,]",name);
       ch_=strstr(name,"->");
       if(ch_) *ch_=0;
       for(i=strlen(name)-1; i>=0 && name[i]==' '; i--) name[i]=0;
       inList[p]=pTabPos(name);
       if(!inList[p]) return -(p+1);
       p++;
       if(ch_) break;     
       ch=strchr(ch,',');
       if(!ch) break; else ch++;
     }  
     inList[p]=0; 
     ch=strstr(txt,"->")+2; 
   }  

   for(p=0;ch; )
   { sscanf(ch," %[^,]",name);
     for(i=strlen(name)-1;i>=0 && name[i]==' '; i--) name[i]=0;
     outList[p]=pTabPos(name);
     if(!outList[p]) return p+1;
     p++;    
     ch=strchr(ch,',');
     if(!ch) break;
     ch++;
   }
   outList[p]=0;
   return 0;
}


 
void massFilter(double M, txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  double Msum=0,dM;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; ch; ch=strchr(ch,','))
     { char buff[10];
       int n;
       char *nm;
        
       ch++;
       sscanf(ch,"%[^,]",buff); 
       
       n=pTabPos(buff);
       nm=ModelPrtcls[abs(n)-1].mass;
       if(nm[0]=='0')   
       {  dM=0;
          switch(abs(ModelPrtcls[abs(n)-1].NPDG))
          {
            case 1: case 2: dM=0.07; break;
            case 3: dM=0.3; break;
            case 4: dM=1.5; break;
            case 5: dM=5. ; break;
          }
       }  
       else  dM = fabs(*(varAddress(nm))); 
       Msum+=dM;
     } 
     lnext=lold->next;
     if(M>Msum) {lold->next=lnew; lnew=lold;}
     else {free(lold->txt); free(lold);}
     lold=lnext;
 } 
 *List=lnew;
}

void gammaGluFilter(txtList * List)
{
  txtList lold=*List, lnew=NULL,lnext;

  while(lold)
  {  int del=0,code;
     char *ch;
     ch=strstr(lold->txt,"->")+2;
     for( ; !del && ch; ch=strchr(ch,','))
     { char buff[10];
       ch++;
       sscanf(ch,"%[^,]",buff); 
       code=pNum(buff);
       if(code==22 || code ==21) { del=1;}
     } 
     lnext=lold->next;
     if(del) {free(lold->txt); free(lold);}
     else    {lold->next=lnew; lnew=lold;}
  
     lold=lnext;
 } 
 *List=lnew;
}


int process2Lib(char * process,char * lib)
{ 
  char * ch, *ch1;
  char bufflib[20];
  char *process_;
  int err=0,nX=0,pos=1;
  process_=malloc(strlen(process)+1);
  strcpy(process_,process);

  ch= strstr(process_,"->");
  ch[0]=' '; ch[1]=','; ch=ch+2; for(;*ch==' ';ch++);

  lib[0]=0;
  ch1=strtok(process_," ,"); 
  for(;ch1 && err==0; pos++)
  { if(ch1==ch) strcat(lib,"_");
    if( strcmp(ch1+1,"*x") &&  strcmp(ch1+1,"*X"))   
    { 
      err=pname2lib(ch1,bufflib);
      if(err) return pos;
      strcat(lib,bufflib);
    }  
    else if(1!=sscanf(ch1,"%d",&nX)) return pos;
    ch1=strtok(NULL," ,");
  }
  if(nX) sprintf(lib+strlen(lib),"x%d",nX);
  free(process_);
  return 0;
}


void process2Mass(char * process,double * mass)
{ 
  char * ch, *ch1;
  char *process_;
  int i;  

  process_=malloc(strlen(process)+1);
  strcpy(process_,process);
    
  ch= strstr(process_,"->");
  ch[0]=' '; ch[1]=','; ch=ch+2; for(;*ch==' ';ch++);

  ch1=strtok(process_," ,"); 
  for(i=0;ch1;i++,ch1=strtok(NULL," ,"))
  { /*if(ch1==ch) strcat(lib,"_");*/
    mass[i]=pMass(ch1);
  }
  free(process_);
}


/*======================  1->2 decay ==================*/

double pWidth2(numout * cc, int nsub)
{
  REAL pvect[12];
  double width=0.;
  REAL m1,m2,m3; 
  int i,ntot,nin,nout;
  double GG;
  procInfo1(cc,&ntot,&nin,&nout);
  if(nsub<1 ||  nsub>ntot|| nin!=1||nout !=2)  return 0;
       
  if(passParameters(cc)) return -1;
  
  cc->interface->pinf(nsub,1,&m1,NULL);
  cc->interface->pinf(nsub,2,&m2,NULL); 
  cc->interface->pinf(nsub,3,&m3,NULL);
  if(cc->SC) GG=*(cc->SC); else GG=sqrt(4*M_PI*alpha_2(m1));
  
  if(m1 >m2 + m3)
  {   int i,err_code=0; 
      double md=m2-m3;
      double ms=m2+m3;
      double pRestOut=sqrt((m1*m1 - ms*ms)*(m1*m1-md*md))/(2*m1);
      double totcoef= pRestOut/(8. * M_PI * m1*m1);
           
      for(i=1;i<12;i++) pvect[i]=0;
      pvect[0]=m1;
      pvect[7]=pRestOut;
      pvect[4]=sqrt(pRestOut*pRestOut+m2*m2);
      pvect[11]=-pRestOut;
      pvect[8]=sqrt(pRestOut*pRestOut+m3*m3);
      width = totcoef * (cc->interface->sqme)(nsub,GG,pvect,&err_code);
  }
  return width;
}

 
double decay2Info(char * pname, FILE* f)
{ int i,j,ntot;
  numout * cc;
  double wtot;
  char pname2[20],process[20],plib[20];
  char * dname[8];

  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=getMEcode(0,ForceUG,process,NULL,"",plib);
  if(!cc) return -1; 
  procInfo1(cc,&ntot,NULL,NULL); 
  if(f) fprintf(f,"\n Partial width for %s->2x decays in GeV\n",pname2); 
  for(wtot=0,i=1;i<=ntot;i++)
  { double w;
    procInfo2(cc,i,dname,NULL);
    w=pWidth2(cc,i);
    if(w!=0)
    { wtot+=w;
      if(f) fprintf(f,"%3.3s %3.3s  %.2E\n",dname[1],dname[2],w); 
    }
  }
  if(f) fprintf(f," Total width %.2E GeV\n",wtot);
  return  wtot;
}

static int chOpen(numout*cc, int k)
{  double m[3],s; 
   int pdg[3],j;
   char*name[3]; 
   for(j=0;j<3;j++) name[j]=cc->interface->pinf(k,j+1,NULL,pdg+j);
   s=pMass(name[0]);
   for(j=1;j<3;j++) s-=pMass(name[j]);

   if( pdg[0]!=23 &&  abs(pdg[0])!=24)
   for(j=1;j<3;j++) if(pdg[j]==23 && VZdecay ) s-=6; else if(abs(pdg[j])==24 && VWdecay) s-=5;
   if(s>0) return 1; else return 0;
}   


numout* xVtoxll(int Nin,int Nout,char**name,int *pdg, int lV, double *wV,  double *br)
{
  int i,err;   
  char* e_=NULL,*E_=NULL,*ne_=NULL,*Ne_=NULL,*m_=NULL,*M_=NULL,*nm_=NULL,*Nm_=NULL,*W_=NULL,*Z_=NULL;
  char processX3[50],plib13[50],exclude[20];
  int lV_=2*Nin+1-lV;
  char *c;
  numout*ccx3;
  
  double ww,wz,wBrE,wBrM,zBrEn,zBrMn;
   
  
  if(pdg[lV]!=23 &&  abs(pdg[lV])!=24) return NULL;
  
//if(Nin==2)printf("%s %s -> %s %s\n", name[0], name[1],name[2],name[3]);  

  for(i=0;i<nModelParticles &&!(e_&&ne_&&m_&&nm_&&e_&&Ne_&&M_&&Nm_&&W_&&Z_ )   ;i++) 
  switch(ModelPrtcls[i].NPDG)
  {  
       case  11: e_ =ModelPrtcls[i].name;  E_=ModelPrtcls[i].aname; break;
       case -11: e_ =ModelPrtcls[i].aname; E_=ModelPrtcls[i].name;  break;
       case  12: ne_=ModelPrtcls[i].name; Ne_=ModelPrtcls[i].aname; break;
       case -12: ne_=ModelPrtcls[i].aname;Ne_=ModelPrtcls[i].name;  break;
       case  13: m_ =ModelPrtcls[i].name;  M_=ModelPrtcls[i].aname; break;
       case -13: m_ =ModelPrtcls[i].aname; M_=ModelPrtcls[i].name;  break;
       case  14: nm_=ModelPrtcls[i].name; Nm_=ModelPrtcls[i].aname; break;
       case -14: nm_=ModelPrtcls[i].aname;Nm_=ModelPrtcls[i].name;  break;
       case  24: W_=ModelPrtcls[i].name;                            break;
       case -24:                           W_=ModelPrtcls[i].aname; break;
       case  23: Z_=ModelPrtcls[i].name;                            break;       
  }
    
  if(!(e_&&ne_&&m_&&nm_&&e_&&Ne_&&M_&&Nm_&&W_&&Z_)) return NULL;
   
  {  txtList wDlist,zDlist;
     double ww_exp=2.085;
     double wz_exp=2.4952; 
     char txt[20];
     ww=pWidth(W_,&wDlist);  
     wz=pWidth(Z_,&zDlist);
     sprintf(txt,"%s,%s",E_,ne_);    wBrE=findBr(wDlist,txt);
     sprintf(txt,"%s,%s",M_,nm_);    wBrM=findBr(wDlist,txt);
     sprintf(txt,"%s,%s",ne_,Ne_);  zBrEn=findBr(zDlist,txt);         
     sprintf(txt,"%s,%s",nm_,Nm_);  zBrMn=findBr(zDlist,txt);

//     wBrE *=ww/ww_exp;  wBrM *=ww/ww_exp;  ww=ww_exp;
//     zBrEn*=wz/wz_exp;  zBrMn*=wz/wz_exp;  wz=wz_exp;
  } 
       
  sprintf(processX3,"%s",name[0]); 
  if(Nin>1) sprintf(processX3+strlen(processX3),",%s",name[1]);
  sprintf(processX3+strlen(processX3),"->%s,",name[lV_]);
  c=processX3+strlen(processX3);
             
  if(abs(pdg[lV_])==11 || abs(pdg[lV_])==12) switch(pdg[lV])
  {   case -24: sprintf(c,"%s,%s",m_,Nm_); *wV=ww; *br=wBrM; break;
      case  24: sprintf(c,"%s,%s",M_,nm_); *wV=ww; *br=wBrM; break;
      case  23: sprintf(c,"%s,%s",nm_,Nm_);  *wV=wz; *br=zBrMn; break;
  } else 
  {  
      switch(pdg[lV])
      { case -24: sprintf(c,"%s,%s",e_,Ne_); *wV=ww; *br=wBrE; break;
        case  24: sprintf(c,"%s,%s",E_,ne_); *wV=ww; *br=wBrE; break;
        case  23: sprintf(c,"%s,%s",ne_,Ne_);  *wV=wz; *br=zBrEn; break; 
       }
  }
  process2Lib(processX3,plib13); 
  if(Nin==2) sprintf(exclude,"%s","%Z+W<1"); else
  { if(pdg[lV]==23) sprintf(exclude,"%s<1",Z_); else sprintf(exclude,"%s<1",W_);}
  strcat(plib13,"V");    
  ccx3=getMEcode(0,ForceUG,processX3,exclude,"",plib13);
  if(ccx3)passParameters(ccx3);
  return ccx3;
}


static double decay22List(char * pname, txtList *LL)
{ int i,j,ntot,no22;
  numout * cc;
  double wtot,w;
  char pname2[20],process[20],plib[120];
  char * dname[8];
  txtList L=NULL,L_;
  char buff[100];
  int pN=pNum(pname);

  if(pN==0) return 0;
  for(i=0,j=0;pname[i];i++)
  if(pname[i]!=' ') pname2[j++]=pname[i];
  pname2[j]=0;
  strcpy(plib,"2width_");
  pname2lib(pname2,plib+7);
  sprintf(process,"%s->2*x",pname2);
  cc=getMEcode(0,ForceUG,process,NULL,"",plib);
  if(!cc) { if(LL) *LL=NULL; return -1;} 
  procInfo1(cc,&ntot,NULL,NULL);
  for(wtot=0,i=1;i<=ntot;i++)  if(chOpen(cc,i))
  {     
    w=pWidth2(cc,i);
    if(w!=0)
    {  
      procInfo2(cc,i,dname,NULL);    
      if(LL)
       { L_=malloc(sizeof(txtListStr));
         L_->next=L;
         L=L_;
         sprintf(buff,"%E  %s -> %s,%s",w,pname2,dname[1],dname[2]);
         L_->txt=malloc(20+strlen(buff));
         strcpy(L_->txt,buff);
       } 
       wtot+=w; 
    }
  }

  no22= L?0:1;
  
  if(pN!=23  && abs(pN)!=24 )
  { int k,l;
    REAL m[5];  
    int pdg[5]; 
    char*name[5];
       
    for(k=1;k<=ntot;k++) if(!chOpen(cc,k))
    { double w1=0,w2=0;
      int vd[3]={0,0,0};
      for(i=0;i<3;i++) name[i]=cc->interface->pinf(k,i+1,m+i,pdg+i);
      if(no22 && m[0]<= m[1]+m[2]) continue;
      
      if(pdg[0]==23 ||  abs(pdg[0])==24) continue;
            
      for(i=1;i<3;i++) vd[i]= (abs(pdg[i])==24 && VWdecay) || (pdg[i]==23 && VZdecay); 
      
      l=0;
      if((vd[1]||vd[2]) && m[0]+VVmassGap > m[1]+m[2])
      for(l=1;l<3;l++) if(vd[l]) break;
      if(l>0 && l<3)      
      {  int nW,iW;
         numout * cc13;
         int err;
         double wV,brV; 
         int l_=3-l; 
         if( vd[l_] && m[l]<m[l_]) {l=l_; l_=3-l;}
         cc13=xVtoxll(1,2,name,pdg, l, &wV, &brV);
         *(cc13->interface->BWrange)=20;
         if(cc13)         
         { double Mmax,C;
           *(cc13->interface->BWrange)=20;              
           passParameters(cc13);
           for(i3_=1;i3_<4;i3_++) if(strcmp(cc13->interface->pinf(1,i3_+1,NULL,NULL),name[l_])==0)break;
           for(i=0;i<4;i++) cc13->interface->pinf(1,1+i,Pmass+i,NULL);
           sqme=cc13->interface->sqme;
           nsub_stat=1;         
           {  double Mmin,Mmax;
              Mmax=Pmass[0]-Pmass[i3_];
              for(Mmin=0,j=1;j<4;j++) if(j!=i3_) Mmin+=Pmass[j]; 
              w=simpson(dWidthdM, Mmin*1.00001 , Mmax*0.9999,1.E-3);
                  
//             if(findVal(ModelPrtcls[abs(pTabPos(name[l]))-1].width,&wVt)) K=1;else K=wVt/wV; 
              w/=brV;
           }
           
           if( vd[l_] )
           {  double w2= pWidth(name[l_], NULL);  
               w*=decayPcmW(m[0],m[l],m[l_],wV,w2,0)/decayPcmW(m[0],m[l],m[l_],wV,0,0);
               if(pdg[l_]==pdg[l])  w/=2;      
           }
           
           if(w!=0)
           {  if(LL)
              { L_=malloc(sizeof(txtListStr));
                L_->next=L;
                L=L_;
                sprintf(buff,"%E  %s -> %s,%s",w,name[0],name[1],name[2]);
                L_->txt=malloc(20+strlen(buff));
                strcpy(L_->txt,buff);
              }           
              wtot+=w; 
           }
         }    
      }
    }
  }

  if(LL)
  {  for(L_=L;L_;L_=L_->next)
     { 
       sscanf(L_->txt,"%lf %[^\n]",&w,buff);
       sprintf(L_->txt,"%E %s",w/wtot,buff);
     }   
    *LL=L;
  }
  return  wtot;
}

static txtList conBrList(txtList BrList)
{ txtList out=NULL;
  char buff[100];
  double br;
  int inCode[10], outCode[10],i;
  for(;BrList;BrList=BrList->next)
  { txtList new=malloc(sizeof(txtListStr));
    new->next=out;out=new;
    sscanf(BrList->txt,"%lf %[^\n]",&br,buff); 
    decodeProcess(buff,inCode,outCode);
    if(inCode[0]>0) sprintf(buff,"%E  %s -> ",br,ModelPrtcls[inCode[0]-1].aname);   
    else            sprintf(buff,"%E  %s -> ",br,ModelPrtcls[-inCode[0]-1].name);
    for(i=0;outCode[i];i++)
    { if(i) strcat(buff,",");
      if(outCode[i]>0) strcat(buff,ModelPrtcls[outCode[i]-1].aname);   
      else             strcat(buff,ModelPrtcls[-outCode[i]-1].name);  
    }
    new->txt=malloc(strlen(buff)+1);
    strcpy(new->txt,buff);
  }
  return out;  
}


void setQforParticle(REAL *Q,char*pname)
{
  char *nm;
  REAL*ma;  
  int n,i,cdim;
  int pdg;
  if(!Q) return;
 
  n=pTabPos(pname);
  if(!n){printf("Wrong particle name '%s'\n",pname); return ;}
  nm=ModelPrtcls[abs(n)-1].mass;
  if(nm[0]=='0') return ; else ma=varAddress(nm);

  cdim=abs(ModelPrtcls[abs(n)-1].cdim);
  pdg=abs(ModelPrtcls[abs(n)-1].NPDG);
  
  if(cdim==1)
  { int err=calcMainFunc();  *Q=fabs(*ma); err=calcMainFunc(); 
    if(err) printf("Cannot calculate %s\n",err,varNames[err]);   
    return;
  }
  switch(pdg)
  { case 1:case 2:case 3: *Q=1; return;
    case 4: *Q=1.5; break;
    case 5: *Q=5;   break;
    case 6: *Q=175; break;
  }  
  calcMainFunc();
  for(i=0;i<10;i++) 
  { 
    if( fabs(*Q-fabs(*ma)) < 1E-2*(*Q)) break;
    *Q=fabs(*ma);
    calcMainFunc();
  }
} 


double pWidth(char *name, txtList * LL)
{
  txtList L,l,Lout;
  char libName[100];
  double sum=0,width;
  int i,i0,j,j0,nout;
  REAL Qstat;
  REAL*Q=NULL;
  for(i=0;i<nModelParticles;i++)
  { char *pnames[2]={ModelPrtcls[i].name,ModelPrtcls[i].aname}; 
    for(j=0;j<2;j++) if(strcmp(name,pnames[j])==0) 
    { if(decayTable[i].status==1)
      { 
        if(LL) *LL=decayTable[i].pdList[j];
        return decayTable[i].width;
      } else if(decayTable[i].status==-1)
      { if(LL) *LL=NULL;
        return 0;
      }  else break;
    } 
    if(j!=2) break;    
  }    

  i0=i,j0=j;
  if(i0==nModelParticles)
  { printf("%s out of model particles\n",name);
    if(LL) *LL=NULL;
    return 0;
  }  

  {  int pdg,pdg0,Len,decay[10];
     double br;
     pdg0=ModelPrtcls[i0].NPDG;
     if(j0) pdg0=-pdg0;
     for(i=1; allDecays(i,0,&pdg,&Len,decay,&width,&br) ;i++)
     {
        if(abs(pdg)==abs(pdg0))
        {  txtListStr*l,*L=NULL;
           l=malloc(sizeof(txtListStr));
           decayTable[i0].width=width;
           decayTable[i0].status=1;  
           for(j=1; allDecays(i,j,&pdg,&Len,decay,&width,&br) ;j++) if(br>0)
           { int k;  
             l=malloc(sizeof(txtListStr));
             l->txt=malloc(100); 
             l->next=L;
             sprintf(l->txt,"%E  %s -> ",br, pdg2name(pdg));
             sprintf(l->txt+strlen(l->txt),"%s",pdg2name(decay[0]));
             for(k=1;k<Len;k++) sprintf(l->txt+strlen(l->txt),", %s",pdg2name(decay[k]));
             L=l;
           }
           if(pdg0==pdg) 
           {  decayTable[i0].pdList[j0]=L; 
              if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname))
                           decayTable[i0].pdList[1-j0]=conBrList(L);
           } else 
           { decayTable[i0].pdList[1-j0]=L;
             if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname))
                           decayTable[i0].pdList[j0]=conBrList(L);
           }                
           if(LL) *LL=decayTable[i0].pdList[j0];
           return width;
        }
     }          
  }
  decayTable[i0].status=-1;  
  if(Q==NULL) for(i=0;i<nModelVars;i++) if(strcmp(varNames[i],"Q")==0){ Q= varValues+i; break;}
  if(Q) { Qstat=*Q; setQforParticle(Q,name);}
  width=decay22List(name,&L);
  if(L) 
  {
    if(LL) *LL=L;
    decayTable[i0].pdList[j0]=L;
    if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) 
                 decayTable[i0].pdList[1-j0]=conBrList(L); 
    decayTable[i0].width=width;
    decayTable[i0].status=1;
    if(Q) {*Q=Qstat; calcMainFunc();}    
    return width;
  }
  Lout=NULL;
  L= makeDecayList(name,3);
  massFilter(pMass(name),&L);
  gammaGluFilter(&L);
  if(L==NULL) 
  { L= makeDecayList(name,4);  
    massFilter(pMass(name),&L);
    gammaGluFilter(&L);
    nout=4;
  }  else nout=3;
      
  for(sum=0,l=L;l;l=l->next)  
  { numout* cc;
    int err=0;
    txtList newr;
    process2Lib(l->txt ,libName);
    cc=getMEcode(0,ForceUG,l->txt,NULL,"",libName);
    width=0;
    if(cc){if(nout==3) width=width13(cc, 1, &err); else width=width14(cc, &err);}
    if(width >0)
    {
      sum+=width;
      newr=malloc(sizeof(txtListStr));
      newr->next=Lout;
      Lout=newr;
      newr->txt=malloc(strlen(l->txt)+20);
      sprintf(newr->txt,"%E  %s",width,l->txt);
    }
  }
  cleanTxtList(L); 
  if(Lout)
  for(L=Lout;L;L=L->next)
  { char buff[100];
    sscanf(L->txt,"%lf %[^\n]",&width,buff);
    sprintf(L->txt,"%E %s",width/sum,buff);  
  }   
  if(LL) *LL=Lout;
  decayTable[i0].pdList[j0]=Lout;
  if(strcmp(ModelPrtcls[i0].name,ModelPrtcls[i0].aname)) 
               decayTable[i0].pdList[1-j0]=conBrList(Lout);
  decayTable[i0].width=sum;
  decayTable[i0].status=1;
  if(Q) { *Q=Qstat; calcMainFunc();}
  return sum;
}

double aWidth(char *name) { return pWidth(name,NULL);}

static int pListEq(char * txt1, char * txt2)  
{  char buff[100];
   char rd1[10][10];
   char rd2[10][10];
   int n1,n2,i1,i2;
   char *ch;
    
   strcpy(buff,txt1); while((ch=strchr(buff,','))) ch[0]=' ';
   
   n1=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd1[0],rd1[1],rd1[2],rd1[3],rd1[4],rd1[5],rd1[6],rd1[7],rd1[8],rd1[9]); 
   
   strcpy(buff,txt2); while((ch=strchr(buff,','))) ch[0]=' ';
   
   n2=sscanf(buff,"%s %s %s %s %s %s %s %s %s %s",
   rd2[0],rd2[1],rd2[2],rd2[3],rd2[4],rd2[5],rd2[6],rd2[7],rd2[8],rd2[9]); 
   
   if(n1!=n2) return 0;
   for(i1=0;i1<n1;i1++)
   { for(i2=0;i2<n2;i2++) if(strcmp(rd1[i1],rd2[i2])==0){rd2[i2][0]=0; break;}
     if(i2==n2) return 0;
   } 
   return 1;
}      

double findBr(txtList L, char * pattern)
{ char buff[100];
  char *ch;
  double width;
  
  for(;L;L=L->next)
  { 
     sscanf(L->txt,"%lf %[^\n]",&width,buff);
     ch=strstr(buff,"->");
     ch+=2;
     if( pListEq(ch,pattern)) return width;
  }
  return 0;   
}

/* =============  ProcInfo ================ */

int  procInfo1(numout*cc, int *nsub, int * nin, int *nout)
{
  if(nin) *nin=cc->interface->nin;
  if(nout)*nout=cc->interface->nout;
  if(nsub)*nsub=cc->interface->nprc;
  return 0;
}

int procInfo2(numout*cc,int nsub,char**name,REAL*mass)
{
  int i;
  int ntot=cc->interface->nin+cc->interface->nout;
    
  if(nsub<1 || nsub> cc->interface->nprc) return 2;

  if(name)for(i=0;i<ntot ;i++) 
  name[i]=(cc->interface->pinf)(nsub,i+1,NULL,NULL);

  if(mass)
  {  
    if(passParameters(cc)) return 4;
    for(i=0;i<ntot ;i++) cc->interface->pinf(nsub,i+1,mass+i,NULL);     
  }
  return 0;
}

/* ======================decayTable ===================*/
 static int nPrtcls_old=0; 
 void cleanDecayTable(void)
 { int i,j;
//printf("cleanDecayTable\n"); 
   if(decayTable) for(i=0;i<nPrtcls_old;i++) for(j=0;j<2;j++) if(decayTable[i].pdList[j]) 
      cleanTxtList(decayTable[i].pdList[j]);    
   decayTable=realloc(decayTable, nModelParticles*sizeof(decayTableStr));
   nPrtcls_old=nModelParticles;
   for(i=0;i<nModelParticles;i++)
   { for(j=0;j<2;j++) decayTable[i].pdList[j]=NULL;
     decayTable[i].width=0;
     decayTable[i].status=0;
   }
 }

/*============= Export of parameters ============*/
int passParameters(numout*cc)
{
   int i,k;
   for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
   k=cc->interface->calcFunc();
   if(k>0) { printf("cannot  calculate constr %s\n", cc->interface->varName[k] ); return 1;}
   return 0;
}

int slhaDecayPrint(char * name,FILE*f)
{
   double w;
   txtList all;
   int i,dim; 
   long PDG;
           
   PDG=qNumbers(name,NULL,NULL,NULL);
   if(!PDG) return 0;
   w=pWidth(name,&all);
   fprintf(f,"DECAY %d  %E  # %s\n",PDG,w,name);
   for(;all;all=all->next)
   {  
      char pn[20], buff[100], *chB,*chE;
      strcpy(buff,all->txt);
      sscanf(buff,"%s", pn);
      chB=strstr(buff,"->");
      chB+=2;
      for(dim=0,chE=chB ; chE;dim++, chE=strchr(chE+1,',')) continue;
      fprintf(f," %s   %d  ",pn,dim);

      for(i=0;i<dim;i++)
      { 
         chE=strchr(chB,',');
         if(chE)chE[0]=0;
         sscanf(chB,"%s",pn);
         fprintf(f," %d", qNumbers(pn,NULL,NULL,NULL));
         if(chE)chB=chE+1;else break;           
      }
      chB=strstr(all->txt,"->");
      fprintf(f,"  # %s \n",chB+2);
   } 
   fprintf(f,"\n");
   return PDG;
} 
