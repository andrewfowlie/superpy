#include"micromegas.h"
#include"micromegas_aux.h"
#include"micromegas_f.h"

#include <sys/types.h>
#include <unistd.h>
              
int readvar_(char *fname, int len)
{ int err;
  char * cname=malloc(len+1);
  fName2c(fname,cname,len);
  err=readVar(cname);
  free(cname);
  return err;
}



void printvar_(int * Nch) 
{ char fname[20];
  FILE*f;
  
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  printVar(f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname); 
}
                                                                                

void printmasses_(int * Nu, int* sort)
{ 
  char fname[20];
  FILE*f;

  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  printMasses(f,*sort);
  fclose(f);
  fortreread_(Nu,fname,strlen(fname));
  unlink(fname);
}

void printhiggs_(int * Nch)
{
  char fname[20];
  FILE*f;
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  printHiggs(f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname);
}


int sortoddparticles_(char * f_name, int len) 
{ 
  char c_name[20];
  int err=sortOddParticles(c_name);
  cName2f(c_name, f_name,len);
  return err;  
}

void nextodd_(int *n, char * pName, double *pMass,int len)
{ 
  char * pm;
  pm=nextOdd(*n,pMass);
  if(pm==NULL) pm=" ";
  cName2f(pm,pName,len);         
}

void pdg2name_(int * pdg, char * name,int len)
{ 
  char *pn;
  pn=pdg2name(*pdg);
  if(pn==NULL) pn=" ";
  cName2f(pn,name,len);
}  
   
double pmass_(char * pName, int len)
{ char c_name[20];   
  fName2c(pName,c_name,len);
  return pMass(c_name);
}

int qnumbers_(char*pname, int *spin2,int*charge3,int*cdim,int len)
{  char cName[20];
   int pdg; 
   
   fName2c(pname,cName,len); 
   pdg=qNumbers(cName, spin2, charge3, cdim);
   return pdg;
}
                                                                                                   
double darkomega_(double * Xf,int*Fast,double *Beps){return darkOmega(Xf,*Fast,*Beps);}
double darkomegafo_(double*Xf,int*fast,double*Beps){return darkOmegaFO(Xf,*fast,*Beps);}


double printchannels_(double*Xf,double*cut,double*Beps,int* prcnt,int *Nu)
{
  char fname[20];
  FILE*f;
  double res;
                                                                                   
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  res=printChannels(*Xf,*cut,*Beps,* prcnt,f);
  fclose(f);
  if(*Nu) fortreread_(Nu,fname,strlen(fname));
  unlink(fname);
 
  return res;
}
double onechannel_(double *Xf,double *Beps,char*fn1,char*fn2,char*fn3,char*fn4,
int len1,int len2,int len3,int len4)
{ 
  char n1[20], n2[20], n3[20], n4[20];
  fName2c(fn1,n1,len1);
  fName2c(fn2,n2,len2);
  fName2c(fn3,n3,len3);
  fName2c(fn4,n4,len4); 
  return  oneChannel(*Xf,*Beps,n1,n2,n3,n4);
}
  


double decay2info_(char * pname, int *Nch, int len)
{ double res;
  char cname[20]; 
  char fname[20];
  FILE*f;
  
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  fName2c(pname,cname,len);    
  res=decay2Info(cname,f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname); 
  return res;
}

int slhadecayprint_(char * pname,int *Nch,int len)
{ double res;
  char cname[20]; 
  char fname[20];
  FILE*f;
  
  sprintf(fname,"%d.tmptxt",getpid());
  f=fopen(fname,"w");
  fName2c(pname,cname,len);    
  res=slhaDecayPrint(cname,f);
  fclose(f);
  fortreread_(Nch,fname,strlen(fname));
  unlink(fname); 
  return res;
}

double vsigma_(double*T,double*Beps,int*Fast){ return vSigma(*T,*Beps,*Fast);}

static int channels(aChannel* Chann, int *i, double *w, int*pdg, char* txt, int len)
{ int k,j;
  if(!Chann || *i<1)  return 0;
  for(k=0;k<*i;k++) if(Chann[k].weight==0) return 0;
  k--;
  *w=Chann[k].weight;
  
  txt[0]=0; 
  for(j=0;j<2;j++) sprintf(txt+strlen(txt),"%s ",Chann[k].prtcl[j]);
  strcpy(txt+strlen(txt),"-> ");
  for(j=2;j<4;j++) sprintf(txt+strlen(txt),"%s ",Chann[k].prtcl[j]);
  for(j=0;j<4;j++) pdg[j]=pNum(Chann[k].prtcl[j]);        
  if(Chann[k].prtcl[4]) 
  {  sprintf(txt+strlen(txt),"%s ",Chann[k].prtcl[j]);
      pdg[j]=pNum(Chann[k].prtcl[j]);
  } else pdg[j]=0;    
  for(j=strlen(txt);j<len;j++) txt[j]=' '; 
  return 1;
}
 
int omegach_(int*i,double*w,int*pdg,char*txt,int len)   {return channels(omegaCh,  i,w,pdg,txt,len);}

int vsigmatch_(int *i,double*w,int*pdg,char*txt,int len){return channels(vSigmaTCh,i,w,pdg,txt,len);}

int vsigmach_(int *i,double*w,int*pdg,char*txt,int len) {return channels(vSigmaCh, i,w,pdg,txt,len);}


static double(*_fDv)(double*);
static double fDv_(double v){ return (*_fDv)(&v);}


int neutrinoflux_( double(*fDv)(double*), int* forSun, double* nu, double * Nu)
{  
  if(fDv  == maxwell_)  return neutrinoFlux(Maxwell, *forSun, nu, Nu);
  else 
  {  _fDv=fDv;
    return   neutrinoFlux(fDv_, *forSun, nu, Nu);
  }    
}

void muonupward_(double*nu, double*Nu, double*mu) { muonUpward(nu, Nu, mu);}
void muoncontained_(double*nu,double*Nu,double *rho, double*mu) { muonContained(nu,Nu,*rho, mu);}

void cleandecaytable_(void) { cleanDecayTable(); } 
void setvvdecay_(int*vwdecay,int*vzdecay ){ VWdecay=*vwdecay;  VZdecay=*vzdecay;  cleanDecayTable(); }

double captureaux_(int*forSun,double(*fDv)(double*), double*csIp,double*csIn,double*csDp, double*csDn)
{
  _fDv=fDv; 
  captureAux(*forSun,fDv_, *csIp,*csIn,*csDp,*csDn);
}

void setbwrange_(double * newValue) { BWrange=*newValue;}
