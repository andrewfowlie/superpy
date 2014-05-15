#include"SLHAplus.h"

#include"aList.h"
/*
       'bridge'  routines  written in C and simulating Fortran code.
       See main rules of bridge constructions in   Chapter 11 of
http://www.nd.edu/~hpcc/solaris_opt/SUNWspro.forte6u1/WS6U1/lib/locale/C/html/manuals/fortran/prog_guide
*/


/*===============  reading SLHA  ==============*/ 

struct 
{ int ferror;} ferror_;

static void fName2c(char*f_name,char*c_name,int len)
{ int i; for(i=len-1;i>=0 &&f_name[i]==' ';i--);
  c_name[i+1]=0;
  for(;i>=0;i--) c_name[i]=f_name[i];
}

static void cName2f(char*c_name,char*f_name,int len)
{ int i; for(i=0;i<len &&c_name[i];i++) f_name[i]=c_name[i];
         for(   ;i<len            ;i++) f_name[i]=' ';
}
         

struct slhacomment_ { char txt[100]; }  slhacomment_;

static int    _N_;

extern int fortranreadline_(int *, char*,int);

static int readLnC(int size, char*buff)
{   int i;
    if(fortranreadline_(&_N_, buff,size-1))  return -1;
    i=size-1;
    buff[i--]=0;
    for(;i>=0 && buff[i]==' ';i--) buff[i]=0;                                                          
    return 0;
}

int  slhareadstream_(int *N, int *mode, char * end, int len )
{
   char end_[100];
   int err,anydate=0;
   fName2c(end,end_,len);
   _N_=*N;
   err=slhaBasicReader(*mode,readLnC,&anydate,end_);
   if((err==0 || err==-1) && anydate==0) {FError=ferror_.ferror=1; return -3;}
   if(err==-1) err=0; 
   if(err) FError=ferror_.ferror=1;
   return err;
}

int slharead_(char * fname, int * mode, int len)
{ char  c_name[200];
  int err;
  fName2c(fname,c_name,len);
  ferror_.ferror=0;
  err=slhaRead(c_name,*mode);
  ferror_.ferror=FError;
  return err;
}

int slhawrite_(char *fname, int len)
{ char  c_name[200];
  fName2c(fname,c_name,len);
  return slhaWrite(c_name);
}

int slhadecayexists_(int *pNum){ return slhaDecayExists(*pNum);}

double slhawidth_(int *pNum)
{double res=slhaWidth(*pNum); ferror_.ferror=FError; return res;}

double slhabranch_(int*pNum,int*N,int*nCh)
{double  res=slhaBranch(*pNum,*N,nCh); ferror_.ferror=FError; return res;}

double slhaval0_(char * Block, double *Q, int len)
{ char  c_name[200];
  double res;
  fName2c(Block,c_name,len);
  res=slhaVal(c_name, *Q, 0);
  ferror_.ferror=FError;
  return res;
} 

double slhaval1_(char * Block, double *Q, int *k1, int len)
{  char  c_name[200];
   double res;
   fName2c(Block,c_name,len);
   res= slhaVal(c_name, *Q, 1,*k1);
   ferror_.ferror=FError;
   return res;
} 

double slhaval2_(char * Block, double *Q, int *k1,int*k2, int len)
{  char  c_name[200];
   double res;
   fName2c(Block,c_name,len);
   res=slhaVal(c_name, *Q, 2,*k1,*k2); 
   ferror_.ferror=FError;
   return res;       
} 

double slhaval3_(char * Block, double *Q, int *k1,int*k2,int*k3, int len)
{  char  c_name[200];
   double res;
   fName2c(Block,c_name,len);
   res=slhaVal(c_name, *Q, 3,*k1,*k2,*k3);
   ferror_.ferror=FError;
   return res;       
} 

double complex cslhaval0_(char * Block, double *Q, int len)
{ char  c_name[200];
  double complex res;
  fName2c(Block,c_name,len);
  res=cslhaVal(c_name, *Q, 0);
  ferror_.ferror=FError;
  return res;
} 


double complex  cslhaval1_(char * Block, double *Q, int *k1, int len)
{  char  c_name[200];
   double complex res;
   fName2c(Block,c_name,len);
   res= cslhaVal(c_name, *Q, 1,*k1);
   ferror_.ferror=FError;
   return res;
} 

double complex cslhaval2_(char * Block, double *Q, int *k1,int*k2, int len)
{  char  c_name[200];
   double complex res;
   fName2c(Block,c_name,len);
   res=cslhaVal(c_name, *Q, 2,*k1,*k2); 
   ferror_.ferror=FError;
   return res;       
} 

double complex cslhaval3_(char * Block, double *Q, int *k1,int*k2,int*k3, int len)
{  char  c_name[200];
   double complex res;
   fName2c(Block,c_name,len);
   res=cslhaVal(c_name, *Q, 3,*k1,*k2,*k3);
   ferror_.ferror=FError;
   return res;       
} 

 
int slhavalexists0_(char * Block, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(c_name, 0);
}

int slhavalexists1_(char * Block, int*k1, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(c_name, 1,*k1);
}

int slhavalexists2_(char * Block, int*k1, int*k2, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(c_name, 2,*k1,*k2);
}

int slhavalexists3_(char * Block, int*k1, int*k2, int*k3, int len)
{  char  c_name[200];
   fName2c(Block,c_name,len);
   return slhaValExists(c_name, 3,*k1,*k2,*k3);
}
/*
int slhawarnings_(int *Nch)
{
  char fname[20];
  FILE*f;
  int err;
  if(*Nch)
  {
     sprintf(fname,"%d.tmptxt",getpid());
     f=fopen(fname,"w");
     err=slhaWarnings(f);
     fclose(f);
     fortreread_(Nch,fname,strlen(fname));
     unlink(fname);   
     return err;
  } else return err=slhaWarnings(NULL);   
}
*/
/*===============  qNumbers ============================= */

int findqnumbers_(int *pdg,int*eQ3,int*spinDim,int*cDim,int*neutral)
{
   return  findQnumbers(*pdg,eQ3,spinDim,cDim,neutral);
}

int allqnumbers_(int *i, int*pdg,int*eQ3,int*spinDim,int*cDim,int*neutral)
{ int k;
  int ret= allQnumbers(*i, pdg,eQ3,spinDim,cDim,neutral);
  strncpy(slhacomment_.txt, slhaComment,100);
  for(k=strlen(slhaComment);k<100;k++) slhacomment_.txt[k]=' ';

  return ret;
}

int allblocks_(int *i,int*j,char*name,int*Len,int*key, double complex * val,int l)       
{ char c_name[40];
  int k, res= allBlocks(*i,*j,c_name,Len,key,val);

  strncpy(slhacomment_.txt,slhaComment,100);
  for(k=strlen(slhaComment);k<100;k++) slhacomment_.txt[k]=' ';
  
  for(k=0;c_name[k]&&k<l;k++) name[k]=c_name[k];
  for(;k<l;k++) name[k]=' ';
  return res;
}
 
int alldecays_(int*i,int*j,int*pdg,int*Len,int*decay,double*width,double*br)
{  int k,ret;
   ret=allDecays(*i,*j,pdg,Len,decay,width,br);
   strncpy(slhacomment_.txt,slhaComment,100);
   for(k=strlen(slhaComment);k<100;k++) slhacomment_.txt[k]=' ';
   return ret;    
}

/*================  writing SLHA ======================== */

int openappend_(char *fileName,int len)
{  int res;
   char  c_name[200];
   fName2c(fileName,c_name,len);
   res= openAppend(c_name);
   ferror_.ferror=FError;
   return res;
} 


int system1_(char *command, int len)
{
   char  c_name[200];
   int res;
   fName2c(command,c_name,len);
   res=System(c_name);
   ferror_.ferror=FError;
   return res;        
}

int system2_(char *format, char*path, int len1,int len2)
{
   char  c_name1[200],c_name2[200];
   int res;
   fName2c(format,c_name1,len1);
   fName2c(path,c_name2,len2);
   res=System(c_name1,c_name2);
   ferror_.ferror=FError;
   return res;        
}

void setsystimelim_(int *lim, int *quant)
{
    sysTimeLim=*lim;
    sysTimeQuant=*quant;  
}
int aprintf0_(char * format, int len)
{
 char  c_name[1000];
 fName2c(format,c_name,len); 
 return aPrintF(c_name);
}

int aprintf1_(char * format, double * x1,int len)
{
  char  c_name[1000];
  fName2c(format,c_name,len); 
  return aPrintF(c_name,*x1);
}

int aprintf2_(char * format, double * x1, double *x2,int len)
{
  char  c_name[1000];
  fName2c(format,c_name,len); 
  return aPrintF(c_name,*x1,*x2);
}

int aprintf3_(char * format, double*x1,double*x2,double*x3, int len)
{
  char  c_name[1000];
  fName2c(format,c_name,len); 
  return aPrintF(c_name,*x1,*x2,*x3);
}

int aprintf4_(char * format, double*x1,double*x2,double*x3,double*x4,int len)
{
  char  c_name[1000];
  fName2c(format,c_name,len); 
  return aPrintF(c_name,*x1,*x2,*x3,*x4);
}

int aprintf5_(char * format, double*x1,double*x2,double*x3,double*x4,double*x5,int len)
{
  char  c_name[1000];
  fName2c(format,c_name,len); 
  return aPrintF(c_name,*x1,*x2,*x3,*x4,*x5);
}

double slhavalformat_(char * Block, double *Q, char * format,int len1,int len2)
{  char  cBlock[1000], cFormat[1000];
   double res; 
   fName2c(Block, cBlock, len1);
   fName2c(format,cFormat,len2); 
   
   res=slhaValFormat(cBlock, *Q, cFormat);

   return res;
}   


int slhastrformat_(char * Block, char * format, char * fRes, int len1,int len2,int len3)
{  char  cBlock[1000], cFormat[1000], cRes[1000];
   int res; 
   fName2c(Block, cBlock, len1);
   fName2c(format,cFormat,len2); 
   
   res=slhaSTRFormat(cBlock, cFormat,cRes);
   if(res!=0) return res;
   cName2f(cRes,fRes,len3);
   return 0;
}   



/*================ Diagonalizing ==========================*/

int rjacobi_(double*a,int*n,double*d,double*v){ return rJacobi(a,*n,d,v); }

int rjacobia_(double*a,int*n,double*d,double*u,double*v)
     {return rJacobiA(a,*n,d,u,v);}
int cjacobih_(double complex*a,int*n,double*d,double complex*v){return cJacobiH(a,*n,d,v);}
int cjacobis_(double complex*a,int*n,double*d,double complex*v){return cJacobiS(a,*n,d,v);}
int cjacobia_(double complex*a,int*n,double*d,double complex*u,double complex*v)
     { return  cJacobiA(a,*n,d,u,v);}

int initdiagonal_(void) { return initDiagonal();}

int rdiagonal2_(aList3(double*))  { return rDiagonal(2,aList3(*));}
int rdiagonal3_(aList6(double*))  { return rDiagonal(3,aList6(*));}
int rdiagonal4_(aList10(double*)) { return rDiagonal(4,aList10(*));}
int rdiagonal5_(aList15(double*)) { return rDiagonal(5,aList15(*));}

int cdiagonalh2_(aList3(double complex*))  { return cDiagonalH(2,aList3(*));}
int cdiagonalh3_(aList6(double complex*))  { return cDiagonalH(3,aList6(*));}
int cdiagonalh4_(aList10(double complex*)) { return cDiagonalH(4,aList10(*));}
int cdiagonalh5_(aList15(double complex*)) { return cDiagonalH(5,aList15(*));}

int cdiagonals2_(aList3(double complex*))  { return cDiagonalS(2,aList3(*));}
int cdiagonals3_(aList6(double complex*))  { return cDiagonalS(3,aList6(*));}
int cdiagonals4_(aList10(double complex*)) { return cDiagonalS(4,aList10(*));}
int cdiagonals5_(aList15(double complex*)) { return cDiagonalS(5,aList15(*));}

int cdiagonala2_(aList4(double complex*))  { return cDiagonalA(2,aList4(*));}
int cdiagonala3_(aList9(double complex*))  { return cDiagonalA(3,aList9(*));}
int cdiagonala4_(aList16(double complex*)) { return cDiagonalA(4,aList16(*));}
int cdiagonala5_(aList25(double complex*)) { return cDiagonalA(5,aList25(*));}

int rdiagonala2_(aList4(double*))  { return rDiagonalA(2,aList4(*));}
int rdiagonala3_(aList9(double*))  { return rDiagonalA(3,aList9(*));}
int rdiagonala4_(aList16(double*)) { return rDiagonalA(4,aList16(*));}
int rdiagonala5_(aList25(double*)) { return rDiagonalA(5,aList25(*));}

double  massarray_(int * id, int *i){ return   MassArray(*id,  *i);}
double   mixmatrix_ (int*id, int*i,int*j){ return    MixMatrix (*id,*i,*j);}
double   mixmatrixu_(int*id, int*i,int*j){ return    MixMatrixU(*id,*i,*j);}
double remixmatrix_ (int*id, int*i,int*j){return creal(cMixMatrix (*id,*i,*j));}
double remixmatrixu_(int*id, int*i,int*j){return creal(cMixMatrixU(*id,*i,*j));}
double immixmatrix_ (int*id, int*i,int*j){return cimag(cMixMatrix (*id,*i,*j));}
double immixmatrixu_(int*id, int*i,int*j){return cimag(cMixMatrixU(*id,*i,*j));}

/*========================  QCD ===================== */

double  initqcd_(double * alfsMZ, double * McMc, double * MbMb, double * Mtp)
{ return  initQCD(*alfsMZ,*McMc,*MbMb,*Mtp); }
double alphaqcd_(double * q) { return alphaQCD(*q);}
double mcrun_(double *q){ return McRun(*q);}
double mbrun_(double *q){ return MbRun(*q);}
double mtrun_(double *q){ return MtRun(*q);}

double mceff_(double *q){ return McEff(*q);}
double mbeff_(double *q){ return MbEff(*q);}
double mteff_(double *q){ return MtEff(*q);}

#include"delList.h"
