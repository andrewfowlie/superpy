#ifndef __DYNAMIC_ME__
#define __DYNAMIC_ME__

#include<stdio.h>

#ifdef __cplusplus
extern "C" {
#endif 

#ifndef  __MICROMEGAS__
#include "../../../include/num_out.h"
#include "../../../include/nType.h"

typedef struct numout
{
  void * handle;
  REAL ** link;
  REAL *Q,*SC;
  int init;
  CalcHEP_interface * interface; 
} numout;

#endif

#ifndef  __MICROMEGAS__



/*======= Subprocesses ===========*/
  typedef struct txtListStr
  {  struct txtListStr * next;
     char  *txt;
  } txtListStr;

  typedef txtListStr * txtList;

#endif



extern char  * libDir;
extern char  * modelDir;
extern char  * compDir;
extern char  * calchepDir;
extern int   modelNum;

extern int prepareWorkPlace(void);
extern int  checkMtime(char * fname);
extern int  checkWorkPlace(void);
extern int getDynamicVP(void);
extern int setModel(char * modelDisp , int nModel);
extern int assignVal(char * name, double val);
extern int assignValW(char * name, double val);
extern int findVal(char * name, double*val);
extern double findValW(char * name);
extern int readVar(char *fname);

extern int    cleanWorkPlace(void);
extern numout*getMEcode(int twidth,int UG,char*Process,
                          char * excludeVirtual,char*excludeOut,char*lib);
extern numout*newProcess(char*Process);

extern txtList  makeDecayList(char * pname, int nx);
extern txtList  makeProcList(char ** InNames, char** OutNames, int nx);
extern void cleanTxtList(txtList L);
extern void printTxtList(txtList L, FILE *f);


/*===========================================================*/
typedef struct{ double width; txtList pdList[2]; int status;}  decayTableStr;

extern double   (*sqme)(int nsub,double GG,REAL *pvect, int * err_code);
extern double   decayPcm(double am0,  double  am1,  double  am2);
extern double   decayPcmW(double m0,double m1,double m2,double w1,double w2,int N);
extern int      ForceUG;
extern int      procInfo1(numout*cc, int *nsub, int * nin, int *nout);
extern int      procInfo2(numout*cc, int nsub,char**name,  REAL *mass);
extern void     massFilter(double M, txtList * List);
extern void     gammaGluFilter(txtList * List);
extern int      process2Lib(char * process,char * lib);
extern          decayTableStr* decayTable;
extern void     cleanDecayTable(void);
extern int      pname2lib(char*pname, char * libname);
extern double   decay2Info(char * pname, FILE* f);
extern double   pWidth(char *name, txtList * LL);
extern double   aWidth(char *name);
extern double   findBr(txtList L, char * pattern);
extern double   pWidth2(numout * cc, int nsub);

extern void     process2Mass(char * process,double * mass);
extern void     delAllLib(void);

extern int VZdecay;
extern int VWdecay;
 
extern numout* xVtoxll(int Nin,int Nout,char**name,int *pdg, int lV, double *wV,  double *br);
extern int passParameters(numout*cc);

extern int slhaDecayPrint(char * name,FILE*f);
extern void setQforParticle(REAL *Q,char*pname);

extern double BWrange;

#ifdef __cplusplus
}
#endif 


#endif
