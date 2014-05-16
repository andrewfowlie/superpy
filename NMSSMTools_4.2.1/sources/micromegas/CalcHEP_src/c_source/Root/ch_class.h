#ifndef __CalcHEP__
#define __CalcHEP__


extern "C" { 

//#include "num_in.h"
#include "num_out.h"
//#include "VandP.h"
#include "dynamic_cs.h"
//#include "rootDir.h"
}

#include "TObject.h"

class ch_class:public TObject  {

private:


 public:
    txtList L;
    numout*cc;
    int Nin,Nout,Nproc;
    ch_class(void);
    int    SetModel(char * modelDisp , int nModel);
    int    AssignValW(char * name, double val);
    double FindValW(char * name);
    int    CalcMainFunc(void); 
    double PMass(char * pName); 
    double PWidth(char * pName);   
    void   PrintTxtList(FILE*f);
    int    SlhaDecayPrint(char*pname,FILE*f);
    int    PNum(char*pName);
    char*  Pdg2name(int pdg);
    int    QNumbers(char *pname, int *spin2, int * charge3, int * cdim);
    double FindBr(char*pattern);
    int    GetMEcode(int twidth,int UG,char*Process,char * excludeVirtual,char*excludeOut,char*lib);
    int    NewProcess(char*Process);
    int    PassParameters(void);
    char * Pinf(int nsub,int nprtcl,double*pmass,int*num);
    int    ProcInfo2(int nsub, char**pnames,double * pmass);    
    double Sqme(int nsub, double GG, double *pvect,int*err);
    ClassDef(ch_class, 1)
 };

#endif
