#include <iostream>
#include "ch_class.h"

extern "C" { 

#include "num_in.h"
#include "num_out.h"
#include "VandP.h"
#include "dynamic_cs.h"
#include "vp.h"
#include "rootDir.h"

}


using namespace std;

ClassImp(ch_class)

       ch_class::ch_class(void) { L=NULL;cc=NULL;Nin,Nout,Nproc=0;}           
int    ch_class::SetModel(char * modelDisp , int nModel) { return setModel(modelDisp ,nModel); } 
int    ch_class::AssignValW(char * name, double val)     { return assignValW(name,val);}
double ch_class::FindValW(char * name)                   { return findValW(name);}
int    ch_class::CalcMainFunc(void)                      { int err=calcMainFunc(); 
                                                           if(err) return err;
                                                           if(cc)  passParameters(cc);
                                                           return err;  
                                                         }
double ch_class::PMass(char * pName)                     { return pMass(pName);}
double ch_class::PWidth(char * pName)                    { return pWidth(pName,&L);}
void   ch_class::PrintTxtList(FILE*f)                    { printTxtList(L,f);}
int    ch_class::PNum(char*pName)                        { return pNum(pName); }
char*  ch_class::Pdg2name(int pdg)                       { return pdg2name(pdg);}
int    ch_class::QNumbers(char*pname,int*spin2,int*charge3,int*cdim) 
                                                         { return qNumbers(pname, spin2, charge3, cdim);}
double ch_class::FindBr(char*pattern)     { if(L) return findBr(L,pattern); else printf("empty decay list\n"); }
int    ch_class::SlhaDecayPrint(char*pname,FILE*f)       { return slhaDecayPrint(pname,f);}
int    ch_class::GetMEcode(int twidth,int UG,char*Process,char * excludeVirtual,char*excludeOut,char*lib)
{ cc=getMEcode(twidth, UG,Process, excludeVirtual,excludeOut,lib);
  if(cc) 
  { Nin=cc->interface->nin;
    Nout=cc->interface->nout;
    Nproc=cc->interface->nprc;
    return 0;
  }  
  return 1;
} 

int ch_class::NewProcess(char*Process)                   { cc=newProcess(Process);
                                                           if(cc)
                                                           {  Nin=cc->interface->nin;
                                                              Nout=cc->interface->nout;
                                                              Nproc=cc->interface->nprc;
                                                              return 0;
                                                            } else return 1;
                                                         }
                                                          
int    ch_class::PassParameters(void)              {  if(cc) return passParameters(cc); else printf("process is not compiled\n"); }
char*  ch_class::Pinf(int nsub,int nprtcl,double*pmass,int*num)
                                                         { if(cc) return cc->interface->pinf(nsub,nprtcl,pmass,num);
                                                           else printf("process is not compiled\n"); }

int ch_class::ProcInfo2(int nsub, char**pnames,double * pmass)
                                                         { if(cc) return procInfo2(cc, nsub, pnames, pmass); else
                                                           {  printf("process is not compiled\n"); return -1;}
                                                         }
                                                         
double ch_class::Sqme(int nsub,double GG, double *pvect,int*err)
                                                         { if(cc) return cc->interface->sqme(nsub,GG,pvect,err);
                                                           else printf("process is not compiled\n");
                                                         } 