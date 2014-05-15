
#include "../../include/extern.h"
#include "../../include/VandP.h"
#include "dynamic_cs.h" 
#include "interface.h"
#include "num1.h"

int nModelParticles=0;
static ModelPrtclsStr ModelPrtcls_[0]={ };

ModelPrtclsStr *ModelPrtcls=ModelPrtcls_; 
int nModelVars=0;
int nModelFunc=0;
static char*varNames_[0]={};

char**varNames=varNames_;
static REAL varValues_[0]={};
REAL*varValues=varValues_;
int calcMainFunc(void) {return 0;}


double aWidth(char*name) { return 0;}
double pWidth(char *name, txtList * LL){ return 0;}

REAL Helicity[2]={0,0};


char   p_names_[MAXNP][NAMELEN];
int    p_codes_[MAXNP]; 
double p_masses_[MAXNP];


char* pinf_ext(int nsub,int num , REAL * mass,int * N)
{  
  if(nsub>1 || num>nin_int+nout_int) return NULL;
  if(mass) *mass=p_masses_[num-1];
  if(N) *N=p_codes_[num-1];
  return  p_names_[num-1];
}

int VWdecay;
int VZdecay;
