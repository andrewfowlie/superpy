#include<math.h>
#include<stdio.h>
#include"micromegas.h"
#include"micromegas_aux.h"

#define makePdtList   makePdtList_mo
#define delPdtList    delPdtList_mo
#define getPdtData    getPdtData_mo
#define freePdtData   freePdtData_mo
#define interFunc     interFunc_mo
#define interAlpha    interAlpha_mo

#include"../CalcHEP_src/c_source/num/pdt.c"

static  pdtStr *data[12]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

static  pdtStr *data_;
static  double q_;


static double x_integrand(double x)
{  if(x==0) return 0; else return  x*interFunc(x, q_, data_);  }


static int pos(int pNum)
{ switch(pNum) 
  { case   5: case -5: return 1;
    case   4: case -4: return 2;
    case   3: case -3: return 3;
    case  -1:          return 4;
    case  -2:          return 5;
    case  21:          return 6;
    case   2:          return 7;
    case   1:          return 8;
    case  81:          return 9;
    case -81:          return 10;
    case  83:          return 11;
    case -83:          return 12;
    default:           return 0;
  }  
}


static pdtStr * readPdtData(int pc)
{ pdtStr * DT;
  static char*pdtFile=NULL;  
  if(pdtFile==NULL)
  {
     pdtFile=malloc(strlen(calchepDir)+30);
     sprintf(pdtFile,"%s/pdTables/cteq6l.pdt",calchepDir);  
  }  
  DT=malloc(sizeof(pdtStr));
  if(getPdtData(pdtFile, pos(pc), DT )==0) return DT; 
  else
  {    
     printf("error in reading PDT file\n");
     free(DT);
     return NULL;
  }
}

static pdtStr * checkdata(int pc)
{ int pp=pos(pc)-1;
  if(data[pp]) return data[pp]; 
  if(abs(pc) <80) data[pp]= readPdtData(pc);
  else
  { pdtStr * dataAux;
    double a,b;
    int i, ntot;
    data[pp]=readPdtData(3);
    ntot=data[pp]->nq*data[pp]->nx; 
    pc=(pc>0)? pc-80:pc+80;
    dataAux=checkdata(pc);
    if(abs(pc)==1) a=0.221*0.221; else a=1-0.221*0.221;
    b=1-a;
    for(i=0;i<ntot;i++) data[pp]->strfun[i]= a*data[pp]->strfun[i]
                                            +b*dataAux->strfun[i]; 
  }
  return data[pp];
}


double parton_x( int pNum, double  Q)
{
  double x1;
  q_=Q;

  data_=checkdata(pNum); 
  if(data_ == NULL || q_< data_->q_min ) return 0;
  
  x1=simpson(x_integrand,data_->x_min,1.,1.E-4);  
  if(pNum==21) return x1;
  if(abs(pNum)>2) return 2*x1;
  
  data_=checkdata(-pNum);
  if(data_ == NULL) return 0;  
  return x1+simpson(x_integrand,data_->x_min,1.,1.E-4);
}

double parton_alpha(double q){return  interAlpha(q,checkdata(21));}

static pdtStr *data1,*data2;
static double x0_;

static double conv_integrand(double y)
{ double x=exp(-y);
  return interFunc(x, q_, data1)*interFunc(x0_/x, q_, data2);
}

static  pdtStr **cData=NULL;

static  pdtStr *convStrFunAux(int pc1, int pc2)
{ pdtStr *D;
  int i,ix,iq,nc,i1=pos(pc1),i2=pos(pc2);
  if(i1>=i2) nc=i1*(i1-1)/2+i2-1; else nc=i2*(i2-1)/2+i1-1;
  if(cData && cData[nc])  return cData[nc];
  if(cData==NULL) 
  { cData=malloc(78*sizeof(pdtStr*));
    for(i=0;i<78;i++) if(cData[i]) cData[i]=NULL;
  }

  D=readPdtData(21);
   
  data1=checkdata(pc1);
  data2=checkdata(pc2);
  if(D->x_grid[0]==0)
  {  D->x_grid[0]=D->x_grid[1]/2;
     if(D->x_grid_aux) D->x_grid_aux[0]=pow(D->x_grid_aux[0],0.3);
  }   
  for(iq=0;iq<D->nq;iq++)for(ix=0;ix<D->nx;ix++)
  { 
    q_ =exp(D->q_grid[iq]);
    x0_=D->x_grid[ix];
    D->strfun[D->nx*iq+ix] = simpson(conv_integrand,0.,-log(x0_),1.E-3);
  }

  cData[nc]=D;
  return D;  
}

double convStrFun2(double x, double q, int pc1, int pc2, int pp)
{ pdtStr *D1,*D2;
  int pc1c,pc2c;
  double strF;
  
  if(pp>0 || pc2==21) pc2c=pc2; else pc2c=-pc2;
  D1=convStrFunAux(pc1, pc2c);
  strF=interFunc(x, q, D1);
  if(pc1!=pc2)
  { if(pp>0 || pc1==21) pc1c=pc1; else pc1c=-pc1;
    D2=convStrFunAux(pc1c, pc2);
    if(D1==D2) strF*=2; else strF+=interFunc(x, q, D2);
  }
  return strF;
}

double alpha_2(double Q) { return parton_alpha(Q);}
