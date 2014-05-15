/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "interface.h"
#include<math.h>
#include"cut.h"
#include"parser.h"
#include"rd_num.h"
#include"crt_util.h"
#include"phys_val.h"
#include"read_func.h"
#include"sortarr.h"

invcut_ invcut_1[60];
int nCuts=0;
table cutTab={"*** Table ***"," Cuts  ",
                   "!|  Parameter  |> Min bound <|> Max bound <|",NULL,0}; 


static void  cutInit(void)
{ int i;
  for(i=0;i<nCuts;i++) 
  { cleanPVlist(invcut_1[i].pLists); invcut_1[i].pLists=NULL;}
  nCuts=0;
}


static int addcut(int Mirr,char* title,int minonw,int maxonw,double cvminw,double cvmaxw)
{ int ok;    
  if(nCuts >= 60) return 2;
  if(strlen(title)>=50) return 3;
  ok=checkPhysValN(title,invcut_1[nCuts].key,&(invcut_1[nCuts].pLists));
  if(!ok) return 4;
  if(!invcut_1[nCuts].pLists) return 0;
    
  strcpy(invcut_1[nCuts].title,title); 
  invcut_1[nCuts].minon = minonw;
  invcut_1[nCuts].maxon = maxonw;
  invcut_1[nCuts].cvmin = cvminw;
  invcut_1[nCuts].cvmax = cvmaxw;
  invcut_1[nCuts].aux=Mirr;
  nCuts++;
  return 0;
} /* addcut_ */


int wrtcut_(FILE *nchan){ fprintf(nchan,"\n"); writetable0(&cutTab,nchan); return 0;}
int rdrcut_(FILE *nchan){ fscanf(nchan,"\n");  return readtable0(&cutTab,nchan); }

int fillCutArray(void)
{
  linelist ln;
  int lineNum=0;
  char cutStr[STRSIZ], minStr[STRSIZ], maxStr[STRSIZ],auxStr[4];
  int minOn,maxOn,aux;
  char fieldName[50];
  
  cutInit();
  for(ln=cutTab.strings ;ln;ln=ln->next)
  {  
    double min_=0,max_=0;
    int k;
    auxStr[0]=0;
    cutStr[0]=0;
    minStr[0]=0;
    maxStr[0]=0;    
    lineNum++;    
    sscanf(ln->line,"%[^|]%*c%[^|]%*c%[^|]%*c%[^|]",auxStr,cutStr,minStr,maxStr);
    trim(auxStr); trim(cutStr); trim(minStr); trim(maxStr);
    
/*================ ! ============*/
    strcpy(fieldName,"Wrong field '!'");
    switch(auxStr[0])
    { case  0 : aux=1; break;
      case '!': aux=-1; break;
      case '%': aux=0; break;
      default : goto errorExit;
    }  
                                      
/*================ MIN bound ============*/
    strcpy(fieldName,"Wrong field 'Min. bound'");
    minOn= strlen(minStr);
    if(minOn) { if(calcExpression(minStr,rd_num,&min_)) goto errorExit;}

/*============ Parameter ===========*/
    strcpy(fieldName,"Wrong field 'Parameter'");
    if(!strlen(cutStr)) goto errorExit;

/*================== MAX bound ==========*/    
    strcpy(fieldName,"Wrong field 'Max bound'");
    maxOn= strlen(maxStr);
    if(maxOn) {if(calcExpression(maxStr,rd_num,&max_)) goto errorExit;}

/* =========== fill array ==========*/

    if(minOn||maxOn)
    {  k=addcut(aux,cutStr,minOn,maxOn,min_,max_);
       switch(k)
       { case 2: strcpy(fieldName,"Too many cuts ");         goto errorExit;
         case 3: strcpy(fieldName,"Too long cut identifier");goto errorExit; 
         case 4: strcpy(fieldName,"Error in parameter");     goto errorExit;
       }  
    }   
  } 
  return 0;
  errorExit:
  { char buff[300];
    sprintf(buff," Error in Cut table line %d \n%s\n%s",lineNum,fieldName,errorText);
    
    if(blind) { printf("%s\n",buff); sortie(150);} else messanykey(2,10,buff);
    return 1;
  }     
}


int rancor_(double *vmin, double *vmax, double shift, double fmult,int nn)
{
    static double vnew,v;

    if(!nn) return 0;
    if(invcut_1[nn-1].aux!=1) return 0;
    if(invcut_1[nn-1].minon && invcut_1[nn-1].key[1]!='^') 
    {  v=invcut_1[nn-1].cvmin;
       vnew = (v*v - shift) * fmult;
       if(fmult > 0) *vmin=MAX(*vmin,vnew); else *vmax=MIN(*vmax,vnew);
    }
    if (invcut_1[nn-1].maxon && invcut_1[nn-1].key[1]!='_') 
    {  v=invcut_1[nn-1].cvmax;
       vnew = (v*v - shift) * fmult;
       if(fmult > 0) *vmax=MIN(*vmax,vnew); else *vmin=MAX(*vmin,vnew);	
    }
    return 0;
} /* rancor_ */

void  rancor_t(double * cosmax, double hsum , double fmult, double Ecm, 
       double mq, double pcmtilda,  double mp1, double ptMin)
{
    double Eq,mt2,Sinv,csnew;
    Eq=sqrt(pcmtilda*pcmtilda+mq*mq);
    mt2=mp1*mp1 + ptMin*ptMin;
    if(Ecm*Ecm < mt2 ) {*cosmax=-1; return; }
    if(Ecm*mq > Eq*sqrt(mt2)) Ecm=Eq*sqrt(mt2)/mq;

    Sinv= mq*mq + mp1*mp1 -2* (Ecm*Ecm*mq*mq + pcmtilda*pcmtilda*mt2 )/
                 (Ecm*Eq + pcmtilda*sqrt(Ecm*Ecm-mt2));

    csnew=(Sinv- hsum) * fabs(fmult);
    if(csnew<*cosmax) *cosmax=csnew; 
}



double  calcCutFactor(double*V)
{  
  int i;
  double val,valM;

  for(i=0; i<nCuts; i++)
  { physValRec* pList;
    for(pList=invcut_1[i].pLists;pList;pList=pList->next)
    {  
       val=calcPhysVal(invcut_1[i].key[0],pList->pstr,V);
       if(!isfinite(val)) continue;
       switch(invcut_1[i].key[1])
       { case '^': if(i){ if(val>valM) valM=val;}  else valM=val; break;
         case '_': if(i){ if(val<valM) valM=val;}  else valM=val; break;
         default:
         if(invcut_1[i].aux==1)
         {
           if ( invcut_1[i].minon && (val < invcut_1[i].cvmin)) return 0;
           if ( invcut_1[i].maxon && (val > invcut_1[i].cvmax)) return 0;
         }
         if(invcut_1[i].aux==-1)
         { 
           if(invcut_1[i].minon && invcut_1[i].maxon)
           {    if( val>= invcut_1[i].cvmin 
                 && val <= invcut_1[i].cvmax)  return 0;
           }
           else if (invcut_1[i].minon) 
           {    if(val >= invcut_1[i].cvmin) return 0;}
           else if(val <= invcut_1[i].cvmax) return 0;
         }                      
       }
    }
    if(invcut_1[i].key[1])
    {  if(invcut_1[i].aux==1)   
       {
         if ( invcut_1[i].minon && (valM < invcut_1[i].cvmin)) return 0;
         if ( invcut_1[i].maxon && (valM > invcut_1[i].cvmax)) return 0;
       }
       else  if(invcut_1[i].aux==-1)
       { 
           if(invcut_1[i].minon && invcut_1[i].maxon)
           {    if( valM>= invcut_1[i].cvmin 
                 && valM <= invcut_1[i].cvmax)  return 0;
           }
           else if (invcut_1[i].minon) 
           {    if(valM >= invcut_1[i].cvmin) return 0;}
           else if(valM <= invcut_1[i].cvmax) return 0;
       }                       
    } 
  }
  return 1;
}  
