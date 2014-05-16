/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "interface.h"
#include"chep_crt.h"
#include"parser.h"
#include"phys_val.h"
#include"read_func.h"
#include"rd_num.h"
#include"regul.h"
/*
#include"VandP.h"
#include"dynamic_cs.h"
*/

table regTab={"*** Table ***"," Regularization ",
" Momentum    |> Mass  <|> Width <| Power|",NULL,0};

invreg_ invreg_1[200];

static void inireg_(void) { invreg_1[0].lvinvr[0] = 0;} 

static int addreg_(char *lv, double rgmasw, double rgwdtw, int nndeg)
{
     int nreg=0;
     int lastrg=-1;
     
     for( ;invreg_1[nreg].lvinvr[0]; nreg++)
     {  if (eqvect_(lv,invreg_1[nreg].lvinvr)) lastrg = nreg; }
    if (nreg >= 199)  return 0;
    strcpy(invreg_1[nreg].lvinvr,lv);
    invreg_1[nreg].rgmass = rgmasw;
    invreg_1[nreg].rgwdth = rgwdtw;
    invreg_1[nreg].ndeg = nndeg;
    invreg_1[nreg].nextrg = 0;
    if (lastrg >=0 ) invreg_1[lastrg].nextrg = nreg+1;  
    invreg_1[nreg+1].lvinvr[0]=0;
    return 0;
} /* addreg_ */



int fillRegArray(void)
{ char charKey; 
  linelist ln=regTab.strings;
  int lineNum=0;
  double mass,width;
  char invStr[STRSIZ], massStr[STRSIZ], widthStr[STRSIZ];

  inireg_();
  while (ln != NULL)
  {  
    int power=0;
    char lv[PLISTLEN]="";
    invStr[1]=0;
    massStr[0]=0;
    widthStr[0]=0;    
    lineNum++;
        
    sscanf(ln->line,"%[^|]%*c%[^|]%*c%[^|]%*c%d",invStr+1,massStr,widthStr,&power);

/*============ Invariant ===========*/
    trim(invStr+1);
    invStr[0]='S';
    if( !checkPhysVal(invStr,&charKey, lv) )
    {
       sprintf(errorText," Error in  regularization table line %d .\n"
                              " Wrong field 'Momentum' .",lineNum);
       goto errorExit;                      
    }

    coninv_(lv);
/*================ Mass ============*/
    if( calcExpression(massStr,rd_num,&mass) )
    {    sprintf(errorText," Error in  regularization table line %d .\n"
                          " Wrong field 'Mass' .",lineNum);
         goto errorExit;
    }                                         
/*==================Width ==========*/    
    if( calcExpression(widthStr,rd_num,&width) )   
    {    sprintf(errorText," Error in  regularization table line %d .\n"
                          " Wrong field 'Width' .",lineNum);
         goto errorExit;
    }                                         
          
/*============ Power ===============*/     
    if( power<1 ||power>2 ) 
    { 
       sprintf(errorText," Error in  regularization table line %d .\n"
                         " Power is out of range.",lineNum);
       goto errorExit;                      
    }
    addreg_(lv,mass,width,power);     
    ln=ln->next;
  }
  
  
  return 0;
  errorExit: messanykey(2,10,errorText);
           return 1;  
}

int wrtreg_(FILE * nchan) { fprintf(nchan,"\n");writetable0(&regTab,nchan); return 0; }
int rdrreg_(FILE * nchan) { fscanf(nchan,"\n"); return readtable0(&regTab,nchan);}


int getreg_(int *nsing, sing_struct *singar, 
           double shift, double fmult, int  n)
{
   double d__1;

   for(;n && *nsing < 100 ;n = invreg_1[n-1].nextrg)
   {
      ++(*nsing);
      d__1 = invreg_1[n-1].rgmass;   /* Computing 2nd power */
      singar[*nsing - 1].pos   = (d__1 * d__1 - shift) * fmult;
      singar[*nsing - 1].width = (d__1 = invreg_1[n-1].rgmass 
      * invreg_1[n-1].rgwdth * fmult, ABS(d__1));
      singar[*nsing -1].power = invreg_1[n-1].ndeg;
   }
   return 0;
} /* getreg_ */
