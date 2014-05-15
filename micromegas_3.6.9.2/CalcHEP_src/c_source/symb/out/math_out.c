/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include"chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "s_files.h"
#include "procvar.h"
#include "out_serv.h"
#include "getmem.h"
#include "saveres.h"
#include "process.h"
#include "writeF.h"
#include "math_out.h"


static void  writeprocessname(int* prtclNum)
{  int i; 

   writeF("  process  ");
   for(i=1;i<=nin;i++) 
   {
     writeF("%s(p%d)",prtclbase[prtclNum[i]-1].name,i);
     if (i<nin) writeF("+"); else writeF("->");
   }
   
   for(i=nin+1;i<=nin+nout;i++) 
   {  writeF("%s(p%d)",prtclbase[prtclNum[i]-1].name,i);
      if (i<nin+nout) writeF("+");else writeF("\n");
   }
} 


static void  emitprocessname(int * prtclNum)
{
   int i;

   writeF("inParticles = {");
   for(i=1;i<=nin;i++) 
   {
      writeF("\"%s\"",prtclbase[prtclNum[i]-1].name);
      if (i<nin) writeF(","); else writeF("}\n");
   }
   
   writeF("outParticles = {");
   
   for(i=nin+1;i<=nin+nout;i++) 
   {  writeF("\"%s\"",prtclbase[prtclNum[i]-1].name);
      if (i<nin+nout) writeF(",");else writeF("}\n");
   }

} 



static void  emitexpression(catrec*  cr)
{ 
  int       i;

   seekArchiv(cr->factpos);
   readvardef(archiv);
   writeF("totFactor = ((");
   rewritepolynom();
   writeF(")/(");
   rewritepolynom();
   writeF("));"); writeF("\n");
   clearvardef();
   
   seekArchiv(cr->rnumpos);
   
   readvardef(archiv);
   writeF("numerator =(");rewritepolynom();writeF(");\n");
   clearvardef();
   seekArchiv(cr->denompos);
   readDenominators();

   writeF("denominator =");  
   for (i = 0; i < denrno; i++)
   {
      char momStr[20];
      if(i)  writeF("*"); else writeF("(");
      momentToString(denom[i].momStr,momStr);  
      writeF("propDen[%s,%s,%s]",momStr,vararr[denom[i].mass].alias,
                                        vararr[denom[i].width].alias);
      if(denom[i].power >1)  writeF("^%d",denom[i].power);
   }
   if(i) writeF(");\n"); else writeF("1;\n");
}


static int modifyF(char * txt)
{ 
  char * Br;
  int extF=0;
  int i,c;
  for(i=0,c=0;txt[i];i++) if(txt[i]!=' ') txt[c++]=txt[i]; 
  txt[c]=0;
  
  for(Br=strchr(txt,'('); Br ;Br=strchr(Br+1,'('))
  {
     if(Br==txt || !(isalnum(Br[-1]) || Br[-1]=='_')) continue; 
     
     for(i=-1; Br+i>=txt && (isalnum(Br[i]) || Br[i]=='_');i--);
     i++;
     Br[0]=0;
     if(strcmp(Br+i,"sqrt")==0 ||
        strcmp(Br+i,"sin" )==0 ||
        strcmp(Br+i,"cos" )==0 ||
        strcmp(Br+i,"tan" )==0 ||
        strcmp(Br+i,"asin")==0 ||
        strcmp(Br+i,"acos")==0 ||
        strcmp(Br+i,"atan")==0 ||
        strcmp(Br+i,"exp" )==0 ||
        strcmp(Br+i,"log" )==0     ) Br[i]=toupper(Br[i]);
     else if( strcmp(Br+i,"one" )==0) strcpy(Br+i,"1  ");    
     else extF=1;
     if(Br[-1]==' ') /* one ! */
     { Br[0]=' ';    
       for(i=1,c=1;c;i++) { if(Br[i]=='(')c++; else if(Br[i]==')')c--;  Br[i]=' ';}    
     } else     
     { Br[0]='[';    
       for(i=1,c=1;c;i++)  if(Br[i]=='(')c++; else if(Br[i]==')')c--;
       Br[i-1]=']';
     }         
  }

  for(i=0,c=0;txt[i];i++) if(txt[i]!=' ') txt[c++]=txt[i]; 
  txt[c]=0;
  
  return extF; 
}


static void  writeparameters(int nsub)
{ 
   int   k = 0; 
   int   first = 1; 
   char  s[STRSIZ]; 
   char  *lch;
 
   writeF("\n"); 
   writeF("parameters={\n"); 

   for(k=1;k<=nmodelvar;k++) if(vararr[k].used)
   {
     if(modelvars[k].func )
     {  
        sscanf(modelvars[k].func,"%[^%|\n]",s); 
        if(!modifyF(s)) continue;                  
     }else 
     {
        sprintf(s,"%17.11E", vararr[k].tmpvalue);
        lch = strchr(s,'E');
        if (lch) 
        {  int d;
           sscanf(lch+1,"%d",&d);
           sprintf(lch,"*10^(%d)",d);
        }   
     }  
     if(first) { first=0; writeF(" ");} else writeF(",");
     writeF("%s -> %s\n",vararr[k].alias,s); 
   } 

   writeF("           };\n"); 
   writeF("\n"); 
   first = 1; 
   writeF("substitutions={\n");
   
   for(k=nmodelvar;k;k--)
   {   
      if (vararr[k].used &&  modelvars[k].func)
      {  sscanf(modelvars[k].func,"%[^%|\n]",s);
         if(modifyF(s)) continue;
         if(first) { first=0; writeF(" ");} else writeF(",");
         writeF("%s->%s",vararr[k].alias,s);
         writeF("\n");
      }
   }
   writeF("              };\n");
} 


static void startMath(int nsub, int* prtclNum,int ncalc)
{
   outputLanguage='M';
   initvararray(nsub,outputLanguage,3);

   outFileOpen("%sresults%csymb%d.m",pathtouser,f_slash,nsub); 

   writeF("(*\n"); writeLabel(' ');
   writeprocessname(prtclNum);
   writeF("*)\n"); 
   writeparameters(nsub); writeF("\n");
   emitprocessname(prtclNum); 
   writeF("\n"); 
   writeF("SetAttributes[ SC, Orderless ];\n"); 
   writeF("\n"); 
   writeF("SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];\n"); 
   writeF("\n"); 
   writeF("SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]\n"); 
   writeF("\n"); 
   writeF("\n"); 

   emitconvlow(prtclNum);
   
   writeF("\ninitSum;\n");
}


static void  diagramMath(vcsect * vcs,catrec * cr)
{ 
   writeF("\n(*\n"); 
   writeF("  Diagram  %d in subprocess %d\n",cr->ndiagr_,cr->nsub_);               
   if (vcs != NULL)  DiagramToOutFile(vcs,0,' ');  
   writeF("*)\n");
   emitexpression(cr);
   
   writeF("\naddToSum;\n");
}

static void endMath(int * prtclNum)
{
   writeF("\nfinishSum;\n");
   outFileClose();
}

void makeMathOutput(void)
{ makeOutput(startMath,diagramMath,endMath);}
