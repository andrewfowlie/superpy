/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su
*/
#include"chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "s_files.h"
#include "procvar.h"
#include "out_serv.h"
#include "saveres.h"
#include "process.h"
#include "red_out.h"
#include "writeF.h"

static void writeprocessname(int * prtclNum)
{
   int i;
    
   writeF("inParticles:={"); 
   for(i=1; i<=nin; i++)   
     if(i==nin) writeF("\"%s\"}$\n",prtclbase1[prtclNum[i]].name);          
          else  writeF("\"%s\",",prtclbase1[prtclNum[i]].name);
      
   writeF("outParticles:={");
   for(i=nin+1; i<=nin+nout; i++)
     if(i==nin+nout) writeF("\"%s\"}$\n",prtclbase1[prtclNum[i]].name);
               else  writeF("\"%s\",",prtclbase1[prtclNum[i]].name);                        
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
        strcmp(Br+i,"log" )==0     ) Br[0]='(' ;
     else if( strcmp(Br+i,"one" )==0) strcpy(Br+i,"1  ");
     else { extF=1; Br[0]='(' ;}
      if(Br[-1]==' ') /* one ! */
     { Br[0]=' ';    
       for(i=1,c=1;c;i++) { if(Br[i]=='(')c++; else if(Br[i]==')')c--;  Br[i]=' ';}    
     }  
  }

  for(i=0,c=0;txt[i];i++) if(txt[i]!=' ') txt[c++]=txt[i]; 
  txt[c]=0;
  
  return extF; 
}


static void  writeparameters(int nsub)
{
   int        k;
   char        ch;
   char        s[STRSIZ];

   writeF("%%\n");
   writeF("parameters:={\n");

   for(k=1,ch=' ';k<=nmodelvar;k++)if(vararr[k].used)
   { if(modelvars[k].func)
     {   
       sscanf(modelvars[k].func,"%[^%|\n]",s);
       if(!modifyF(s)) continue;
     } else  sprintf(s,"%E",vararr[k].tmpvalue);

     writeF("%c%s=>%s\n",ch,vararr[k].alias,s); ch=',';       
   }
  
   writeF("}$\n");

   writeF("%%\n");
   writeF("substitutions:={\n");
   for(k=nmodelvar,ch = ' '; k ;k--) if(vararr[k].used && modelvars[k].func)
   { 
      sscanf(modelvars[k].func,"%[^%|\n]",s);
      if(modifyF(s)) continue;
      writeF("%c%s=>%s\n",ch,vararr[k].alias,s);
      ch=',';
   }
   writeF("}$\n");
}

static void  emitexpression(catrec* cr)
{
   int  i;

   seekArchiv(cr->factpos);
   readvardef(archiv);
   writeF("totFactor:=(");
   rewritepolynom();
   writeF(")/(");
   rewritepolynom();
   writeF(")$");
   writeF("\n");
   clearvardef();
   seekArchiv(cr->rnumpos);
   readvardef(archiv);
   writeF("numerator:="); rewritepolynom();writeF("$\n");
   clearvardef(); 

   seekArchiv(cr->denompos);
   readDenominators();

   writeF("denominator:=");
   if(denrno)    
     for (i = 0; i < denrno; i++)
     {  char momStr[20]; 
        momentToString(denom[i].momStr,momStr);   
      
        if(i) writeF("*");
        writeF("propDen(%s,%s,%s)",momStr, vararr[denom[i].mass].alias,
        vararr[denom[i].width].alias);
        if(denom[i].power!=1) writeF("^%d",denom[i].power);
     }
   else writeF("1");
   writeF("$\n");
}


static void startReduce(int nsub,int* prtclNum,int ncalc)
{  
   outputLanguage='R';
   initvararray(nsub,outputLanguage,3);
   outFileOpen("%sresults%csymb%d.red",pathtouser,f_slash,nsub);
   writeLabel('%');
   writeprocessname(prtclNum);
   writeparameters(nsub);
   writeF("\n\n");
   writeF("vector p1,p2,p3,p4,p5,p6$\n");
   emitconvlow(prtclNum);

   writeF("\nvector !=p_,!=q_$\n");
   writeF("operator propDen$\n");
   writeF("for all p_,q_,m,w let propDen(0*p_+q_,m,w)=propDen(q_,m,w)$\n");
   writeF("for all p_,m,w such that ordp(p_,-p_) "
                             "let propDen(p_,m,w)=propDen(-p_,m,w);$\n\n");
   writeF("initSum();\n");
}

static void diagramReduce( vcsect * vcs,  catrec * cr )
{
   writeF("\n");
   writeF("DiagrNumber:=\"%d_%d\"$\n",cr->nsub_,cr->ndiagr_);
   writeF("\n");
   if (vcs)  DiagramToOutFile(vcs,0,'%');
   emitexpression(cr);
   writeF("\n\n");
   writeF("addToSum()$\n");
}

static void endReduce(int * prtclNum)
{ 
  writeF("finishSum();\n");
  writeF("End$\n");
  outFileClose();
}

void makeReduceOutput(void)
{  makeOutput(startReduce,diagramReduce,endReduce);}

