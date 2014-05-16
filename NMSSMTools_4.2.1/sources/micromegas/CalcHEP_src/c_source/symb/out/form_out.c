/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include"chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "s_files.h"
#include "procvar.h"
#include "getmem.h"
#include "out_serv.h"
#include "saveres.h"
#include "form_out.h"
#include "writeF.h"

static void  writeparameters(int nsub)
{
 int        k;
 int     first;
 char        ch;


   writeF("\n");
   first = 1;
   ch = ' ';
   writeF("Symbols\n");

   for(k=1;k<=nmodelvar;k++)
   {  
      if ( vararr[k].used && ! modelvars[k].func )
      { 
         writeF( " %c%s\n",ch,vararr[k].alias);
         ch = ',';
      }
   }  
   writeF(" ;\n");
}


static void  writefunctions(int nsub)
{
 int        k;
 char        s[STRSIZ];
 
   writeF("\n");
   for(k=1;k<=nmodelvar;k++)
   {
      if (vararr[k].used && modelvars[k].func)
      {  sscanf(modelvars[k].func,"%[^|]",s);
         trim(s);
         writeF("id %s=%s;\n",vararr[k].alias,s);
      }
   }
}


static void startForm(int nsub, int* prtclNum,int ncalc)  
{  
  outputLanguage='F';
  initvararray(nsub,outputLanguage,3);
  outFileOpen("%sresults%cusr%d.frm",pathtouser, f_slash,nsub);
  
  writeF("\n#-\nCFunction Sqrt;\n#procedure userWork(nnn)\n");        

  writefunctions(nsub);
  
  emitconvlow(prtclNum);

  writeF("#if 'nnn' == %d\n.end\n#endif\n#endprocedure\n",ncalc);
  outFileClose();

  outFileOpen("%sresults%csum%d.frm",pathtouser, f_slash,nsub);                                        
  writeLabel('*');
  writeF("vector p1,p2,p3,p4,p5,p6;\n");         
  writeparameters(nsub);       
  writeF("CFunction den;\n");
  writeF("Symbol factor;\n");
  writeF("#include usr%d.frm\n",nsub);
}


static void  diagramForm(vcsect * vcs, catrec * cr )
{
 int         i;

   writeF("\n*Diagrama number %d-%d;\n",cr->nsub_,cr->ndiagr_); 

   if (vcs != NULL)  DiagramToOutFile(vcs,0,'*');

   seekArchiv(cr->factpos);
   readvardef(archiv);
   writeF("Local FACTOR%d = (",cr->ndiagr_); 
   rewritepolynom();writeF(")/(");rewritepolynom();writeF(");\n");
   clearvardef();
   writeF("Local ANS%d = factor*(",cr->ndiagr_);
   
   seekArchiv(cr->rnumpos);
   readvardef(archiv);
   rewritepolynom();writeF(")");
   clearvardef();
   seekArchiv(cr->denompos);
   readDenominators();

   for (i = 0; i <denrno; i++)
   {
      char momStr[20];
                   
      momentToString(denom[i].momStr,momStr);
      writeF("\n *den(%s,%s,%d)", momStr,
      vararr[denom[i].mass].alias, -denom[i].power);
   }
   writeF(";\n");
   writeF( "\n#call userWork{%d}\n",cr->ndiagr_);
}

static void endForm(int * prtclNum)
{
  writeF("*\n");
  outFileClose();
}  


void  makeFormOutput(void)
{ makeOutput(startForm,diagramForm,endForm); }

