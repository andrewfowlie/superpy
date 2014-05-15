/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include <unistd.h>
#include "chep_crt.h"
#include "getmem.h"
#include "syst2.h"
#include "s_files.h"
#include "physics.h"
#include "sos.h"
#include "c_out.h"      /* consLow */
#include "prepdiag.h"   /* longDouble */
#include "process.h"
#include "cweight.h"
#include "dynamic_cs.h"
#include "VandP.h"
#ifdef pow
#undef pow
#endif

int forceUG=0;
int newCodes=0;
static void init_safe(int * exitlevel)
{  int i;

   strcpy(processch,"e,E -> m,M");
   for(i=0;i<MAXINOUT;i++) 
   { strcpy(hadrons[i].name,"");
     strcpy(hadrons[i].contents,"");
     hadrons[i].pow=0;
   } 
   strcpy(limpch,"");
   strcpy(deloutch,"");
   n_model = 1;
   newCodes=0;
   *exitlevel = 0;
}

void restoreent(int * exitlevel)
{
   FILE * ff=fopen("tmp/safe","r"); 

   if (ff != NULL)
   { int OK=1;
     int i,ntot;
     
      if(1!=fscanf(ff,"#Model %d\n",&n_model) 
      || 1!=fscanf(ff,"#ForceUG %d\n",&forceUG)
      || 1!=fscanf(ff,"#nIn %d\n",&nin)       
      || 1!=fscanf(ff,"#nOut %d\n",&nout)     
      || 1!=fscanf(ff,"#Process%[^\n]\n",processch) 
      || 1!=fscanf(ff,"#nTot %d\n",&ntot)
        ) OK=0;
      if(!OK || ntot>MAXINOUT) { init_safe(exitlevel); return;}     
      
      for(i=0;i<ntot;i++)
      {
        OK=fscanf(ff,"*%[^=]=%[^\n]\n",hadrons[i].name,hadrons[i].contents);   
        if(OK!=2) { init_safe(exitlevel); return;}
      }  
      
      if(1!=fscanf(ff,"#Remove_Virtual%[^\n]\n",limpch) 
       ||1!=fscanf(ff,"#Remove_X%[^\n]\n",deloutch)
       ||1!=fscanf(ff,"#nSubproc(ampl) %d\n", &subproc_f)
       ||1!=fscanf(ff,"#nSubproc(squared) %d\n", &subproc_sq)
/*       ||1!=fscanf(ff,"#ConservationLow %d\n",&consLow)  */
       ||1!=fscanf(ff,"#Nc==inf  %d\n",&NcInfLimit)
       ||1!=fscanf(ff,"#NoColorChains %d\n",&noCChain)
       ||1!=fscanf(ff,"#NoDiagrams  %d\n", &noPict)
       ||1!=fscanf(ff,"#T-widths    %d\n", &tWidths)    
       ||2!=fscanf(ff,"#VVdecays  %d  %d \n", &VWdecay,&VZdecay)
       ||1!=fscanf(ff,"#NewCodes  %d\n", &newCodes) 
       ||1!=fscanf(ff,"#ExitCode %d\n",exitlevel) )  init_safe(exitlevel);

      fclose(ff);
      trim(processch); trim(limpch); trim(deloutch);
      
   } else init_safe(exitlevel);
}


void saveent(int  exitlevel)
{  FILE * ff; 
   int i;
   int ntot=0;
   if(strcmp(outputDir,"results/")) chdir("../"); 
   ff=fopen("tmp/safe","w");
   fprintf(ff,"#Model %d\n",n_model);
   fprintf(ff,"#ForceUG %d\n",forceUG);
   fprintf(ff,"#nIn %d\n",nin);
   fprintf(ff,"#nOut %d\n",nout);
   fprintf(ff,"#Process %s\n",processch);
   for(ntot=0;ntot<MAXINOUT && hadrons[ntot].name[0] ;ntot++)continue;
     fprintf(ff,"#nTot %d\n",ntot);
   for(i=0;i<ntot;i++) 
   fprintf(ff,"*%s=%s\n",hadrons[i].name,hadrons[i].contents);
   fprintf(ff,"#Remove_Virtual %s\n",limpch);
   fprintf(ff,"#Remove_X %s\n",deloutch);
   fprintf(ff,"#nSubproc(ampl) %d\n", subproc_f);
   fprintf(ff,"#nSubproc(squared) %d\n", subproc_sq);
/*   fprintf(ff,"#ConservationLow %d\n",consLow); */
   fprintf(ff,"#Nc==inf  %d\n",NcInfLimit);
   fprintf(ff,"#NoColorChains %d\n",noCChain);
   fprintf(ff,"#NoDiagrams  %d\n", noPict); 
   fprintf(ff,"#T-widths    %d\n", tWidths);
   fprintf(ff,"#VVdecays    %d  %d\n", VWdecay,VZdecay);
   fprintf(ff,"#NewCodes  %d\n", newCodes);   
   fprintf(ff,"#ExitCode %d\n",exitlevel);
//   fprintf(ff,"#Model paramters:\n");
//   { for(i=0;i<nModelVars;i++) fprintf(ff,"%s\n",varNames[i]);}
   fclose(ff);
}

void  save_sos(int ercode)
{
   unsigned         nproc;
   csdiagram    cd;
   marktp mark;

   if (ercode == -2) /* User Break */
   {
      finish();
      sortie(20);  /*  Restart  */
   }
 
   if (ercode == -1) /* Heap is empty */
   { mark.blk_=NULL;
     mark.pos_=0;
     release_(&mark);
   } /* Heap is empty, continue */
     
     /*  TooLargeNumber  */
     /*  TooManyIdentifiers  */
     /*  RangeCheckError  */
   if ((ercode < 0)   || (ercode == 7) ||
       (ercode == 11) || (ercode == 12))   /*  Restart  */
   {
      saveent(10);
      nproc = ftell(diagrq) - sizeof(cd);
      fseek(diagrq,nproc,SEEK_SET);
      FREAD1(cd,diagrq);
      cd.status = -2;
      fseek(diagrq,nproc,SEEK_SET);
      FWRITE1(cd,diagrq);
      finish();
      sortie(20);  /*  Restart  */
   }

/*  not disk space  */
   if ((ercode == 40))
   {
      finish();
      sortie(65);  /*  End of work  */
   }
   if(ercode ==14)   
   {
      messanykey(10,10," Check model !");
      finish();
      sortie(62);  
   }
}
