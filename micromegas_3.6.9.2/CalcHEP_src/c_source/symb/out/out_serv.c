/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <limits.h>

#include "chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "s_files.h"
#include "diaprins.h"
#include "procvar.h"
#include "lnum.h"
#include "pvars.h"
#include "polynom.h"
#include "pre_read.h"
#include "parser.h"
#include "saveres.h"
#include "process.h"
#include "writeF.h"
#include "../../../include/version.h"
#include "out_serv.h"



static  polyvars  varsInfo={0,NULL};



int readSize;
poly readBuff=NULL;

int outputLanguage=0;

void seekArchiv(long n)
{  fseek(archiv,n,SEEK_SET);}


void DiagramToOutFile(vcsect * vcs, int label,char comment)
{
   if (xpos !=1) writeF("\n");
    writeTextDiagram(vcs,label,comment,outFile);
}

void readvardef(FILE*archive)
{
  monom template;
  char * end;
  
  vardef=&varsInfo;  
  fread(&vardef->nvar,sizeof(vardef->nvar),1,archive);

  if (vardef->nvar) 
  { vardef->vars= m_alloc(vardef->nvar*sizeof(varinfo));
    fread(vardef->vars,vardef->nvar*sizeof(varinfo),1,archive);
    end=(char*)&template.power[vardef->vars[vardef->nvar-1].wordpos];
  } else 
  { end=(char*)&template.power[0];
    vardef->vars=NULL;
  }
  readSize = end- (char*)&template.num;
  if(readBuff) free(readBuff);
  readBuff = m_alloc(end- (char*)&template);
}

void clearvardef(void)
{
  clearVars(vardef);
  if(readBuff) {free(readBuff); readBuff=NULL;}
}


static int  readmonom(char* txt)
{NUM_TYPE     l;

 int       deg,n;
 
   fread(&readBuff->num,readSize,1,archiv);
   if(!txt) return readBuff->num;
   
   l=readBuff->num;
   
   if(!l) return 0;
         
   if( l == 1 || l == -1) txt[0]=0; 
   else 
   { 
#ifdef NUM_DOUBLE   
     if(ABS(l)>1.E15) sprintf(txt,"%+"NUM_STR".",l);
     else       
#endif
                      sprintf(txt,"%+"NUM_STR,l);
   }

   for (n = 0; n < vardef->nvar; n++)
   {
       deg = (readBuff->power[vardef->vars[n].wordpos-1] /
       vardef->vars[n].zerodeg) %
       vardef->vars[n].maxdeg;
       if(deg)
       { sprintf(txt+strlen(txt),"*%s",vararr[vardef->vars[n].num].alias);
         if(deg>1)sprintf(txt+strlen(txt),"^%d",deg);
       }
   }

   if(!txt[0]) 
   {
#ifdef NUM_DOUBLE   
     if(ABS(l)>1.E15) sprintf(txt,"%+"NUM_STR".",l);
     else       
#endif
                      sprintf(txt,"%+"NUM_STR,l);
   }
   else
   {       if (l ==  1) txt[0]='+';
      else if (l == -1) txt[0]='-';
   }
   return 1;
}





void readDenominators(void) 
{
   int i,m;

   FREAD1(denrno,archiv);   /*  number of demominatirs  */
   for (i = 0; i < denrno; i++)
   {
       FREAD1(denom[i].power,archiv);   /*  power  1 or 2  */
       FREAD1(denom[i].mass,archiv);
       FREAD1(denom[i].width,archiv);
       m=0;  
       do FREAD1(denom[i].momStr[m],archiv); while(denom[i].momStr[m++]);
   }
}



void momentToString(char * momStr, char * outstr)
{
  int m=0;

  strcpy(outstr,"");
   while(momStr[m]) 
   {   
   if(momStr[m]<=nin) strcat(outstr,"-");
      else if(m)   strcat(outstr,"+");
      sprintf(outstr+strlen(outstr),"p%d",momStr[m++]);
   }    
}


void  rewritepolynom(void)
{ 
  char  monomtxt[STRSIZ];
  if(!readmonom(monomtxt)) {wrt_0("0");return;}
  if (monomtxt[0] == '+') wrt_0(monomtxt+1); else wrt_0(monomtxt);                                            
  while(readmonom(monomtxt)) wrt_0(monomtxt);
}


void findPrtclNum (char * procName, int * prtclNum)
{
 int	k;
 char * pname,txt1[STRSIZ]; 
 int pnum;
 
   strcpy(txt1,procName);
   memcpy(strstr(txt1,"->"),", ",2);
   pname=strtok(txt1,",");
   k=1;
   while ( pname != NULL)
   { 
     trim(pname);locateinbase(pname,&pnum);
     prtclNum[k]=pnum; k++;       
     pname=strtok(NULL,",");
   } 
}

void  emitconvlow(int* prtclNum)
{
 int	i, j, ntot;
 char  	pmass[MAXINOUT+1][11]  ; /* for locateinbase */
 char   s[MAXINOUT+1]; 
 
   ntot = nin + nout;
   for(i=1;i<=ntot;i++) strcpy(pmass[i],prtclbase[prtclNum[i]-1].massidnt);

   for (i = 1; i <= nin; i++) s[i] = 1;
   for (i = nin + 1; i <= ntot; i++) s[i] = -1;
   
   switch (outputLanguage) 
   {case 'R': writeF("\nlet p%d = ",ntot); break;
    case 'F': writeF("\nid p%d = ",ntot); break;
    case 'M': writeF("\np%d = ",ntot); break;
   }
   for (i = 1; i <= nin ; i++) writeF("+p%d",i);
   for (i = 1; i <= nout - 1; i++) writeF("-p%d",i+nin);
   writeF(";\n");

   for (i = 1; i <= ntot - 1; i++)
   switch (outputLanguage)     
   {case 'R': writeF("mass p%d  = %s; Mshell p%d;\n",i,pmass[i],i); break;
    case 'F': writeF("id p%d.p%d  = %s^2;\n",i,i,pmass[i]); break;
    case 'M': writeF("p%d/: SC[p%d,p%d] =%s^2;\n",i,i,i,pmass[i]); break; 
   }
              
   switch (outputLanguage)     
   {case 'R':writeF("let p%d.p%d = ",ntot - 2,ntot - 1); break;
    case 'F':writeF("id p%d.p%d = ",ntot - 2,ntot - 1);break;
    case 'M':writeF("p%d/: SC[p%d,p%d] = ",ntot-2,ntot - 2,ntot - 1);break;
   }           

   writeF("%d*(%s^2",s[ntot-1]*s[ntot-2],pmass[ntot]);
   for (i = 1; i <= ntot - 1; i++) writeF("-%s^2",pmass[i]);
   for (i = 2; i <= ntot - 1; i++)
   for (j = 1; j <= i - 1; j++)
   if (j < (ntot - 2))    switch (outputLanguage)     
   {case 'R':
    case 'F':writeF("%+d*p%d.p%d",-2*s[j]*s[i],j,i);break;
    case 'M':writeF("%+d*SC[p%d,p%d]",-2*s[j]*s[i],j,i);break; 
   }
   writeF(")/2;\n");
}

void writeLabel(char comment)
{
 if (xpos !=1) writeF("\n");
 fprintf(outFile,"%c    ==============================\n",comment );
 fprintf(outFile,"%c    *  %s *\n",comment,VERSION_);
 fprintf(outFile,"%c    ==============================\n",comment);         
}



void  makeOutput(  void (*startOutput)(int,int*,int),
                   void (*diagramOutput)(vcsect*, catrec* ),
                   void (*endOutput)(int*)            
                )
{  catrec    cr;
   int  ndel, ncalc, nrest;
   long recpos;
   long   count, ntot;
   int graphOn;
   shortstr  txt;
   int prtclNum[MAXINOUT+1];
   csdiagram   csdiagr;
   vcsect      vcs; 
 
   informline(0,1);
   catalog=fopen(CATALOG_NAME,"rb");
   menuq=fopen(MENUQ_NAME,"rb");
   count = 0;
   graphOn=1;
 
   if (graphOn) diagrq=fopen(DIAGRQ_NAME,"rb");

   for(nsub=1,ntot=0;nsub<=subproc_sq;nsub++)   
   { rd_menu(2,nsub,txt,&ndel,&ncalc,&nrest,&recpos);
     ntot+=ncalc;
   }
 
   for(nsub=1;nsub<=subproc_sq;nsub++) 
   {
      rd_menu(2,nsub,txt,&ndel,&ncalc,&nrest,&recpos);
      findPrtclNum (txt, prtclNum);
      if(ncalc != 0)
      { 
         startOutput(nsub,prtclNum,ncalc);  
         fseek(catalog,0,SEEK_SET);
         while (FREAD1(cr,catalog))
         {
            if (cr.nsub_ == nsub)
            {  
               whichArchive(cr.nFile,'r');      
               if (graphOn)
               {  fseek(diagrq,(cr.ndiagr_+recpos-1)*sizeof(csdiagram),SEEK_SET);
                  FREAD1(csdiagr,diagrq);
                  transfdiagr(&csdiagr,&vcs);
                  diagramOutput(&vcs,&cr);
               }
                else diagramOutput(NULL,&cr);
                ++(count); 
	       if(informline(count,ntot)) goto escexit;
            }
         }
escexit: 
         endOutput(prtclNum);
      }
   }
   fclose(catalog);
/*   if(ArcNum)fclose(archiv);*/
   whichArchive(0,0);
   fclose(menuq);
   if(graphOn) fclose(diagrq);
   informline(ntot,ntot);
}
