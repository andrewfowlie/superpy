/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <limits.h>

#include"chep_crt.h"
#include "syst2.h"
#include "crt.h"
#include "physics.h"
#include "pvars.h"
#include "sos.h"
#include "s_files.h"
#include "ghosts.h"
#include "rfactor.h"
#include "polynom.h"
#include "screen.h"

#include "saveres.h"

/*
 int     den_power[2 * maxvert - 2];
 short   den_width[2 * maxvert - 2];
 short   den_mass[2 * maxvert - 2];
*/

 denom_struct   denom[2 * maxvert - 2];
 
 int     denrno;

 char    denStr[2 * maxvert - 2][MAXINOUT];

static void wAbort(void)
{
      saveent(menulevel);
      messanykey(5,20,"Error in writing on the disk. \n"
                      "Check the existence of the \n" 
                      "'tmp' and 'results' directories \n"
                      "or the existence of free disk space");
      finish();
      sortie(65);  /*  End of work  */
}



static void savevardef(void)
{ 
  FWRITE1(vardef->nvar,archiv);
  if(vardef->nvar && 
        fwrite(vardef->vars, (vardef->nvar)*sizeof(varinfo) ,1,archiv)!=1)
   wAbort();
}

static void  savepoly(poly p)
{

 int width;
 char * b, *e;
 poly zero=plusone();

 zero->num=0;    
 
 b=(char*)&(zero->num);
 if(vardef->nvar) e=(char*)&(zero->power[vardef->vars[vardef->nvar-1].wordpos]); 
       else       e=(char*)&(zero->power[0]); 
 width=e-b;
 
 while(p) { if(fwrite(&(p->num),width,1,archiv)!=1) wAbort(); p=p->next; }

 if(fwrite(&(zero->num),width,1,archiv)!=1) wAbort();  

}


void  saveanaliticresult(poly rnum,poly factn,poly factd, vcsect vcs, int nFile)
{catrec      cr;
 int        i;
 int m;

   diskerror = wAbort;

   cr.nsub_ = nsub;
   cr.ndiagr_ = ndiagr;
       
   cr.factpos = ftell(archiv);
   cr.nFile=nFile;
   vardef++;
   savevardef();   
   savepoly(factn);
   savepoly(factd);
   cr.rnumpos = ftell(archiv);

   vardef--;
   savevardef();
   savepoly(rnum);

   cr.denompos = ftell(archiv);
   calcdenominators(vcs );   

   FWRITE1(denrno,archiv);   /*  number of demominatirs  */
 
   for (i = 0; i < denrno; i++)
   {
       FWRITE1(denom[i].power,archiv);   /*  power  1 or 2  */
       FWRITE1(denom[i].mass,archiv);
       FWRITE1(denom[i].width,archiv);
       m=0;  do FWRITE1(denom[i].momStr[m],archiv); while(denom[i].momStr[m++]);
   }
   FWRITE1(cr,catalog);
   diskerror = NULL;
}
