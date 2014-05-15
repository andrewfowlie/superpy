#ifndef  __S_FILES__ 
#define __S_FILES__

#include"files.h"
#include "syst2.h"

typedef struct catrec
   {
      int      nsub_, ndiagr_,nFile;
      long      factpos, rnumpos, denompos;
   }      catrec;

    extern FILE *  menup;
    extern FILE *  menuq;
    extern FILE * diagrp;   /* file of Adiagram; */
    extern FILE * diagrq;   /* file of CSdiagram; */

    extern FILE * archiv;
    extern FILE * catalog;
     
    extern char* mdFls[5];
    extern shortstr  pathtouser;
    
#define MENUQ_NAME   "./tmp/menuq.ch"  
#define MENUP_NAME   "./tmp/menup.ch"               
#define DIAGRP_NAME  "./tmp/proces.tp"
#define DIAGRQ_NAME  "./tmp/csproces.tp"
#define ARCHIV_NAME  "./tmp/archive_%d.bt"
#define CATALOG_NAME "./tmp/catalog.tp"


extern void  wrt_menu(int menutype, int k,
              char*txt,int ndel,int ncalc,int nrest,long recpos);
extern int  rd_menu(int   menutype, int k,
              char*txt,int*ndel,int*ncalc,int*nrest,long*recpos);

extern int whichArchive(int nFile, int rw);
#define  FREAD1(d,f)   fread(&(d),sizeof(d),1,f)
#define  FWRITE1(d,f)  f_write(&(d),sizeof(d),1,f)
#define  MAXARCHIVE 100000000L


#endif
