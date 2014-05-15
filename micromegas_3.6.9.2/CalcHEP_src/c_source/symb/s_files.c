/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include"syst.h"
#include"syst2.h" 
#include"s_files.h"

 FILE * menup;
 FILE * menuq;
 FILE * diagrp;   /* file of Adiagram; */
 FILE * diagrq;   /* file of CSdiagram; */

 FILE * catalog;
 FILE * archiv;


char*mdFls[5] = {"vars","func","prtcls","lgrng","extlib"};
shortstr  pathtouser;

void wrt_menu(int menutype,int k,char*txt,int ndel,int ncalc,int nrest,long recpos)
{
   if (menutype == 1)
   {
      fseek(menup,(k - 1)*71 + 2,SEEK_SET);
      fprintf(menup,"%4d| %-44.44s|%5d|%5d|%-8d",k,txt,ndel,nrest,recpos);
   }
   else
   {
      fseek(menuq,(k - 1)*77 + 2,SEEK_SET);
      fprintf(menuq,"%4d| %-44.44s|%5d|%5d|%5d|%-8d",k,txt,ndel,ncalc,nrest,recpos);
   }
}


int rd_menu(int menutype,int k,char*txt,int*ndel,int*ncalc,int*nrest,long*recpos)
{
   if (menutype == 1)
   {
      fseek(menup,(k - 1)*71 + 2,SEEK_SET);
      if(5!=fscanf(menup,"%d| %[^|]%*c%d|%d|%ld",&k,txt,ndel,nrest,recpos)) return 0;
      *ncalc=0;
   }
   else
   {
      fseek(menuq,(k - 1)*77 + 2,SEEK_SET);
      if(6!=fscanf(menuq,"%4d| %[^|]%*c%d|%d|%d|%ld",&k,txt,ndel,ncalc,nrest,recpos))
           return 0;
   }
   return 1;
}

int whichArchive(int nFile,int rw)
{ static int ArchNum=0, rw_=0;
  char archivName[40];  
  if(nFile==0)
  { if(ArchNum!=0) {ArchNum=0; fclose(archiv);}
    return 0;
  }
  
  if(ArchNum==nFile &&  rw==rw_ )
  {  if(rw=='w' && ftell(archiv) >= MAXARCHIVE) 
     { fclose(archiv);
       ArchNum++;
       sprintf(archivName,ARCHIV_NAME,ArchNum);
       archiv=fopen(archivName,"w");
     } 
     return ArchNum;
  } 
  
  if(ArchNum) fclose(archiv); 
  ArchNum=nFile;  rw_=rw;
  sprintf(archivName,ARCHIV_NAME,ArchNum);
  if(rw=='w') 
  { if(access(archivName,F_OK)==0)  archiv=fopen(archivName,"a");
    else archiv=fopen(archivName,"w");
  }  
  else  archiv=fopen(archivName,"r");
   
  return ArchNum;
} 

