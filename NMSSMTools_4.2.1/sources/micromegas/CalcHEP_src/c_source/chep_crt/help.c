#include "syst.h"
#include "crt_util.h"
#include "file_scr.h"

#include "help.h"

char pathtohelp[256]="";

int show_help(char * fname)
{
   char f_name[STRSIZ];
   FILE* f;
   int z[4];
   sprintf(f_name,"%s%s%s",pathtohelp,fname,".txt");

   f=fopen(f_name,"r"); 
   if (f == NULL) return 0;
   fgets(f_name,STRSIZ,f);
   sscanf(f_name,"%d %d %d",&z[0],&z[1],&z[2]);
   z[2]=z[2]+z[0]+1;

   showtext (z[0],z[1],z[2],z[1],/*" Help "*/fname,f);
   fclose(f);
   return 1;
}

