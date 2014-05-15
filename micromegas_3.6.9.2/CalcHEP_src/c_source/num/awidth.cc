#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
              
#include "../../include/num_out.h"
#include "../dynamicME/dynamic_cs.h"

double aWidth_stat(char * pName)
{  int dim;
   txtList LL;
printf("label 0 aWidth_stat\n");   
   if(!compDir)
   { int size=100;
     char *buff=NULL,*err;
     for(;;)
     {  compDir=realloc(compDir,size); 
        if(getcwd(compDir,size)) break; else size*=2;
     }
     libDir=malloc(strlen(compDir)+20);
     sprintf(libDir,"%s/so_generated",compDir); 
     modelNum=1;
     calchepDir=getenv("CALCHEP");
     if(!calchepDir) calchepDir=interface_ext.CALCHEP;
     ForceUG=interface_ext.forceUG;
   }

printf(" label 1 aWidth_stat\n");

   { double r=  pWidth(pName,&LL,&dim);
     printf("r=%E\n",r);
     return r;
   }  
}

double (*aWidth)(char *)=aWidth_stat;