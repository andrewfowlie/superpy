#include <dirent.h>
#include"chep_crt.h"

#include "viewdir.h"


static void  readtext(char* fname)
{
   FILE * txt;
   trim(fname);
   txt=fopen(fname,"r");
   if (txt == NULL)
   {
      messanykey(10,10," File not found");
	return;
   }
   showtext (1, 1, 80,1,"",txt);
   fclose(txt);
}



int viewDir(char * dirName)
{
  int   i,k;
  void *  pscr2 = NULL;
  char f_name[STRSIZ];
  char  menustr[2020];
  DIR *dirPtr=opendir(dirName);
  struct dirent * dp;
  menustr[0]=16;
  
  if(!dirPtr) return 1;

  k=1;
  
  while((dp=readdir(dirPtr)) && k<=2000)
  if(strcmp(dp->d_name,".")&&strcmp(dp->d_name,"..")&& strcmp(dp->d_name,"aux") )
  {
    for(i=0; (i <strlen(dp->d_name))&&(i<16);i++) menustr[k++]=dp->d_name[i];
    for(; i <16; i++) menustr[k++]=' ';
  }
  closedir(dirPtr);
  menustr[k]=0;
  if (menustr[1] == 0 )
  {
    messanykey(10,15,"directory  is empty");
    return 0;
  }     
        
  for(k=1;k;)
  {
    menu1(10,10,"",menustr,"",&pscr2,&k);
    if(k > 0)
    { sprintf(f_name,"%s/%.16s",dirName,menustr+k*16-15);
      readtext(f_name);                           
    }                                             
  }
  return 0;                                               
}                                                 


int checkDir(char * dirName)
{
  DIR *dirPtr=opendir(dirName);
  struct dirent * dp;
  int res=0;
  
  if(!dirPtr) 
  { char buff[200]; 
    sprintf(buff,"mkdir %s",dirName);
    res=system(buff);
    if(res) 
    { sprintf(buff,"Fatal error. Directory\n%s\n"
       "is absent and can not be created",dirName);
      sortie(99);
    }
    return 0;    
  }
  while((dp=readdir(dirPtr)))
  if(strcmp(dp->d_name,".")&&strcmp(dp->d_name,"..")&& strcmp(dp->d_name,"aux")  )
  {   closedir(dirPtr); return 1;
  }
  closedir(dirPtr); return 0;
}
