/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#include <dirent.h>
        
#include "chep_crt.h"
#include "files.h"


 char * outputDir = "";

 char  pathtocalchep[STRSIZ];


char * unixPath(char *name)
{  static struct passwd *passwdPrt;
   char * cslash, *beg, *res;
   int l;
   trim(name);
   if(name[0]!='$' && name[0] != '~') return name;
   
   cslash=strchr(name,'/');
   if(cslash) cslash[0]=0;
   
   beg=NULL;
   if(name[0]=='$') beg=getenv(name+1);
   if(strcmp(name,"~")==0) beg=getenv("HOME");
   if(!beg)
   { passwdPrt=getpwnam(name+1);  
     if(passwdPrt) beg=passwdPrt->pw_dir;
   }
   if(cslash) cslash[0]='/'; else cslash="";      
   if(beg)
   { res=malloc(strlen(beg)+strlen(cslash)+2);
     sprintf(res,"%s%s",beg,cslash);
     return res;
   }else return NULL;  
}


void copyfile(char* namefrom,char* nameto)
{ 

  FILE * filefrom;
  FILE * fileto;
  char  s [STRSIZ];
  filefrom= fopen(namefrom,"r");
  fileto  = fopen(nameto,"w");
  if ( (filefrom==NULL) || (fileto==NULL) ) return;

  for(;;)
  if (fgets(s,STRSIZ,filefrom) != NULL) f_printf(fileto,s) ; else
  {
     fclose(fileto);
     fclose(filefrom);
     return;
  }
}

void nextFileName(char* f_name,char * firstname,char * ext)
{  int tabnum;
   for(tabnum=1;;tabnum++)
   {  sprintf(f_name,"%s%s%d%s",outputDir,firstname,tabnum,ext);
      if(access(f_name,F_OK)) return;
   } 
}

#define W 50
static char* view_Dir(char * dirName)
{
  int   i,k,k1;
  void *  pscr2 = NULL;
  static char f_name[W+1];
  char  *menustr=NULL,*menustr1=NULL;
  DIR *dirPtr=opendir(dirName);
  struct dirent * dp;
  menustr=malloc(10);
  menustr1=malloc(10);
  menustr[0]=W; 
  
  if(!dirPtr) { messanykey(10,10,"Wrong directory name");return NULL;}

  k=1;
  k1=0;

  while((dp=readdir(dirPtr)))
  if(strcmp(dp->d_name,".") && strcmp(dp->d_name,"..") )
  { int err;
    char * longname=malloc(strlen(dirName)+strlen(dp->d_name)+2);
    struct stat buf;
    sprintf(longname,"%s/%s",dirName,dp->d_name); 
    if(stat(longname,&buf)==0)
    {
      if(S_ISDIR(buf.st_mode) )
      {  menustr=realloc(menustr,k+W+5);
         for(i=0; (i <strlen(dp->d_name))&&(i<W-1);i++) menustr[k++]=dp->d_name[i];
	 menustr[k++]='/'; 
         for(; i <W-1; i++) menustr[k++]=' ';
      } else
      {  menustr1=realloc(menustr1,k1+W+5);
         for(i=0; (i <strlen(dp->d_name))&&(i<W);i++) menustr1[k1++]=dp->d_name[i];
         for(; i <W; i++) menustr1[k1++]=' ';
      }
    }
    free(longname);           
  }
  closedir(dirPtr);
  menustr[k]=0; menustr1[k1]=0;
  menustr=realloc(menustr, 1+strlen(menustr)+strlen(menustr1));
  strcat(menustr,menustr1);
  free(menustr1);
  
  if (menustr[1] == 0 )
  {
    messanykey(10,15,"directory is empty");
    return NULL;
  }

  for(k=1;k;)
  {
    menu1(3,4,"Contents of directory",menustr,"",NULL,&k);
    if(k > 0)
    { strncpy(f_name,menustr+(k-1)*W+1,W);  
      f_name[W]=0;
      trim(f_name); 
      { free(menustr);  return f_name;}
    }else { free(menustr);  return NULL;}                                           
  }
}                                                 

int findCalcHEPfile(char * name) 
{ int key;
  static char fName[200]="$CALCHEP/*"; 
  char *fName_=NULL;
  void * pscr;
  int X=3,Y=2,xx;
  
  get_text(X,Y,maxCol(),Y,&pscr);
  scrcolor(White,Black);
  goto_xy(X,Y); print("File Search:");   
  xx=where_x()+1;
  for(;;)
  {  char *file;
     goto_xy(xx,Y);
     key  = str_redact(fName,strlen(fName),200-W);     
     if(key==KB_ESC) { put_text(&pscr);return 0;}
     if(key==KB_F1) show_help("findFile");
     if(key==KB_ENTER)
     { char *buff; 
       fName_=unixPath(fName);
       if(!fName_)
       { messanykey(10,17, "Can not interpret beginning of path");  continue; }
       if(strlen(fName)>1 && strcmp(fName+strlen(fName)-2,"/*")) break;
       buff=malloc(1+strlen(fName_));
       strcpy(buff, fName_);
       buff[strlen(buff)-2]=0;
       file=view_Dir(buff);
       free(buff);
       if(file)
       {
	  sprintf(fName+strlen(fName)-1,"%s",file);
	  if(fName[strlen(fName)-1]=='/') strcat(fName,"*");
       }
       if(fName_!=fName) free(fName_);
     }       
  }
  strcpy(name,fName_);
  if(fName_!=fName) free(fName_); 
  put_text(&pscr);

  return 1;
}
#undef W
