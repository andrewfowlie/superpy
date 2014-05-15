#define _FILE_OFFSET_BITS 64 
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <errno.h>
#include "event2pyth.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


//extern FILE *fopen64(const char *filename, const char *type);

static FILE * eF=NULL;
static double cs;

static void fName2c(char*f_name,char*c_name,int len)
{ int i; for(i=len-1;i>=0 &&f_name[i]==' ';i--);
  c_name[i+1]=0;
  for(;i>=0;i--) c_name[i]=f_name[i];
}

int openeventfile_(char *fname, int len)
{
  char * cname;
 
  if(len<0) {eF=stdin; return 0;}
  cname=malloc(len+1);
  fName2c(fname,cname,len);
/*  channel=open(cname,O_RDONLY);*/
 
  printf("cname=|%s|\n",cname);
  eF=fopen(cname,"r");
  printf("eF=%p\n",eF);
/*  printf("errno=%d\n",errno); */ 
  free(cname);
  if(eF) return 0; else return 1;  
}


void closeeventfile_(void){ if(eF){ fclose(eF); eF=NULL; }}

int readeventheader_(void)
{ char buff[200];
  if(!eF) return -1;
  fseek(eF,SEEK_SET,0);
  for(;;) 
  { int n=fscanf(eF,"%s",buff);
    if(n!=1) return -2;
    if(strcmp(buff,"<init>")==0)
    {
       n=fscanf(eF," %d %d %lf %lf %d %d %d %d %d %d %lf %lf %lf %d",
         R_.IDBMUP,R_.IDBMUP+1,R_.EBMUP, R_.EBMUP+1,R_.PDFGUP,R_.PDFGUP+1,
         R_.PDFSUP,R_.PDFSUP+1,&R_.IDWTUP,&R_.NPRUP,R_.XSECUP,R_.XERRUP,R_.XMAXUP, 
         R_.LPRUP);
      if(n==14)  return 0; else return 2; 
    }
  }
  return 0;
}      

int  readevent_(void)
{ char buff[200]; 
  if(!eF) return -1;
  for(;;) 
  { int n=fscanf(eF,"%s",buff);
    if(n!=1) return -2;
    if(strcmp(buff,"<event>")==0)
    { int I,J;
      n=fscanf(eF,"%d %d %lf %lf %lf %lf\n",
          &E_.NUP,&E_.IDPRUP,&E_.XWGTUP,&E_.SCALUP,&E_.AQEDUP,&E_.AQCDUP);
      if(n!=6) return 1;     
      for(I=0;I<E_.NUP;I++)
      {
        n+=fscanf(eF," %8d %4d",  E_.IDUP+I,E_.ISTUP+I);
        for(J=0;J<2;J++)  n+=fscanf(eF," %d", E_.MOTHUP[I]+J);
        for(J=0;J<2;J++)  n+=fscanf(eF," %d", E_.ICOLUP[I]+J);
        for(J=0;J<5;J++)  n+=fscanf(eF," %lf",E_.PUP[I]+J);
        n+=fscanf(eF," %lf %lf\n",E_.VTIMUP+I,E_.SPINUP+I);
      }
      if(n==6+13*E_.NUP) return 0; else return 2;
    }
  }
  return 0;  
}
