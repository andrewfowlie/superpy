#include<stdio.h>
#include <stdarg.h>
#include<string.h>

#include"syst.h"
#include"writeF.h"

#define xmax 76

FILE * outFile=NULL;

int    xpos=1;


static void  wrtln(void)
{
   if(!outFile) outFile=stdout;
   fprintf(outFile,"\n");
   xpos = 1;
}


void  wrt_0(char* s)
{
   
   int  l= strlen(s);
   if(!outFile) outFile=stdout;
   if(!l) return;
   if(xpos>xmax) 
   {  fprintf(outFile,"\n ");
      xpos=2;
      wrt_0(s);
      return;
   }

   if ((xpos -1 + l) <= xmax)
   {
      fprintf(outFile,"%s",s);
      xpos += l;
      return;
   }
  
   l = xmax - xpos + 1;
   
/*   if(strchr("*=",s[l]) && strchr("*-/+",s[l-1])) l--;  */

   while(l && (!strchr("*+-)(^=<> ",s[l-1]) || strchr("*+-)^=<>;",s[l])) ) l--;
     
      
   if(l==0)
   { 
      if (xpos>2) 
      {  fprintf(outFile,"\n "); 
         xpos=2; 
         wrt_0(s);
         return;
      } else
      {      
        l = xmax - xpos +1;
        if(strchr("*=",s[l]) && strchr("*-/+",s[l-1])) l++; 
        while(s[l] && strchr("*+-)(^= ",s[l]) == NULL) l++;               
      }
   }
      
   fprintf(outFile,"%.*s\n ",l,s); xpos=2;
   wrt_0(s+l);
}



void writeF(char * format,...)
{
  va_list args;
  char dump[2*STRSIZ] , *beg, *nn;

  va_start(args,format);   
  vsprintf(dump,format,args);
  va_end(args);


  beg=dump;   
  while(1)
  {  nn=strchr(beg,'\n');
     if(nn==NULL){(*wrt_0)(beg);return;}
     nn[0]=0;
     (*wrt_0)(beg);
     wrtln();
     beg=nn+1;
  }
}

void  outFileOpen(char * format, ...)
{
  va_list args;
  char dump[1024];

  va_start(args,format);   
  vsprintf(dump,format,args);
  va_end(args);
  outFile= fopen(dump,"w" );xpos=1;
  if(!outFile) outFile=stdout;
}

void outFileClose(void) { fclose(outFile); outFile=stdout; xpos=1; }
