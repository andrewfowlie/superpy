#include"micromegas_aux.h"

void fName2c(char*f_name,char*c_name,int len)
{ int i; for(i=len-1;i>=0 &&f_name[i]==' ';i--);
  c_name[i+1]=0;
  for(;i>=0;i--) c_name[i]=f_name[i];
}

void cName2f(char*c_name,char*f_name,int len)
{ int i; for(i=0;i<len &&c_name[i];i++) f_name[i]=c_name[i];
         for(   ;i<len            ;i++) f_name[i]=' ';
}
