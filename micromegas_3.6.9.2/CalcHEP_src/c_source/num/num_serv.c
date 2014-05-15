/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "chep_crt.h"
#include "plot.h"
#include "param.h"
#include "interface.h"
#include "err_code.h"

#include "num_serv.h"

void  paramdependence(r_func  ff, char*  procname, char*  resultname)
{
   double      minprm, maxprm;
   int         npoints;
   double      stepprm;
   unsigned    count;
   double f[201];
   int Esc=0,mPos=1;
   double prmval; 
   char txt[100];
   char name[20];
   REAL *vPos, memprm;   
   if(!selectParam(54,11,"Choose parameter",NULL,nin_int==2,0,1,0,&vPos,name,&mPos)) return; 

   memprm=*vPos;
 
   minprm = memprm; 
   maxprm = minprm; 
    
label1:
   sprintf(txt,"'%s' min=",name); 
   if (!correctDouble(55,14,txt,&minprm,0))  return;

label2:
   sprintf(txt,"'%s' max=",name);    
   if (!correctDouble(55,15,txt,&maxprm,0)) 
   { 
      goto_xy(55,14);
      clr_eol(); 
      goto label1;
   } 
   if (maxprm <= minprm) 
   { 
      messanykey(55,17,"Range check error"); 
      goto_xy(55,15);
      clr_eol();
      goto label2;
   } 

label4: npoints = 101;
   if (correctInt(55,16,"Number of points= ",&npoints,0)) 
   { 
      if (npoints < 3) 
      { 
         messanykey(55,17,"Too few points!"); 
         goto label4;
      } 
      if (npoints > 201) 
      { 
          messanykey(55,17,"Too many points!"); 
          goto label4;
      } 
   } 
   else
   { 
      goto_xy(55,15);
      clr_eol(); 
      goto label2;
   }

   goto_xy(55,14); clr_eol();
   goto_xy(55,15); clr_eol();
   goto_xy(55,16); clr_eol();

   stepprm = (maxprm - minprm)/(npoints - 1); 

   informline(0,npoints);

   stepprm = (maxprm - minprm) / (npoints - 1); 
   prmval=minprm;
   err_code = 0; 

   for(count = 1; count <= npoints; count++)
   { 
      *vPos=prmval;
      err_code=checkParam();
      if(err_code>1) break; 

      f[count-1] = ff();
      if(err_code>1) break;
      
      Esc=informline(count,npoints);
      if(Esc) break;       
      prmval += stepprm;
   } 
 
   if(err_code) errormessage();
   strcpy(txt,name);
   if (err_code <=1 && Esc==0)  plot_1(minprm,maxprm ,npoints,f,NULL,procname,txt,resultname); 
   *vPos=memprm; 
   checkParam();
} 
