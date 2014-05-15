/*
 Copyright (C) 1997,2003  Alexander Pukhov 
*/

#include "polynom.h"
#include "syst2.h"
#include "pvars.h"
#include "tensor.h"
#include "symb_wrt.h"
#include "writeF.h"


void  writepoly(poly p)
{  char txt[STRSIZ], numtxt[STRSIZ];
   int  i, deg;
   int  plus, first;
   int  wpos;
   unsigned long   wpower;

   if(!p){wrt_0("0"); return;}
   if(p->next) plus=0; else plus=1;
 

   for(;p;p=p->next)
   {
      strcpy(txt,"");
      i = 1;
      first = 1;
      wpos = 1;
      wpower = p->power[wpos-1];
      for (i = 1; i <= vardef->nvar; i++)
      {
         deg = (p->power[vardef->vars[i-1].wordpos-1] /
                vardef->vars[i-1].zerodeg) %
                vardef->vars[i-1].maxdeg;

         if(deg > 0)
         {
            if(first)  first = 0; else        strcat(txt,"*");
            strcat(txt,vardef->vars[i-1].name);
            if(deg > 1)  sprintf(txt+strlen(txt),"^%d",deg);
         }
      }

      sprintf(numtxt,"%"NUM_STR, p->num);
       
      if(plus && numtxt[0]!='-') wrt_0("+"); 

      if(strlen(txt))
      {  if(strcmp(numtxt,"1"))
         {
            if(strcmp(numtxt,"-1")==0) wrt_0("-");
            else {wrt_0(numtxt);wrt_0("*");}
         }    
         wrt_0(txt); 
      }
      else wrt_0(numtxt);
      plus=1;
   }
}  


void  writetens(tensor p)
{  char  txt[STRSIZ];
   int   i,s;
   int   plus, first;
   int   isMono;

   if(!p) {wrt_0("0"); return;}

   if(p->next) plus=0; else plus=1;    

   for(; p; p=p->next)
   {
      strcpy(txt,"");
      first = 1;
      for (i = 1; i <= maxIndex; i++)
      {
         s = p->tens[i-1];
         if(s < i && s)
         {
            if(first) first = 0; else strcat(txt,"*"); 
            if(s > 0) sprintf(txt+strlen(txt),"m%d.m%d",i,s); 
            else      sprintf(txt+strlen(txt),"p%d.m%d",-s,i); 
         } 
      }
      
      isMono=0;
      if(p->re && p->re->next==NULL && p->im==NULL) isMono=1;
      if(p->im && p->im->next==NULL && p->re==NULL) isMono=1;

      if(!isMono){if(plus) wrt_0("+("); else wrt_0("(");}
      plus=1;

      if(p->re) writepoly(p->re); 
      if(p->im)
      { if(p->im->next) 
        {  wrt_0("+i*(");
           writepoly(p->im);
           wrt_0(")");
        } else {writepoly(p->im);wrt_0("*i");}
      }      

      if(!isMono) wrt_0(")");
       
      if(strlen(txt)) { wrt_0("*"); wrt_0(txt);} 
   } 
} 


void  writeEtens(Etens p)
{  char  txt[20]; 
   int   i,l; 
   int plus=0;
   if(!p) {wrt_0("0"); return;} 
   for(;p;p=p->next) 
   { 
      if((p->tcoef->next||((p->tcoef->re&&p->tcoef->im))) && p->eps[0]!=X_MARK)
      {  if(plus) wrt_0("+(");else wrt_0("(");   
         writetens(p->tcoef);   wrt_0(")");
      }
      else   writetens(p->tcoef);
      plus=1;
   
      if(p->eps[0]!=X_MARK)
      {  sprintf(txt,"*eps(");
         for(i=0;i<4;i++)
         { int  c=p->eps[i];
           l=strlen(txt);
           if(c>0)sprintf(txt+l,"m%d,",c);else sprintf(txt+l,"p%d,",-c);
         }
         txt[strlen(txt)-1]=')';
         wrt_0(txt);
      }     
   }
}

void  writespinor(SpinTensor p)
{  char  txt[200]; 
   int   i,l; 
   int plus=0;
   if(!p) {wrt_0("0"); return;} 
   for(;p;p=p->next) 
   { 
      if(p->tcoef->next && (p->g5 || p->l))
      { { if(plus) wrt_0("+(");else wrt_0("(");}  writetens(p->tcoef); wrt_0(")");}
      else   writetens(p->tcoef);
         
      if(p->g5 || p->l)
      { 
         sprintf(txt,"*g(ln");
         if(p->g5) strcat(txt,",a");
         for(i=0;i<p->l;i++)
         { int  c=p->g[i];
           l=strlen(txt);
           if(c>0)sprintf(txt+l,",m%d",c);else sprintf(txt+l,",p%d",-c);
         }
         strcat(txt,")");
         wrt_0(txt);
      }     
   }
}

void  symb_print(char* txt, symb_data  m)
{  char s[5];
   int n;
   tensor p;

   wrt_0(txt);
   switch(m.type)  
   { case numbertp:
     case polytp : writepoly(m.expr.p);  break;
     case tenstp : writetens(m.expr.t);  break;
     case spintp : writespinor(m.expr.s);break;
     case etenstp: writeEtens(m.expr.et);break;
     default:
        if(m.type == vectortp)				       
        {								       
           p=m.expr.t;
           if(p == NULL) wrt_0("0");					       
           else							       
           for(;p;p=p->next)							       
           {  n=-p->tens[0];
	      sprintf(s,"p%d*(",n);
              wrt_0(s);						       
              if(p->re) writepoly(p->re);	       
              if(p->im)					       
              { wrt_0("i*("); writepoly(p->re);wrt_0(")");}	       
              wrt_0(")");						       
              if(p->next) wrt_0("+");				       
           }   				       
        }								       
        else								       
        {	
           n = m.expr.t->tens[0];	 
           if(!n) n = 1;						      
           sprintf(s,"l%d",n);
           wrt_0(s);
           writetens(m.expr.t);
           
        }							       
   }
   writeF(";\n");
}
