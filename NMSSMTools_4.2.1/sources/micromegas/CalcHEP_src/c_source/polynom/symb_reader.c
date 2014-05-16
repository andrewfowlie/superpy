/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <ctype.h>
#include "syst2.h"
#include "pvars.h"
#include "polynom.h"
#include "tensor.h"
#include "spinor.h"
#include "parser.h"
#include "sos.h"
#include "reader0.h"
#include "symb_reader.h"
#include "symb_tot.h"

/*#define STRACE*/
#ifdef STRACE
#include "symb_wrt.h"
#endif



symb_all rd_symb(char*s)
{  
   symb_all p;
   long L=0;   
   if(s[0]=='p' || s[0]=='m')
   { int i;
     for(i=1;isdigit(s[i]);i++);
     if(s[i]==0 && i>1) sscanf(s+1,"%ld",&L);
   }
   
   NewUnit(p);
   if(isdigit(s[0]))
   { sscanf(s,"%ld",&L);
     p->type=numbertp;
     p->expr.p=plusone();
     multpolyint(&(p->expr.p),L);
   } 
   else if(!strcmp(s,"G5")) {p->type=spintp;p->expr.s=spin1(); p->expr.s->g5=1;}
   else if(!strcmp(s,"i")) 
   { p->type=tenstp;
     p->expr.t=newtensor1();  
     p->expr.t->im= p->expr.t->re; 
     p->expr.t->re=NULL;
   } 
   else if(s[0]=='p' && L>0)
   { 
     p->type=vectortp;
     p->expr.t=newtensor1();
     p->expr.t->tens[0] = -L;
   }                                    
   else if(s[0]=='m' &&  L>0 )     
   {      
     p->type=indextp;
     p->expr.t=newtensor1();
     if(L!=1){ p->expr.t->tens[0]=L; p->expr.t->tens[L-1]=1;}
   }    
   else
   { int i;
     for(i=0;i<vardef->nvar;i++) if(strcmp(vardef->vars[i].name,s)==0)
     { p->expr.p=plusone();
       p->expr.p->power[vardef->vars[i].wordpos-1] =vardef->vars[i].zerodeg;
       p->type=polytp;
       break;
     } 
     if(i==vardef->nvar)
     { p->expr.p=NULL;
       p->type=errortp;
       printf("can not find |%s|\n",s); 
     }  
   }
#ifdef STRACE
	printf(" \n rd_ (%s) ->  (type=%d) ",s,p->type);
	 symb_print("?",*p);
	printf(";\n"); 
#endif
   return p;
}


static symb_all uact_(char * ch, symb_all mm)
{ 

#ifdef STRACE
	printf("\n uact_ (%s)\n",ch);
	 symb_print("?",*mm);
	printf(" -> \n\n\n");
#endif
   if(strcmp(ch,"-")==0) { *mm=symb_imult(*mm,1,-1);
#ifdef STRACE
         symb_print("?",*mm);
#endif   
    return mm;}
   else if( (strcmp(ch,"g")==0 || strcmp(ch,"G")==0)
          &&(mm->type==vectortp|| mm->type==indextp) )  
   { SpinTensor sum=NULL;
     tensor t=mm->expr.t;
     for(;t;t=t->next)
     {  SpinTensor s=spin1();
        int np=t->tens[0];
        s->l=1;s->g[0]=np?np:1;
        delpoly(&(s->tcoef->re)); s->tcoef->re=copypoly(t->re);
        delpoly(&(s->tcoef->im)); s->tcoef->im=copypoly(t->im);
        addSpin(&sum,s);
     }
     mm->type=spintp;
     deltensor(&(mm->expr.t)); 
      mm->expr.s=sum; 
   } else 
   { mm->type=errortp;
     mm->expr.p=NULL;
   }  
   return mm;
}


static symb_all  bact_(char ch,symb_all mm1,symb_all mm2)
{
 symb_all mm3;
 NewUnit(mm3); mm3->type= errortp; mm3->expr.p=NULL;
 
#ifdef STRACE
	printf("\n bact_ (%d)(%c)(%d)\n",mm1->type,ch,mm2->type);
	 symb_print("?",*mm1);printf(" \n|%c|\n ",ch);
	 symb_print("?",*mm2);printf("\n\n ->\n\n ");
#endif

   switch (ch)
   {
      case '/': printf("unexpected division\n"); sortie(60);
      case '+': *mm3=symb_sum(*mm1,0,*mm2,0); break; 
      case '*': *mm3=symb_mult(*mm1,0,*mm2,0); break;
      case '^':
         { int i, d = mm2->expr.p->num;
           mm3->type=numbertp;
           mm3->expr.p=plusone(); 
           for(i=0; i < d; i++) *mm3=symb_mult(*mm3,1,*mm1,0);
         } break;

      case '.':
         mm3->expr.t=multtwotens(mm1->expr.t,mm2->expr.t);
         mm3->type=tenstp;       
   }
#ifdef STRACE
	 symb_print("?",*mm3);
#endif
   return mm3;
}

static symb_all act_symb(char * name, int n, symb_all * args)
{  
  if(n==1) return uact_(name,args[0]);
  if(name[0]=='-') { uact_(name,args[0]); strcpy(name,"+");} 
  if(n==2) return bact_(name[0],args[0],args[1]);
  if(n==4 && !strcmp(name,"eps"))
  {  tensor parg[4];
     int i;
     Etens ans=NULL;
     
     for(i=0;i<4;i++) parg[i]=args[i]->expr.t;
     
     for(;parg[0];parg[0]=parg[0]->next)
     for(;parg[1];parg[1]=parg[1]->next)
     for(;parg[2];parg[2]=parg[2]->next)
     for(;parg[3];parg[3]=parg[3]->next)
     { int l[4], sgn;
       Etens eps=etens1();
       for(i=0;i<4;i++) 
       {  tensor t;
          void * Next=parg[i]->next;
          parg[i]->next=NULL; t=copytens(parg[i]); parg[i]->next=Next;
          l[i]=t->tens[0];
          if(l[i]==0) l[i]=1;
          t->tens[0]=0;
          if(l[i]>0) t->tens[l[i]-1]=0;
          multEtensTens(&eps,t);
       }
       
       for(sgn=1, i=0; i<3 && sgn; )
       {
                if (l[i+1] < l[i]) i++;
          else  if (l[i+1] > l[i])
          { int c=l[i];
             l[i]=l[i+1]; l[i+1]=c;
             sgn*=-1;
             if (i) i--;  else  i++;
          } else  sgn = 0;
       }
       
       for(i=0;i<4;i++) eps->eps[i]=l[i];
       multEtensInt(&eps,sgn);
       addEtens(&ans,eps);
       { symb_all p;
         NewUnit(p) 
         p->type=etenstp;
         p->expr.et=ans;
         return p;
       } 
     }      
  }
  return NULL;
}

symb_data symb_read(char * s)
{  symb_data ret, *m;
   
   m=(symb_all)readExpression(s,(rdelement)rd_symb,(operation)act_symb,NULL);
   if(m) { ret=*m; delunit(m);} else {m->type=errortp; m->expr.p=NULL;}
   return ret;
}
