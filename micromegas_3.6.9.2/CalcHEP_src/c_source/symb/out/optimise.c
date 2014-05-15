/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include <limits.h>

#include "syst2.h"
#include "s_files.h"
#include "getmem.h"
#include "physics.h"
#include "out_serv.h"
#include "pvars.h"
#include "procvar.h"

#include "optimise.h"

 infoptr info;
 int    firstVar;

static infoptr infoone;

static int  int_strlen(int * s)
{
  int i=0;
  while( s[i]) i++;
  return i;
}

static void int_strcpy(int * to  , int * from)
{
  int i=0;
  
  while( from[i]) { to[i]=from[i];i++;}
  to[i]=0;
} 

static int  int_strcmp(int * s1, int * s2)
{
   int i=0;
   while( s1[i]&&s2[i]&&s1[i]==s2[i] ) i++;
   return s1[i]-s2[i];
}


static int* int_strchr(int * str, int s)
{
  int i=0;
  while( str[i] && str[i]!=s) i++;
  if( str[i]==s) return (str+i); else return NULL;
}


void  initinfo(void)
{
   info = (infoptr)getmem_(sizeof(struct inforec));
   info->next = NULL;
   strcpy(info->name,"1");
   info->ival = 1;
   info->consttype = numb;
   infoone=info;
}

int      equalexpr(varptr v1,varptr v2)
{
   while (v1 != NULL || v2 != NULL)
   {
      if (v1 == NULL || v2 == NULL) return 0;
      if (v1->sgn == v2->sgn && v1->coef == v2->coef &&
      int_strcmp(v1->vars,v2->vars) == 0)
      {  v1 = v1->next;
         v2 = v2->next;
      }
      else return 0;
   }
   return 1;
}

static void  readmonom(int * varstr,int * conststr,NUM_TYPE * numc)
{  
   int iv=0,ic=0;
   int       deg,n,k,pos;

   fread(&readBuff->num,readSize,1,archiv);
   *numc = readBuff->num;
   if(!*numc ) return ;

   for (n = 0; n < vardef->nvar; n++)
   {
       deg = (readBuff->power[vardef->vars[n].wordpos-1] /
       vardef->vars[n].zerodeg) %
       vardef->vars[n].maxdeg;
 
       pos=vardef->vars[n].num;
       for(k=0;k<deg;k++)
       { 
         if(pos >= firstVar)  varstr[iv++] = pos; else conststr[ic++] = pos;
       }


   }

   varstr[iv] = 0;
   conststr[ic] = 0;

}   /* ReadMonom */

static void addnum(NUM_TYPE n,char* signum,infoptr* ans)
{
   if (n>0) *signum='+'; else { *signum='-';n= -n;}
   *ans=info;
   while (*ans != NULL)
   {  if ( ((*ans)->consttype == numb) && ((*ans)->ival == n) ) return;
      *ans= (*ans)->next;
   }
   *ans=(infoptr) getmem_( sizeof(struct inforec));
   (*ans)->next=info;
   (*ans)->consttype=numb;
   (*ans)->ival=n;
   if(ABS(n) >= LONG_MAX) sprintf((*ans)->name,"%"NUM_STR".",n);
        else              sprintf((*ans)->name,"%"NUM_STR   ,n);
   info = *ans;
}

static int addtmpconst(varptr tmpconst,int* s,infoptr* coeff )
{
   varptr     c;
   int count;
 
   revers((void**)&tmpconst);
   if (tmpconst->sgn =='-')
   {  *s='-';
      for(c = tmpconst;c;c = c->next)
      { if (c->sgn=='-' ) c->sgn='+'; else c->sgn='-'; }
   }
   else  *s='+';

   if( (tmpconst->next == NULL ) && (tmpconst->vars[0] == '\0') )
   {  *coeff = tmpconst->coef; return  0; }

   *coeff = info;
   count=0;
   for( *coeff = info,count=0; *coeff; *coeff = (*coeff)->next)
   if ((*coeff)->consttype == expr)
   { count++;
     if(equalexpr((*coeff)->const_,tmpconst))return 0; 
   } 
   
   *coeff = (infoptr)getmem_(sizeof(struct inforec));
   (*coeff)->next = info;
   (*coeff)->const_ = tmpconst;
   (*coeff)->consttype = expr;
   sprintf((*coeff)->name,"C[%d]",count);

   info = *coeff;
   return 1;
} /*  AddTmpConst */


void  readpolynom(varptr* expr_)
{      
   int  varstr[STRSIZ], conststr[STRSIZ];
   NUM_TYPE      n;
   void*      pntr;
   varptr       tmpconst;
   char         s;
   infoptr      coeff;
   marktp  tmpmark;

   readmonom(varstr,conststr,&n);
   if(!n) { *expr_ = NULL;  return;}

   *expr_ = (varptr)getmem_(minvarrec + sizeof(int)*int_strlen(varstr));
   (*expr_)->next = NULL;
   int_strcpy((*expr_)->vars,varstr);

	addnum(n,&s,&coeff);
	mark_(&tmpmark);
   tmpconst = (varptr)getmem_(minvarrec + sizeof(int)*int_strlen(conststr));
   tmpconst->next = NULL;
   int_strcpy(tmpconst->vars,conststr);
   tmpconst->sgn=s;
   tmpconst->coef=coeff;

   while(1)
   {
      readmonom(varstr,conststr,&n);
      if(!n) break;
      if (int_strcmp( varstr,(*expr_)->vars) != 0)
      {  if (! addtmpconst(tmpconst,&((*expr_)->sgn),&((*expr_)->coef) ))
                                   release_(&tmpmark); 
         pntr = (void*)(*expr_);
         *expr_ = (varptr)getmem_(minvarrec + sizeof(int)*int_strlen(varstr));
         (*expr_)->next = (varptr)pntr;
         int_strcpy((*expr_)->vars,varstr);
         pntr = NULL;
         addnum(n,&s,&coeff);
         mark_(&tmpmark);
      }
      else
      {
          pntr = (void*)tmpconst;
          addnum(n,&s,&coeff);
      }
      tmpconst = (varptr)getmem_(minvarrec + sizeof(int)*int_strlen(conststr));
      tmpconst->next = (varptr)pntr;
      int_strcpy(tmpconst->vars,conststr);
      tmpconst->sgn = s;
      tmpconst->coef =coeff;
   }
	if (! addtmpconst(tmpconst,&((*expr_)->sgn),&((*expr_)->coef)) )
	  release_(&tmpmark);
}   /* ReadPolynom */


static void  findmaxvar(varptr ex,unsigned* n,int* ch,int * power)
{ int   *   nterms;
  int   *   minpower;
  int       k, bt, d, nv;

  nterms=  (int*) m_alloc(sizeof(int)*nProcessVar);
  minpower=(int*) m_alloc(sizeof(int)*nProcessVar);
  
   for (k = firstVar; k < nProcessVar; k++)
   {  nterms[k] = 0;
      minpower[k] = 0;
   }

   while (ex != NULL)
   {
      if (int_strlen(ex->vars))
      {
         d = 1;
         bt = ex->vars[0];
         for (k = 1; k <int_strlen(ex->vars); k++)
         {
            nv = ex->vars[k];
            if (nv != bt)
            {
               minpower[bt] = minpower[bt] ?  MIN(d,minpower[bt]):d;
               nterms[bt]++;
               d = 1;
               bt = nv;
            }
            else ++(d);
         }
         minpower[bt] = minpower[bt] ?  MIN(d,minpower[bt]):d; 
         nterms[bt]++;
      }
      ex = ex->next;
   }

   bt = firstVar;
   *n = nterms[firstVar];
   *power = minpower[firstVar];
   
   for (k = firstVar+1; k <nProcessVar; k++)
   if (*n < nterms[k] || (*n == nterms[k] && *power <  minpower[k]))
   {
      *n =  nterms[k];
      bt = k;
      *power = minpower[k];
   }
   *ch =  bt;
   free(nterms); free(minpower);
}


static void  findmaxcoef(varptr ex,infoptr* i,unsigned* n)
{
 varptr jj;
   jj=ex;  while (jj != NULL) { (jj->coef)->count=0;jj=jj->next;}
   jj=ex;  while (jj != NULL) { (jj->coef)->count++;jj=jj->next;}

   *n=0;
   infoone=info;
   while (infoone->next != NULL ) infoone=infoone->next;
   (*i)=infoone;

   (*i)->count=0;
   jj=ex;
   while(jj != NULL)
   {  if (*n < (jj->coef)->count)

      {
         *i =  jj->coef;
         *n =  (*i)->count;
      }
      jj = jj->next;
   }
}

static void  clipvar(varptr ex,int ch,int power,varptr* ex1,varptr* ex2)
{ 
 varptr       exnext;

 var_rec     ex1rec;
 var_rec     ex2rec;
 varptr      ex1_,ex2_;
 int * u;
 
  ex1_ =  & ex1rec;
  ex2_ =  & ex2rec;
   while (ex != NULL)
   {
      u = int_strchr(ex->vars,ch);
      exnext = ex->next;
      if (u)
      { int i=0;
        do u[i]=u[i+power]; while(u[i++]); 
        ex1_->next = ex;
        ex1_ = ex;
      }
      else
      {
         ex2_->next = ex;
         ex2_ = ex;
      }
      ex = exnext;
   }


   ex1_->next = NULL;  *ex1=ex1rec.next;
   ex2_->next = NULL;  *ex2=ex2rec.next;
}


static void  clipconst(varptr ex,infoptr i,varptr* ex1,varptr* ex2)
{

 varptr       exnext;
 infoptr      one;

   var_rec     ex1rec, ex2rec;
   varptr      ex1_,ex2_;

        one = info;
        while (one->next != NULL) one = one->next;

   ex1_ =  & ex1rec;
   ex2_ =  & ex2rec;

   while (ex != NULL)
   {
      exnext = ex->next;
      if (ex->coef != i)
      {
         ex2_->next = ex;
         ex2_ = ex;
      }
      else
      {
         ex->coef = one;
         ex1_->next = ex;
         ex1_ = ex;
      }
      ex = exnext;
   }


   ex1_->next = NULL;  *ex1=ex1rec.next;
   ex2_->next = NULL;  *ex2=ex2rec.next;

}

void*      emitexpr(varptr ex,smplemit smplemitfun,vfact vfactfun,
                               cfact cfactfun)
{unsigned     nv, nc;
 int        ch;
 infoptr      i;
 varptr       ex1, ex2;
 void         *pmult, *psum;
 int          deg;

   findmaxcoef(ex,&i,&nc);
   findmaxvar(ex,&nv,&ch,&deg);
   if (nc < 2 && nv < 2) return smplemitfun(ex);
   if (nv >= nc)
   {
      clipvar(ex,ch,deg,&ex1,&ex2);
      pmult = emitexpr(ex1,smplemitfun,vfactfun,cfactfun);
      psum  = emitexpr(ex2,smplemitfun,vfactfun,cfactfun);
      return  vfactfun(ch,deg,pmult,psum);
   }
   else
   {
      clipconst(ex,i,&ex1,&ex2);
      pmult = emitexpr(ex1,smplemitfun,vfactfun,cfactfun);
      psum = emitexpr(ex2,smplemitfun,vfactfun,cfactfun);
      return cfactfun(i,pmult,psum);
   }
}

