/*
 Copyright (C) 1997, Alexander Pukhov, e-mail pukhov@theory.npi.msu.su
*/
#include "syst.h"
#include "syst2.h"
#include "pvars.h"
#include "parser.h"

#include "pre_read.h"

preres pregarbage=NULL;

void  clearpregarbage(void)
{
   while (pregarbage)
   { preres  p = pregarbage->next;
       if(pregarbage->varsdeg) free(pregarbage->varsdeg);
       free(pregarbage);
       pregarbage = p;
   }

}


static void  rangecheck(preres q)
{int    i;

  if( q->degp > 255  || q->maxg > 20) {rderrcode = rangecheckerror; return;}
  for (i = 0; i < vardef->nvar; i++)
  if (q->varsdeg[i] > 127) {rderrcode = rangecheckerror; return; }
}


static void  newrecord(preres* q)
{  int i;
   *q = pregarbage;
   while (*q != NULL && !(*q)->free) *q = (*q)->next;
   if (*q == NULL)
   {
      *q = (preres)m_alloc(sizeof(struct preresrecord));
      (*q)->next = pregarbage;
      pregarbage = *q;
      (*q)->nvar=0;
      (*q)->varsdeg=NULL;
   }
   else for(i=0; i<(*q)->nvar;i++) (*q)->varsdeg[i]=0;

   (*q)->free = 0;
   (*q)->num = 1;
   (*q)->maxp = 0;
   (*q)->degp = 0;
   (*q)->g5 = 0;
   (*q)->maxg = 0;
   (*q)->indlist = 0;
}


void * rd_pre(char* s)
{
 int       i,k;
 preres     m;
 long    li;
 char rest[100];
 
   if(s[0]=='"')
   {  rderrcode=unexpectedcharacter;
      return NULL;
   }else if(isdigit(s[0]))
   {  
      if(1!=sscanf(s,"%ld%s",&li,rest)    ) rderrcode = toolargenumber; else
      {  newrecord(&m);
         m->num = labs(li);
         m->tp = numbertp;
      }
   }
   else
   {
      if (strlen(s) >= VAR_NAME_SIZE)  rderrcode = toolongidentifier;
      else
      {
         newrecord(&m);
         m->tp = polytp;

         if (strlen(s) == 2 && isdigit(s[1]) && s[1] != '0')
         {
            switch (s[0])
            {
               case 'p':
               case 'P':
                  m->tp = vectortp;
                  m->maxp =s[1] - '0';
                  m->degp = 1;
               break;

               case 'm':
                  m->tp = indextp;
                  m->indlist = 1<< (s[1] - '1');
                  break;
               case 'M':
                  m->tp = indextp;
                  m->indlist = 1<< (s[1] - '1' +4);                  
            } 

            if (strcmp(s,"G5") == 0)
            {
               m->tp = spintp;
               m->g5 = 1;
            }
         }

         if (m->tp == polytp)
         {
            i = 0;
            while (i < vardef->nvar && strcmp(s,vardef->vars[i].name)) i++;
            if (i == vardef->nvar)
            {  
               increaseVars(vardef);
               strcpy(vardef->vars[i].name,s);
            }
            m->nvar=ALIG(vardef->nvar);
            m->varsdeg=re_alloc(m->varsdeg,m->nvar*sizeof(unsigned));
            for(k=0; k< m->nvar;k++)  m->varsdeg[k] = 0;
            m->varsdeg[i] = 1;
         }
      }
   }
   return (void *)m;
}


static void * uact(char* ch,void * mm)
{preres   m;

   m = (preres)mm;
   if (strcmp(ch,"G") == 0 || strcmp(ch,"g") == 0)
   {
      if(m->tp!=vectortp && m->tp!=indextp)
      {
         rderrcode = typemismatch;
         return mm;
      }
      m->tp = spintp;
      m->maxg = 1;
   }
   else if (strcmp(ch,"-") == 0)
   {
      if (m->tp == indextp)
      {
         rderrcode = typemismatch;
         return mm;
      }
      m->num = -m->num;
   } else rderrcode=unexpectedoperation;
   return mm;
}


static void *  bact(char ch,void * mm1,void * mm2)
{  preres    m1, m2, m3;
  int i,d,k;
  
  m1 = (preres)mm1;
  m2 = (preres)mm2;

  if(m1->nvar < vardef->nvar)
  {  int newsize=ALIG(vardef->nvar);
     m1->varsdeg=re_alloc(m1->varsdeg,newsize *sizeof(unsigned)); 
     for (i=m1->nvar; i<newsize ;i++) m1->varsdeg[i]=0;
     m1->nvar=newsize;
  }
 
  if(m2->nvar < vardef->nvar)
  {  int newsize=ALIG(vardef->nvar);
     m2->varsdeg=re_alloc(m2->varsdeg, newsize*sizeof(unsigned)); 
     for (i=m2->nvar; i<newsize;i++) m2->varsdeg[i]=0;
     m2->nvar=newsize;
  }

     
   switch (ch)
   {
      case '+':
/*         if(!sum_perm[m1->tp][m2->tp])
         { rderrcode=typemismatch; return NULL; }
*/
         if(m1->indlist != m2->indlist)
         { rderrcode = indexuncompatibility; return NULL; }

         if(m1->tp < m2->tp) { m3 = m1; m1 = m2; m2 = m3; }

/*
         if (m1->tp == indextp ||
             (m1->tp == vectortp && m2->tp != vectortp) ||
             (m2->tp > polytp    && m1->tp != m2->tp))
         {
            rderrcode = typemismatch;
            return NULL;
         }
*/
         m1->num += m2->num;
         m1->maxp = MAX(m1->maxp,m2->maxp);
         m1->degp = MAX(m1->degp,m2->degp);
         m1->g5 = m1->g5 || m2->g5;
         m1->maxg = MAX(m1->maxg,m2->maxg);
         if (m1->tp == rationtp)
            for(i=0; i<vardef->nvar; i++) m1->varsdeg[i] += m2->varsdeg[i];
         else
            for(i=0; i<vardef->nvar; i++) m1->varsdeg[i] = 
                                     MAX(m1->varsdeg[i],m2->varsdeg[i]);
      break;



      case '*':
/*         if(!mult_perm[m1->tp][m2->tp])
         {  rderrcode = typemismatch; return NULL; }
*/
         if((m1->indlist & m2->indlist))
         {  rderrcode = indexuncompatibility; return NULL; }

         
         if (m1->tp < m2->tp) { m3 = m1; m1 = m2; m2 = m3; }
/*
         if (m1->tp == indextp || (m2->tp > polytp && m1->tp != m2->tp))
         {
            rderrcode = typemismatch;
            return NULL;
         }
*/
         if(m1->tp > m2->tp &&  m2->tp ==rationtp) 
         {  rderrcode = unexpectedoperation;
            return 0;
         }

         m1->num *= m2->num;
         m1->maxp = MAX(m1->maxp,m2->maxp);
         m1->degp += m2->degp;
         m1->g5 = m1->g5 || m2->g5;
         m1->maxg += m2->maxg;
         for (i = 0; i < vardef->nvar; i++) m1->varsdeg[i] += m2->varsdeg[i];
         m1->indlist |= m2->indlist;
      break;

      case '/':
         if(m2->tp > rationtp || m1->tp > rationtp || m1->tp == indextp)
         {
            rderrcode = typemismatch;
            return NULL;
         }

         m1->maxp = MAX(m1->maxp,m2->maxp);
         m1->degp += m2->degp;
         m1->g5 = m1->g5 || m2->g5;
         for(i=0; i<vardef->nvar; i++) m1->varsdeg[i] += m2->varsdeg[i];

         m1->tp = rationtp;
      break;

      case '^':
         if (m2->tp != numbertp )
         {
            rderrcode =  unexpectedoperation;
            return NULL;
         }
      
         d = m2->num;
         if( d <= 0 || d > 255)
         {
            rderrcode = rangecheckerror;
            return NULL;
         }

         if (m1->tp > rationtp)
         {
            rderrcode = unexpectedoperation;
            return NULL;
         }

         k = m1->num;
         for (i = 1; i <= d-1; i++) m1->num *= k;
         m1->degp *= d;
         for (i = 1; i <= vardef->nvar; i++) m1->varsdeg[i-1] *= d;
      break;

      case '.':
         if (m1->tp < vectortp || m2->tp < vectortp)
         {
            rderrcode = typemismatch;
            return NULL;
         }
         m1->tp = m1->tp == vectortp && m2->tp == vectortp ?
            polytp : tenstp;
         if (m1->indlist && m1->indlist == m2->indlist)
         {
            rderrcode = indexuncompatibility;
            return NULL;
         }
         else   m1->indlist |= m2->indlist;
         m1->num *= m2->num;
         m1->maxp = MAX(m1->maxp,m2->maxp);
         m1->degp += m2->degp;
         m1->g5 = m1->g5 || m2->g5;
         m1->maxg += m2->maxg;
         for (i = 0; i < vardef->nvar; i++) m1->varsdeg[i] += m2->varsdeg[i];
         break;
     default:
         rderrcode=unexpectedoperation;
         m1->free=1;
         m2->free=1;
         return NULL;
   } /*  Case  */

   m2->free = 1;

   rangecheck(m1);
   return (void *)m1;
}

void * act_pre(char * name, int n, void ** args)
{ void* args2[2];
  if(n==2 && !strcmp(name,"-")) strcpy(name,"+");
  if(n==1) return uact(name,args[0]);
  if(n==2) return bact(name[0],args[0],args[1]);
  if(n==4 && !strcmp(name,"eps") )
  { 
     if(!(args2[0]=bact('.',args[0],args[1]))) return NULL;
     if(!(args2[1]=bact('.',args[2],args[3]))) return NULL;
     return bact('*',args2[0],args2[1]);
  }
  rderrcode=unexpectedoperation;
  return NULL;  
}

void * act_preF(char *name,int n,void ** args)
{
  if (name[0]=='+'||(n==2 && name[0]=='-')||name[0]=='.' )
  {  rderrcode = unexpectedoperation;
	  return NULL;
  }
  return act_pre(name,n,args);
}
