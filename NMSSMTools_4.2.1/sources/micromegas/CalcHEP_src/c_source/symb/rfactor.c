/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "syst2.h"
#include "physics.h"
#include "parser.h"
#include "ghosts.h"
#include "prepdiag.h"

#include "rfactor.h"


void  copysmonom(s_monom src,s_monom* dst)
{vmrec     rec;
 vmptr     vs, vd;

   dst->c = src.c;
   vs = src.v;
   vd = &rec;
   while (vs != NULL)
   {
      vd->next = (vmptr) m_alloc(sizeof(struct vmrec));
      vd = vd->next;
      strcpy(vd->name,vs->name); 
      vd->deg = vs->deg;
      vs = vs->next;
   }
   vd->next = NULL;
   dst->v = rec.next;
}


static void  reducec(NUM_TYPE* l1,NUM_TYPE * l2)
{ NUM_TYPE    c, i1, i2;

   i1 = ABS(*l1); i2 = ABS(*l2);
   if (i2 > i1) { c = i1; i1 = i2; i2 = c; }
   while (i2 != 0) { c = i2; i2 = REST(i1,i2); i1 = c; }
   (*l1) = DIV(*l1,i1);
   (*l2) = DIV(*l2,i1);

}


void  clrvm(vmptr c)
{vmptr  p;

   while (c != NULL)
   {
      p = c;
      c = c->next;
      free(p);
   }
}


static void  comparenames(char* s1,char* s2,int* eq,int* gt)
{int    i, l1, l2, n1, n2;

   l1 =  strlen(s1); 
   l2 =  strlen(s2); 
   *eq = 0; 
   *gt = 0; 
   if (l1 > l2) 
   {  *gt = 1;
      return;
   } 
   else 
      if (l1 < l2) 
         return;
      else 
         for (i = 1; i <= l1; i++) 
      { 
         n1 = s1[i-1];
         n2 = s2[i-1];
         if (n1 > n2) 
         {  *gt = 1; 
            return;
         } 
         else 
            if (n1 < n2) return;
      }
   *eq = 1; 
} 


void  sew_vm(vmptr* p1,vmptr p2,int mlt)
{vmptr    m, mm, m1, m2; 
 int  gt, eq;
 vmrec    mrec; 

   if (p2 == NULL) return;
   if (*p1 == NULL) { *p1 = p2; return; } 
   m1 = *p1; 
   mrec.next = *p1; 
   m = &mrec;
   *p1 = m; 
   m2 = p2; 

label_1: 
   comparenames(m1->name,m2->name,&eq,&gt); 

label_2: 
   while (gt) 
   { 
      m = m1;
      m1 = m1->next; 
      if (m1 == NULL) 
      {  m->next = m2; 
         goto label_3;
      } 
      comparenames(m1->name,m2->name,&eq,&gt); 
   }
   if (eq) 
   { 
      if (mlt) 
         m1->deg += m2->deg; 
      else 
         m1->deg = MAX(m1->deg,m2->deg); 
      mm = m2; 
      m2 = m2->next; 
      free(mm); 
      if (m2 == NULL) goto label_3;
      goto label_1;
   } 
   mm = m1; 
   m1 = m2; 
   m2 = mm; 
   m->next = m1; 
   gt = 1; 
   goto label_2;

label_3: 
   *p1 = (*p1)->next; 
} 


static void  reducev(vmptr* p1,vmptr* p2)
{vmptr       m, mm1, mm2, m1, m2; 
 int     gt, eq; 
 vmrec       mrec1, mrec2; 
 unsigned        d;

   if (*p1 == NULL || *p2 == NULL) return;
   m1 = *p1; 
   mrec1.next = *p1; 
   mm1 = &mrec1; 
   *p1 = mm1; 

   m2 = *p2; 
   mrec2.next = *p2; 
   mm2 = &mrec2; 
   *p2 = mm2; 
   while (1) 
   { 
      comparenames(m1->name,m2->name,&eq,&gt); 
      if (eq) 
      { 
         d = MIN(m1->deg,m2->deg); 
         m1->deg -= d; 
         m2->deg -= d; 
         if (m1->deg == 0) 
         { 
            m = m1; 
            mm1->next = m1->next; 
            free(m); 
         } 
         else 
            mm1 = m1; 
         if (m2->deg == 0) 
         { 
            m = m2; 
            mm2->next = m2->next; 
            free(m); 
         } 
         else 
            mm2 = m2; 
         m1 = mm1->next; 
         if (m1 == NULL) goto exi;
         m2 = mm2->next; 
         if (m2 == NULL) goto exi;
      } 
      else 
      { 
         if (gt)
         {
            mm1 = m1;
            m1 = m1->next;
            if (m1 == NULL) goto exi;
         }
         else
         {
            mm2 = m2;
            m2 = m2->next;
            if (m2 == NULL) goto exi;
         }
      }
   }

exi:
   *p1 = (*p1)->next;
   *p2 = (*p2)->next;
}


static void  delsqrt2(s_monom* s)
{vmrec       rec;
 vmptr       m, m1;

   rec.next = s->v;
   m = &rec;
   m1 = s->v;
   while (m1 != NULL)
   if (strcmp(m1->name,"Sqrt2") == 0 && m1->deg != 1)
   {
      while (m1->deg > 1)
      {
			s->c *= 2;
         m1->deg -= 2;
      }
      if (m1->deg == 0)
      {
         m->next = m1->next;
         free(m1);
         m1 = m->next;
      }
   }
   else
   {
      m = m1;
      m1 = m1->next;
   }
   s->v = rec.next;
}


static void  del_i(s_monom* s)
{vmrec       rec;
 vmptr       m, m1;

   rec.next = s->v;
   m = &rec;
   m1 = s->v;
   while (m1 != NULL)
	if (strcmp(m1->name,"i") == 0 && m1->deg != 1)
   {
      while (m1->deg > 1)
      {
			s->c = -s->c;
         m1->deg -= 2;
      }
      if (m1->deg == 0)
      {
         m->next = m1->next;
         free(m1);
         m1 = m->next;
      }
   }
   else
   {
      m = m1;
      m1 = m1->next;
   }
   s->v = rec.next;
}







void  reduce_s(s_monom* s1,s_monom* s2)
{
   reducec(&(s1->c),&(s2->c));
   reducev(&s1->v,&s2->v);
}


void  mult_s(s_monom* s1,s_monom* s2)
{
   s1->c *= s2->c;
   sew_vm(&s1->v,s2->v,1);
	delsqrt2(s1);
	del_i(s1);
}

void  mult_rptr(rmptr* m1,rmptr* m2)
{
   reduce_s(&(*m1)->n,&(*m2)->d);
   reduce_s(&(*m1)->d,&(*m2)->n);

   mult_s(&(*m1)->n,&(*m2)->n);
   mult_s(&(*m1)->d,&(*m2)->d);

   free(*m2);
   delsqrt2(&(*m1)->n);
   del_i(&(*m1)->n);

   reducec(&(*m1)->n.c, &(*m1)->d.c);

   if ((*m1)->d.c < 0)
   { (*m1)->d.c = - (*m1)->d.c;
     (*m1)->d.c = - (*m1)->d.c;
   }
}


static void  revol(rmptr r)
{s_monom     s;
 vmrec       rec;
 vmptr       m, m1;
 int i;

	s = r->n;
	r->n = r->d;
	r->d = s;

	rec.next = (r->d).v;
   m = &rec;
	m1 = m->next;
	while (m1 != NULL)
	{
		if (strcmp(m1->name,"i") == 0 )
		{
			m->next = m1->next;
			if ( (m1->deg & 1) == 0)  s.c= 1 ; else s.c= -1;
			s.v=m1;
			m1->next=NULL;
			mult_s(&(r->n),&s);
			m1 = m->next;
		} else
		if (strcmp(m1->name,"Sqrt2") == 0 )
		{
			m->next = m1->next;
		        for(i=0; i<m1->deg; i++) (r->d).c = (r->d).c*2;
			s.c=1;
			s.v=m1;
			m1->next=NULL;
			mult_s(&(r->n),&s);
			m1 = m->next;
		} else
		{
			m = m1;
			m1 = m1->next;
		}
	}
	(r->d).v = rec.next;

        reducec(&r->n.c,& r->d.c);
        
	if (r->d.c < 0)
	{ r->d.c = - r->d.c;
	  r->d.c = - r->d.c;
	}

}


static void *  rd_rat(char* s)
{ 
  rmptr       m;
  NUM_TYPE    li;

   if (isdigit(s[0]))
   {
      sscanf(s,"%"STR_NUM,&li);

      m = (rmptr) m_alloc(sizeof(struct r_monom));
      m->n.c = li;
      m->d.c = NUM_ONE;
      m->n.v = NULL;
      m->d.v = NULL;
   }
   else
   {
//      if (strlen(s) > 6)
//         rderrcode = toolongidentifier;
 //     else
      {
         m = (rmptr) m_alloc(sizeof(struct r_monom));
         m->n.c = 1;
         m->d.c = 1;
         m->d.v = NULL;
         m->n.v = (vmptr) m_alloc(sizeof(struct vmrec));
         strcpy(m->n.v->name,s);
         m->n.v->next = NULL;
         m->n.v->deg = 1;
      }
   }
   return (void *) m;
}


static void *  uact_r(char* ch,void * mm)
{rmptr      m;

   if (strcmp(ch,"-") == 0)
   {
      m = (rmptr) mm;
      m->n.c = -m->n.c;
   }
   else
      rderrcode = unexpectedoperation;
   return mm;
}


static void *  bact_r(char ch,void * mm1,void * mm2)
{rmptr     m1, m2;
 int    i;
 NUM_TYPE   ln, ld,d;
 vmptr     p;

   m1 = (rmptr) mm1;
   m2 = (rmptr) mm2;
   
   switch (ch)
   {
      case '+':
      case '.':
         rderrcode = unexpectedoperation;
      break;

      case '*':
         mult_rptr(&m1,&m2);
      break;

      case '/':
         revol(m2);
         mult_rptr(&m1,&m2);
      break;

      case '^':
         if (m2->n.v == NULL && m2->d.v == NULL && m2->d.c == 1)
         {
            d = m2->n.c;
            free(m2);
            if (d < 0)
            {
               revol(m1);
               d = -d;
            }
            ln = 1;
            ld = 1;
            for (i = 1; i <= d; i++)
            {
               ln *= m1->n.c;
               ld *= m1->d.c;
            }
            m1->n.c = ln;
            m1->d.c = ld;

            p = m1->n.v;
            while (p != NULL)
            {
               p->deg *= d;
               p = p->next;
            }
	    delsqrt2(&m1->n);
	    del_i(&m1->n);
            p = m1->d.v;
            while (p != NULL)
            {
               p->deg *= d;
               p = p->next;
            }
            reducec(& m1->n.c ,& m1->d.c);
         }
      break;
   }  /*  Case  */
   
   
   return (void *) m1;
}

static void * act_rat(char * name,int n, void ** args)
{
  if(n==1) return uact_r(name,args[0]);
  if(n==2 && name[0]=='-') {rderrcode=unexpectedoperation; return NULL;}
  return bact_r(name[0],args[0],args[1]); 
}

void *  read_rmonom(char* txt)
{ return  readExpression(txt,rd_rat, act_rat, NULL);}

char *  smonomtxt(s_monom s)
{vmptr        p;
 int      first;
 static char  ss[STRSIZ];

   first = 1;
   if (s.c != 1) 
   { 
      sprintf(ss,"%"NUM_STR,s.c); 

      first = 0; 
   } 
   else  strcpy(ss,"");
    
   p = s.v; 
   while (p != NULL) 
   { 
      if (first)  first = 0; else  strcat(ss,"*"); 
      sprintf(ss+strlen(ss),p->name); 
      if (p->deg > 1) sprintf(ss+strlen(ss),"^%d",p->deg); 
      p = p->next; 
   } 


   if(strcmp(ss,"") == 0 )  strcpy(ss,"1");
   return ss;    
} 


char  * rmonomtxt(r_monom r)
{static char  snum[STRSIZ]; 
 char         sden[STRSIZ]; 

   strcpy(snum,smonomtxt(r.n));
   strcpy(sden,smonomtxt(r.d));
   if (strcmp(sden,"1") != 0) sprintf(snum+strlen(snum),"/(%s)",sden);
   return snum;
}


void  diagramsrfactors(hlpcsptr gst,s_listptr* s,rmptr* totf)
{s_listptr    s1, s2, ss, dl;
 vcsect       vcs_copy;
 int          i;
 rmptr        rcoef, rrcoef;
 rmptr        r;
 s_monom      stmp, stmp2; 
 int      first; 


   s1 = NULL; 
   s2 = NULL; 

   first = 1; 
   vcs_copy = vcs; 
   while (gst != NULL) 
   { 
      ss = (s_listptr) m_alloc(sizeof(struct s_listrec));
      ss->next = s1; s1 = ss; 
      ss = (s_listptr) m_alloc(sizeof(struct s_listrec));
      ss->next = s2; s2 = ss; 
      coloringvcs(gst); 
      attachvertexes(); 
      rcoef = (rmptr) readExpression(vertexes[0].lgrnptr->comcoef,
                                 rd_rat,act_rat,NULL); 
      for (i = 2; i <= vcs.sizet; i++) 
      { 
         rrcoef = (rmptr) readExpression(vertexes[i-1].lgrnptr->comcoef,
                                     rd_rat,act_rat,NULL); 
         mult_rptr(&rcoef,&rrcoef); 
      } 
      if (first) 
      { 
         r = rcoef; 
         s1->monom.c = 1; 
         s1->monom.v = NULL; 
         s2->monom.c = 1; 
         s2->monom.v = NULL; 
         first = 0; 
      } 
      else 
      { 
         copysmonom(r->n,&s1->monom); 
         s2->monom = rcoef->n; 
         reduce_s(&s1->monom,&s2->monom); 
         copysmonom(s1->monom,&stmp); 
         reduce_s(&r->n,&stmp); 
         mult_s(&r->n,&stmp);   /*  for case Stmp = -1  */ 
         copysmonom(r->d,&stmp); 
         reduce_s(&stmp,&rcoef->d); 
         copysmonom(rcoef->d,&stmp2); 
         mult_s(&s1->monom,&rcoef->d); 
         mult_s(&s2->monom,&stmp); 
         mult_s(&r->d,&stmp2); 
         free(rcoef); 
      } 
      vcs = vcs_copy; 
      gst = gst->next; 
   } 

   vcs = vcs_copy; 
   ss = s2; 
   stmp.v = NULL; 
   stmp.c = 1; 
   while (s2 != NULL) 
   { 
      copysmonom(stmp,&stmp2); 
      mult_s(&s2->monom,&stmp2); 
      mult_s(&stmp,&s1->monom); 
      dl = s1; 
      s1 = s1->next; 
      free(dl); 
      s2 = s2->next; 
   } 
   clrvm(stmp.v); 
   revers((void **)&ss); 
   *s = ss; 
   *totf = r; 
   
} 


void  eraseslist(s_listptr s)
{s_listptr   sdel; 

   while (s != NULL) 
   { 
      clrvm(s->monom.v); 
      sdel = s; 
      s = s->next; 
      free(sdel);
   } 
} 
