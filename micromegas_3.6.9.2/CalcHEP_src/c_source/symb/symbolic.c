/*
 Copyright (C) 1997,2006, Alexander Pukhov
*/

#include "chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "ghosts.h"
#include "cweight.h"
#include "prepdiag.h"
#include "s_files.h"
#include "pvars.h"
#include "chess.h"
#include "polynom.h"
#include "tensor.h"
#include "spinor.h"
#include "saveres.h"
#include "parser.h"
#include "pre_read.h"
#include "reader0.h" 
#include "symb_reader.h"
#include "rfactor.h"
#include "getmem.h"
#include "process.h"
#include "symbolic.h"
#include "sos.h"
#include "sets.h"
#include "symb_tot.h"

/*#define STRACE */
#ifdef STRACE
#include "writeF.h"
#include "symb_wrt.h"
#endif

static int  gamma_map[2 * maxvert];
static int  maxmomdeg;
static int  calcdiag_sq;

static symb_data  vert[2*MAXINOUT];
static symb_data  block[6*MAXINOUT],fermres[maxvert]; 

static symb_data rnum; 
static poly factn, factd;


static poly  polyfactor(poly p)
{ int i;
  poly m;
  poly q;
  NUM_TYPE c,num1,num2,numS;
 
  NewUnit(m);
  m->next=NULL;
  if(!p) for(i=0, m->num=1; i<monomLength; i++) m->power[i] = 0;
  else
  {
     num1=p->num;
     if(num1<0) {num1=-num1; numS=-1;} else  numS=1;
               
     q=p->next;
     while(q)
     {
        num2=q->num >0 ? q->num: -q->num;
        if (num2 > num1)  { c = num1; num1 = num2; num2 = c; }
        while (num2 != 0) { c = num2; num2 = REST(num1,num2); num1 = c; }
        q=q->next;
     }
     for(i=0; i<monomLength; i++) m->power[i]=0;

     for (i=0; i<vardef->nvar; i++)
     if( vardef->vars[i].num<=nmodelvar)
     {  int wrd=vardef->vars[i].wordpos-1; 
        unsigned long z_d=vardef->vars[i].zerodeg;
        unsigned long m_d=vardef->vars[i].maxdeg;
        int deg = (p->power[wrd]/z_d) %m_d;
         
        q=p->next;
        while(q && deg)
        { int deg2= (q->power[wrd]/z_d) %m_d; 
          if(deg2<deg) deg=deg2;      
          q=q->next;  
        }
        if(deg)  m->power[wrd]+=deg*z_d;
     }
     num1 *=numS;
     while(p)
     {
        for(i=0; i<monomLength; i++) p->power[i] -= m->power[i];
        p->num=DIV(p->num,num1);
        p=p->next;
     }
     m->num=num1;      
  }
  return m;             
} 


static void memoryInfo_(int used)
{ 
  goto_xy(14,16);  
  print("%d Kb    ",used >> 10);
  if (escpressed()) save_sos(-2);
}

      
static void  wrtoperat(char* s)
{
  scrcolor(Blue,BGmain);
  goto_xy(14,17); clr_eol();
  print(s);
}


static void  firstvertexreading(void)
{  preres m, m_; 
   int    v, n; 

   clearVars(vardef);
 
   for (v = 0; v < vcs.sizet; v++)
   {  void * agrs[2];
      m_ = (preres) readExpression(vertexes[v].lgrnptr->description,
                          rd_pre,act_pre,NULL);
                          
      if (rderrcode) save_sos(rderrcode);
      if (m_->g5) fermloops[fermmap[v]-1].g5 = 1;
      gamma_map[v] = m_->maxg;
      m_->g5 = 0;
      m_->maxg = 0;
      m_->indlist = 0;
      m_->tp = polytp;
      if(v){agrs[0]=m; agrs[1]=m_; m=(preres)act_pre("*",2,agrs);} else m=m_;
      if (rderrcode) save_sos(rderrcode);
   }

   for (n=0;n<vardef->nvar;n++) vardef->vars[n].maxdeg=m->varsdeg[n]+1;
   maxmomdeg = m->degp;
}

static int   findMaxIndex(void)
{  int  v, l, vln, maxI;

   for(maxI = 0 ,v = 0; v < vcs.sizet; v++)
   {
      vln = vcs.valence[v];
      for(l=0;l<vln;l++)  maxI=MAX(maxI,vcs.vertlist[v][l].lorentz);
   }
   if(!maxI) maxI=1;

#ifdef STRACE
printf("maxI=%d\n",maxI);
#endif

   return maxI;
}


static int  findSpinLength(void)
{  int  n, currentlen, nl, sLen;

   sLen = 0;
   for (nl = 0; nl < nloop; nl++)
   {  fermloopstp * with1 = &fermloops[nl];
      currentlen = 0;
      for (n = 0; n < with1->len; n++)
      {
         currentlen += gamma_map[with1->vv[n]-1] + with1->spin[n];
         if (sLen < currentlen) sLen = currentlen;
      }
   }
#ifdef STRACE
   printf("Max Spin Length=%d\n", sLen);
#endif      
   return sLen;

/*   spinLength =1+ (spinLength +1)/sizeof(long);*/
}



static void  propagatorsfirstreading(void)
{  int  ninout, v, d, n, l;
   char *name;

   ninout=nin + nout;
   for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++)
   {  edgeinvert *L= &vcs.vertlist[v][l];
      if(L->moment>0)
      {
         d=prtclbase1[L->partcl].spin;
         if(d==2 && !set_in(L->lorentz,setmassindex)
                 && !(PLR_PRTCL&L->prop)) d=0;
	 if(d==4 && tolower(prtclbase1[L->partcl].hlp)=='t') d=2;
         maxmomdeg +=d;
         if(d && L->moment>ninout)
         {
            name=prtclbase1[L->partcl].massidnt;
            if(strcmp(name,"0") != 0)
            {
               for(n=0; n<ninout && strcmp(name,inoutmasses[n]);n++);
               if(n>=ninout) addvar(name,d);
               addvar(name,d);
            }
         }
      }
   } 
}


static void  addinoutmasses(void)
{  int  n, m;
   char * massname;

   for(n=0; n<nin + nout; n++)
   {
      massname=inoutmasses[n];
      if(strcmp(massname,"0"))
      {
         for(m=0;strcmp(massname,inoutmasses[m]);m++);
         if(n==m) addvar(massname,maxmomdeg);
      }
   }
}

/*1==========================================================*/
static void  addscmult( int nSpin )
{  int n, m, i,j, px,py;
   symb_data ee;
   char     name[8];

   for(i=0,px=0,py=0;i<3*maxvert; i++)
   {
      if (momdep[i][0]==1) px++;  /* For momdep this is a length. V.E. */
      if (momdep[i][0] > 0) py++;
   }

   maxmomdeg /= 2;
   for(n=1; n<px; n++) for(m=0; m<n; m++)
   { sprintf(name,"p%d.p%d",m+1,n+1); addvar(name,maxmomdeg);}

   addvar("Helicity1",1);  
   addvar("Helicity2",1);
   addvar("HelicityN1",1);  
   addvar("HelicityN2",1);
   
   sortvar();
   symb_start(vardef->nvar, vardef->vars, nSpin, findMaxIndex(),py);

   for(n=1; n<=nin+nout; n++)
   if(strcmp(inoutmasses[n-1],"0")==0) assignsclmult(-n,-n,NULL);
   else
   {
      ee=symb_read(inoutmasses[n-1]);
      assignsclmult(-n,-n,multtwopoly(ee.expr.p,ee.expr.p));
      delpoly(&(ee.expr.p)); 
   }

   for(n=2; n<=px; n++) for(m=1; m<=n-1; m++)
   {  symb_all mm;
      sprintf(name,"p%d.p%d",m,n);
      mm=rd_symb(name);
      assignsclmult(-m,-n,mm->expr.p);
      delunit(mm);      
   }

   for(n=px; n<py; n++) for(m=0; m<=n; m++) if((n!=m || n>= nin+nout))
   {  poly qq=NULL;
      for(i=1; i<=momdep[n][0]; i++)
      { int p1=momdep[n][i];
        for(j=1; j<=momdep[m][0]; j++)
        {  int  p2=momdep[m][j];
           poly pp=copypoly(scalarmult(-abs(p1),-abs(p2)));
           if(p1*p2<0) multpolyint(&pp,-1);
           sewpoly(&qq,&pp);
        }
      }
      assignsclmult(-m-1,-n-1,qq);
   }
}
/*2==================================================*/

static void  secondVertexReading(void)
{  
   int v;

#ifdef STRACE
	writeF("\n vertex reading \n");
	writeF("tensLength=%d\n", tensLength);
#endif

   for(v=0; v<vcs.sizet; v++)
   {  char *pstr;
     pstr=parseVertex(v,0);
#ifdef STRACE
printf("Parced Vert=%s\n",pstr);
#endif
      vert[v]=symb_read(pstr);
#ifdef STRACE 
   writeF("\n VERTEX %d  %s\n type=%d \n",v,vertexes[v].lgrnptr->description,vert[v].type); 
   symb_print("=>",vert[v]);
#endif
      free(pstr);
   }

}

/*3============================================================*/
static symb_data  masscalc(int v,int l,int deg)
{  char  *mass=prtclbase1[vcs.vertlist[v][l].partcl].massidnt;
   char masstxt[20];
   
   if(strcmp(mass,"0") == 0) return symb_read("0");
   if(deg!=1) {sprintf(masstxt,"%s^%d",mass,deg); mass=masstxt;}
   return symb_read(mass);
}


static symb_data  kmorgcalc(int v,int l)
{  int k=vcs.vertlist[v][l].moment, m=vcs.vertlist[v][l].lorentz;
   char txt[20];

   if(k>0) sprintf(txt,"p%d.m%d",k,m); else sprintf(txt,"-p%d.m%d",-k,m);
   return symb_read(txt);
}


static symb_data  kmsubcalc(int v,int l)
{  int     ll,k,m;
   symb_data p=symb_read("0"),q; 
   m = vcs.vertlist[v][l].lorentz;
  
   for(ll=0; ll<vcs.valence[v]; ll++) if(l!=ll)
   {  char txt[20];
      
      k=vcs.vertlist[v][ll].moment;
      if(k>0) sprintf(txt,"-p%d.m%d",k,m); else sprintf(txt,"p%d.m%d",-k,m);
      q=symb_read(txt);
      p=symb_sum(p,1,q,1);
   }
   return p;
}

static symb_data fermpropag(int v,int l,int spin)
{ symb_data  m;
  char *proptxt;

  if(spin==1) proptxt=fermPropagTxt(v,l,0);
  else        proptxt=spin3_2_propagator(v,l,0);   
  m=symb_read(proptxt);
  free(proptxt);
  return m;
}

/*4=============================================================*/
static symb_data tPropagator(int v,int l)
{ 
  char*txt=tPropagator6mass4(v,l,NULL);
  symb_data tmp=symb_read(txt);
  free(txt);
  return tmp;
}

static void  calcFermLoops(void)
{  int n,lpcount;

   symb_data frmprpg,fctmp1,fctmp2,fctmp3;
 
   for(lpcount=0; lpcount<nloop; lpcount++) 
   {  fermloopstp *Lp = &fermloops[lpcount];

      vert[Lp->vv[0]-1]=symb_typeUp(vert[Lp->vv[0]-1],1,spintp);

      frmprpg=fermpropag(Lp->vv[0]-1,Lp->ll[0]-1,Lp->spin[0]);
      fctmp1=symb_mult(vert[Lp->vv[0]-1],1,frmprpg,1);
      for(n=1; n<Lp->len; n++)
      {  symb_data V;
         int i,l,v,s;
         V=vert[Lp->vv[n]-1];
         for(i=0;i<Lp->nint[n];i++)
         { edgeinvert*L;
           l=Lp->intln[n][i]-1;
           v=Lp->vv[n]-1;
           L=&(vcs.vertlist[v][l]);
           s=prtclbase1[L->partcl].spin;
           if(s==2 && set_in(L->lorentz,setmassindex))
           {  symb_data Vk=symb_copy(V), m2=masscalc(v,l,2), km=kmsubcalc(v,l);
              V=symb_mult(V,1,m2,0); 
              Vk=symb_mult(Vk,1,km,1);
              km=kmsubcalc(L->link.vno, L->link.edno);
              Vk=symb_mult(Vk,1,km,1); 
              V=symb_sum(V,1,Vk,1);
           } else if(s==4 && tolower(prtclbase1[L->partcl].hlp)!='t')
                   V=symb_mult(V,1,tPropagator(v,l),1);
         }

         frmprpg=fermpropag(Lp->vv[n]-1,Lp->ll[n]-1,Lp->spin[n]);
/*symb_print("NewPropag:",frmprpg);  */       
         fctmp2 = symb_mult(V,1,frmprpg,1);
/*
symb_print("first mult",fctmp2);
printf("address:%p %p\n", fctmp1.expr.s, fctmp2.expr.s);
symb_print("first coeff",fctmp1);          
*/
         fctmp3 = symb_mult(fctmp1,1,fctmp2,1);
/*         
symb_print("second mult",fctmp3);     
*/    
         fctmp1=fctmp3;
      }
#ifdef STRACE
   symb_print("Fermion loop: ",fctmp1);
#endif
      fermres[lpcount]=symb_spur(fctmp1,1);
#ifdef STRACE
   symb_print("Its trace: ",fermres[lpcount]);
#endif
         
   }
}


static void  multblocks(int v1,int v2,int saveEps)
{  
   symb_data  t1, t2, sum, q, mult_1, mult_2,sub[15],mass2[15];
   int      i, v, l, vmin, nmassind, n;
   set      msind, ind1, ind2;
   int      lastn;
   char     messtxt[80];

   t1 = block[v1-1];
   t2 = block[v2-1];

   ind1=vertinfo[v1-1].ind;
   ind2=vertinfo[v2-1].ind;

   strcpy(messtxt,"Indices contraction   ");
   for (i = 1; i <= 15; i++) if (set_in(i,ind1) && set_in(i,ind2)) 
           sprintf(messtxt+strlen(messtxt)," L%d",i);
   strcat(messtxt,"           ");
   wrtoperat(messtxt);

   msind=set_and(ind1,ind2);
   msind=set_and(msind,setmassindex);

   if(set_eq0(msind)) sum=symb_mult(t1,1,t2,1);
   else
   { 
      nmassind = 0; 
      for(i=0; i<15; i++)  if(set_in(i+1,msind)) 
      { 
         v = massindpos[i].vrt1; 
         l = massindpos[i].ln1;
         sub[nmassind]=kmorgcalc(v-1,l-1);
         mass2[nmassind]=masscalc(v-1,l-1,2);
         nmassind++;
      } 
      sum.type=polytp; sum.expr.p=NULL;
      lastn = (1 << nmassind) - 1; 
      
      for(n=0;n<=lastn; n++) 
      { 
         if(n==lastn) { mult_1 = t1; mult_2 = t2; } 
         else  { mult_1 = symb_copy(t1); mult_2 = symb_copy(t2);}
         for(i=0; i<nmassind; i++)  if(((1 << i) & n))
         { 
            mult_1=symb_mult(mult_1,1,sub[i],0);
            mult_2=symb_mult(mult_2,1,sub[i],0);
            mult_2=symb_imult(mult_2,0,-1);
         }  else  mult_1= symb_mult(mult_1,1,mass2[i],0);       
         q=symb_mult(mult_2,0,mult_1,0);
         sum=symb_sum(sum,1,q,1);
      } 
      for (i = 0; i<nmassind; i++) {symb_clean(sub[i]); symb_clean(mass2[i]);} 
   }
 
   vmin = MIN(v1,v2);
   block[vmin-1] = sum;
   vertinfo[vmin-1].g5 = vertinfo[v1-1].g5 + vertinfo[v2-1].g5; 
   vertinfo[vmin-1].ind=set_aun(set_or(ind1,ind2),set_and(ind1,ind2));
} 


static void  del_pp(poly* p, poly * fact, long * del)
{
   int      d1, n, i, j;
   unsigned long  z_d, m_d;
   poly     psub, fact1, ans,q,*pp_powers;
   long dmax;   

   if(*p == NULL) { *fact=plusone(); *del=1; return; }
   n = nin + nout;
   psub = scalarmult(-n,-n);
   multpolyint(&psub,-1);
   for(i=1; i<n; i++) { q=scalarmult(-i,-i); sewpoly(&psub,&q);}

   for(i=2; i<n; i++) for(j=1; j<i; j++) if(j != n - 2)
   {
      q = scalarmult(-i,-j);
      if(i>nin && j<=nin)  multpolyint(&q,-2); else multpolyint(&q,2);
      sewpoly(&psub,&q);
   }

   if (nout != 2) multpolyint(&psub,-1);

   z_d = vardef->vars[0].zerodeg;
   m_d = vardef->vars[0].maxdeg;
   fact1=polyfactor(*p);
   dmax=((*p)->power[0] / z_d) % m_d;

   pp_powers=m_alloc((dmax+1)*sizeof(poly));
   pp_powers[0]=plusone();
   for(i=1;i<=dmax;i++) pp_powers[i]=multtwopoly(pp_powers[i-1],psub);


   ans=NULL;
   q=*p;
   while(q)
   {  poly mon=q;
      poly tmp; 
      q=q->next;

      mon->next=NULL;
      d1 = (mon->power[0] / z_d) % m_d;
      mon->power[0] -= d1 * z_d;
      tmp=multtwopoly(pp_powers[d1], mon);
      multpolyint(&tmp, 1<< (dmax-d1));
      sewpoly(&ans, &tmp);     
      delpoly(&mon);
   }
   q=polyfactor(ans);
   (*fact)=multtwopoly(fact1,q);
   delpoly(&fact1);
    delpoly(&q);
   *del= 1<<dmax;   
   (*p)=ans;
   for(i=0;i<=dmax;i++) delpoly(pp_powers +i); 
   free(pp_powers);
}


static void  transformfactor(rmptr* t_fact,poly mon,long del)
{  int     i, deg;
   rmptr   vard;
   int     v, l, n, s, pnum;
   rmptr   tf_add;
   int     c;
   long    factnum, factdenum;
   char    factortxt[STRSIZ];
   preres  m;
   symb_data mm;

   sprintf(factortxt,"%"NUM_STR,mon->num);

   for(i=0; i<vardef->nvar; i++)
   {
      deg=(mon->power[vardef->vars[i].wordpos-1] /
             vardef->vars[i].zerodeg) %  vardef->vars[i].maxdeg;
      if(deg)
      {
         sprintf(factortxt+strlen(factortxt),"*%s",vardef->vars[i].name);
         if (deg > 1) sprintf(factortxt+strlen(factortxt),"^%d",deg);
      }
   }


   vard = (rmptr) read_rmonom(factortxt);
   
   mult_rptr(t_fact,&vard);
/* ------- Symmetry and Color  Factors ------- */
   factnum = vcs.symnum * vcs.clrnum;
   factdenum = vcs.symdenum * vcs.clrdenum * del;

/* -----  average factor  --------- */
   for(c=1,v=0; v<vcs.sizel; v++) for(l=0; l<vcs.valence[v]; l++)
   if(IN_PRTCL & vcs.vertlist[v][l].prop)
   {
      pnum=vcs.vertlist[v][l].partcl;
      s=prtclbase1[pnum].spin;
      switch (s)
      {
        case 1: if(!strchr("LR",prtclbase1[pnum].hlp)) c*=2; break;
        case 2: if (zeromass(pnum))c*=2;else c*=3; break;
        case 3: if (zeromass(pnum))c*=2;else c*=4; break;
        case 4: if(zeromass(pnum))  c*=2;else c*=5;
      }
      c*= abs(prtclbase1[pnum].cdim);
   }
   factdenum *= c;
   
/* ----- Fermion factor  --------- */

   for(c=1,v=0; v<vcs.sizel; v++) for(l=0; l<vcs.valence[v]; l++)
   if(IN_PRTCL & vcs.vertlist[v][l].prop
   &&(prtclbase1[vcs.vertlist[v][l].partcl].spin&1)) c*=-1;
   for(v=0; v<nloop; v++) c*=-4;
   factnum*=c;
   
/* ----- Vector/left spinor  factor  --------- */
   for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++)
   if(vcs.vertlist[v][l].moment > 0
   && strchr("LR",prtclbase1[vcs.vertlist[v][l].partcl].hlp)) factdenum *= 2;


/* ----------- end of numeric factors collecting--------------- */
   sprintf(factortxt,"%d",factnum);

   tf_add = (rmptr)read_rmonom(factortxt);
   mult_rptr(t_fact,&tf_add);
   sprintf(factortxt,"1/(%d",factdenum);

   for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++)
   {  edgeinvert *L = &(vcs.vertlist[v][l]);
      char * mass=prtclbase1[L->partcl].massidnt;
      if(L->moment>0 && strcmp(mass,"0"))
      {
        if(pseudop(L->partcl)||(L->lorentz&& set_in(L->lorentz,setmassindex0)))
          sprintf(factortxt+strlen(factortxt),"*%s^2",mass);
        else if(prtclbase1[L->partcl].spin==3)
          sprintf(factortxt+strlen(factortxt),"*3*%s^2",mass); 
        else if(prtclbase1[L->partcl].spin==4)
          sprintf(factortxt+strlen(factortxt),"*6*%s^4",mass);
      }
   }

   strcat(factortxt,")");    
   tf_add = (rmptr) read_rmonom(factortxt);
   mult_rptr(t_fact,&tf_add);

   vardef++;
   clearVars(vardef);


   m = (preres)readExpression(rmonomtxt(**t_fact),rd_pre,act_pre,NULL);

   if(rderrcode) save_sos(rderrcode);
   for(n=0; n<vardef->nvar; n++) vardef->vars[n].maxdeg=m->varsdeg[n]+1;
   sortvar();
   
   symb_start(vardef->nvar,vardef->vars, 0,0,0);

   mm = symb_read(smonomtxt((**t_fact).n));
   factn = mm.expr.p;
   mm = symb_read(smonomtxt((**t_fact).d));
   factd = mm.expr.p;
   
   clrvm((*t_fact)->n.v);
   clrvm((*t_fact)->d.v);
   free(*t_fact);
}


static void  formBlocks(void)
{  int  v, vv, l, lori, v1, v2, count; 
   int  vertmap[2 * maxvert];

   for(v=0; v<2*maxvert; v++)  vertmap[v]=fermmap[v]; 
   for(count=0; count<nloop; count++)  block[count]=fermres[count];

   count=nloop; 
   for(v=0; v<vcs.sizet; v++)  if(!vertmap[v]) 
   { symb_data tmp=symb_typeUp(vert[v], 1,etenstp);
     block[count]=tmp;
     vertmap[v]=++count;
   }
   for(v=0; v<count; v++) 
   {
      vertinfo[v].vlnc=0; 
      vertinfo[v].weight=1;
      vertinfo[v].ind=set_constr(_E);
      vertinfo[v].g5=0;
   } 
   for(v=0; v<vcs.sizet; v++) 
   { 
      v1=vertmap[v]-1;
      if(v1<nloop) 
      { 
        vertinfo[v1].weight += 2;
        if(fermloops[v1].g5) { vertinfo[v].g5=1; vertinfo[v1].weight++;} 
      } 

      for(l=0; l<vcs.valence[v]; l++)
      {  int np=vcs.vertlist[v][l].partcl;
         vv=vcs.vertlist[v][l].link.vno;
         v2=vertmap[vv]-1;
         lori=vcs.vertlist[v][l].lorentz; 
         if(lori && (prtclbase1[np].spin!=4||tolower(prtclbase1[np].hlp)=='t'))
         {
            if(v1!=v2)
            { 
               vertinfo[v1].link[vertinfo[v1].vlnc++]=v2+1;
               set_add1(&(vertinfo[v1].ind),lori);
               if(tolower(prtclbase1[np].hlp)=='t') set_add1(&(vertinfo[v1].ind),lori-1);
            }
            vertinfo[v1].weight += 2;
            if(set_in(lori,setmassindex)) vertinfo[v1].weight++;
         }
      }
   }
   for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++)
   {  int np=vcs.vertlist[v][l].partcl;
      if(prtclbase1[np].spin==4 && tolower(prtclbase1[np].hlp)!='t')
      { vv=vcs.vertlist[v][l].link.vno;
        if(v<vv && vertmap[v]!=vertmap[vv] )
        {  int ll=vcs.vertlist[v][l].link.edno;
           int m2=vcs.vertlist[v][l].lorentz,   m1=m2-1;
           int n2=vcs.vertlist[vv][ll].lorentz, n1=n2-1;
           symb_data m=tPropagator(v,l);             
           v1 = vertmap[v]-1;
           v2 = vertmap[vv]-1;

           vertinfo[count].vlnc = 1;
           vertinfo[count].link[0]=v1+1;
           vertinfo[v1].link[vertinfo[v1].vlnc++]=count+1;
           if(v1!=v2)
           {
              vertinfo[count].vlnc = 2;
              vertinfo[count].link[1]=v2+1;
              vertinfo[v2].link[vertinfo[v2].vlnc++]=count+1;
           }
           vertinfo[count].weight = 2;
           vertinfo[count].ind=set_constr(m1,m2,n1,n2,_E);
           vertinfo[v1].ind=set_or(vertinfo[v1].ind,set_constr(m1,m2,_E));
           vertinfo[v2].ind=set_or(vertinfo[v2].ind,set_constr(n1,n2,_E));
           vertinfo[count].g5 = 0;
           block[count]=m; 
           count++;
        }
      }
      else  if(prtclbase1[np].spin==2 && PLR_PRTCL&vcs.vertlist[v][l].prop )
      { int P=vcs.vertlist[v][l].moment;
        vv=vcs.vertlist[v][l].link.vno;
        if(P>0 /*&& vertmap[v]!=vertmap[vv]*/ )
        {  int ll=vcs.vertlist[v][l].link.edno;
           int m1=vcs.vertlist[v][l].lorentz;
           int m2=vcs.vertlist[vv][ll].lorentz;
           symb_data m;
           char txt[150];
           int Paux;
           if(P==1) Paux=2; else Paux=1;
           sprintf(txt,"m%d.m%d+i*HelicityN%d*eps(p%d,p%d,m%d,m%d)",
           m1,m2, P,P,Paux,m1,m2);           
            
                                                      
           m=symb_read(txt);

           v1 = vertmap[v]-1;
           v2 = vertmap[vv]-1;

           vertinfo[count].vlnc = 1;
           vertinfo[count].link[0]=v1+1;
           vertinfo[v1].link[vertinfo[v1].vlnc++]=count+1;
           if(v1!=v2)
           {
              vertinfo[count].vlnc = 2;
              vertinfo[count].link[1]=v2+1;
              vertinfo[v2].link[vertinfo[v2].vlnc++]=count+1;
           }
           vertinfo[count].weight = 1;
           vertinfo[count].ind=set_constr(m1,m2,_E);
           vertinfo[v1].ind=set_or(vertinfo[v1].ind,set_constr(m1,_E));
           vertinfo[v2].ind=set_or(vertinfo[v2].ind,set_constr(m2,_E));
           vertinfo[count].g5 = 1;
           block[count]=m; 
           count++;
        }
      }    
   }      
   n_vrt = count; 

   if (MEMORY_OPTIM)
   {
      for (v = 0; v < n_vrt; v++) vertinfo[v].weight = 0;
      for (v = 0; v < vcs.sizet; v++) vertinfo[vertmap[v]-1].weight++;
   }
}


static int  delImageryOne(rmptr   t_factor)
{
  vmrec       rec;
  vmptr       m, m1;

  rec.next =(t_factor->n).v;
  m = &rec;
  m1 = m->next;
  while (m1 != NULL && (strcmp(m1->name,"i") != 0 ) )
  {
     m = m1;
     m1 = m1->next;
  }

  if(m1)
  {  
     m->next = m1->next;
     free(m1);
     m1 = m->next;
     (t_factor->n).v=rec.next;
     return 1; 
  } else return 0;
        
}



static int  symbcalc(hlpcsptr ghst)
{  vcsect       vcs_copy;
   hlpcsptr     gstcopy;
   int          maxmom_s;
   int          spinl_s;
   int          first;
   long         del;
   int          i,j;
   s_listptr    d_facts, df;
   rmptr        t_fact;
   vmptr        coefvar;
   poly          mon;
   polyvars *   vardef_s;

   goto_xy(14,15); clr_eol();
   wrtoperat("Factors normalization");
   diagramsrfactors(ghst,&d_facts,&t_fact);
   wrtoperat("Preparing for calculation");
   vcs_copy = vcs;
   first = 1;

   gstcopy = ghst;
   df = d_facts;
   do
   {  
      coloringvcs(ghst);
      attachvertexes();
      firstvertexreading();
      propagatorsfirstreading();
      coefvar = df->monom.v;
      while (coefvar != NULL)
      {
         addvar(coefvar->name,coefvar->deg);
         coefvar = coefvar->next;
      }
      if (first)
      {
         vardef_s = vardef;
         vardef++;
         maxmom_s = maxmomdeg;
         spinl_s = findSpinLength();
         first = 0;
      }
      else
      {  int sl=findSpinLength();
         unite_vardef(vardef_s,vardef);
         spinl_s = MAX(sl,spinl_s);
         maxmom_s = MAX(maxmomdeg,maxmom_s);
      }

      vcs = vcs_copy;
      df = df->next;
      ghst = ghst->next;
   }  while (ghst);
   vardef = vardef_s;
   maxmomdeg = maxmom_s;

   addinoutmasses();  
   addscmult(spinl_s);
   clearpregarbage();  
   ghst = gstcopy;
   df = d_facts;
   rnum=symb_read("0");

   do
   {   
      goto_xy(14,15);
      print("%d(of %d)  ",ghst->num,ghst->maxnum);
      coloringvcs(ghst);
      attachvertexes();
      findReversVert();
      secondVertexReading();
      wrtoperat("Fermion loops calculation  ");

      calcFermLoops();
#ifdef STRACE
printf("calcFermLoops() OK\n");
#endif

      formBlocks();

#ifdef STRACE
printf("formBlocks() OK\n");
#endif
      makeprgcode();

#ifdef STRACE
printf("makeprgcode(); OK\n");      
#endif
#ifdef STRACE
writeF("n_vrt= %d",n_vrt);
for (i = 0;i<=n_vrt-1;i++)
{
   writeF("\n vrt%d:=",i);
   symb_print("",block[i]);
}
#endif
      for (i = n_vrt - 2; i >= 0; i--)
      {

#ifdef STRACE
 writeF("\n multiplication level %d\n",i);
 symb_print("T1:=",block[prgcode[i][0]-1]);
 symb_print("T2:=",block[prgcode[i][1]-1]);
#endif
          multblocks(prgcode[i][0],prgcode[i][1],i);
#ifdef STRACE
 symb_print("T1*T2:=",block[MIN(prgcode[i][0],prgcode[i][1])-1]);
#endif
      }

      {int mom, np;
       symb_data pp,mm;
          for (i = 0; i < vcs.sizet; i++)
          for (j = 0; j < vcs.valence[i]; j++)
          {
	     mom=vcs.vertlist[i][j].moment;
	     np=vcs.vertlist[i][j].partcl;
             if ((mom >0)&& (tolower(prtclbase1[np].hlp) == 't') )
             {   mm=masscalc(i,j,2);
                 if (prtclbase1[ghostmother(np)].hlp=='*') 
                         block[0]=symb_mult(block[0],1,mm,1);
                 else
                 {  char txt[20];
                    sprintf(txt,"-p%d.p%d",mom,mom);
                    pp = symb_read(txt);
                    pp=symb_sum(pp,1,mm,1);
                    block[0]=symb_mult(block[0],1,pp,1); 
		 }
		 if(set_in(vcs.vertlist[i][j].lorentz,setmassindex0))
	         {
	            mm=masscalc(i,j,2);
		    block[0]=symb_mult(block[0],1,mm,1);
		 }
             }
          }
      }
      {  int k;
         int sgn=ghst->sgn;
         for(k=0;k<vcs.sizet;k++) sgn*=vertexes[k].lgrnptr->factor;        
         if(sgn!=1) block[0]=symb_imult(block[0],1,sgn);
      }
      block[0]=symb_mult(block[0],1,symb_read(smonomtxt(df->monom)),1); 
      block[0]=symb_delEps(block[0],1);
      rnum=symb_sum(rnum,1, block[0],1);

      vcs = vcs_copy;
      df = df->next;
      ghst = ghst->next;
   }  while (ghst);
   eraseslist(d_facts);

   if(symb_iszero(rnum)) return 2;
     
   if(delImageryOne(t_fact)) rnum=symb_mult(rnum,1, symb_read("i"),1);

      
   rnum=symb_real(rnum,1);

   if(symb_iszero(rnum)) return 2;

   if (/*consLow*/ nin+nout<=5) del_pp(&(rnum.expr.t->re), &mon,&del);
       else  { mon= polyfactor(rnum.expr.t->re); del=1;}

   if(rnum.expr.t->re==NULL) return 2;

   wrtoperat("Total factor calculation");
   transformfactor(&t_fact,mon,del);
   wrtoperat("Denominator calculation");
   if(factn) return 1; else return 2;
}


static void  calcproc(csdiagram* csdiagr)
{  hlpcsptr    gstlist;

   transfdiagr(csdiagr,&vcs);
   cwtarg(&vcs);
   if (vcs.clrnum == 0)  csdiagr->status = 2;
   else
   {
      generateghosts(&vcs,&gstlist);
      if (gstlist == NULL) csdiagr->status = 2; else
      {
         preperdiagram();
         csdiagr->status = symbcalc(gstlist);
      }
      eraseghosts(gstlist);
   }
}


static void  writestatistic(unsigned noutmemtot,char* txt)
{
   goto_xy(1,12);  print("%u\n",calcdiag_sq);
   if(noutmemtot) {goto_xy(1,13);  print("%u Out of memory   ",noutmemtot);}
   goto_xy(17,14); print("%4d",ndiagr);
   goto_xy(42,14); print("%s      ",txt);
   goto_xy(14,17); clr_eol();
}


static void heap_is_empty(void) { save_sos(-1);}

void  calcallproc(void)
{  int        ndel, ncalc, nrest;
   long       nrecord;
   csdiagram  csd;
   unsigned   noutmemtot;
   shortstr   txt;
   marktp     heap_beg;
   int        ArcNum;

   polyvars  varsInfo[3]={ {0,NULL}, {0,NULL}, {0,NULL} };
   
   memerror=heap_is_empty;
   memoryInfo=memoryInfo_;

   goto_xy(1,14);
   scrcolor(FGmain,BGmain);

//   print("0     Out of memory\n");
   
   print("current diagram          in (Sub)process \n");
   print("Subdiagram  :\n");
   print("Used memory :%3d Kb       \n",(int)(usedmemory/1000));
   print("Operation   :\n");
   scrcolor(Yellow,Blue);
   print("\n");
   print(" Press Esc to halt calculations \n");
   
   scrcolor(Blue,BGmain);

   
   catalog = fopen(CATALOG_NAME,"rb");
   if(catalog==NULL)  ArcNum=1; else
   { long pos=fseek(catalog,0, SEEK_END);
     if(pos==0)  ArcNum=1; else 
     { catrec cr;
       fseek(catalog,-sizeof(catrec), SEEK_END);
       fread(&cr, sizeof(cr), 1, catalog);
       ArcNum=cr.nFile;  
     }
     fclose(catalog);
   }     
                 
   catalog = fopen(CATALOG_NAME,"ab");
  
   calcdiag_sq = 0;
   noutmemtot = 0;
   diagrq=fopen(DIAGRQ_NAME,"rb");
   while(FREAD1(csd,diagrq))
   {
      switch (csd.status)
      {
         case 1:  ++(calcdiag_sq);  break;
         case -2: ++(noutmemtot);    break;
         case 2:  ++(calcdiag_sq);
      }
   }
   fclose(diagrq);
   diagrq=fopen(DIAGRQ_NAME,"r+b");
   menuq=fopen(MENUQ_NAME,"r+b");
   for(nsub=1;nsub<=subproc_sq;nsub++)
   {  int naux;  
      rd_menu(2,nsub,txt,&ndel,&ncalc,&nrest,&nrecord);
      naux = ndel + ncalc + nrest;
      for (ndiagr = 1; ndiagr <= naux; ndiagr++)
      {
         writestatistic(noutmemtot,txt);
         fseek(diagrq,sizeof(csd)*(nrecord+ndiagr-1),SEEK_SET)  ;
         FREAD1(csd,diagrq);
         if (csd.status == 0)
         {  
            vardef=&varsInfo[0];
            mark_(&heap_beg);
            calcproc(&csd);
            fseek(diagrq,sizeof(csd)*(nrecord+ndiagr-1),SEEK_SET);
            FWRITE1(csd,diagrq);
            if(csd.status == 1)
            {
               wrtoperat("Writing result             ");
               vardef=&(varsInfo[0]);     
               ArcNum=whichArchive(ArcNum,'w'); 
               saveanaliticresult(rnum.expr.t->re,factn,factd, vcs,ArcNum); 
               newCodes=1;
            }
            release_(&heap_beg);
            {int i; for(i=0;i<3;i++) clearVars(varsInfo+i);}                     
            ncalc++;
            nrest--;
            calcdiag_sq++;            
            wrt_menu(2,nsub,txt,ndel,ncalc,nrest,nrecord);
            if (escpressed()) goto exi;
         }         
      }
   }
exi:

   fclose(diagrq);
   fclose(menuq);
   
   fclose(catalog);
   whichArchive(0,0);      
   scrcolor(FGmain,BGmain);
   clrbox(1,14,70,20);
   memerror=NULL;
}
