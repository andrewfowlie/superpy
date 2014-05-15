/*
 Copyright (C) 1997, Alexander Pukhov, e-mail: pukhov@theory.npi.msu.su
*/

#include "chep_crt.h"
#include "syst2.h"
#include "physics.h"
#include "parser.h"
#include "ghosts.h"
#include "cweight.h"
#include "prepdiag.h"
#include "reader0.h"
#include "s_files.h"
#include "diaprins.h"
#include "chess.h"
#include "out_serv.h"
#include "rfactor.h"
#include "denominators.h"
#include "process.h"
#include "writeF.h"

#include "r_code.h"

static int    vertmap[3 * maxvert];
static int    gammaflag;
static char * parsedverts[2*maxvert];
static int    numberg5;


static char * pexpr(char p)
{  static char  snum[8];
   if(p>=0) sprintf(snum,"P%d",p);else sprintf(snum,"(-P%d)",-p);
   return snum;
}


static char * iexpr(int i)
{  static char  snum[6];
   sprintf(snum,"m%d",i);
   return snum;
} 

 
static void  head(void)
{  int  l;

   writeF("%% ----------- VARIABLES ------------------ \n");
   writeF(" vector  A,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ZERO_;\n");
   writeF(" vector  m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16;\n");
   writeF("%%\n");
   writeF("%%--------- Mass shell declarations -----------\n");
   for (l = 1; l <= nin + nout; l++)
      writeF(" MASS  P%d = %s$  MSHELL P%d$\n",l,inoutmasses[l-1],l);
   writeF("%%\n");
   writeF("operator propDen$\n");
}


static void  writesubst(void)
{ int i,j,l,c;

  writeF("%%-------- Moment substitutions --------\n");
  for (i = 0; i < 3 * maxvert; i++)
  {  l = momdep[i][0]; /* For momdep this is a length V.E. */
     if (l > 1)
     {
        writeF(" Let  p%d = ",i+1);
        for (j = 1; j <= l; j++)
        {
           c =momdep[i][j];
           if (c > 0)  writeF("+p%d",c); else writeF("-p%d",-c);
        }
        writeF("$\n");
     }
  }
  writeF(" Let Sqrt2=sqrt(2)$\n");
  writeF("%%\n");
}

static void  emitfactors(void)
{
   int	 c, i, j, s, vln, pnum;

   writeF("%%---------- Factors ---------------\n");
   /* ------- Symmetry Factor ------- */
   writeF(" SymmFact:=%d/%d$    %% Diagram symmerty factor\n",
                    vcs.symnum,vcs.symdenum);
   /* -----  average factor  --------- */
   for (c=1,i=0; i < vcs.sizel; i++)
   {
      vln = vcs.valence[i];
      for (j = 0; j < vln; j++)
      if (IN_PRTCL & vcs.vertlist[i][j].prop)
      {
         pnum = vcs.vertlist[i][j].partcl;
         s = prtclbase1[pnum].spin;
         switch (s)
         {
            case 1: if(!strchr("LR",prtclbase1[pnum].hlp)) c*=2; break;
            case 2: if(zeromass(pnum)) c*=2; else c*=3; break;
            case 3: if(zeromass(pnum)) c*=2; else c*=4; break;
            case 4: if(zeromass(pnum)) c*=2; else c*=5;                   
         } 
         c *= abs(prtclbase1[pnum].cdim);
      }
   }
   writeF(" AverFact:=1/%d$       %% Normalization factor of polarization average\n",c);
   /* ----- Fermion factor  --------- */
   for (c=1,i= 0; i < vcs.sizel; i++)
   for (j=0; j <vcs.valence[i]; j++)
   {
     pnum=vcs.vertlist[i][j].partcl;
     if(prtclbase1[pnum].spin&1 && (IN_PRTCL & vcs.vertlist[i][j].prop))c*=-1;
   }
   
   writeF(" FermFact:=%d$      %% (-1)**(number of in-fermion particles)\n",c);
   /* ----- Color factor  --------- */
   writeF(" ColorFact:=%d/%d$    %%  QCD color weight of diagram  \n%%\n",
         vcs.clrnum,vcs.clrdenum);
}


 /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 | Program for emitting denominators of propagators in CrossSection |
 |                        diagramm.                                 |
 |                September 19, 1989                                |
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static void  emitdenoms(void)
{int   j; /* 1..2 * maxvert */
 int   k; /* 1..4 */
 int denrno;

   /*  EmitDenoms -- main  */
   
   for (j = 1; j <= vcs.sizet; j++)
   for (k = 1; k <= vcs.valence[j-1]; k++)
   { edgeinvert *ln = &vcs.vertlist[j-1][k-1];
     if(ln->moment>0)
     { char*mass= prtclbase1[ln->partcl].massidnt;
        if(!(ln->prop&(IN_PRTCL|OUT_PRTCL)) && pseudop(ln->partcl))
            writeF("totFactor_:=totFactor_/%s^2$\n",mass);
        else 
        { int spin=prtclbase1[ln->partcl].spin;
          if(spin==3)writeF("totFactor_:=totFactor_/(3*%s^2)$\n",mass);
          else if(spin==4)writeF("totFactor_:=totFactor_/(6*%s^4)$\n",mass);
        }  
     }                    
   } 
   
   writeF("denominator_:=");
   denrno = 0;
   for (j = 1; j <= vcs.sizet; j++)
   for (k = 1; k <= vcs.valence[j-1]; k++)
   {edgeinvert *ln = &vcs.vertlist[j-1][k-1];

     if (ln->moment > 0 && !(ln->prop & (IN_PRTCL|OUT_PRTCL)) )
     {
        if (!pseudop(ln->partcl))
        {  char width[VAR_NAME_SIZE]="0";
           if(denrno) writeF("*");
           denrno++;
/*           if(!ttypepropag(j,k))*/ strcpy(width,prtclbase1[ln->partcl].imassidnt);
           writeF("propDen(p%d,%s,%s)",ln->moment,
                            prtclbase1[ln->partcl].massidnt,width);
        }
     }
   } 
   if(!denrno) writeF("1");
   writeF("$\n");

}  /*   EmitDenoms   */

/*  Programm for constructing vertex-oriented representation for !
! amplitude diagram.                                             !
! A.Taranov 20.07.89                                            */


static void  emitindex(set *index)
{  int n,l;

   l=set_first(index,1);

   if(!l) return;
   writeF(" Index ");
   for(n=0;(n<4)&&l;n++)
   {
      if(n) writeF(",");
      writeF("m%d",l);
      set_del1(index,l);
      l=set_first(index,1);
   }
   writeF("$\n");
}

static void reducemult(char* nameres,char* name1,char* name2,set index)
{  int l, n;
   set  index1=index,index2=index;

   emitindex(&index1);
   writeF(" %s:=%s*%s$\n",nameres,name1,name2);
   while (!set_eq0(index1))
   {
      emitindex(&index1);
      writeF(" %s:=%s$\n",nameres,nameres);
   }
   if(!set_eq0(index2))
   {
      writeF(" RemInd  ");
      for(l=1,n=1;!set_eq0(index2);l++)
      {
         if(set_in(l,index2))
         {
            if (n != 1) writeF(",");
            n++;
	    writeF("m%d",l);
         }
         set_del1(&index2,l);
      }
      writeF("$\n");
   }
}


static int  fermvrt(int v)
{int  i;

   for (i = 0; i < vcs.valence[v-1]; i++)
      if (fermionp(vcs.vertlist[v-1][i].partcl)) return 1;
   return 0;
}

static char * subst(int ind,int arg)
{  int         v, l, i;
   static shortstr ans;

   v = massindpos[ind-1].vrt1;
   if (vertmap[v-1] == arg) l = massindpos[ind-1].ln1;
   else
   {
      v = massindpos[ind-1].vrt2;
      l = massindpos[ind-1].ln2;
   }

   sprintf(ans,"%s=>(",iexpr(ind));
   if (gammaflag && fermvrt(v))
   {
      for (i = 0; i < vcs.valence[v-1]; i++) if (i != l-1)
        sprintf(ans+strlen(ans),"+%s",pexpr(vcs.vertlist[v-1][i].moment));
   }
   else sprintf(ans+strlen(ans),"%s",pexpr(-vcs.vertlist[v-1][l-1].moment));

   sprintf(ans+strlen(ans),")/%s",prtclbase1[vcs.vertlist[v-1][l-1].partcl].massidnt);
   return ans;
}


static void  r_mult(int arg1,int arg2,set index1)
{  char      name0[8],name1[8],name2[8];
   set       mind;
   unsigned  m,mm,maxmm;
   int       pstn[15];
   int       n,l;
   set       index;

   index=index1;
   sprintf(name1,"Vrt_%d",arg1);
   sprintf(name2,"Vrt_%d",arg2);
   if(arg1 < arg2) strcpy(name0,name1); else strcpy(name0,name2);
   mind=set_and(index,setmassindex);
   if(set_eq0(mind)) reducemult(name0,name1,name2,index); else
   {
      reducemult("Vrt_0",name1,name2,index);
      l = 1; n = 0;
      while (!set_eq0(mind))
      {
         if (set_in(l,mind)) pstn[n++] = l;
         set_del1(&mind,l++);
      }
      maxmm =  (1 << n) - 1 ;
      for (mm = 1;mm <= maxmm ; mm++)
      {
         writeF(" Vrt_L:=%s$   Vrt_R:=%s$\n",name1,name2);
         mind=set_constr(_E);
         m = mm; n = 1;
         while (m != 0)
         {
            if ((m & 1) != 0 && (l = pstn[n-1]) != 0)
            {
               writeF(" Vrt_L:=(Vrt_L where %s)$\n",subst(l,arg1));
               writeF(" Vrt_R:=(Vrt_R where %s)$\n",subst(l,arg2));
               set_add1(&mind,l);
            }
            ++(n);
            m >>= 1;
         }
         reducemult("Vrt_0","Vrt_0 + Vrt_L","Vrt_R",set_aun(index,mind));
         writeF(" Clear Vrt_L $    Clear Vrt_R $\n");
      }
      writeF(" %s:=Vrt_0 $    Clear Vrt_0 $\n",name0);
   }
   if (arg1 < arg2)
   {
      writeF(" Clear %s$\n",name2);
      for (n = 1; n <= vcs.sizet; n++)
         if (vertmap[n-1] == arg2) vertmap[n-1] = arg1;
   }
   if (arg2 < arg1)
   {
      writeF(" Clear %s$\n",name1);
      for (n = 1; n <= vcs.sizet; n++)
      if (vertmap[n-1] == arg1) vertmap[n-1] = arg2;
   }
}  /*  R_mult  */


static void  fermprgemit(void)
{  int  v, v1, l, lpcount, nn;
   set  indexs;
   char  * vertStr, * propStr;

   if (nloop == 0) { gammaflag = 0; return; }
   gammaflag = 1;

   for(v=0; v<nloop; v++)
   for(l=0;!fermloops[v].g5 && l<fermloops[v].len;l++) 
   if(strstr(parsedverts[fermloops[v].vv[l]-1],"G(ln,A)")) fermloops[v].g5=1;
   
   for(v=0,numberg5=0; v<nloop; v++) if(fermloops[v].g5) numberg5++;

for(v=0; v<2*maxvert; v++) vertmap[v]=0;

   for(lpcount=0; lpcount<nloop; lpcount++)
   {  fermloopstp * FL=&(fermloops[lpcount]);

      writeF("%%------- Fermion loop calculation ------- \n");
      writeF(" NoSpur ln$\n");
      
      for(nn=0; nn<FL->len; nn++)
      {  int v=FL->vv[nn]-1, l=FL->ll[nn]-1;
      
         if(vertexes[v].r_vert)
	 writeF("%% reversed vertex\n");

	 vertStr=parsedverts[v];
         if(FL->spin[nn]==1) propStr=fermPropagTxt(v,l,1);
         else                propStr=spin3_2_propagator(v,l,1);
         if( FL->ind[nn][0]) writeF("Index m%d;\n",FL->ind[nn][0]); 
	 writeF(" Vrt_%d:=(%s)*(%s)$\n",1+nn+lpcount,vertStr,propStr);
         if( FL->ind[nn][0]) writeF("RemInd m%d;\n",FL->ind[nn][0]);	 
	 free(propStr);
 vertmap[FL->vv[nn]-1] = nn +1+ lpcount;
      }
      
      for (v1=1; v1<FL->len; v1++)
      {  int l,v,s,i;

         indexs=set_constr(_E);
                  
         for(i=0;i<FL->nint[v1];i++)
         { l=FL->intln[v1][i]-1; 
           v=FL->vv[v1]-1;
           s=prtclbase1[vcs.vertlist[v][l].partcl].spin; 
           if(s==2) set_add1(&indexs, vcs.vertlist[v][l].lorentz);
         }  
         if(FL->ind[v1-1][1]) set_add1(&indexs,FL->ind[v1-1][1]);
         if(v1==FL->len-1 && FL->ind[v1][1]) set_add1(&indexs,FL->ind[v1][1]);
         r_mult(lpcount+1,lpcount+1 + v1,indexs);
                  
         for(i=0;i<FL->nint[v1];i++)
         { l=FL->intln[v1][i]-1; 
           v=FL->vv[v1]-1;
           s=prtclbase1[vcs.vertlist[v][l].partcl].spin;
           if(s==4 && prtclbase1[vcs.vertlist[v][l].partcl].hlp!='t')
           {  char vert[10];
              char*txt=tPropagator6mass4(v,l,&indexs);               
              writeF("tProp:=(%s)/(6*%s^4)$\n",txt,prtclbase1[vcs.vertlist[v][l].partcl].massidnt);
              free(txt); 
              sprintf(vert,"Vrt_%d",lpcount+1);

              reducemult(vert,vert,"tProp",indexs);
              writeF("Clear tProp$\n");
           }
         }           
      }


/*===========================================================*/

      if(numberg5 == 1 && FL->g5)
      {
         writeF(
                   " Vrt_%d:=Sub(A=0*A,Vrt_%d)$\n"
	        /* " Vrt_%d:=(Vrt_%d where A=>0*ZERO_)$\n"*/ 
	                                           ,lpcount+1,lpcount+1);
         numberg5 = 0;
         FL->g5 = 0;
      }
/*
      for (k = 1; k <= strlen(FL->invrt); k++)
      {
         v1 = FL->invrt[k-1];
         indexs=set_constr(_E);
         for (l = 1; l <= vcs.valence[v1-1]; l++)
         {
            i = vcs.vertlist[v1-1][l-1].lorentz;
            if(i) set_add1(&indexs,i);
         }
         vertmap[v1-1] = lpcount + 2;
         writeF("Vrt_%d:=%s$\n",lpcount + 2,parsedverts[v1-1]);
         r_mult(lpcount+1,lpcount + 2,indexs);
      }
*/
      writeF(" Spur ln $\n");
      writeF(" Vrt_%d:=-4*Vrt_%d$\n",lpcount+1,lpcount+1);
      writeF("%%\n");
      writeF(" Fl%d:=Vrt_%d$\n",lpcount+1,lpcount+1);
      writeF("%%\n");
   }

   for (v1 = 0; v1 < 2 * maxvert; v1++)  fermmap[v1] = vertmap[v1];
   writeF("%%\n");
   gammaflag = 0;
}

static void  formblocks(int* vrtsize)
{  int   v, l, lori, v1, vv, v2, count;

   for (v = 0; v < 2 * maxvert; v++)  vertmap[v] = fermmap[v]; 
   for (count = 1; count <= nloop; count++) /*  ferm loop vertex  */
   {
      writeF(" Vrt_%d:=Fl%d$",count,count);
      writeF(" Clear Fl%d$\n",count);
   }

   count = nloop;
   for (v = 0; v < vcs.sizet; v++) if(!vertmap[v])
   {  vertmap[v] = ++count;
      writeF(" Vrt_%d:=%s$\n",count,parsedverts[v]);
   } 

   for (v = 0; v < count; v++)
   {  vertinfo[v].vlnc = 0;
      vertinfo[v].weight = 1;
      vertinfo[v].ind=set_constr(_E);
      vertinfo[v].g5 = 0;
   }

   for(v=0; v<vcs.sizet; v++)
   {
      v1 = vertmap[v]-1;
      if(v1 < nloop)
      {
         vertinfo[v1].weight += 2;
         if(fermloops[v1].g5) vertinfo[v1].g5 = 1;  
      }
      for(l=0; l<vcs.valence[v]; l++)
      {  int np=vcs.vertlist[v][l].partcl;
         vv=vcs.vertlist[v][l].link.vno;
         v2=vertmap[vv]-1;
         lori=vcs.vertlist[v][l].lorentz;
         if(lori && (prtclbase1[np].spin!=4 || prtclbase1[np].hlp=='t'))
         { 
            if(v1!=v2)
            {  
               vertinfo[v1].link[vertinfo[v1].vlnc++]=v2+1;
               set_add1(&(vertinfo[v1].ind),lori);
               if(prtclbase1[np].hlp=='t') set_add1(&(vertinfo[v1].ind),lori-1);
            }
            vertinfo[v1].weight += 2;
            if(set_in(lori,setmassindex))  ++(vertinfo[v1].weight);
         }
      }
   }
   
   for(v=0; v<vcs.sizet; v++) for(l=0; l<vcs.valence[v]; l++)
   {  int np=vcs.vertlist[v][l].partcl;
      if(prtclbase1[np].spin==4 && prtclbase1[np].hlp!='t')   
      { vv=vcs.vertlist[v][l].link.vno;
        if(v<vv && vertmap[v]!=vertmap[vv] )  
        {  int ll=vcs.vertlist[v][l].link.edno;
           int m2=vcs.vertlist[v][l].lorentz,   m1=m2-1;
           int n2=vcs.vertlist[vv][ll].lorentz, n1=n2-1;

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
           count++;
           { 
             int P=abs(vcs.vertlist[v][l].moment);
             char *mass=prtclbase1[np].massidnt;
              
             writeF(" Vrt_%d:=(m%d.m%d*m%d.m%d+m%d.m%d*m%d.m%d-m%d.m%d*m%d.m%d)/2",
             count, m1,n1, m2,n2,  m1,n2,m2,n1,  m1,m2,n1,n2);
             writeF("-(m%d.m%d*P%d.m%d*P%d.m%d+m%d.m%d*P%d.m%d*P%d.m%d+"
              "m%d.m%d*P%d.m%d*P%d.m%d+m%d.m%d*P%d.m%d*P%d.m%d )/(2*%s^2)",
             m1,n1,P,m2,P,n2, m2,n2,P,m1,P,n1, m1,n2,P,m2,P,n1, m2,n1,P,m1,P,n2,mass);
             
             writeF("+(m%d.m%d+2*P%d.m%d*P%d.m%d/%s^2)*(m%d.m%d+2*P%d.m%d*P%d.m%d/%s^2)/6$\n",
             m1,m2, P,m1,P,m2,mass, n1,n2, P,n1,P,n2,mass);
           }
        }
      }

      else if(prtclbase1[np].spin==2 && PLR_PRTCL&vcs.vertlist[v][l].prop )   
      { int  P=vcs.vertlist[v][l].moment;
        vv=vcs.vertlist[v][l].link.vno;
        if(P>0 /*&& vertmap[v]!=vertmap[vv]*/ )  
        {  int ll=vcs.vertlist[v][l].link.edno;
           int m1=vcs.vertlist[v][l].lorentz;  
           int m2=vcs.vertlist[vv][ll].lorentz;
           
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
           count++;
           { int Paux;
             if(P==1) Paux=2; else Paux=1;
             writeF(" Vrt_%d:=m%d.m%d*p%d.p%d - p%d.m%d*p%d.m%d - p%d.m%d*p%d.m%d"
                             " -i*lPolar%d*eps(p%d,p%d,m%d,m%d)$",
             count, m1,m2, P,Paux,  P,m1,Paux,m2,  P,m2,Paux,m1,P,P,Paux,m1,m2);
           }
        }
      }  
   }  
   *vrtsize = count;
} 



static void  emitreducecode(void)
{
   int  i, j, v1, v2;
   set  ind1, ind2;

   for(i=0; i<vcs.sizet; i++) parsedverts[i]=parseVertex(i,1);
   fermprgemit();
   formblocks(&n_vrt);
   makeprgcode();
   for(i=0; i<vcs.sizet; i++) free(parsedverts[i]);

   for (i = n_vrt - 1; i >= 1; i--)
   {
      v1 = prgcode[i-1][0];
      v2 = prgcode[i-1][1];
      ind1=vertinfo[v1-1].ind;
      ind2=vertinfo[v2-1].ind;

      r_mult(v1,v2,set_and(ind1,ind2));
      vertinfo[MIN(v1,v2)-1].g5 = vertinfo[v1-1].g5 + vertinfo[v2-1].g5;
      vertinfo[MIN(v1,v2)-1].ind=
      set_aun(set_or(ind1,ind2),set_and(ind1,ind2));
   }
   
   
   for (i = 0; i < vcs.sizet; i++)
   for (j = 0; j <vcs.valence[i]; j++)
   { int mom,np;     
     mom=vcs.vertlist[i][j].moment;
     np =vcs.vertlist[i][j].partcl;
     if (mom>0 && prtclbase1[np].hlp == 't')
     {  writeF("Vrt_1:=Vrt_1*(%s**2",prtclbase1[np].massidnt);
       if(prtclbase1[ghostmother(np)].hlp=='*') writeF(")$\n");
        else                                    writeF("-p%d.p%d)$\n",mom,mom); 
     }              
   }   
   writeF("%%\n");
}


void  mk_reduceprograms(void)
{
   int          ndel, ncalc, nrest, i;
   long         nrecord, naxu;
   csdiagram    csd;
   unsigned     ncalctot;
   shortstr     txt;
   hlpcsptr     gstlist, c;
   vcsect       vcs_copy;
   s_listptr    d_facts, df;
   rmptr        t_fact;

   goto_xy(1,21); scrcolor(Yellow,Blue);
   print("  REDUCE code generation \n");
   scrcolor(Red,BGmain);
   print(" Generated........\n");
   print(" current diagram :\n");
   scrcolor(Yellow,Blue);
   print(" Press Esc to halt REDUCE codes generation ");
   scrcolor(FGmain,BGmain);
   diagrq=fopen(DIAGRQ_NAME,"rb");
   ncalctot = 0;
   menuq=fopen(MENUQ_NAME,"rb");

   for(nsub=1;nsub<=subproc_sq;nsub++) 
   {
      rd_menu(2,nsub,txt,&ndel,&ncalc,&nrest,&nrecord);
      fseek(diagrq,nrecord*sizeof(csdiagram),SEEK_SET);
      naxu = ndel + ncalc + nrest;
      for (ndiagr = 1; ndiagr <= naxu; ndiagr++)
      {
	 goto_xy(20,22); print("%u",ncalctot);
	 goto_xy(20,23); print("%u",ndiagr); clr_eol();
         FREAD1(csd,diagrq);
         if (csd.status != -1)
         { 
            outFileOpen("%sresults%cp%d_%d.red",pathtouser,f_slash,nsub,ndiagr);
            writeLabel('%');
            writeF("%%\n");
            
            transfdiagr(&csd,&vcs);

            cwtarg(&vcs);
            if (vcs.clrnum == 0)
            {
               writeF(
                  "%%-------  Zero color factor --------\n");
               writeF("totFactor_:=0$\n");
               writeF("numerator_:=0$\n");
               writeF("denominator_:=1$\n");
            }
            else
            {
               generateghosts(&vcs,&gstlist);
               if (gstlist == NULL)
               {
                  writeF( "%%-------  non-existent diagram  --------\n");
                  writeF("totFactor_:=0$\n");
                  writeF("numerator_:=0$\n");
                  writeF("denominator_:=1$\n");
               }
               else
               {
		  goto_xy(40,23);
		  print("(%% %4d subdiagrams)",gstlist->maxnum);
		  writeF("%% The total number of diagrams %d\n",gstlist->maxnum);
                  preperdiagram();
                  head();
                  emitfactors();
                  diagramsrfactors(gstlist,&d_facts,&t_fact);
                  writeF("totFactor_:=%s$\n",rmonomtxt(*t_fact));

                 writeF("totFactor_:="
                       "totFactor_*SymmFact*AverFact*FermFact*ColorFact$\n");
                 
                  clrvm(t_fact->n.v);
                  clrvm(t_fact->d.v);
                  free(t_fact);

                  writesubst();
                  writeF("numerator_:=0$\n");

                  c = gstlist;
                  df = d_facts;
                  vcs_copy = vcs;
                  while (c != NULL)
                  {
                     coloringvcs(c);
                     writeF("%%  diagram  number =   %d\n", c->num);
                     DiagramToOutFile(&vcs,1,'%');

                     {int k; int sgn=c->sgn;
                      for(k=0;k<vcs.sizet;k++) sgn*=vertexes[k].lgrnptr->factor;
                      writeF("  GhostFact:=%d$\n",sgn);
                     }
                     
		     findReversVert();
                     attachvertexes();

                     emitreducecode();

                     writeF(" numerator_:=numerator_ +(%s)*GhostFact*Vrt_1 $\n",
                        smonomtxt(df->monom));
                     writeF(" Clear Vrt_1,GhostFact$\n");
                     writeF("%%\n");

                     vcs = vcs_copy;
                     c = c->next;
                     df = df->next;
                  }

                  eraseslist(d_facts);
                  eraseghosts(gstlist);

                  vcs = vcs_copy;
                  emitdenoms();
	          writeF(" Clear p%d",nin + nout + 1);
                  for (i = nin + nout + 2; i <= 12; i++)
	          writeF(",p%d",i);
                  writeF("$\n");
                  writeF("%%\n");

               }
            }
            writeF("End$\n");
            outFileClose();
            --(nrest);
            ++(ncalctot);
            if (escpressed()) goto exi;
         }
      }
   }

exi: 
   fclose(diagrq); fclose(menuq);
   clrbox(1,21,70,24);

}

void  makeghostdiagr(int dnum,char * fname)
{csdiagram    csd;
 hlpcsptr     gstlist, c;
 vcsect       vcs_copy;
 FILE * diagrq;
   diagrq=fopen(DIAGRQ_NAME,"rb");
   fseek(diagrq,dnum*sizeof(csdiagram),SEEK_SET);
   FREAD1(csd,diagrq);
   {
      outFileOpen(fname);
      writeLabel('%');
      transfdiagr(&csd,&vcs);

      cwtarg(&vcs);
      if (vcs.clrnum == 0) writeF("Color weight equal zero \n"); else
      {
         generateghosts(&vcs,&gstlist);
         if (gstlist == NULL)
			{
            writeF(
               "Diagrams of this type are absent\n");
            csd.status = 2;
            fseek(diagrq,dnum*sizeof(csdiagram),SEEK_SET);
            FWRITE1(csd,diagrq);
         }
         else
         {
            emitfactors();
            c = gstlist;
            vcs_copy = vcs;
            writeF(" total number diagrams of this type is  %u\n",c->maxnum);
            while (c != NULL)
            {
               coloringvcs(c);
               {
                   writeF("  diagrams number =   %u\n",c->num);
                   DiagramToOutFile(&vcs,0,' ');                     
               }               
	       writeF("  GhostFact:=%d$\n",c->sgn);
               vcs = vcs_copy;
               c = c->next;
            }

         }
         eraseghosts(gstlist);
         vcs = vcs_copy;
      }
      writeF("End$\n");
      outFileClose();
      if (escpressed()) goto exi;
   }

exi:
   fclose(diagrq);
}
