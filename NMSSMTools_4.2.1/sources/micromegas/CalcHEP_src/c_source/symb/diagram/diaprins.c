/*
 Copyright (C) 1997, Victor Edneral 
*/
/****************************************
*                                       *          Written by V.Edneral         *
*              SCL of MSU               *
*               25.02.90                *
*                                       *
****************************************/
#include "physics.h"
#include "chep_crt.h"
#include "diagrams.h"
#include "sets.h"
#include "diaprins.h"

#define xstep 6
#define ystep 4  
#define chsp ' '  
#define chp  'P'  
#define koder (maxvert + 1)  /* was 20 ? */
#define xstart 0  
#define xend (xstep * koder - 1)
#define ystart 2  
#define yend (ystep * (koder + 1) + 1)
#define scrwidth 65  

#ifdef BLOCK_GRAPHICS

#define vertsign 0xDB  
#define chrC0    0xC0
#define chrD9    0xD9
#define chrBF    0xBF
#define chrDA    0xDA
#define chrC8    0xC8
#define chrBC    0xBC
#define chrBB    0xBB
#define chrC9    0xC9
#define chrC4    0xC4
#define chrCD    0xCD
#define chrB3    0xB3
#define chrBA    0xBA
#define chrC5    0xC5
#define chrD8    0xD8
#define chrD7    0xD7
#define chrCE    0xCE
#else

#define vertsign '@'  
#define chrC0    '\\'
#define chrD9    '/'
#define chrBF    '\\'
#define chrDA    '/'
#define chrC8    '\\'
#define chrBC    '/'
#define chrBB    '\\'
#define chrC9    '/'
#define chrC4    '-'
#define chrCD    '='
#define chrB3    '|'
#define chrBA    '|'
#define chrC5    '+'
#define chrD8    '+'
#define chrD7    '+'
#define chrCE    '+'
#endif

typedef struct intsect 
   { 
      int   y, xst, sp, selfvert, selfslot, 
             aliasvert, aliasslot; 
   }  intsect; 

typedef struct dscreen 
   { 
      int         scr[xend][yend]; 
      intsect      ints[maxvert + 1]; 
      int         mver, mhor, mxy; 
   }  dscreen; 

typedef int interscreen[koder][yend + ystep]; 


#define pnt  struct spool * 

typedef struct spool 
   { 
      pnt          add; 
      char         stn[100]; 
   }  spool;

#undef pnt

typedef spool* pnt; 



static  dscreen    * res;          /* from diaprin */
static  int      hcontr;       /* from diaprin */
static  int         isect, jsect; /* from isector */
static  int      ssect;        /* from isector */
static  int      left;         /* from diaprt  */
static  dscreen      resl, resr;   /* from diaprt  */
static  interscreen  is;           /* from diaprt  */
static  vcsect       diag;         /* from writepict  */
static  int      debug;        /* from writepict  */


static void cr(int n,int x,int y,int ivr,int slt)
{
   res->scr[x-1][y-1] =
      prtclbase[diag.vertlist[ivr-1][slt-1].partcl-1].spin % 2 == 0 ?
         ((n == 1 && left) || (n == 2 && !left)) ? chrC0 :
         ((n == 2 && left) || (n == 1 && !left)) ? chrD9 :
         ((n == 3 && left) || (n == 4 && !left)) ? chrBF :
                                                   chrDA            :
         ((n == 1 && left) || (n == 2 && !left)) ? chrC8 :
         ((n == 2 && left) || (n == 1 && !left)) ? chrBC :
         ((n == 3 && left) || (n == 4 && !left)) ? chrBB :
                                                   chrC9            ;
}   /* The end of Cr's body */



static void hline(int xs,int xe,int y,int ivr,int slt)
{int    linech, arrow; 
 int i, j, k, l; 
 char    nm[P_NAME_SIZE]; 

   k = hcontr ? (xs + xe) / 2 : res->mhor - (xstep / 2) + 1; 
   if (y > res->mxy) res->mxy = y;
   
   {edgeinvert * with1 = &diag.vertlist[ivr-1][slt-1];     
      linech = (prtclbase[with1->partcl-1].spin % 2) == 0 ?
               chrC4 : chrCD; 
      l = left ? prtclbase[with1->partcl-1].anti : with1->partcl; 
      arrow = l > prtclbase[l-1].anti ? '<' : '>'; 
      if (!(((IN_PRTCL|OUT_PRTCL) & with1->prop) || l < prtclbase[l-1].anti)) 
         l = prtclbase[l-1].anti; 
      strcpy(nm,prtclbase[l-1].name); 
      l = strlen(nm); 
      for (i = xs; i <= xe; i++)
         res->scr[i-1][y-1] = linech; 
      if (left) 
      { 
         j = l > 1 ? k + 1 - l : k - 1; 
         for (i = 1; i <= l; i++)  res->scr[i + j-1][y - 1-1] = nm[i-1]; 
         if (((IN_PRTCL|OUT_PRTCL) & with1->prop)  || debug) 
         { 
            if ((left ? 1 : 0) == 
                (with1->moment > 0 ? 1 : 0)) 
               res->scr[k - 1-1][y + 1-1] = '-'; 
            res->scr[k-1][y + 1-1] = chp; 
            res->scr[k + 1-1][y + 1-1] = 
               abs(with1->moment) > 9        ? 
               abs(with1->moment) + '7' :
               (abs(with1->moment) + '0') ; 
            if (debug && with1->lorentz > 0) 
               res->scr[k - 1-1][y-1] =
                  with1->lorentz > 9         ? 
                  (with1->lorentz + '7')  : 
                  (with1->lorentz + '0')  ; 
         } 
      } 
      else 
      { 
         j = l > 1 ? k + l - 1 : k + 1; 
         for (i = 1; i <= l; i++) 
            res->scr[j - i-1][y - 1-1] = nm[i-1]; 
         if (((IN_PRTCL|OUT_PRTCL) & with1->prop) || debug) 
         {  if ((left ? 1 : 0) == 
                (with1->moment > 0 ? 1 : 0)) 
               res->scr[k + 1-1][y + 1-1] = '-'; 
            res->scr[k-1][y + 1-1] = chp; 
            res->scr[k - 1-1][y + 1-1]        = 
               abs(with1->moment) > 9         ?
               (abs(with1->moment) + '7')  : 
               (abs(with1->moment) + '0')  ; 
            if (debug && with1->lorentz > 0) 
               res->scr[k - 1-1][y-1]         = 
               with1->lorentz > 9             ? 
               (with1->lorentz + '7')      :
               (with1->lorentz + '0')      ;  
         } 
      } 
      if (with1->partcl != prtclbase[with1->partcl-1].anti) 
         res->scr[k-1][y-1] = arrow;
   }
}   /* The end of Hline's body */ 


static void vline(int x,int ys,int ye,int ivr,int slt,  int txt) 
{int     linech /*, arrow */; 
 int  i, j, k, l; 
 char     nm[P_NAME_SIZE]; 

   if (ys > ye) 
   {  --(ys);
      ++(ye);
   } 
   else 
   {  ++(ys);
      --(ye);
   } 
   k = (ys + ye) / 2;
    
   {edgeinvert *with1 = &diag.vertlist[ivr-1][slt-1]; 
      
      linech = (prtclbase[with1->partcl-1].spin % 2) == 0 ?
               chrB3 : chrBA; 
      l = left ? prtclbase[with1->partcl-1].anti : with1->partcl; 

      if (!(((IN_PRTCL|OUT_PRTCL)&with1->prop)   ||  l < prtclbase[l-1].anti)) 
         l = prtclbase[l-1].anti; 
      strcpy(nm,prtclbase[l-1].name); 
      l = strlen(nm); 
      if (ye > ys) 
         for (i = ys; i <= ye; i++) res->scr[x-1][i-1] = linech; 
      else 
         for (i = ye; i <= ys; i++) res->scr[x-1][i-1] = linech; 
      if (txt) 
      { 
         if (left) 
         { 
            j = x - 1 - l; 
            for (i = 1; i <= l; i++) res->scr[i + j-1][k-1] = nm[i-1]; 
            if (debug) 
            { 
               if ((left ? 1 : 0) == 
                   (with1->moment > 0 ? 1 : 0)) 
               { 
                  res->scr[x + 1-1][k-1] = '-'; 
                  j = 2; 
               } 
               else 
                  j = 1; 
               res->scr[j + x-1][k-1] = chp;
               res->scr[j + x + 1-1][k-1]       = 
                  abs(with1->moment) > 9        ? 
                  (abs(with1->moment) + '7') : 
                  (abs(with1->moment) + '0') ; 
               if (with1->lorentz > 0)
                  res->scr[x-1][k + 1-1]        = 
                     with1->lorentz > 9         ?
                     (with1->lorentz + '7')  :  
                     (with1->lorentz + '0')  ; 
            } 
         } 
         else 
         { 
            j = x + l + 1; 
            for (i = 1; i <= l; i++) res->scr[j - i-1][k-1] = nm[i-1]; 
            if (debug) 
            { 
               if ((left ? 1 : 0) == 
                   (with1->moment > 0) ? 1 : 0) 
               { 
                  res->scr[x - 1-1][k-1] = '-'; 
                  j = 2; 
               } 
               else 
                  j = 1; 
               res->scr[x - j-1][k-1] = chp; 
               res->scr[x - j - 1-1][k-1]       = 
                  abs(with1->moment) > 9        ? 
                  (abs(with1->moment) + '7') : 
                  (abs(with1->moment) + '0') ; 
               if (with1->lorentz > 0) 
                  res->scr[x-1][k + 1-1]        = 
                     (with1->lorentz > 9)       ?
                     (with1->lorentz + '7')  :
                     (with1->lorentz + '0')  ; 
            } 
         } 
/*       if (partcl != prtclbase[partcl-1].anti) scr[x-1,k-1] = arrow; */ 
      }        
   } 
}   /* The end of Vline's body */ 


static void mline(int xs,int ys,int xe,int ye,
                 int ivr,int slt,int hor)
{int lg; 

   if (xe == xend)    
   {intsect *with1 = &res->ints[res->mver + 1-1]; 

      ++(res->mver);        
   /*  Range check error for Z->A,A,A,e1,E1 process  !  */ 
      with1->y = ye; 
      with1->xst = xs + 1; 
      with1->sp = prtclbase[diag.vertlist[ivr-1][slt-1].partcl-1].spin; 
      with1->selfvert = ivr; 
      with1->selfslot = slt; 
      {vertlink *with2 = &diag.vertlist[ivr-1][slt-1].nextvert;
         with1->aliasvert = with2->vno; 
         with1->aliasslot = with2->edno; 
      } 
   } 
         
   lg = xs == xe ? 1 : 0; 
   if (lg || !hor) 
      vline(xs,ys,ye,ivr,slt,lg); 
   else 
      if ((ys == ye || hor) && xe != xend) 
         hline(xs + 1,xe - 1,ys,ivr,slt); 
   if (ys != ye && !lg) 
   { 
      if (hor) 
      { 
         cr(3,xe,ys,ivr,slt); 
         vline(xe,ys,ye,ivr,slt,0); 
      } 
      else 
      { 
         if (ye > ys) 
            cr(1,xs,ye,ivr,slt); 
         else 
            cr(4,xs,ye,ivr,slt); 
         if (xe != xend) 
            hline(xs + 1,xe - 1,ye,ivr,slt); 
      }
   }    
}  /* The end of Mline's body */ 


static int howin(int i)
{int k, r=0; 

   for (k = 1; k <= diag.valence[i-1]; k++) 
      if (IN_PRTCL & diag.vertlist[i-1][k-1].prop) ++(r); 
   return r; 
}   /* The end of Howin's body */ 


static int wherein(int i)
{int  k; 

   k = diag.valence[i-1]; 
   while (k > 0 && !(IN_PRTCL & diag.vertlist[i-1][k-1].prop)) --(k); 
   return k; 
}  /* The end of Wherein's body */ 


static void lot(int* ylo,int* yhi,int ivr,int slt,int hor)
{int    iv,/* sl, */ tn, tp, yl, yh, l; 
 int flg; 

   *ylo = 0; 
   *yhi = 0; 
   {edgeinvert *with1 = &diag.vertlist[ivr-1][slt-1]; 
      if ((IN_PRTCL|OUT_PRTCL)&with1->prop)  return;
      iv = with1->nextvert.vno; 
   } 
      
   tp = howin(iv); 
   flg = 1; 
   tn = diag.valence[iv-1]; 
   if (tp == 0  && (tn == 4 || (hor && tn == 3  && 
       (OUT_PRTCL&diag.vertlist[iv-1][1].prop) && 
       (OUT_PRTCL&diag.vertlist[iv-1][2].prop)))) 
   { 
      lot(&yl,&yh,iv,tn,1); 
      *ylo += yl + yh + 1; 
      --(tn); 
   } 
   for (l = tn; l >= 2; l--) 
   if (!(IN_PRTCL&diag.vertlist[iv-1][l-1].prop)) 
   { 
      lot(&yl,&yh,iv,l,flg); 
      if (flg) 
      { 
         flg = 0; 
         *ylo += yl; 
         *yhi += yh; 
      } 
      else 
         *yhi += yl + yh + 1; 
   } 
}   /* The end of Lot's body */ 

//static int deep=0;

static void painter(int xc,int yc,int ivr,int slt,int hor)
{int    tp, l, xcnt, ycnt, xx, yy, hh, 
         ylo, yhi, start, ks, kk, pp; 
 int flag; 
/*
char buff[10];
for(l=0;l<deep;l++) buff[l]=' ';
buff[l]=0; deep++;

printf("New vertex vrt=%d slt=%d \n",ivr,slt);
*/
   flag = 0; 
   start = 2; 
   xcnt = xc; 
   ycnt = yc; 
   tp = howin(ivr); 
   ks = wherein(ivr); 
   kk = ks == diag.valence[ivr-1] ?  
        diag.valence[ivr-1] - 1   : 
        diag.valence[ivr-1]       ; 
   lot(&ylo,&yhi,ivr,kk,hor); 
   
//printf("%s vert=%d  kk=%d hor=%d ylo=%d yhi=%d\n",buff,ivr,kk,hor,ylo,yhi);
   
   pp = ks == kk - 1 ? kk - 2 : kk - 1; 
   if (tp == 2) 
   { 
      if (ylo < 2) 
         ycnt += ystep; 
      else 
      {  ycnt += ylo * ystep; 
         yc += (ylo - 1) * ystep;
      } 
      mline(xstart,yc,xcnt,ycnt,ivr,slt,1); 
      mline(xstart,ycnt,xcnt,ycnt,ivr,ks,1); 
   } 
   else 
      if (tp == 1) 
      { 
         if (yc == ystart) ycnt += ylo * ystep; 
         if (hor && diag.valence[ivr-1] == 4) 
         { 
            lot(&hh,&yy,ivr,pp,hor); 
//printf("%s tp=1, ivr=%d pp=%d hh=%d yy=%d\n",buff, ivr,pp,hh,yy);            
            ycnt += (yhi + hh + 1) * ystep; 
            start = 1; 
         } 
         mline(xstart,ycnt,xcnt,ycnt,ivr,ks,1); 
      } 
      else 
         if (diag.valence[ivr-1] == 4) 
         { 
            start = 1; 
            if (!hor) cr(1,xcnt++,ycnt,ivr,slt); 
         } 
         else 
            if (hor && (OUT_PRTCL&diag.vertlist[ivr-1][1].prop) 
                    && (OUT_PRTCL&diag.vertlist[ivr-1][2].prop)) 
               start = 1; 

//printf("%s first start=%d\n",buff,start);

   res->scr[xcnt-1][ycnt-1] = vertsign; 
   if (xcnt > res->mhor) 
      res->mhor = xcnt; 
   for (l = kk; l >= 2; l--) 
   {edgeinvert *with1 = &diag.vertlist[ivr-1][l-1]; 
      if (!(IN_PRTCL & with1->prop)) 
      { 
         if (flag) 
{
//            lot(&ylo,&yhi,ivr,l,start != 3); 

             
                         lot(&ylo,&yhi,ivr,l,start==2);
//printf("%s start=%d  vrt=%d l=%d   hor=%d ylo=%d yhi=%d\n",buff, start, ivr,l, start==2, ylo,yhi);  
}          
         else 
            flag = 1; 
         xx = (OUT_PRTCL&with1->prop) ? xend : 
              start == 3 ? xcnt : xcnt + xstep; 
         {vertlink *with2 = &with1->nextvert; 
            if (start == 1) 
            { 
               
//               lot(&hh,&yy,ivr,pp,1); 
                 lot(&hh,&yy,ivr,l-1,1);
//printf("%s statr1 vrt=%d l-1=%d hh=%d yy=%d\n",buff, ivr, l-1,hh,yy);
               yy = ycnt - (hh +1) * ystep;
               mline(xcnt,ycnt,xx,yy,with2->vno,with2->edno,0); 
               if (!(OUT_PRTCL & with1->prop))  painter(xx,yy,with2->vno,with2->edno,1); 
            } 
            else 
               if (start == 2) 
               { 
                  lot(&hh,&yhi,ivr,l,1);
                  mline(xcnt,ycnt,xx,ycnt,with2->vno,with2->edno,1); 
                  if (!(OUT_PRTCL&with1->prop)) painter(xx,ycnt,with2->vno,with2->edno,1); 
                  yy = yhi;
//                  yy=ylo; 
               } 
               else 
               { 
//printf("%s start=? %d yy=%d  ylo=%d\n",buff,start,yy,ylo);
                  yy += ylo + 1;
//                    yy+=yhi+1;
//printf("%s down %d\n",buff,yy);
                  yy = ycnt + yy * ystep; 
                  
                  mline(xcnt,ycnt,xx,yy,with2->vno,with2->edno,0); 
                  if (!(OUT_PRTCL&with1->prop)) painter(xx,yy,with2->vno,with2->edno,0); 
               } 
         } 
         ++(start); 
      } 
   } 
}   /* The end of Painter's body */ 


static void diaprin(dscreen* res1)
{  int i, j; 

   res = res1; 
   for (i = 0; i < xend; i++) for (j = 0; j < yend; j++)  res->scr[i][j] = chsp; 
   res->mver = 0; 
   res->mhor = 0; 
   res->mxy  = 0; 

   hcontr = 1; 
   if (left)  painter(xstart + xstep,ystart,1,1,1); 
        else  painter(xstart + xstep,ystart,diag.sizel + 1,1,1); 
   res->mhor += xstep - 1; 
   hcontr = 0; 

// outgoing particles 
   for (i = 0; i < res->mver; i++)     
   {  intsect *with1 = &res->ints[i]; 
      hline(with1->xst,res->mhor,with1->y,with1->selfvert,with1->selfslot); 
   }

}   /* The end of Diaprin's body */ 


static int lh(int y)
{   
   if (is[jsect-1][y-1] == chsp) 
      is[jsect-1][y-1]  = ssect ? chrC4 : chrCD; 
   else 
      if (is[jsect-1][y-1] == chrB3) 
         is[jsect-1][y-1]  = ssect ? chrC5 : chrD8; 
      else 
         if (is[jsect-1][y-1] == chrBA) 
            is[jsect-1][y-1]  = ssect ? chrD7 : chrCE; 
         else 
            return 1;
   return 0; 
}  /* The end of Lh's body */ 


static void lv(void)
{int c; 

   c = is[isect-1][jsect-1];
   
   is[isect-1][jsect-1]           = 
      c == chsp                   ? 
         (ssect ? chrB3 : chrBA)  : 
      c == chrC4                  ? 
         (ssect ? chrC5 : chrD7)  : 
      c == chrCD                  ? 
         (ssect ? chrD8 : chrCE)  : c; 
}   /* The end of Lv's body */ 


static int isector(void)
{int  k, y1, y2, vv, ss, pp=0, tt=0, sh=0; 
 int count=0;

   /* Nested function: lh */ 
   /* Nested function: lv */ 

   k = resl.mver; 
   for (isect = 1; isect <= k; isect++) 
   { 
      {intsect *with1 = &resl.ints[isect-1]; 
         y1 = with1->y; 
         vv = with1->aliasvert; 
         ss = with1->aliasslot;    
      } 
      jsect = 1; 
      while (resr.ints[jsect-1].selfvert != vv || 
             resr.ints[jsect-1].selfslot != ss) 
         ++(jsect); 
      {intsect *with1 = &resr.ints[jsect-1]; 
         if (y1 == with1->y) 
            ++(sh); 
         else 
            if (y1 == with1->y + ystep) 
               ++(pp); 
            else 
               if (y1 == with1->y - ystep) 
                  ++(tt); 
      } 
   }
   sh = tt > sh ? - ystep : pp > sh ?  ystep : 0; 
   goto b;

a: count++;
    if(sh==0) sh=-1;   else {if(sh<0) sh=-sh; else sh=-sh-1;}
//   if (sh < 0) 
//      ++(sh); 
//   else 
//      --(sh); 

b: pp = sh < 0 ? resl.mxy - sh : resr.mxy + sh; 
   if (resr.mxy > pp) pp = resr.mxy; 
   if (resl.mxy > pp) pp = resl.mxy; 
   ++(pp); 
   
   for (isect = 1; isect <= k; isect++) 
      for (jsect = 1; jsect <= pp; jsect++) 
         is[isect-1][jsect-1] = chsp; 

   for (isect = 1; isect <= k; isect++) 
   {  
       
      {intsect *with1 = &resl.ints[isect-1]; 
         y1 = with1->y; 
         ssect = with1->sp % 2 == 0 ? 1 : 0; 
         vv = with1->aliasvert; 
         ss = with1->aliasslot;    
      }    
      jsect = 1; 
      while (resr.ints[jsect-1].selfvert != vv || 
             resr.ints[jsect-1].selfslot != ss) 
         ++(jsect); 
          
      y2 = resr.ints[jsect-1].y; 
      if (sh < 0) 
         y1 -= sh; 
      else 
         y2 += sh; 
      if (y1 == y2) 
         for (jsect = 1; jsect <= k; jsect++) 
            tt = lh(y1); 
      else 
         if (y1 < y2) 
         { 
            for (jsect = 1; jsect <= isect - 1; jsect++) 
               if (lh(y1) == 1) 
                  {if(count<5)goto a;}
            if (is[isect-1][y1-1] != chsp) 
               {if(count<5)goto a;}
            is[isect-1][y1-1] = ssect ? chrBF : chrBB; 
            for (jsect = y1 + 1; jsect <= y2 - 1; jsect++) 
               lv(); 
            if (is[isect-1][y2-1] != chsp) 
               {if(count<5)goto a;}
            is[isect-1][y2-1] = ssect ? chrC0 : chrC8; 
            for (jsect = isect + 1; jsect <= k; jsect++) 
               if (lh(y2) == 1) 
                  {if(count<5)goto a;}
         } 
         else 
         { 
            for (jsect = 1; jsect <= isect - 1; jsect++) 
               if (lh(y1) == 1) 
                  {if(count<5)goto a;}
            if (is[isect-1][y1-1] != chsp) 
               {if(count<5)goto a;}
            is[isect-1][y1-1] = ssect ? chrD9 : chrBC; 
            for (jsect = y2 + 1; jsect <= y1 - 1; jsect++) 
               lv(); 
            if (is[isect-1][y2-1] != chsp) 
               {if(count<5)goto a;}
            is[isect-1][y2-1] = ssect ? chrDA : chrC9; 
            for (jsect = isect + 1; jsect <= k; jsect++) 
               if (lh(y2) == 1) 
                  {if(count<5)goto a;}
         } 
   } 
   return sh; 
}  /* The end of Isector's body */ 


static void diaprt(pnt* fout)
{int  i, j, r; 
 int  sh, p, k, kk; 
 pnt      po, pn = NULL; 

   left = 1; 
   diaprin(&resl);
   left = 0; 
   diaprin(&resr); 


//   sh = isector(); 
sh=0;   
//printf("isector ok  sh=%d  \n",sh);
   kk = scrwidth - resl.mhor - resl.mver - resr.mhor - 2; 
   if (kk < 0) kk = 0; 
   kk = kk / 2; 
   r = sh < 0 ? resl.mxy - sh : resr.mxy + sh; 
   if (resr.mxy > r) r = resr.mxy; 
   if (resl.mxy > r) r = resl.mxy; 

//for(j=0;j<r;j++) for(i=0;i<resl.mver;i++)  is[i][j]=' ';

   
   p = 68 - 1 - kk - resl.mhor - resl.mver - resr.mhor;   
/* Pukhov   80 -> 68  */ 
   ++(r); 
   for (j = 1; j <= r; j++) 
   { 
      po = pn; 
      pn=(pnt)m_alloc(sizeof(spool));
      if (j == 1) 
         *fout = pn; 
      else 
         po->add = pn;
      k = kk; 
      if (k > 0) memset(pn->stn,chsp,k);
      if (sh < 0) 
         if (j <= -sh)
         {
            memset(pn->stn + k,chsp,resl.mhor);
            k += resl.mhor;
         }
         else 
         {
            for (i = 0; i < resl.mhor; i++) 
               pn->stn[k + i] = resl.scr[i][j + sh-1]; 
            k += i;
         }
      else 
         if (j <= r)
         { 
            for (i = 0; i < resl.mhor; i++) 
               pn->stn[k + i] = resl.scr[i][j-1];
            k += i;
         }
         else
         {  
            memset(pn->stn + k,chsp,resl.mhor);
            k += resl.mhor;
         }
      for (i = 0; i < resl.mver; i++) 
         pn->stn[k + i] = ' ';  //   is[i][j-1];
      k += i;
/*      
#ifdef BLOCK_GRAPHICS
      pn->stn[k++] = '|';
#else
      pn->stn[k++] = '!';
#endif
*/
      pn->stn[k++] = ' ';
      if (sh > 0) 
         if (j <= sh)
         {
            memset(pn->stn + k,chsp,resr.mhor);
            k += resr.mhor;
         }
         else
         {
            for (i = 0; i < resr.mhor; i++) 
               pn->stn[k + i] = resr.scr[resr.mhor - i-1][j - sh-1];
            k += i;
         }
      else 
         if (j <= r) 
         {
            for (i = 0; i < resr.mhor; i++) 
               pn->stn[k + i] = resr.scr[resr.mhor - i-1][j-1];
            k += i;
         }
         else 
         {
            memset(pn->stn + k,chsp,resr.mhor);
            k += resr.mhor;
         }
      if (p > 0) 
      {  memset(pn->stn + k,chsp,p);
         k += p;
      }
       pn->stn[k] = '\0';
   } 
   pn->add = NULL; 
}  /* The end of Diaprt's body */ 



void writeAmDiagram(vampl* diag1,int label,char comment,FILE* outfile)
{  int i,j, i2,j2;

   diag.sizel=diag1->size;
   diag.sizet=diag.sizel+ (diag1->outno-1)/MAXVALENCE+1;
   for(i=0;i<diag.sizel;i++)
   { diag.valence[i]=diag1->valence[i];
     for(j=0;j< diag.valence[i]; j++) diag.vertlist[i][j]=diag1->vertlist[i][j];
   }

   
   for(i=0,i2=diag1->size,j2=0;i<diag1->outno;i++)
   { 
     int i1=diag1->outer[i].vno;
     int j1=diag1->outer[i].edno;

     if(j2>=MAXVALENCE) {i2++;j2=0;}
     
     diag.valence[i2]=j2+1;     
     diag.vertlist[i2][j2]=diag.vertlist[i1][j1];

     diag.vertlist[i2][j2].moment=-diag.vertlist[i2][j2].moment;
     diag.vertlist[i2][j2].partcl= prtclbase[diag.vertlist[i2][j2].partcl-1].anti;
     diag.vertlist[i1][j1].nextvert.vno=i2+1;
     diag.vertlist[i1][j1].nextvert.edno=j2+1;
     diag.vertlist[i2][j2].nextvert.vno=i2+1;
     diag.vertlist[i2][j2].nextvert.edno=j2+1;
     j2++;
   }
   
   debug=label;
   left = 1;
   diaprin(&resl); 
   for(i=0;i<yend;i++)
   {char ss[xend+1];
       
      for(j=0;j<xend;j++) ss[j]=resl.scr[j][i];
      ss[xend]=0;
      for(j=xend-1; j>=0 &&(ss[j]==' ');j--) ss[j]=0;
    
      if(ss[0]) fprintf(outfile,"%c%s\n",comment,ss);
   }
} 


void writeTextDiagram(vcsect* diag1,int label,char comment,FILE* outfile)
{  pnt  pp, pict; 
//   outfile=stdout;
// deep=0;
   diag= *diag1;
   debug=label;
   diaprt(&pp); 
   while (pp != NULL) 
   { 
     pict = pp; 
     pp = pict->add; 
     fprintf(outfile,"%c%s\n",comment,pict->stn); 
     free(pict);
   } 
} 
