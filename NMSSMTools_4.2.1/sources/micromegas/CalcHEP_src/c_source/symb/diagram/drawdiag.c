/*
 Copyright (C) 1997, Victor Edneral 
*/
#include "chep_crt.h"
#include <math.h>
#include "syst2.h"
#include "physics.h"
#include "process.h"
#include "drawdiag.h"

   typedef struct intsect {
           int      y1, y2;
           int      nline;
                          } intsect;

   typedef struct knot {
           int      ty, pt, fr, e1, e2, e3, n1, n2, n3;
                       } knot;

static   int      xn1, xn2; 
static   int      upr, ur, ni;
static   int      ys, yh, bh, bv;

static   int      x[4];
static   int      y[9];

static   int      x1,x2;  /* From mline */
static   int      np;     /* From mline */
static   knot     kn[4];  /* From picture */
static   int      kt, kl, kk, ks, nk, se, se1, ku,
                  tpc, pext; /* From picture */
static   intsect  intar[5]; /* For mline */
static   int      sflag = 0;
static   int cHeight, cWidth;

static void mtg_line(int x1,int y1,int x2,int y2,int arrowflag)
{
   switch (arrowflag)
   { case  0: tg_line(x1,y1,x2,y2);      break;
     case  1: tg_arrowline(x1,y1,x2,y2); break;
     case -1: tg_arrowline(x2,y2,x1,y1);
   }   
}

static void  mouttextxy(int xm,int y)
{
  if(texflag) tg_outtextxy(xm,y, prtclbase[np-1].latex); else
/*  if(prtclbase[np-1].name[0]=='~')
  {   tg_outtextxy(xm,y,prtclbase[np-1].name+1);
      if(!prtclbase[np-1].name[2]) tg_outtextxy(xm,y-tg_textheight("H")/2,"~"); 
      else                         tg_outtextxy(xm,y-tg_textheight("H")/2,"~~");
  } else 
*/ tg_outtextxy(xm,y,prtclbase[np-1].name);
}  

static int rightPart=0;
static void  mline(int x1_,int y1,int x2_,int y2,
        int r,int np_,int mm)
{
   int      h, y, xc, yc;
   int         sp, lt=SolidLn, uk;
   int arrowflag;
   int arrowshift=3;
   int nameOut=0;

   x1=x1_; x2=x2_; np=np_;
   
   sp = prtclbase[np-1].spin;
   switch(sp)
   { case 0: lt = DottedLn;   break;
     case 1: lt = SolidLn;    break;
     case 2: lt = DashedLn;   break;
     case 3: lt = SolidLn;   break;
     case 4: lt = DashedLn;   break;
   }
   tg_setlinestyle(lt,NormWidth);
   
   arrowflag =0 ;  /*Shichanin-insert*/

   if (r >= 0)
   {
      xc = (x1 + x2) / 2;
      yc = (y1 + y2) / 2;
      uk = (prtclbase[np-1].anti < np) ;
      if (prtclbase[np-1].anti != np ) { if(uk) arrowflag= -1; else arrowflag=1;}
            
      if (uk && (x1 == x2 || r == 1 || r == 3))
         np = prtclbase[np-1].anti;
      if (ur < 0)
      {
         h = x1; x1 = x2; x2 = h;
         y = y1; y1 = y2; y2 = y;
         if (r < 3) r = 2 - r;
      }
      if (r == 3)  /* unknown case */
      { 
         y = yc - 2;
         if (ur > 0) h = x1; else h = x2;
         tg_settextjustify(1 - ur,BottomText);
      }
      else
      {
         if (x1 != x2)
         switch(r)
         {  
            case 0: /*left end text */ 
               y = y1;
               h=x1-1;
               tg_settextjustify(RightText,CenterText);
nameOut=1;   
              break;
            case 1: /* horizontal line ? */  
                h = xc; y = yc-1;
                if(arrowflag) y=y-arrowshift;
                tg_settextjustify(CenterText,BottomText);
                break;
            case 2: /* right end text */
                 y=y2;
                  h=x2+1;
                  tg_settextjustify(LeftText,CenterText);
         } 
         else /* vertical line */
         { 
            switch(r)
            { case 1: /* unknown case */  
                   h = xc; 
                   y = yc + 4;
                   break;
                   
              case 0: /* line with left side label */
                   h = xc-1; 
                   if (arrowflag) h-=arrowshift;
                   y=yc;
                   tg_settextjustify(RightText,CenterText);
                   break;
                   
              case 2: /* line with right side label */
                   h = xc+1;
                   if (arrowflag) h+=arrowshift; 
                   y=yc;
                   tg_settextjustify(LeftText,CenterText);
                   break; 
            }
         }
      }
      if(texflag || !nameOut || !rightPart)mouttextxy(bh + h,bv + y);
   }

   mtg_line(bh + x1,bv + y1,bh + x2,bv + y2,arrowflag);
   if (x1 == x[3])
   {
      intar[mm-1].y2 = bv + y1;
      if (ur == 1) intar[mm-1].nline = lt;
   }
   else
      if (x2 == x[3])
      {
         intar[mm-1].y2 = bv + y2;
         if (ur == 1) intar[mm-1].nline = lt;
      }
}

static int  elong(int m)
{int     el;

   if (m < 0)
      el = 1;
   else
   {  knot *with1 = &kn[m];
      el = elong(-with1->e1) + elong(-with1->e2);
   }
   return el;
}

static int     elongl(int m)
{ int     k, el;

   if (m == se || m == se1)
      el = 0;
   else
      if (m < 0)
         el = 100;
      else
      {  knot *with1 = &kn[m];
         k = MIN(elongl(-with1->e1),elongl(-with1->e2));
         if ((with1->e3 < 0)) k = MIN(k,elongl(-with1->e3));
         el = k + 1;
      }
   return el;
}

static void  order(int* e1,int* e2,int* m1,int* m2)
{ int     e;

      if ((*e2 > 0 && *e1 < 0) || (*e1 < 0 && *e2 < 0 &&
          elong(-(*e1)) > elong(-(*e2))))
      { e = *e1; *e1 = *e2; *e2 = e; e = *m1; *m1 = *m2; *m2 = e; }
}

static int  empt(int l)
{ int     k, em;

   if (l < 0)
      em = 0;
   else
   {  knot *with1 = &kn[l];
      k = empt(-with1->e1) + empt(-with1->e2);
      if (with1->pt == 0)
         em = k + 1;
      else
         em = k;
   }
   return em;
}


static void  triplet(int l)
{ int  l1, l2, l3, ll, m1, m2, m3, k;

   if (kn[l].pt == 0)
   {
      k = kn[l].fr;
      l1 = kn[k].e1; l2 = kn[l].e1; l3 = kn[l].e2;
      m1 = kn[k].n1; m2 = kn[l].n1; m3 = kn[l].n2;
      if (l1 == -l)
      {  l1 = kn[k].e2;
         m1 = kn[k].n2;
      }
      else
         if (tpc == 2)
         {   if (se == k)
               se1 = l;
            else
               if (se == l)
               {  ll = l1;
                  l1 = l2;
                  l2 = ll;
                  ll = m1;
                  m1 = m2;
                  m2 = ll;
                  se1 = k;
               }
         }    
      order(&l2,&l3,&m2,&m3);
      if (tpc == 1 || (se != l && se != k))
      {  order(&l1,&l2,&m1,&m2);
         order(&l2,&l3,&m2,&m3);
      }

      if (l3 > 0)
      {
         ll = 3;
         if ((se != l) && (se != k)) ks++;
      }
      else
         if (l2 > 0)
            ll = 2;
         else
            if (l1 > 0)
               ll = 1;
            else
               ll = 0;
      {  knot *with1 = &kn[k];
         with1->e1 = l1; with1->e2 = l2; with1->e3 = l3;
         with1->n1 = m1; with1->n2 = m2; with1->n3 = m3;
         with1->ty = ll;
      }

      {  knot *with1 = &kn[l];
         with1->e1 = l1; with1->e2 = l2; with1->e3 = l3;
         with1->n1 = m1; with1->n2 = m2; with1->n3 = m3;
         with1->ty = ll;
      }

      kn[l].fr = kn[k].fr; kn[l].pt = kn[k].pt;
   }
   else
   {  knot *with1 = &kn[l];
      l1 = with1->e1;
      l2 = with1->e2;
      l3 = with1->e3;
   }

   if (l1 < 0) triplet(-l1);
   if (l2 < 0) triplet(-l2);
   if (l3 < 0) triplet(-l3);
}


static void  ladd(int* l)
{
   (*l)++;
   if ((*l < 11) && (nk < 4)) return;
   print("  There is a bad diagram:");
   print("  Type ENTER, please:");
   getchar();
}


static void  filler(decayDiagram ar)
{int     nf=0, k, l=1;

   ni = 0; nk = 0; se = 100;

nn:if (nk > 3) ladd(&l);
   {  knot *with1 = &kn[nk];

      with1->pt = abs(ar[l-1]);
      ladd(&l); with1->fr = nf;
      with1->ty = 0; with1->e3 = 0;

ll:   if (ar[l-1] > 0)
      {
         if (tpc == 1 || se != 100) ni++;
         if (se == 100) se = nk;
         with1->ty++;
         if (with1->ty == 2)
         {  with1->e2 = ar[l-1];
            with1->n2 = ni;
            ladd(&l); goto pp;
         }
         with1->e1 = ar[l-1];
         with1->n1 = ni;
         ladd(&l);
         goto ll;
      }
      nf = nk;
      with1->e2 = -(++nk);
      goto nn;
   }

pp:for (k = nk; k >= 0; k--)
   {  knot *with1 = &kn[k];
      if (with1->ty == 0)
      {
         if (ar[l-1] > 0)
         {
            with1->e1 = ar[l-1];
            if (tpc == 1 || se != 100)
               with1->n1 = ++ni;
            if (se == 100) se = k;
            ladd(&l);
            with1->ty = 1;
         }
         else
         {  nf = k;
            with1->ty = -1;
            with1->e1 = -(++nk);
            goto nn;
         }
       }  
   }
   for (k = 0; k <= nk; k++)
   {  knot *with1 = &kn[k];
      if (with1->ty < 0) with1->ty = 0;
      if (tpc == 1 || k != se)
         order(&with1->e1,&with1->e2,&with1->n1,&with1->n2);
   }

   kl = elong(0) + 1;
   kt = empt(0);
   kk = kl - 2 - kt;
   if (tpc == 2) kn[se].e1 = prtclbase[kn[se].e1-1].anti;

   se1 = se;
   ks = 0;
   if (kt > 0) triplet(0);
}

static void  knot1(int l,int* yl,int* yyl)
{ int     m, m1, m2, m3, p, p1, p2, p3;
  int      z;

   if (l > 0) { m1 = l; m2 = 0; m3 = 0; p1 = pext; goto bb; }

aa:{ knot *with1 = &kn[-l];
      if ((tpc > 1) && ((l == -se) || (l == -se1)))
      {
         if (tpc == 3) mline(x[0],*yl,x[2],*yl,0,with1->e1,0);
         m1 = with1->e2; m2 = with1->e3; m3 = 0;
         p1 = with1->n2; p2 = with1->n3; 
      } 
      else 
      {  m1 = with1->e1; m2 = with1->e2; m3 = with1->e3; 
         p1 = with1->n1; p2 = with1->n2; p3 = with1->n3;
      } 
   } 
   if (ku == 1)
      {  ku = 0; 
         if (m1 < 0) { m = m1; m1 = m2; p = p1; p1 = p2; } 
         else        { m = m2; p = p2; }
         m2 = m3;
         p2 = p3; 
         if (m1 < 0) m3 = m1;
         if ((tpc == 3) && (elongl(-m3) > elongl(-m)))
            {
               if (m1 < 0) { m1 = m; p1 = p; }
               else        { m2 = m; p2 = p; }
               m = m3;
               p = p3; 
            }
         m3 = 0;

         { knot *with1 = &kn[-m]; 

               mline(x[2],*yl,x[2],*yl - ys,0,with1->pt,0);
               mline(x[2],*yl - ys,x[3],*yyl,2,with1->e1,with1->n1);
               *yyl += ys; 
               mline(x[2],*yl - ys,x[3],*yyl,2,with1->e2,with1->n2);
               *yyl += ys; 
         } 
      }


bb: if ((m1 > 0))
      {
         mline(x[2],*yl,x[3],*yyl,2,m1,p1);
         *yyl += ys;
         l = m2;
      }
   else l = m1; 

   if ((m2 > 0)) 
      { 
         mline(x[2],*yl,x[3],*yyl,2,m2,p2);
         *yyl += ys; 
         l = m3; 
      } 

   if ((m3 > 0)) 
         { 
            mline(x[2],*yl,x[3],*yyl,2,m3,p3);
            *yyl += ys;
         }

   if ((l >= 0)) return;

   { knot *with1 = &kn[-l];

         z = *yl + ys;
         if (tpc == 3 && (l == -se || l == -se1))  m = with1->ty - 1; 
         else m = with1->ty;
         if ((m > 1) && (z < *yyl)) z += ys; 
         mline(x[2],*yl,x[2],z,0,with1->pt,0);
         *yl = z;
   } 
      
   goto aa;
} 


static void  knot2(int* yl,int* yyl)
{ int     l, m, m1, m2, p1, p2;
  int      yz;

   l = 0;
   yz = *yl;
   { knot *with1 = &kn[-l];
      if (tpc == 2) { m1 = with1->e2; m2 = with1->e3;
                      p1 = with1->n2; p2 = with1->n3; } 
      else { m1 = with1->e1; m2 = with1->e2;
             p1 = with1->n1; p2 = with1->n2; } 
   }

aa: if (tpc == 3) 
      if ((m1 < 0) && (m2 < 0) && (elongl(-m1) < elongl(-m2))) 
         { m = m1; m1 = m2; m2 = m; }

   if ((tpc == 3) && ((l == -se) || (l == -se1))) 
      {
         mline(x[0],*yl,x[1],*yl,0,m1,0);
         m1 = m2; m2 = 0; p1 = p2; 
      } 

   if (m1 > 0) mline(x[1],yz,x[2],*yl,-1,m1,0);
   else mline(x[1],yz,x[2],yz,1,kn[-m1].pt,0);
   pext = p1;
   knot1(m1,&yz,yyl);

   if (m2 == 0) return;

   l = m2;
   yz += ys; 
   if ((yz + ys < *yyl) || ((yz < *yyl) && ((tpc == 1) ||
       (((l == -se) || (l == -se1)) && (kn[-l].e2 < 0))))) 
      yz += ys; 
   { knot *with1 = &kn[-l];
      
         if ((tpc < 3) && (with1->e2 > 0)) 
              mline(x[1],*yl,x[1],yz,-1,with1->pt,0);
         else mline(x[1],*yl,x[1],yz,0,with1->pt,0);
         *yl = yz; 
         m1 = with1->e1; p1 = with1->n1; 
         m2 = with1->e2; p2 = with1->n2; 
         if ((tpc < 3) && (m2 > 0))
            {  mline(x[1],yz,x[2],yz,1,with1->pt,0);
               knot1(l,yl,yyl);
               return;
            }
   }

   goto aa;
}


static void   draw(void)
{ 
   int      yl, yyl; 
   int     l; 

   if ((tpc == 2) && (se != 0) && (se1 != 0)) tpc = 3; 

   ku = 0; 

   if ((tpc == 2) && (kl == 3)) yyl = y[7];  else yyl = y[8];

   if (((tpc == 3) && (kl == elongl(0) + 3))) yl = y[8];
   else
   if (((tpc < 3) || ((kt > 0) && (ks == 0))) &&
       ((kn[0].ty == 0) || ((kn[0].e2 < 0) &&
       (kn[0].e3 < 0)) || ((tpc == 2) &&
       (kn[0].e2 < 0) && (kn[-(kn[0].e2)].ty == 0)))) 
      { yl = y[5]; ku = 1; }
   else yl = y[7]; 


   if ((tpc == 2) && (kl == 4) && (kk == 1))  { yyl = y[6]; yl = y[5]; } 

   if ((tpc == 2)) 
      if (((tpc == 2) && (kl == 4) && (kk == 1)) || 
          ((kn[0].e2 < 0) && (kn[0].e3 < 0)) ) 
         { 
            mline(x[0],y[7],x[2],yl,0,kn[0].pt,0);
            mline(x[0],y[3],x[2],yl,0,kn[0].e1,0);
            knot1(0,&yl,&yyl); 
         }
      else
         {
            mline(x[0],yl - yh,x[1],yl,0,kn[0].pt,0);
            mline(x[0],yl + yh,x[1],yl,0,kn[0].e1,0);
            l = kn[0].e2;
            if (l < 0)
            { 
               if(ku == 1) mline(x[1],yl,x[2],yl,3,kn[-l].pt,0);
               else        mline(x[1],yl,x[2],yl,1,kn[-l].pt,0);
               knot1(l,&yl,&yyl);
            } 
            else knot2(&yl,&yyl); 
         } 
   else 
   if (((tpc == 3) && ((kt == 0) || (ks != 0)) &&
       (kl != elongl(0) + 3)) || ((tpc == 1) && 
       (kn[0].e2 < 0) && (kn[-(kn[0].e2)].ty == 0))) 
      {
         mline(x[0],yl,x[1],yl,0,kn[0].pt,0);
         knot2(&yl,&yyl); 
      }
   else
      {
         mline(x[0],yl,x[2],yl,0,kn[0].pt,0);
         knot1(0,&yl,&yyl);
      }
}


static void  reverse(void)
{ int      xx;

   xx = x[0]; x[0] = x[3]; x[3] = xx;
   xx = x[1]; x[1] = x[2]; x[2] = xx; ur = -ur;
}


void  picture(int squared,void * buff, int x, int y)
{ int tpp;
  int q;
  bh=x;
  bv=y;

   ur=1; 
   tpc =nin;
   
   if (squared) {  upr=2; filler(((csdiagram*)buff)->dgrm1);}
   else         {  upr=1; filler((( adiagram*)buff)->dgrm0);}   
 
   if (tpc == 2 && kl == 4 && kk == 1)  tpp = 0; else tpp = tpc + 4;
 
   if (kl < tpp) {bv += ys; sflag = 1;}
   rightPart=0;
   draw();
 
   if (kl < tpp) {bv -= ys; sflag = 0;}
 
   if (squared)
   {  int dWidth= texflag? 0: cWidth; 
      tg_settextjustify(RightText,BottomText);
      for (q = 0; q < ni; q++) intar[q].y1 = intar[q].y2;
      tpc = nin;
      filler(((csdiagram*)buff)->dgrm2);
      if (tpc == 2 && kl == 4 && kk == 1)   tpp = 0; else tpp = tpc+4;
      if (kl < tpp) bv += ys;
      bh += xn2    -2*dWidth;
      reverse();
      rightPart=1;
      draw();
      if (kl < tpp) bv -= ys;
      reverse();
      bh -= xn2;
      for(q = 0; q < ni; q++)
      {
         tg_setlinestyle(intar[q].nline,NormWidth);
	 tg_line(bh+xn1+1 +6*dWidth ,intar[q].y1,
	           bh+xn2 +4*dWidth,intar[((csdiagram*)buff)->lnk[q]-1].y2);
      }
   }
}

void setPictureScale(int squared, int * xn,int *ynu)
{  int quant;
   int i;
   int len;
   int dWidth;
   int nout_;

   if(nout==1) nout_=2; else nout_=nout;
   cHeight=tg_textheight("H");
   cWidth =  tg_textwidth("H");
   
   if(texflag) dWidth=0; else dWidth=cWidth;   
      
   len=2*cWidth;
   quant=7*len/4;

   for(i=0;i<=3;i++) x[i]=4+len+i*quant +2*dWidth;

   if(squared)
   {  *xn=12+4*len+8*quant              +2*dWidth ;
      xn1=6+2*len + 3*quant;
      xn2=xn1+2*quant;
   }else  *xn=8+2*len+3*quant           +4*dWidth;
       
   quant=quant/2;

/*   *ynu=5+2*cHeight+2*MAX(3,nout)*quant; */
   *ynu=5+2*cHeight+2*nout_*quant;
   y[0]=*ynu-3-cHeight; 
   if(nout_==2) y[0]+=2*quant;
   for(i=1;i<=8;i++) y[i]=y[i-1] - quant;
      
   ys = 2*quant;
   yh = quant;    
}
