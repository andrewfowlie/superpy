/*
 Copyright (C) 1997,2006,  Alexander Pukhov 
*/

#include "chep_crt.h"
#include "plot.h"
#include "num_serv.h"
#include "simpson.h"
#include "cs_22.h"
#include "interface.h"
#include "subproc.h"
#include "mc_menu.h"
#include "err_code.h"
#include "mc_menu.h"
//#include "cut.h"
#include "kinaux.h"
#include "param.h"
#include "alphas2.h"
#include "usrfun.h"
#include "nType.h"
#include "rw_sess.h"

static REAL pvect4[16];

static double         totcoef;
static double         cos1, cos2;
static double         eps=0.001;
static int            recalc;
static char           procname[STRSIZ];

static REAL pmass[4];
static REAL pRestIn,pRestOut;


static void  infotext(void)
{
	goto_xy(1,4); 
   scrcolor(Red,BGmain);
	print(" %-53s\n","P(c.m.s.)    :");
	
	print(" Cos(p1,p3): min=                   max=            \n");
	print(" %-53s","Cross Section: ");
   scrcolor(FGmain,BGmain);
}


static void  writeinformation(void)
{  double Pcm=Pcm22;
   scrcolor(FGmain,BGmain);
   goto_xy(18,4); print("%12f [GeV]    ",Pcm);   /*  Energy  */
   goto_xy(18,5); print("%8.6f",cos1);
   goto_xy(42,5); print("%8.6f",cos2);
}


static void  calccoef(void)
{int i;
   REAL  lambda12, lambda34, s_, ms, mdiff, sqrt_S;
   
   pRestIn=Pcm22;
   err_code = 0;
   
   for(i=0;i<4;i++) pinf_int(Nsub,i+1,pmass+i,NULL);
   
   sqrt_S=sqrt(pmass[0]*pmass[0]+pRestIn*pRestIn)
         +sqrt(pmass[1]*pmass[1]+pRestIn*pRestIn);    

   s_=sqrt_S*sqrt_S; 
   
   lambda12=2*sqrt_S*pRestIn;
   
   ms = pmass[2] + pmass[3];
   if (ms >= sqrt_S) goto errorexit;
   mdiff=pmass[2] - pmass[3];
   lambda34 = sqrt((s_ - ms*ms) * (s_ - mdiff*mdiff));

   pRestOut=lambda34/(2*sqrt_S);
   
   totcoef = 3.8937966E8 * lambda34 /(32.0 * M_PI * lambda12 * s_);

   for(i=0;i<16;i++)pvect4[i]=0;
   
   pvect4[3] = pRestIn; 
   pvect4[7] =-pRestIn;
   pvect4[0] = sqrt(pRestIn*pRestIn + pmass[0]*pmass[0]);
   pvect4[4] = sqrt(pRestIn*pRestIn + pmass[1]*pmass[1]);
   pvect4[8] = sqrt(pRestOut*pRestOut + pmass[2]*pmass[2]);
   pvect4[12]= sqrt(pRestOut*pRestOut + pmass[3]*pmass[3]);

   err_code = 0;
   return;

errorexit:
   if (err_code == 0)  err_code = 4;
}


static void  calcscalars(double  cos_f)
{
   REAL sin_f=sqrt(fabs((1-cos_f)*(1+cos_f)));
   pvect4[11]=pRestOut*cos_f;
   pvect4[15]=-pvect4[11];
   pvect4[10]=pRestOut*sin_f;
   pvect4[14]=-pvect4[10];
} 


static double  cross_section(double  x)
{ double  r,GG,qF,qR,pvect4_[16];
  int i; 
  
  calcscalars(x);
  for(i=0;i<16;i++) pvect4_[i]=pvect4[i];
 
  Scale(pvect4_,&qF,&qR);
  GG=sqrt(4*M_PI*alpha_2(qR));
      
  r = sqme_int( Nsub,GG,pvect4,&err_code)*usrFF(2,2,pvect4_,p_names,p_codes);
  if (err_code != 0) return 0; 
  return r * totcoef; 
} 


static int fillseq(int  n,double * f)
{

 int        i;
 double      step;
 
   err_code = 0;
   step = (cos2 - cos1) / (n - 1);
   for (i = 0; i < n; i++)
   {
      f[i]=cross_section(cos1+i*step);
      if (err_code > 1)  return 0;
   }
   return 1; 
}


static void  drawgraph(void)
{
 int         n=101;
 double  f[202]; 
   
  calccoef();
  if(err_code) { messanykey(10,10,"Error in kinematics"); return; }
          
  do
  {  if (correctInt(56,8,"Number of points=",&n,1))
     {
        if (n < 3) messanykey(56,8,"Too few points");
        if (n > 201) messanykey(56,8,"Too many points");
     }
      else return;
  }  while (n < 3 || n > 201 );

           
  if( !fillseq(n,f)) 
  {  messanykey(10,10,"Error in calculation"); 
     return;
  }

  plot_1(cos1,cos2,n,f,NULL,procname,"cos(p1,p3)", "Diff. cross section [pb]");

}


static double  totcs(void)
{  double  int_val = 0.0;
   calccoef();
   if (err_code == 0) int_val=gauss345(cross_section,cos1,cos2,eps,&err_code);
   return int_val;
}


static void  total_cs(void)
{  double  totcs;

   goto_xy(18,6); 
   scrcolor(FGmain,BGmain);
   print("?                            "); 
   goto_xy(18,6);
   refresh_scr();
   calccoef(); 
   if (err_code ) print("incorrect"); else 
   {
      totcs= gauss345(cross_section, cos1,cos2,eps,&err_code); 
      if (err_code<=0 ) { print("%-G [pb]",totcs);} 
      if(err_code==1)  print("?");
   }
}


static double vtotcs(void)
{
  double cs=totcs();
  double m1q=pmass[0]*pmass[0];
  double m2q=pmass[1]*pmass[1];
  double pq=pRestIn*pRestIn;
        
  double E1=sqrt(m1q+pq);
  double E2=sqrt(m2q+pq);
  double pij=E1*E2+pq;
  double v=sqrt(pij*pij -m1q*m2q)/(E1*E2);
               
  return v*cs;
}
                  


int  cs_numcalc(double Pcm)
{
   int  k,l;
   void * pscr0=NULL;
   void * pscr = NULL;
    
   get_text(1,3,60,11,&pscr0); 
   k=Nsub;   
   sprintf(procname,"%s,%s ->%s,%s",pinf_int(k,1,NULL,NULL),
           pinf_int(k,2,NULL,NULL),pinf_int(k,3,NULL,NULL),pinf_int(k,4,NULL,NULL));              

   Pcm22=Pcm;
   
    cos1=-0.999;
    cos2= 0.999;

   infotext();
   writeinformation();
   k = 1;
   l = 1;
   
   recalc = 1;
   do
   {  char menuTxt[]="\030"
         " Change parameter       "
         " Set precision          "
         " Cos13(min) = cosmin    "
         " Cos13(max) = cosmax    "
         " Angular dependence     "
         " Parameter dependence   "
         " sigma*v plots          ";

      if (recalc)
      {  
         total_cs();
         recalc = 0;
        if (err_code) errormessage();
        if(err_code==3) err_code=0;                  
      }

      improveStr(menuTxt,"cosmin","%.6f",cos1);
      improveStr(menuTxt,"cosmax","%.6f",cos2);
      
      menu1(54,4,"",menuTxt,"n_22_*",&pscr,&k);
         
      switch (k)
      {
         case 0:  break;
         case 1: if(change_parameter(54,5,1)) recalc=1; break;
         case 2: 
            do {   /* Precision */
               recalc = correctDouble(1,23," Enter precision : ",&eps,1);
               if (eps < 1.E-10 || eps > 0.0011)
                  messanykey(10,12,"Range check error");
            }  while (!(eps >= 1.E-10 && eps <= 0.03)); 
         break;
         case 3: recalc=correctDouble(15,10,"Min[cos(p1,p3)]=",&cos1,1); break;
         case 4: recalc=correctDouble(15,10,"Max[cos(p1,p3)]=",&cos2,1); break;
         case 5: if(err_code>1)  errormessage(); else drawgraph();    break;
         case 6: paramdependence(totcs,procname,"Cross Section [pb]"); break;
         case 7: paramdependence(vtotcs,procname,"v*sigma[pb]"); break;
      }  /*  switch  */
      if (k > 0) writeinformation();
   }  while (k != 0);
   put_text(&pscr0);
   return 0;
}
