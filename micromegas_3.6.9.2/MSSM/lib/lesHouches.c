#include<math.h>
#include"pmodel.h"
#include"pmodel_aux.h"
#include"pmodel_f.h"
#include"../../sources/micromegas.h"

#include<ctype.h>
#include<stdio.h>
#include<stdarg.h>

#define alfEMZ  0.00781806

static double findX0(int NX,...)
{
   double X[20],Y[20],Xav,dX;
   int pow[12],maxPow,pos;
   int i,j,N;
   va_list ap;
   
   va_start(ap,NX);
    for(i=0;i<NX;i++){ X[i]=va_arg(ap, double); Y[i]=0;}
   va_end(ap);
                           
   for(i=0,Xav=0;i<NX;i++) {pow[i]=0; Xav+=X[i];}
   
   Xav/=NX;
   for(i=0,N=0;i<NX;i++)
   {
     for(j=0;j<N;j++) if(Y[j]==X[i]) { pow[j]++; break;}
     if(j==N) {Y[N]=X[i];pow[N]=1;N++;}  
   }
   
   pos=0;maxPow=pow[0];dX=fabs(Xav-Y[0]);
   
   for(i=1;i<N;i++)
   if(pow[i]>maxPow){ pos=i;maxPow=pow[i]; dX=fabs(Xav-Y[i]);}
   else if(pow[i]==maxPow && dX>fabs(Xav-Y[i])) {pos=i; dX=fabs(Xav-Y[i]);}
   
   return Y[pos];
}



int sugraLesH(char *fname,  double tb, double gMG1,double gMG2,double gMG3,
    double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
    double gMl2,double gMl3,double gMr2,double gMr3,
    double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3) 
{  double m0,mh,a0, m0_,mh_,a0_;
   double c=1.E-8;
   FILE*  f=fopen(fname,"w");
   if(f==NULL) return -1;
   
   fprintf(f,"Block MODSEL                 # Select model\n"   
           " 1    1                      # sugra\n"
           "Block SMINPUTS               # Standard Model inputs\n");
   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/alfEMZ);
   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5); 
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f," 7   %.8E       #  Mtau     \n",      findValW("Ml"));
   fprintf(f,"Block MINPAR                 # Input parameters\n");
   
   a0=findX0(3,gAl, gAt, gAb);
   mh=findX0(3,gMG1,gMG2,gMG3); 
   m0=findX0(17,gMHu,gMHd,gMl2,gMl2,gMl3,gMr2,gMr2,gMr3,gMq2,gMq2,gMq3,gMu2,gMu2,gMu3,gMd2,gMd2,gMd3);
   if(m0<0) m0=-m0;
   
   a0_=fabs(a0)*c;
   m0_=fabs(m0)*c;
   mh_=fabs(mh)*c;    

   fprintf(f," 1   %.8E       # m0\n",       m0); 
   fprintf(f," 2   %.8E       # m1/2\n",     mh);   
   fprintf(f," 3   %.8E       # tanb\n",     tb);
   fprintf(f," 4   %.0f     # sign(mu)\n", sgn);
   fprintf(f," 5   %.8E       # A0\n",       a0);


   fprintf(f,"Block EXTPAR\n");


   if(fabs(mh-gMG1)>mh_)  fprintf(f," 1   %.8E      # MG1\n",  gMG1);
   if(fabs(mh-gMG2)>mh_)  fprintf(f," 2   %.8E      # MG2\n",  gMG2);
   if(fabs(mh-gMG3)>mh_)  fprintf(f," 3   %.8E      # MG3\n",  gMG3);

   if(fabs(a0-gAt)>a0_)   fprintf(f," 11  %.8E      # At \n",  gAt);
   if(fabs(a0-gAb)>a0_)   fprintf(f," 12  %.8E      # Ab \n",  gAb);
   if(fabs(a0-gAl)>a0_)   fprintf(f," 13  %.8E      # Al\n", gAl);

   if(fabs(m0-gMHd)>m0_) fprintf(f," 21  %.8E      # MHd^2\n",gMHd*fabs(gMHd));
   if(fabs(m0-gMHu)>m0_) fprintf(f," 22  %.8E      # MHu^2\n",gMHu*fabs(gMHu));

   if(fabs(m0-gMl2)>m0_) {fprintf(f," 31  %.8E      # Ml1\n",  gMl2);
                          fprintf(f," 32  %.8E      # Ml2\n",  gMl2);}
   if(fabs(m0-gMl3)>m0_)  fprintf(f," 33  %.8E      # Ml3\n",  gMl3);
   if(fabs(m0-gMr2)>m0_) {fprintf(f," 34  %.8E      # MR1\n",  gMr2);
                          fprintf(f," 35  %.8E      # MR2\n",  gMr2);}
   if(fabs(m0-gMr3)>m0_)  fprintf(f," 36  %.8E      # MR3\n",  gMr3);

   if(fabs(m0-gMq2)>m0_) {fprintf(f," 41  %.8E      # Mq1\n",  gMq2);
                          fprintf(f," 42  %.8E      # Mq2\n",  gMq2);}
   if(fabs(m0-gMq3)>m0_)  fprintf(f," 43  %.8E      # Mq3\n",  gMq3);
   if(fabs(m0-gMu2)>m0_) {fprintf(f," 44  %.8E      # Mu1\n",  gMu2);
                          fprintf(f," 45  %.8E      # Mu2\n",  gMu2);}
   if(fabs(m0-gMu3)>m0_)  fprintf(f," 46  %.8E      # Mu3\n",  gMu3);
   if(fabs(m0-gMd2)>m0_) {fprintf(f," 47  %.8E      # Md1\n",  gMd2);
                          fprintf(f," 48  %.8E      # Md2\n",  gMd2);}
   if(fabs(m0-gMd3)>m0_)  fprintf(f," 49  %.8E      # Md3\n",  gMd3);
   fclose(f);
   return 0;
}



int amsbLesH(char * fname, double m0,double m32, double  tb, int  sgn)
{
   FILE*  f=fopen(fname,"w");
   if(f==NULL) return -1;

   fprintf(f,"Block MODSEL                 # Select model\n"
             " 1    3                      # amsb\n"
             "Block SMINPUTS               # Standard Model inputs\n");
   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/alfEMZ);
   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5);
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f," 7   %.8E       #  Mtau     \n",      findValW("Ml"));

   fprintf(f,"Block MINPAR                 # Input parameters\n");
   fprintf(f," 1   %.8E       # m0\n",       m0);    
   fprintf(f," 2   %.8E       # m3/2\n",     m32);
   fprintf(f," 3   %.8E       # tanb\n",     tb);
   fprintf(f," 4   %d         # sign(mu)\n", sgn);
   fclose(f);
   return 0;
}

int gmsbLesH(char *fname, double L, double Mmess, double tb, int sgn,int  N5, double cGrav)
{
   FILE*  f=fopen(fname,"w");
   if(f==NULL) return -1;
   fprintf(f,"Block MODSEL                 # Select model\n"
             " 1    3                      # gmsb\n"
             "Block SMINPUTS               # Standard Model inputs\n");
   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/alfEMZ);
   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5);
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f,"Block MINPAR                 # Input parameters\n");
   fprintf(f," 1   %.8E       # Scale \n",   L);    
   fprintf(f," 2   %.8E       # Mmess\n",    Mmess);
   fprintf(f," 3   %.8E       # tanb\n",     tb);
   fprintf(f," 4   %d         # sign(mu)\n", sgn);
   fprintf(f," 5   %d         # index\n",    N5);
   fprintf(f," 6   %.8E       # cGrav\n",    cGrav);
   fclose(f);
   return 0;
}

int EWSBLesH(char * fname, double tb, double MG1, double MG2, double MG3, double Al, double At, double Ab, 
            double mu, double MH3, double Ml1, double Ml2, double Ml3, double Mr1, double Mr2, double Mr3, 
            double Mq1, double Mq2, double Mq3, double Mu1, double Mu2, double Mu3, 
            double Md1, double Md2, double Md3)
{  
   int i;
   FILE*  f=fopen(fname,"w");
   if(f==NULL) return -1;
                                    
   fprintf(f,"Block MODSEL                 # Select model\n"   
           " 1    0                        # Low energy MSSM\n"
           "Block SMINPUTS               # Standard Model inputs\n");
   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/alfEMZ);
   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5); 
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f," 7   %.8E       #  Mtau     \n",      findValW("Ml"));

   fprintf(f,"Block MINPAR                 # Input parameters\n");
   fprintf(f," 3   %.8E       # tanb\n",  tb);

   fprintf(f,"Block EXTPAR\n");
   fprintf(f," 0   -1        # EWSB\n"); 
   fprintf(f," 1   %.8E      # MG1\n", MG1);
   fprintf(f," 2   %.8E      # MG2\n",MG2);
   fprintf(f," 3   %.8E      # MG3\n",MG3);

   fprintf(f," 11  %.8E      # At \n",At);
   fprintf(f," 12  %.8E      # Ab \n",Ab);
   fprintf(f," 13  %.8E      # Atau\n",Al);
         
   fprintf(f," 23  %.8E      # mu\n",mu);
/*   fprintf(f," 24  %.8E      # MA\n",MH3);   */
   fprintf(f," 26  %.8E      # MA\n",MH3);
   fprintf(f," 31  %.8E      # Ml1\n",Ml1);
   fprintf(f," 32  %.8E      # Ml2\n",Ml2); 
   fprintf(f," 33  %.8E      # Ml3\n",Ml3);
   fprintf(f," 34  %.8E      # MR2\n",Mr1);
   fprintf(f," 35  %.8E      # MR2\n",Mr2); 
   fprintf(f," 36  %.8E      # MR3\n",Mr3);

   fprintf(f," 41  %.8E      # Mq1\n",Mq1);
   fprintf(f," 42  %.8E      # Mq2\n",Mq2); 
   fprintf(f," 43  %.8E      # Mq3\n",Mq3);
   fprintf(f," 44  %.8E      # Mu1\n",Mu1);
   fprintf(f," 45  %.8E      # Mu2\n",Mu2); 
   fprintf(f," 46  %.8E      # Mu3\n",Mu3);
   fprintf(f," 47  %.8E      # Md1\n",Md1);
   fprintf(f," 48  %.8E      # Md2\n",Md2); 
   fprintf(f," 49  %.8E      # Md3\n",Md3);
   fclose(f);
   return 0;
}

int sugraHiggsLesH(char *fname,  double tb, double gMG1,double gMG2,double gMG3,
    double gAl, double gAt, double gAb, double gMl2,double gMl3,double gMr2,double gMr3,
    double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3,double mu,double MA) 
{  double m0,mh,a0, m0_,mh_,a0_; 
   double c=1.E-8;
   FILE*  f=fopen(fname,"w");
   if(f==NULL) return -1;
   
   fprintf(f,"Block MODSEL                 # Select model\n"   
           " 1    1                      # sugra\n"
           "Block SMINPUTS               # Standard Model inputs\n");
   fprintf(f," 1   %.8E       # alpha_em^(-1)(MZ) SM MSbar\n",1/alfEMZ);
   fprintf(f," 2   %.8E       # G_Fermi \n",1.16637E-5); 
   fprintf(f," 3   %.8E       # alpha_s(MZ) SM MSbar\n",findValW("alfSMZ"));
   fprintf(f," 5   %.8E       # mb(mb) SM MSbar\n", findValW("MbMb"));
   fprintf(f," 6   %.8E       # mtop(pole)\n",      findValW("Mtp"));
   fprintf(f," 7   %.8E       #  Mtau     \n",      findValW("Ml"));
   fprintf(f,"Block MINPAR                 # Input parameters\n");
   
   a0=findX0(3,gAl, gAt, gAb);
   mh=findX0(3,gMG1,gMG2,gMG3); 
   m0=findX0(15,gMl2,gMl2,gMl3,gMr2,gMr2,gMr3,gMq2,gMq2,gMq3,gMu2,gMu2,gMu3,gMd2,gMd2,gMd3);
   if(m0<0) m0=-m0;
   
   a0_=fabs(a0)*c;
   m0_=fabs(m0)*c;
   mh_=fabs(mh)*c;    
   
   fprintf(f," 1   %.8E       # m0\n",       m0); 
   fprintf(f," 2   %.8E       # m1/2\n",     mh);   
   fprintf(f," 3   %.8E       # tanb\n",     tb);
   fprintf(f," 4   %d     # sign(mu)\n", mu>0? 1:-1);
   fprintf(f," 5   %.8E       # A0\n",       a0);


   fprintf(f,"Block EXTPAR\n");

   fprintf(f,"# GUT parameters\n");
   if(fabs(mh-gMG1)>mh_)  fprintf(f," 1   %.8E      # MG1\n",  gMG1);
   if(fabs(mh-gMG2)>mh_)  fprintf(f," 2   %.8E      # MG2\n",  gMG2);
   if(fabs(mh-gMG3)>mh_)  fprintf(f," 3   %.8E      # MG3\n",  gMG3);

   if(fabs(a0-gAt)>a0_)   fprintf(f," 11  %.8E      # At \n",  gAt);
   if(fabs(a0-gAb)>a0_)   fprintf(f," 12  %.8E      # Ab \n",  gAb);
   if(fabs(a0-gAl)>a0_)   fprintf(f," 13  %.8E      # Al\n",   gAl);
   if(fabs(m0-gMl2)>m0_) {fprintf(f," 31  %.8E      # Ml1\n",  gMl2);
                          fprintf(f," 32  %.8E      # Ml2\n",  gMl2);}
   if(fabs(m0-gMl3)>m0_)  fprintf(f," 33  %.8E      # Ml3\n",  gMl3);
   if(fabs(m0-gMr2)>m0_) {fprintf(f," 34  %.8E      # MR1\n",  gMr2);
                          fprintf(f," 35  %.8E      # MR2\n",  gMr2);}
   if(fabs(m0-gMr3)>m0_)  fprintf(f," 36  %.8E      # MR3\n",  gMr3);

   if(fabs(m0-gMq2)>m0_) {fprintf(f," 41  %.8E      # Mq1\n",  gMq2);
                          fprintf(f," 42  %.8E      # Mq2\n",  gMq2);}
   if(fabs(m0-gMq3)>m0_)  fprintf(f," 43  %.8E      # Mq3\n",  gMq3);
   if(fabs(m0-gMu2)>m0_) {fprintf(f," 44  %.8E      # Mu1\n",  gMu2);
                          fprintf(f," 45  %.8E      # Mu2\n",  gMu2);}
   if(fabs(m0-gMu3)>m0_)  fprintf(f," 46  %.8E      # Mu3\n",  gMu3);
   if(fabs(m0-gMd2)>m0_) {fprintf(f," 47  %.8E      # Md1\n",  gMd2);
                          fprintf(f," 48  %.8E      # Md2\n",  gMd2);}
   if(fabs(m0-gMd3)>m0_)  fprintf(f," 49  %.8E      # Md3\n",  gMd3);
   fprintf(f,"# EWSB  parameters\n");
   fprintf(f," 23  %.8E      # mu\n",  mu);
   fprintf(f," 26  %.8E      # MA\n",  MA);
   
   fprintf(f,"Block QEXTPAR\n");
   fprintf(f," 23  -1      # mu\n");

   fclose(f); 
   return 0;
}
