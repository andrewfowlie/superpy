
#ifndef __MSSM_AUX__
#define __MSSM_AUX__


extern void FillVal(int mode);


/* ============================================
    function for CalcHEP interaction session
===============================================*/ 


extern double  ewsbMSSMc(
 double tb, double MG1,double MG2,double MG3,
 double Am, double Al, double At, double Ab,
 double MH3,double mu,
 double Ml1,double Ml3,
 double Mr1,double Mr3,double Mq1,double Mq3,
 double Mu1,double Mu3,double Md1,double Md3,
         double lCon);

extern double suspectSUGRAc(
  double tb, double gMG1,double gMG2,double gMG3,
  double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
  double gMl2,double gMl3,double gMr2,double gMr3,
  double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3);

extern double isajetSUGRAc(
  double tb, double gMG1,double gMG2,double gMG3,
  double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
  double gMl2,double gMl3,double gMr2,double gMr3,
  double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3);

extern double softSusySUGRAc(
  double tb, double gMG1,double gMG2,double gMG3,
  double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
  double gMl2,double gMl3,double gMr2,double gMr3,
  double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3);   

extern double sphenoSUGRAc(
  double tb, double gMG1,double gMG2,double gMG3,
  double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
  double gMl2,double gMl3,double gMr2,double gMr3,
  double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3);   


/*=============================================
   Evaluation of mass spectrum in case of general MSSM via suspect
===============================================*/
extern int MSSMspect(int LCon);


typedef struct read_param_tag
{ double U[3][3], V[3][3], T[3][3], mc[3], mst[3], msq1, mglu, At, Ab;
 double B[3][3], msb[3], mss[3], mn[5], N[5][5],msnmu, tb;
  double m1, m2, mu,MW, MZ, Mt, Mb, Mhc, Vts,Vtb,Vcb,ee, sw, alphS_MZ;
}  read_param_tag;

extern void calc_eps(struct read_param_tag* param, double *eps_b, double *eps_bp, double *eps_tps);
extern int read_prm(struct read_param_tag* param);

extern int sugraLesH(char *fname,  double tb, double gMG1,double gMG2,double gMG3,
    double gAl, double gAt, double gAb, double sgn, double gMHu, double gMHd,
    double gMl2,double gMl3,double gMr2,double gMr3,
    double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3);

extern int gmsbLesH(char *fname, double L, double Mmess, double  tb, int  sgn, int N5, double cGrav);
extern int amsbLesH(char *fname, double m0,double m32,   double  tb, int  sgn);
extern int EWSBLesH(char *fname, double tb,  double MG1, double MG2, double MG3, double Al, double At,  double Ab, 
                    double mu,   double MH3, double Ml1, double Ml2, double Ml3, double Mr1,double Mr2, double Mr3, 
                    double Mq1,  double Mq2, double Mq3, double Mu1, double Mu2, double Mu3,double Md1, double Md2, double Md3);

int sugraHiggsLesH(char *fname,  double tb, double gMG1,double gMG2,double gMG3,
    double gAl, double gAt, double gAb, double gMl2,double gMl3,double gMr2,double gMr3,
    double gMq2,double gMq3,double gMu2,double gMu3,double gMd2,double gMd3,double mu,double MA); 
                                    

extern void CheckNCsector(double *zero_, double*m1_, double*m2_,double*mu_, double*tb_,
double *mz_, double* sw_);

extern double * VarDump(void);
extern void VarRest(double *b);

extern double width2(char * pName, int *first);

extern int  tree2LesH(void);

#endif
