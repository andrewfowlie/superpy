extern double xxxSUGRAnuhc( double tb, double  gMG1, double gMG2, double gMG3, 
    double gAl, double  gAt, double  gAb, double gMl2, double gMl3, double gMr2, double gMr3, 
    double gMq2, double gMq3, double gMu2, double gMd2, double gMu3, double gMd3,double mu,double MA);


int xxxSUGRAnuh(double tb, double gMG1,double gMG2,double gMG3,
   double gAl, double gAt, double gAb, double gMl2,double gMl3,double gMr2,double gMr3,
   double gMq2,double gMq3,double gMu2,double gMd2,double gMu3,double gMd3, double mu,double MA)
{ int err;
   err= xxxSUGRAnuhc( tb,  gMG1, gMG2, gMG3, gAl,  gAt,  gAb,
    gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMd2, gMu3, gMd3,mu,MA);
    assignValW("tb",tb);
   if(err==0) FillVal(1);
   return err;
}

int  xxxsugranuh_(double* tb, double* gMG1,double* gMG2,double* gMG3,
   double* gAl, double* gAt, double* gAb, double* gMl2,double* gMl3,double* gMr2,double* gMr3,
   double* gMq2,double* gMq3,double* gMu2,double* gMd2,double* gMu3,double* gMd3,double*mu,double*MA)
{  return  xxxSUGRAnuh(*tb, *gMG1,*gMG2,*gMG3, *gAl, *gAt, *gAb,
   *gMl2,*gMl3,*gMr2,*gMr3, *gMq2,*gMq3,*gMu2,*gMd2,*gMu3,*gMd3,*mu,*MA);
}
