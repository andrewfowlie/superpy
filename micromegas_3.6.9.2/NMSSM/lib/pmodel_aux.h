
#ifndef __NMSSM_AUX__
#define __NMSSM_AUX__


extern   double ewsbNMSSM(double tb,  double MG1, double MG2, double MG3,   
  double Ml2, double Ml3, double Mr2, double Mr3, double Mq2, double Mq3,   
  double Mu2, double Mu3, double Md2, double Md3, 
  double At,  double Ab,  double Al, double Am, double mu,     
  double LambdQ,double KappaQ,double aLmbdQ,double aKappQ);

extern   double  ewsb_n_MSSM(double tb,  double MG1, double MG2, double MG3,   
  double Ml2, double Ml3, double Mr2, double Mr3, double Mq2, double Mq3,   
  double Mu2, double Mu3, double Md2, double Md3, 
  double At,  double Ab,  double Al, double Am, double mu,     
  double LambdQ,double KappaQ,double aLmbdQ,double aKappQ,double MA,double MP);


extern double sugraNMSSM( double m0, double mhf, double a0, double tb, double sgn, 
       double Lambda,double aLambda, double aKappa);


#endif
