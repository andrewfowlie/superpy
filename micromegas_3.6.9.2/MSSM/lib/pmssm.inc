
extern double xxxEwsbMSSMc(double tb, double  MG1, double  MG2, double  MG3, 
double  Al, double  At, double  Ab, double   mu, double  MH3, 
double  Ml1, double  Ml2, double  Ml3, double  Mr1, double  Mr2, double  Mr3, 
double  Mq1, double  Mq2, double  Mq3, double  Mu1, double  Mu2, double  Mu3, 
double  Md1, double  Md2, double  Md3);

int xxxEwsbMSSM(void)
{
     int err;
     double tb  = findValW("tb"  );
     double MG1 = findValW("MG1" );
     double MG2 = findValW("MG2" );
     double MG3 = findValW("MG3" );
/*     double Am  = findValW("Am"  );*/
     double Al  = findValW("Al"  );
     double At  = findValW("At"  );
     double Ab  = findValW("Ab"  );
     double MH3 = findValW("MH3" );
     double mu  = findValW("mu"  );
     double Ml1 = findValW("Ml1" );
     double Ml2 = findValW("Ml2" );  
     double Ml3 = findValW("Ml3" );  
     double Mr1 = findValW("Mr1" );
     double Mr2 = findValW("Mr2" );  
     double Mr3 = findValW("Mr3" );
     double Mq1 = findValW("Mq1" ); 
     double Mq2 = findValW("Mq2" );
     double Mq3 = findValW("Mq3" );
     double Md1 = findValW("Md1" );
     double Md2 = findValW("Md2" );
     double Md3 = findValW("Md3" );
     double Mu1 = findValW("Mu1" );
     double Mu2 = findValW("Mu2" );
     double Mu3 = findValW("Mu3" );

     err=xxxEwsbMSSMc(tb, MG1, MG2, MG3, Al, At, Ab, 
            mu, MH3, Ml1, Ml2, Ml3, Mr1, Mr2, Mr3, 
            Mq1, Mq2, Mq3, Mu1, Mu2, Mu3, 
            Md1, Md2, Md3);
//     if(err==0) FillVal(0);
     return err;
}

int  xxxewsbmssm_(void){ return  xxxEwsbMSSM();}
