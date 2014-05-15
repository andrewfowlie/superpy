                                                                                 
extern void  nnaaini_(int*h1,int*h2,int*h3,int*h4,double *v);
extern void  differentialcs_(double *result,double * costh);

double vcsnngg(double v)
{
   int h1[4]={0,0,1,1};
   int h2[4]={1,1,-1,-1};
   double vcs1,vcs2;
   double res[2];
   double coss=0;
                                                                                 
    nnaaini_(h1,h1+1,h1+2,h1+3,&v);
    differentialcs_(res,&coss);
    vcs1=res[1];
    /*nnaaini_(h2,h2+1,h2+2,h2+3,&v);
    differentialcs_(res,&coss);
    vcs1+=res[1];*/
    vcs1*=v*4.0;
    return vcs1;
}

double vcsnngg_(double *err){ return vcsnngg(*err);}

