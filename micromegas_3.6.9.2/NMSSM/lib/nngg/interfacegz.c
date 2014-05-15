
extern int  nnazini_(int*h1,int*h2,int*h3,int*h4,double *v);
extern void  azdifferentialcs_(double *result,double * costh);

double vcsnngz(double v)
{
   int h1[4]={1,1,1,1};
   int h2[4]={1,1,-1,-1};
   double vcs1,vcs2;
   double res[2];
   double coss=0;
   int ok=nnazini_(h1,h1+1,h1+2,h1+3,&v);
   if(!ok) return 0.;
    gzdifferentialcs_(res,&coss);
    vcs1=res[1];
    nnazini_(h2,h2+1,h2+2,h2+3,&v);
    gzdifferentialcs_(res,&coss);
    vcs1+=res[1];
    vcs1*=v;
    return vcs1;
}

double vcsnngz_(double *v){ return vcsnngz(*v);}
