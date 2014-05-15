#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"alpha.h"

#define NX 49
#define NQ 37

int  main(int npar, char ** parch)
{
  int i,j;
  int NfMx=5;
  double q[NQ]= { 1.25,1.5,2.,2.5,3.2,4.,5.,6.4,8.,10.,12.,18.,
                  26.,40.,64.,100.,160.,240.,400.,640.,1e3,1800.,3200.,5600.,1e4,
                  1.8e4,3.2e4,5.6e4,1e5,1.8e5,3.2e5,5.6e5,1e6,1.8e6,3.2e6,5.6e6,1e7};
  double x[NX]={ 1e-5,2e-5,4e-5,6e-5,8e-5,1e-4,2e-4,4e-4,6e-4,
                 8e-4,.001,.002,.004,.006,.008,.01,.014,.02,.03,.04,.06,.08,.1,
                 .125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.375,.4,.425,.45,.475,
                 .5,.525,.55,.575,.6,.65,.7,.75,.8,.9,1. };

  char * version;
  double Mc=sqrt(2.045);
  double Mb=sqrt(18.5);
  double Mt=175;
  
  char names[4][10]={"(6 -6)", "(5 -5)", "(4 -4)", "(3 -3)"};
  double f1[49*37],f2[49*37],f3[49*37],f4[49*37], f5[49*37],f6[49*37],f7[49*37],f8[49*37];

  if( npar!=2 && npar!=5 ) 
  { fprintf(stderr,"This routine needs 1 parameters: identifier of the set, or 4 parameters:\n"
                   "1. the identifier; 2. nf; 3. order (lo, nlo,nnlo); 4. alpha_QCD(MZ)\n");  
    return 1;
  }
  
  version=parch[1];
  
  printf("#distribution \"%s(proton)\"  2212 =>    ",version);
  for(i= 6-NfMx; i<4;i++) printf(" %s",names[i]);
  printf(" -1 -2 21 2 1 \n");

  printf("#distribution \"%s(anti-proton)\" -2212 => ",version);
  for(i=6-NfMx; i<4;i++) printf(" %s",names[i]);
  printf(" 1 2 21 -2 -1\n");

  printf("\n#Q_grid\n"); 
  for(i=0;i<NQ;i++) q[i]=sqrt(q[i]); 
  for(i=0,j=1;i<NQ;i++,j++) 
  {  printf(" %.9E",q[i]); 
     if(j==10) {printf("\n"); j=0;}
  }

  if(npar==5)
  {
     int ordr, nf, nfMx;
     double lambda, alphMZ;

     if(sscanf(parch[2],"%d",  &nf)!=1 || nf>6 || nf<4)
     { fprintf(stderr,"Second parameter should be a hole number between 4 and 6\n");
       exit(1);
     }
          if(strcmp(parch[3],"lo")==0)    ordr=1;
     else if(strcmp(parch[3],"nlo")==0)   ordr=2;
     else if(strcmp(parch[3],"nnlo")==0)  ordr=3;
     else 
     { fprintf(stderr,"Third parameter should be 'lo', 'nlo', or 'nnlo' \n");
       exit(1);
     }        
     if(sscanf(parch[4],"%lf", &alphMZ)!=1 || alphMZ<0.09 || alphMZ>0.15)
     { fprintf(stderr,"Fourth parameter should be a  number between 0.09 and 0.15\n");
       exit(1);
     } 
     if(nf==6){nf=5;nfMx=6;} else nfMx=nf;
     lambda=findLambda(nf,ordr,alphMZ,91.187);
     writeAlpha(stdout,nf,ordr,lambda,nfMx,Mc,Mb,Mt,NQ,q);
  }  

  printf("\n#X_grid\n");
  for(i=0,j=1;i<NX;i++,j++) 
  {
     printf(" %.5E",x[i]); 
     if(j==10) {printf("\n"); j=0;}
  }     
  printf("\n");

  printf("\n#Interpolation MRST2001\n"); 

  for(i=0; i<NX-1; ++i) for(j=0; j<NQ;++j) 
  {  int k=i+j*NX;
     scanf("%lf %lf %lf %lf %lf %lf %lf %lf",f1+k,f2+k,f3+k,f4+k,f5+k,f7+k,f6+k,f8+k);
  }
 

   /* notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea */

/*(b,B) (c,C) (s,S) D U G u d */
           
    printf("\n#q_threshold %.7E\n",Mb);
    printf("\n#1-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.8E",f7[j*NX+i]/x[i]);
      printf(" 0.\n");
    }

    printf("\n#q_threshold %.7E\n",Mc);
    printf("\n#2-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",f5[j*NX+i]/x[i]);
      printf(" 0.\n");
    }

    printf("\n#q_threshold 0\n");
    printf("\n#3-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",f6[j*NX+i]/x[i]);
      printf(" 0.\n");
    }

    printf("\n#4-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",f8[j*NX+i]/x[i]);
      printf(" 0.\n");
    }

    printf("\n#5-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",f4[j*NX+i]/x[i]);
      printf(" 0.\n");
    }

    printf("\n#6-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",f3[j*NX+i]/x[i]);
      printf(" 0.\n");
    }

    printf("\n#7-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",(f1[j*NX+i]+f4[j*NX+i])/x[i]);
      printf(" 0.\n");
    }

    printf("\n#8-parton\n");
    for(j=0;j<NQ;j++)
    { for(i=0;i<NX-1;i++)  printf(" %.6E",(f2[j*NX+i]+f8[j*NX+i])/x[i]);
      printf(" 0.\n");
    }
    return 0;
}
