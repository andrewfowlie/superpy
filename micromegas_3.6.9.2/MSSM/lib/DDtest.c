#include<math.h>
#include<stdlib.h>
#include"../../sources/micromegas.h"
#include"../../sources/micromegas_aux.h"
#include"pmodel.h"
#define SQ(x)  ((x)*(x))

extern double GluWeight(void);

int MSSMDDtest(int loop, double*pA0,double*pA5,double*nA0,double*nA5)
{
  double  ApB[2][7],AmB[2][7],MI[2][2][7],NL[5],T3Q[2],EQ[2],mq[7],mqSM[7],
          msq[2][7],Aq[7];
  double mh,mH,ca,sa,mu;
  double o1o1h,o1o1H,SQM[2][2][2][2],capb,sapb,w2s3,w4s3;
  char buffName[10];
  char * ZqNames[6]={"Zdd" ,"Zuu" ,"Zss", "Zcc", "Zb",  "Zt"};
  char * MS1mass[6]={"MSdL","MSuL","MSsL","MScL","MSb1","MSt1"};
  char * MS2mass[6]={"MSdR","MSuR","MSsR","MScR","MSb2","MSt2"};
  char * AqNames[6]={"Ad","Au", "Ad","Au","Ab","At"};
  char mess[10];
  double  ALPE,SW,CW,MW,MZ,E,G,mne,beta,sb,cb,s;
  int i,II,IQ,i1,i2;
  double Ampl0,Ampl2;
  double MN=0.939;
  double MqPole[7]={0,0,0,0,1.67,4.78,173.};
  double qcdNLO,qcdNLOs;
  
  double wS0P__[6],wS0N__[6]; /*scalar */
  double wV5P__[3],wV5N__[3]; /*pseudo-vector*/
  
   for(i=0;i<3;i++)
   { wS0P__[i]= *(&(ScalarFFPd)+i);
     wS0N__[i]= *(&(ScalarFFNd)+i);
     wV5P__[i]= *(&(pVectorFFPd)+i);
     wV5N__[i]= *(&(pVectorFFNd)+i);
  }
  
  for(s=0,i=0;i<3;i++) s+= wS0P__[i];
  for(s=2./27.*(1-s),i=3;i<6;i++)wS0P__[i]=s;

  for(s=0,i=0;i<3;i++) s+= wS0N__[i];
  for(s=2./27.*(1-s),i=3;i<6;i++)wS0N__[i]=s;
    
  *pA0=0,*pA5=0,*nA0=0,*nA5=0;  

  if(sortOddParticles(mess)) return 0;
  if(strcmp(mess,"~o1")!=0) 
  { printf("qbox returns 0 because WINP is not ~o1\n"); return 0;} 
/*ccccccccccccccccccc CONSTANTS ccccccc*/   

  ALPE=1/127.994;
  SW=findValW("SW");
  CW=sqrt(1.-SW*SW);
  MZ=findValW("MZ");
  MW=MZ*CW;  
  E=sqrt(4*M_PI*ALPE);
  G =E/SW;
  mne=fabs(findValW("MNE1"));
/*=======*/
  beta=atan(findValW("tB"));
  sb=sin(beta);
  cb=cos(beta);
  mu=findValW("mu");
/*========  Quark,SQUARK masses and mixing  ======*/
  for(IQ=1;IQ<=6;IQ++)
  { 
    mqSM[IQ]= pMass(pdg2name(IQ));
    mq[IQ]= mqSM[IQ];
        
    msq[0][IQ]=findValW(MS1mass[IQ-1]);
    msq[1][IQ]=findValW(MS2mass[IQ-1]);
    for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++)
    {
       sprintf(buffName,"%s%d%d", ZqNames[IQ-1],i1+1,i2+1);
       MI[i1][i2][IQ]=findValW(buffName);
    }
    Aq[IQ]=findValW(AqNames[IQ-1]);
  }
  
  mq[1]/=1+deltaMd();
  mq[3]/=1+deltaMd();
  mq[5]/=1+deltaMb();
  
  for(i=1;i<=4;i++){sprintf(buffName,"Zn1%d",i); NL[i]=findValW(buffName);}

  T3Q[0]=0.5; EQ[0]=2/3.; T3Q[1]=-0.5; EQ[1]=-1/3.;
  
  for(IQ=1;IQ<=6;IQ++) for( II=0;II<2;II++)
  {  double X,Y,Z,A,B;
     X=-(T3Q[IQ&1]*NL[2]+SW/CW*NL[1]/6);
     Y=SW/CW*EQ[IQ&1]*NL[1];
     Z= -0.5*( (IQ&1)? mq[IQ]*NL[3]/cb: mq[IQ]*NL[4]/sb )/MW;
     A=   G*(MI[II][0][IQ]*(X+Z)+MI[II][1][IQ]*(Y+Z));
     B=   G*(MI[II][0][IQ]*(X-Z)+MI[II][1][IQ]*(-Y+Z));

     ApB[II][IQ]  =(A*A+B*B)/2;        
     AmB[II][IQ]  =(A-B)*(A+B)/2;  /* Normalized like in D&N */ 
  }
/* Higgs sector */ 
  mh=findValW("Mh");
  mH=findValW("MHH");
  sa=findValW("sa");
  ca=findValW("ca");
  o1o1h= E*(ca*NL[4]+sa*NL[3])*(CW*NL[2]-SW*NL[1])/CW/SW;
  o1o1H=-E*(ca*NL[3]-sa*NL[4])*(CW*NL[2]-SW*NL[1])/CW/SW;

/*====================================================================== */ 
/*========= STARTING OF SUMMATION OF DIFFERENT CONTRIBUTIONS =========== */

/*================= THE SD AMPLITUDED ================= */
/******  light squarks SD contribution  */

  for(IQ=1;IQ<=3;IQ++) for(II=0;II<2;II++)
  {  double D=SQ(msq[II][IQ])-SQ(mne)-SQ(mqSM[IQ]);
     double D2=D*D-SQ(2*mne*mqSM[IQ]);
     double f=sqrt(3.)*0.25*ApB[II][IQ]*D/D2;
     *pA5+=f*wV5P__[IQ-1];
     *nA5+=f*wV5N__[IQ-1];
  }

/******  Z  SD contribution */
  for(IQ=1;IQ<=3;IQ++)
  { double f=sqrt(3.)*0.5*(SQ(NL[4])-SQ(NL[3]))*SQ(G/2/MW)*T3Q[IQ&1];
    *pA5+=f*wV5P__[IQ-1];  
    *nA5+=f*wV5N__[IQ-1];
  }
  *pA5/=sqrt(3);
  *nA5/=sqrt(3);

/*================= THE SI AMPLITUDED =================*/

/****** light quarks-squarks SI contribution (plus heavy squarks at tree level)*/

  for(IQ=1;IQ<=(loop? 3:6);IQ++) for(II=0;II<2;II++)
  {
    double g,f,D,D2;
    qcdNLO=qcdNLOs=1;
    if(QCDcorrections)
    { double alphaMq;
    
      if(IQ>3)
      { double alphaMq;
        switch(IQ)
        { case 4: alphaMq=0.39;break;
          case 5: alphaMq=0.22;break;
          default:alphaMq=parton_alpha(mqSM[IQ]);
        } 
        qcdNLO=1+(11./4.-16./9.)*alphaMq/M_PI;
      } 
      qcdNLOs=1+(25./6.-16./9.)*parton_alpha(msq[II][IQ])/M_PI;
    }
/* q,~o1 reaction */
    D=SQ(msq[II][IQ])-SQ(mne)-SQ(mqSM[IQ]);
    D2=D*D-SQ(2*mne*mqSM[IQ]);
    g=-0.25*ApB[II][IQ]/D2;
    f=-0.25*AmB[II][IQ]*D/D2;

    if(!Twist2On) g*=4;    
    *pA0+= (f/mqSM[IQ]-g*mne/2)* MN*wS0P__[IQ-1]*qcdNLO;
    *nA0+= (f/mqSM[IQ]-g*mne/2)* MN*wS0N__[IQ-1]*qcdNLO; 

/******  squarks from nucleon   (~q,~o1 reaction)*/
    D=-SQ(msq[II][IQ])-SQ(mne)+SQ(mqSM[IQ]),
    D2=D*D-SQ(2*mne*msq[II][IQ]);

    g=mqSM[IQ]*AmB[II][IQ]*D/D2;  
    f=mne*ApB[II][IQ]*(SQ(msq[II][IQ])-SQ(mne)+SQ(mqSM[IQ]))/D2;
    *pA0+=(f+g)/SQ(msq[II][IQ])*MN*wS0P__[5]/8*qcdNLOs ;
    *nA0+=(f+g)/SQ(msq[II][IQ])*MN*wS0N__[5]/8*qcdNLOs ;

    if(Twist2On && IQ!=6)
    { double D,g;
      int IQn;
      qcdNLO=1;

      switch(IQ)
      { case 1: IQn=2;break;
        case 2: IQn=1;break;
        default: IQn=IQ;
      } 
    
      D=SQ(msq[II][IQ])-SQ(mne)-SQ(mqSM[IQ]);                                                                               
      g=-0.25*ApB[II][IQ]/(D*D-4*mne*mne*mqSM[IQ]*mqSM[IQ]); 
      *pA0-=1.5*g*mne*MN*parton_x(IQ, msq[II][IQ]-mne);
      *nA0-=1.5*g*mne*MN*parton_x(IQn,msq[II][IQ]-mne);
    }  
  }

/****** Heavy squarks in case of loops   */ 
  if(loop)for(IQ=4;IQ<=6;IQ++) for(II=0;II<2;II++)
  {  double f,g,bd,b1d,bs,b1s,b2s; 
     if(QCDcorrections )
     { double alphaMq;
       switch(IQ)
       { case 4: alphaMq=0.39;break;
         case 5: alphaMq=0.22;break;
         default:alphaMq=parton_alpha(mqSM[IQ]);
       }
       qcdNLO=1+(11./4.-16./9.)*alphaMq/M_PI;
     }  
     else qcdNLO=1;
                               
     bd  = AmB[II][IQ]*mqSM[IQ]*LintIk(1,msq[II][IQ],mqSM[IQ],mne)*3/8.;
     b1d = AmB[II][IQ]*mqSM[IQ]*LintIk(3,msq[II][IQ],mqSM[IQ],mne);      
     bs  =ApB[II][IQ]*mne     *LintIk(2,msq[II][IQ],mqSM[IQ],mne)*3/8.;
     b1s =ApB[II][IQ]*mne     *LintIk(4,msq[II][IQ],mqSM[IQ],mne);	
     b2s =ApB[II][IQ]         *LintIk(5,msq[II][IQ],MqPole[IQ],mne)/4.;
     f=-(bd+bs-mne*b2s/2-mne*mne*(b1d+b1s)/4); 
     *pA0+=f*MN*wS0P__[IQ-1]*qcdNLO;
     *nA0+=f*MN*wS0N__[IQ-1]*qcdNLO;

     if(Twist2On) 
     { double Ampl2;
       Ampl2=parton_alpha(mqSM[IQ])/(12*M_PI)*(b2s+mne*(b1s+b1d)/2)*parton_x(21,MqPole[IQ]);   
      *pA0+=1.5*Ampl2*mne*MN;
      *nA0+=1.5*Ampl2*mne*MN;             
     }     
  } 

/******  higgs-quark-anitiquark */
  for(IQ=1;IQ<=6;IQ++)
  { double fh,fH;
    if(QCDcorrections && IQ>3) 
    { double alphaMq;
      switch(IQ)
      { case 4: alphaMq=0.39;break;
        case 5: alphaMq=0.22;break;
        default:alphaMq=parton_alpha(mqSM[IQ]);
      } 
      qcdNLO=1+(11/4.-16./9.)*alphaMq/M_PI;
    }
    else qcdNLO=1;

    if(IQ&1)
    { 
       double dMq=mqSM[IQ]/mq[IQ]-1;
       fh=o1o1h*E*sa*mq[IQ]*(1-dMq*ca*cb/sa/sb)/(2*MW*cb*SW);   
       fH=-o1o1H*E*ca*mq[IQ]*(1+dMq*sa*cb/ca/sb)/(2*MW*cb*SW);
    }
    else 
    {
       fh=-o1o1h*E*ca*mq[IQ]/(2*MW*sb*SW);   
       fH=-o1o1H*E*sa*mq[IQ]/(2*MW*sb*SW);
    }   
    *pA0+=0.5*(fh/(mh*mh)+fH/(mH*mH))/mqSM[IQ]*MN*wS0P__[IQ-1]*qcdNLO;
    *nA0+=0.5*(fh/(mh*mh)+fH/(mH*mH))/mqSM[IQ]*MN*wS0N__[IQ-1]*qcdNLO;
  } 

/******  higgs squark-antisquark */
  capb=ca*cb-sa*sb;
  sapb=sa*cb+ca*sb;
  w4s3=4*SW*SW-3;
  w2s3=2*SW*SW-3;

#define h 0
#define H 1
#define L 0
#define R 1
#define U 0
#define D 1

  SQM[h][L][L][U]= sapb/SW*(-w4s3/2);
  SQM[h][L][L][D]= sapb/SW*( w2s3/2);
  SQM[h][R][R][U]= sapb*SW*( 2);
  SQM[h][R][R][D]= sapb*SW*(-1);
  SQM[H][L][L][U]= capb/SW*( w4s3/2);
  SQM[H][L][L][D]= capb/SW*(-w2s3/2);
  SQM[H][R][R][U]= capb*SW*(-2);
  SQM[H][R][R][D]= capb*SW*( 1);
  
  SQM[h][L][R][U]=SQM[h][R][L][U]=0;
  SQM[h][L][R][D]=SQM[h][R][L][D]=0;
  SQM[H][L][R][U]=SQM[H][R][L][U]=0;
  SQM[H][L][R][D]=SQM[H][R][L][D]=0;
 
  for(IQ=1;IQ<=6;IQ++)for(II=0;II<2;II++)
  { double fh,fH;
    int i,j;

    if(QCDcorrections)qcdNLOs=1+(25./6.-16./9.)*parton_alpha(msq[II][IQ])/M_PI;
    else qcdNLOs=1;

    for(fh=0,fH=0,i=0;i<2;i++)for(j=0;j<2;j++)
    {  double dSQMh,dSQMH;
       double b=-T3Q[IQ&1]+0.5,
              t= T3Q[IQ&1]+0.5;   
       if(i==j) 
       { dSQMh=( ca*t - sa*b )*3*SQ(mq[IQ]*CW/MW)/SW;
         dSQMH=( sa*t + ca*b )*3*SQ(mq[IQ]*CW/MW)/SW;
       }  
       else
       {
         dSQMh=( (ca*Aq[IQ]+mu*sa)*t - (sa*Aq[IQ]+mu*ca)*b )*1.5*mq[IQ]*SQ(CW/MW)/SW;  
         dSQMH=( (sa*Aq[IQ]-mu*ca)*t + (ca*Aq[IQ]-mu*sa)*b )*1.5*mq[IQ]*SQ(CW/MW)/SW;
       }
       
       if(IQ&1) {dSQMh/=cb;dSQMH/=cb;} else {dSQMh/=sb;dSQMH/=sb;} 
       
       fh+=(SQM[h][i][j][IQ&1]-dSQMh)*MI[II][i][IQ]*MI[II][j][IQ];
       fH+=(SQM[H][i][j][IQ&1]-dSQMH)*MI[II][i][IQ]*MI[II][j][IQ];  
    }   
    fh*=o1o1h*E*MW/(3*CW*CW);
    fH*=o1o1H*E*MW/(3*CW*CW);     
    *pA0+=(fh/(mh*mh)+fH/(mH*mH))/(2*msq[II][IQ]*msq[II][IQ])*MN*wS0P__[5]/8*qcdNLOs;
    *nA0+=(fh/(mh*mh)+fH/(mH*mH))/(2*msq[II][IQ]*msq[II][IQ])*MN*wS0N__[5]/8*qcdNLOs;
  } 
}
