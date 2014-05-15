#include "pmodel.h"
#include "pmodel_aux.h"
#include"../../sources/micromegas_aux.h"


static double neuz[4][4], neur[4];

static double symbdiag4(double MG1, double MG2, double sw, double MZ, double sb,
	double cb, double mu)
{       
	double cw=sqrt(1.0-sw*sw),  s2b=2.*sb*cb;
	double  Znx0 = (MG1+MG2)/4.,
          Znxc2 = MG1*MG2 - MZ*MZ - mu*mu - 3./8.*(MG1+MG2)*(MG1+MG2),
          Znxc3a = -1./8.*(MG1+MG2)*(MG1+MG2)*(MG1+MG2)
			+1./2.*(MG1+MG2)*(MG1*MG2-MZ*MZ-mu*mu),
          Znxc3 = Znxc3a + (MG1+MG2)*mu*mu + (MG1*cw*cw+MG2*sw*sw)*MZ*MZ
                   - mu*MZ*MZ*s2b,
          Znxc4b = (MG1*cw*cw+MG2*sw*sw)*mu*MZ*MZ*s2b - MG1*MG2*mu*mu,
          Znxc4a = Znxc4b
                   + Znx0*Znx0*(MG1*MG2-MZ*MZ-mu*mu)- 3./256.*
			(MG1+MG2)*(MG1+MG2)*(MG1+MG2)*(MG1+MG2),
		   
          Znxc4= Znxc4a + Znx0*( (MG1+MG2)*mu*mu+(MG1*cw*cw
	                      +MG2*sw*sw)*MZ*MZ  -mu*MZ*MZ*s2b),
	  
	  Znxs = - Znxc3*Znxc3 - 2./27.*Znxc2*Znxc2*Znxc2 + 8./3.*Znxc2*Znxc4,
          Znxu = - 1./3.*Znxc2*Znxc2 - 4.*Znxc4;
	  
      double  Zncxd =  -4.*Znxu*Znxu*Znxu - 27.*Znxs*Znxs,
	  
          Zncua=atan2(sqrt(Zncxd/27.)/2.,-Znxs/2.)/3.,
	  Zncxa=pow(Zncxd/108.+Znxs*Znxs/4.,1./6.)*cos(Zncua),
	
	  Zncxb = 8.*Zncxa - 8./3.*Znxc2,
	  
	  Zncx1 =  Zncxa/2. - Znxc2/6.,
          Zncx2 = -Zncxa/2. - Znxc2/3.,
          Zncx3 = Znxc3/sqrt(Zncxb);

int i;

       neur[0] = Znx0 - sqrt(Zncx1) + sqrt(Zncx2+Zncx3),
       neur[1] = Znx0 + sqrt(Zncx1) - sqrt(Zncx2-Zncx3),
       neur[2] = Znx0 - sqrt(Zncx1) - sqrt(Zncx2+Zncx3),
       neur[3] = Znx0 + sqrt(Zncx1) + sqrt(Zncx2-Zncx3);
 
  for(i=0;i<3;)
  { if(fabs(neur[i])>fabs(neur[i+1]))  
    { double c=neur[i];
      neur[i]=neur[i+1];
      neur[i+1]=c; 
      if(i) i--; else i++;
    } else i++;
  }

  for(i=0;i<4;i++)
  {  double Lmbd,x1,x2,x3,x4,n;
     Lmbd=neur[i];
     x1=1;
     x4=-(MG1-Lmbd)*(Lmbd*sb+mu*cb)/((Lmbd+2*mu*cb*sb)*MZ*sw);          
     x3=-x4*(mu*sb+Lmbd*cb)/(Lmbd*sb+mu*cb);
     x2=(MZ*sb*sw-mu*x3-Lmbd*x4)/(MZ*sb*cw);
     n=1/sqrt(1+x2*x2+x3*x3+x4*x4);
     neuz[i][0]=x1*n; neuz[i][1]=x2*n; neuz[i][2]=x3*n; neuz[i][3]=x4*n;
  }
/*
  { double a[4][4];
    int k1,k2,j;
     a[0][0]=MG1;       a[0][1]=0;         a[0][2]=-MZ*cb*sw; a[0][3]= MZ*sb*sw;
     a[1][0]=0  ;       a[1][1]=MG2;       a[1][2]= MZ*cb*cw; a[1][3]=-MZ*sb*cw;
     a[2][0]=-MZ*cb*sw; a[2][1]= MZ*cb*cw; a[2][2]=0;         a[2][3]=-mu;	
     a[3][0]= MZ*sb*sw; a[3][1]=-MZ*sb*cw; a[3][2]=-mu;       a[3][3]=0;

    for(i=0;i<4;i++) for(j=0;j<4;j++)
    { double sum=0;
      for(k1=0;k1<4;k1++) for(k2=0;k2<4;k2++)
      sum+=neuz[i][k1]*a[k1][k2]*neuz[j][k2];
      if(fabs(sum)<1.E-10) sum=0;  
      printf("Test of symbdiag4: i=%d j=%d MassMatrix=%e\n",i,j,sum);
    }
  }
*/
  	return 1.0;
}


static double ne4mixm(double di, double dj, double ok)
{
  int i=di+0.01,j=dj+0.01; 
  return neuz[i-1][j-1];
}

static double ne4mass(double di, double ok)
{ int i=di+0.01;
  return  neur[i-1];
}


static double MC1,MC2,Zu11,Zu21,Zu12,Zu22,Zv11,Zv21,Zv12,Zv22;

static double chargDiag(double MG2,double MW,double mu,double sb,double cb)
{
  double muQ=mu*mu;
  double MWQ=MW*MW;
  double MG2Q=MG2*MG2;
  double c2b=cb*cb-sb*sb;
  double Sqrt2=sqrt(2.);

  double Zctx=(MG2Q-muQ)*(MG2Q-muQ)/4+MWQ*MWQ*c2b*c2b+MWQ*(MG2Q+muQ+2*mu*MG2*2*sb*cb); 
  double Zcc2u=-(MG2Q-muQ-2*MWQ*c2b)/(2*sqrt(Zctx));
  double Zcc2v=-(MG2Q-muQ+2*MWQ*c2b)/(2*sqrt(Zctx));

  double Zcsigu= MG2*cb+mu*sb>0 ? -1 :1 ;

  double Zcsigv= MG2*sb+mu*cb>0 ? -1 :1 ;

  double Zccu=sqrt((1+Zcc2u)/2);
  double Zcsu=sqrt((1-Zcc2u)/2)*Zcsigu;
  double Zccv=sqrt((1+Zcc2v)/2);
  double Zcsv=sqrt((1-Zcc2v)/2)*Zcsigv;

  double det= mu*MG2-2*MWQ*sb*cb>0 ? 1 : -1;
  double MC1_=MW*Sqrt2*(cb*Zcsu*Zccv+sb*Zccu*Zcsv)+mu*Zcsv*Zcsu+MG2*Zccu*Zccv;
  double tsgn= MC1_>0? det :-det;
  double MC2_=-MW*Sqrt2*(cb*Zccu*Zcsv+sb*Zcsu*Zccv)+mu*Zccv*Zccu+MG2*Zcsu*Zcsv;
/*
 printf("chargDiag  MG2=%f, MW=%f, mu=%f, sb=%f, cb=%f\n", MG2, MW, mu, sb, cb);
*/
  MC1=fabs(MC1_);
  MC2=fabs(MC2_);
  Zu11=Zccu*det*tsgn;
  Zu21=-Zcsu*tsgn;
  Zu12=Zcsu*det*tsgn;
  Zu22=Zccu*tsgn;
  Zv11=Zccv;
  Zv21=-Zcsv;
  Zv12=Zcsv;
  Zv22=Zccv;

  return 1;
}

static double chargMass(double i_, double ok)
{ int i=i_+0.001;
  if(i==1) return MC1;
  else if(i==2) return MC2; else return 0;
}


static double chargZu(double i_,double j_,double ok)
{ int i=i_+0.001, j=j_+0.001;

  switch(i) 
  { case 1 :if(j==1) return Zu11; else if(j==2) return Zu12; else return 0;
    case 2 :if(j==1) return Zu21; else if(j==2) return Zu22; else return 0;
    default:return 0;
  }
}

static double chargZv(double i_,double j_,double ok)
{ int i=i_+0.001, j=j_+0.001;
 
  switch(i) 
  { case 1 :if(j==1) return Zv11; else if(j==2) return Zv12; else return 0;
    case 2 :if(j==1) return Zv21; else if(j==2) return Zv22; else return 0;
    default:return 0;
  }
}
  
int  tree2LesH(void)
{  
  double  SW = findValW("SW"),        SW_2=SW*SW;
  double  MZ = findValW("MZ"),        MZ_2=MZ*MZ;
  double  Ml = findValW("Ml"),        Ml_2=Ml*Ml;
/*  double  Mtp =findValW("Mtp"); */
  double  tb = findValW("tb"),        sb=tb/sqrt(1+tb*tb),
                                      cb=sqrt(1-sb*sb),
                                      c2b=(cb-sb)*(cb+sb);
  double  mu = findValW("mu");
  double  MG1= findValW("MG1");
  double  MG2= findValW("MG2");
  double  MG3= findValW("MG3");
/*  double  MH3= findValW("MH3");  */
  double  Ml1= findValW("Ml1"),       Ml1_2=Ml1*Ml1;
  double  Ml2= findValW("Ml2"),       Ml2_2=Ml2*Ml2;
  double  Ml3= findValW("Ml3"),       Ml3_2=Ml3*Ml3;
  double  Mr1= findValW("Mr1"),       Mr1_2=Mr1*Mr1;
  double  Mr2= findValW("Mr2"),       Mr2_2=Mr2*Mr2;
  double  Mr3= findValW("Mr3"),       Mr3_2=Mr3*Mr3;
  double  Mq1= findValW("Mq1"),       Mq1_2=Mq1*Mq1;
  double  Mq2= findValW("Mq2"),       Mq2_2=Mq2*Mq2;
  double  Mq3= findValW("Mq3"),       Mq3_2=Mq3*Mq3;
  double  Mu1= findValW("Mu1"),       Mu1_2=Mu1*Mu1;
  double  Mu2= findValW("Mu2"),       Mu2_2=Mu2*Mu2;
  double  Mu3= findValW("Mu3"),       Mu3_2=Mu3*Mu3;
  double  Md1= findValW("Md1"),       Md1_2=Md1*Md1;
  double  Md2= findValW("Md2"),       Md2_2=Md2*Md2;
  double  Md3= findValW("Md3"),       Md3_2=Md3*Md3;
  double  Al = findValW("Al");
  double  Ab = findValW("Ab");
  double  At = findValW("At");
  
  double  Mb,Mt,Mb_2,Mt_2,ok,Q,mll,mrr,mlr,det,v;
  double  MSlth, MSbth, MStth;
  double MbMb;
  int i,j;
  char name[10];   

  ok=symbdiag4(MG1, MG2, SW, MZ, sb, cb, mu);
  if(ok<0) return 2;
 
  for(i=1;i<=4;i++){sprintf(name,"MNE%d",i);  assignValW(name, ne4mass(i, ok)); }
  
  for(i=1;i<=4;i++) 
  for(j=1;j<=4;j++){sprintf(name,"Zn%d%d",i,j);assignValW(name,ne4mixm(i,j,ok));}

  ok=chargDiag(MG2,MZ*sqrt(1-SW_2),mu,sb,cb); if(ok<0) return 2;

  for(i=1;i<=2;i++) { sprintf(name,"MC%d",i);  assignValW(name, chargMass(i, ok)); } 
  for(i=1;i<=2;i++) for(j=1;j<=2;j++)
  { sprintf(name,"Zu%d%d",i,j);assignValW(name,chargZu(i,j,ok)); 
    sprintf(name,"Zv%d%d",i,j);assignValW(name,chargZv(i,j,ok));
  }
  assignValW("MSG",MG3);
  MbMb=findValW("MbMb");
  Q=initQCD(findValW("alfSMZ"),1.4,MbMb,findValW("Mtp"));
  if(Q<0 || Q>=MbMb) return  3;

  Q=fabs(findValW("MNE1"));

  Mt=MtEff(2*Q); Mt_2=Mt*Mt;
  Mb=MbEff(2*Q); Mb_2=Mb*Mb;


  v=c2b*MZ_2/2+Ml1_2;              if(v<0)return 2; assignValW("MSne",sqrt(v));
  v=c2b*MZ_2/2+Ml2_2;              if(v<0)return 2; assignValW("MSnm",sqrt(v));
  v=c2b*MZ_2/2+Ml3_2;              if(v<0)return 2; assignValW("MSnl",sqrt(v));

  v=Ml1_2-MZ_2*(0.5-SW_2)*c2b;     if(v<0)return 2; assignValW("MSeL",sqrt(v));
  v=Mr1_2-MZ_2*SW_2*c2b;           if(v<0)return 2; assignValW("MSeR",sqrt(v));

  v=Ml2_2-MZ_2*(0.5-SW_2)*c2b;     if(v<0)return 2; assignValW("MSmL",sqrt(v));
  v=Mr2_2-MZ_2*SW_2*c2b;           if(v<0)return 2; assignValW("MSmR",sqrt(v));

  mll=Ml3_2+Ml_2-MZ_2*(0.5-SW_2)*c2b;
  mlr=Ml*(Al-mu*tb);                             
  mrr=Mr3_2+Ml_2-MZ_2*SW_2*c2b; 
  det=sqrt((mll-mrr)*(mll-mrr)+4*mlr*mlr);

  v=(mll+mrr-det)/2;               if(v<0)return 2; assignValW("MSl1",sqrt(v)); 
  v=(mll+mrr+det)/2;               if(v<0)return 2; assignValW("MSl2",sqrt(v)); 
  MSlth=atan2(-2*mlr,mrr-mll)/2;

  v=Mq1_2+MZ_2*(0.5-2/3.*SW_2)*c2b;if(v<0)return 2; assignValW("MSuL",sqrt(v));
  v=Mu1_2+MZ_2*(2/3.*SW_2)*c2b;    if(v<0)return 2; assignValW("MSuR",sqrt(v));
  
  v=Mq2_2+MZ_2*(0.5-2/3.*SW_2)*c2b;if(v<0)return 2; assignValW("MScL",sqrt(v));
  v=Mu2_2+2/3.*MZ_2*SW_2*c2b;      if(v<0)return 2; assignValW("MScR",sqrt(v));

  mll =Mq3_2+Mt_2+MZ_2*(0.5-2/3.*SW_2)*c2b;
  mlr =Mt*(At-mu/tb);
  mrr =Mu3_2+Mt_2+2/3.*MZ_2*SW_2*c2b;
  det=sqrt((mll-mrr)*(mll-mrr)+4*mlr*mlr);

  v=(mll+mrr-det)/2;               if(v<0)return 2; assignValW("MSt1",sqrt(v));
  v=(mll+mrr+det)/2;               if(v<0)return 2; assignValW("MSt2",sqrt(v));
  MStth=atan2(-2*mlr,mrr-mll)/2;


  v=Mq1_2-MZ_2*(0.5-1/3.*SW_2)*c2b;if(v<0)return 2; assignValW("MSdL",sqrt(v));
  v=Md1_2-MZ_2*(1/3.*SW_2)*c2b;    if(v<0)return 2; assignValW("MSdR",sqrt(v));

  v=Mq2_2-MZ_2*(0.5-1/3.*SW_2)*c2b;if(v<0)return 2; assignValW("MSsL",sqrt(v));
  v=Md2_2-MZ_2*(1/3.*SW_2)*c2b;    if(v<0)return 2; assignValW("MSsR",sqrt(v));

  mll =Mq3_2+Mb_2-MZ_2*(0.5-1/3.*SW_2)*c2b;
  mlr =Mb*(Ab-mu*tb);
  mrr =Md3_2+Mb_2-1/3.*MZ_2*SW_2*c2b;
  det=sqrt((mll-mrr)*(mll-mrr)+4*mlr*mlr);

  v=(mll+mrr-det)/2;               if(v<0)return 2; assignValW("MSb1",sqrt(v));
  v=(mll+mrr+det)/2;               if(v<0)return 2; assignValW("MSb2",sqrt(v));
  MSbth=atan2(-2*mlr,mrr-mll)/2;
  assignValW("Zl11", cos(MSlth)); 
  assignValW("Zl21",-sin(MSlth));
  assignValW("Zl12", sin(MSlth));
  assignValW("Zl22", cos(MSlth));
  assignValW("Zt11", cos(MStth));
  assignValW("Zt21",-sin(MStth));
  assignValW("Zt12", sin(MStth));
  assignValW("Zt22", cos(MStth));
  assignValW("Zb11", cos(MSbth));
  assignValW("Zb21",-sin(MSbth));
  assignValW("Zb12", sin(MSbth));
  assignValW("Zb22", cos(MSbth));
  return 0;
}


#ifdef OLD
int readEwsb(char *fname)
{
  double val;
  int n,i;
  char name[80];
  char * vName[25]={"tb","MG1","MG2","MG3","Am","Al","At","Ab","MH3","mu",
"Ml1","Ml2","Ml3","Mr1","Mr2","Mr3","Mq1","Mq2","Mq3","Mu1","Mu2","Mu3","Md1","Md2","Md3"};  
 

  FILE * f=fopen(fname,"r");

  if(f==NULL) return -1;
  for(n=1;;n++)
  { if(fscanf(f,"%s",name)!=1) { n=0; break;}
    if(name[0]=='#') { fscanf(f,"%*[^\n]"); continue;}
    if(fscanf(f,"%lf",&val)!=1) break;
    fscanf(f,"%*[^\n]");
    { int err=assignVal(name,val);
      if(err==1) break;
      for(i=0;i<25;i++) if(vName[i]&&strcmp(vName[i],name)==0)vName[i]=NULL;
    }
  }
  fclose(f);
  if(n) return n;
  for(i=0; i<25;i++) if(vName[i]) 
  {  printf("readEwsb: '%s' is not defined\n",vName[i]);  
     n=-2;
  } 
  return n;
}

int  readewsb_(char *f_name, int len)
{
   char *name=malloc(len+10);
   int err;
   fName2c(f_name,name,len);
   err= readEwsb(name);
   free(name);
   return err;
}

#endif
