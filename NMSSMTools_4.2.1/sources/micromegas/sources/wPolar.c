#include "micromegas.h"
#include "micromegas_aux.h"

int vPolar( int out1,int out2,int out3, double*left,double*right,double*lng)
{ double pvect[20],pcm1,pcm2,ms,md,chY,shY;
  int i,err_code;
  int iW,ie,in;
  double m[5];
  int code[5];
  double massMin=-1;
  int oId;
  char n1[10],n2[10],n3[10];
  
  numout * cc;
  char *Wp=NULL,*el=NULL,*Ne=NULL,*o1=NULL,*O1=NULL;
  double r[3];
  double GG=sqrt(4*M_PI*parton_alpha(2*Mcdm));
  
  for(i=0;i<nModelParticles;i++)
  { 
    if(ModelPrtcls[i].NPDG== out1) { Wp=ModelPrtcls[i].name; sprintf(n1,"p%d",i+1);}
    if(ModelPrtcls[i].NPDG==-out1) { Wp=ModelPrtcls[i].aname;sprintf(n1,"m%d",i+1);}

    if(ModelPrtcls[i].NPDG== out2) { el=ModelPrtcls[i].name; sprintf(n2,"p%d",i+1);}
    if(ModelPrtcls[i].NPDG==-out2) { el=ModelPrtcls[i].aname;sprintf(n2,"m%d",i+1);}

    if(ModelPrtcls[i].NPDG== out3) { Ne=ModelPrtcls[i].name;sprintf(n3,"p%d",i+1);}
    if(ModelPrtcls[i].NPDG==-out3) { Ne=ModelPrtcls[i].aname; sprintf(n3,"m%d",i+1);}

    if(ModelPrtcls[i].name[0]=='~')
    { double  mass=fabs(findValW(ModelPrtcls[i].mass));
      if(mass<massMin || massMin<0) 
      {
        o1=ModelPrtcls[i].name; O1=ModelPrtcls[i].aname;
        massMin=mass;
        oId=i+1;
      }
    }  
 }

 if(!o1 || !O1 || !Wp || !el || !Ne) return 1;

  { char lib[40], process[30], cond[20];
    sprintf(lib,"p%d_Polar_%s%s%s",oId,n1,n2,n3);
    sprintf(process,"%s,%s->%s,%s,%s",o1,O1,Wp,el,Ne);
    sprintf(cond,"%s!=2",Wp);  
    cc=getMEcode(0,1,process,cond,NULL,lib); 
    if(!cc) { printf("can not generate\n");return 2;} 
  }
   
  for(i=1;i<=cc->interface->nvar;i++) if(cc->link[i]) cc->interface->va[i]=*(cc->link[i]);
  if( cc->interface->calcFunc()>0 ) { return -1;}  
  for(i=0;i<5;i++) cc->interface->pinf(1,i+1,m+i,code+i);
  for(i=0;i<20;i++) pvect[i]=0;

  pvect[0]=m[0];
  pvect[4]=m[1];
  
  iW=ie=in=-1;
  for(i=2;i<5;i++) 
  { if(code[i]==out1) iW=i;
    if(code[i]==out2) ie=i; 
    if(code[i]==out3) in=i;
  } 
  if(iW<0 || ie <0 || in<0) return 1; 
  if(m[0]+m[1]<=2*m[iW]) { printf("Mcdm too low\n"); return 3;}
  pcm1=sqrt((m[0]+m[1])*(m[0]+m[1]) - 4*m[iW]*m[iW])/2;
  ms=m[ie]+m[in];
  md=m[ie]-m[in];
  pcm2=sqrt((m[iW]*m[iW] - ms*ms)*(m[iW]*m[iW]-md*md))/(2*m[iW]);

  for(i=0;i<3;i++)
  {  double csfi=i-1;
     pvect[4*iW]= sqrt(m[iW]*m[iW]+pcm1*pcm1);
     pvect[4*iW+3] = -pcm1;    
     pvect[4*ie]=sqrt(m[ie]*m[ie]+pcm2*pcm2);
     pvect[4*ie+3]=pcm2*csfi;
     pvect[4*ie+2]=pcm2*sqrt(1-csfi*csfi);
  
     pvect[4*in]=sqrt(m[in]*m[in]+pcm2*pcm2);
     pvect[4*in+3]=-pcm2*csfi;
     pvect[4*in+2]=-pcm2*sqrt(1-csfi*csfi);  

     chY=sqrt(1+pcm1*pcm1/m[iW]/m[iW]);
     shY=sqrt(pcm1*pcm1/m[iW]/m[iW]);  
  
     { double p0=pvect[4*ie], p3=pvect[4*ie+3];
       pvect[4*ie]=  chY*p0 + shY*p3;
       pvect[4*ie+3]=shY*p0 + chY*p3;

       p0=pvect[4*in]; p3=pvect[4*in+3];
       pvect[4*in]=  chY*p0 + shY*p3;
       pvect[4*in+3]=shY*p0 + chY*p3;
     }
     r[i]=(cc->interface->sqme)(1,GG,pvect,&err_code);
  }

  { double s;
  
    r[2]/=4;
    r[0]/=4;
    r[1]=(r[1]-r[0]-r[2])/2;
    
    s=r[0]+r[1]+r[2];
                  
    s=r[0]+r[1]+r[2];
    if(left)  *left=r[2]/s;
     
    if(right) *right=r[0]/s;
    if(lng)   *lng=r[1]/s;
  }
  return 0;
}  
