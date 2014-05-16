#include<stdio.h>
#include<math.h>
#include <ctype.h>
#include <string.h>
#include"alpha.h"

static void skipline(int l)
{  char c; 
   for(;l;l--)  { do scanf("%c",&c); while(c!='\n'); }
}

int  main(int npar, char ** parch)
{
  double buff;
  int i,j,k;
  int nx, nt, NfMx;
  char c;
  char *cteq="CTEQ";
  double *q,*x,*d;
  char version[30];
  int ordr, nf,nf6,nf5,nf4,nf3;
  double ordrf, nff;
  double Mc,Mb,Mt;
  double lambda,L6,L5,L4,L3;
  int Ncteq;
  char names[4][10]={"(6 -6)", "(5 -5)", "(4 -4)", "(3 -3)"};
  char * p;
  
  do scanf("%c",&c); while(c!=':');
  scanf("%s",version);
  p=strstr(version,".tbl");
  if(p) p[0]=0; 
  for(i=0;i<4;i++) 
  if(toupper(version[i]) !=cteq[i]) {fprintf(stderr,"Non expected file format\n"); return 1;}  
  Ncteq=version[4]-'0';
  switch(Ncteq)
  { case 5: Ncteq=4; break;
    case 4:
    case 6: break;
    default:{fprintf(stderr,"Unknown version CTEQ%d \n",Ncteq); return 2;}
  }
                    
  skipline(2);
  if(6!=scanf("%lf %lf %lf %*lf %*lf %*lf %lf %lf %lf", &ordrf, &nff,&lambda,&Mc,&Mb,&Mt )) goto errExit; 
  nf=nff; ordr=ordrf;
  skipline(2);
  if(3!=scanf("%d %d %d",&nx,&nt,&NfMx))goto errExit;
  printf("#distribution \"%s(proton)\"    2212 =>    ",version);
  for(i= 6-NfMx; i<4;i++) printf(" %s",names[i]);
  printf(" -1 -2 21 2 1 \n");

  printf("#distribution \"%s(anti-proton)\" -2212 =>  ",version);
  for(i=6-NfMx; i<4;i++) printf(" %s",names[i]);
  printf(" 1 2 21 -2 -1\n");
    
  skipline(2);
  scanf("%lf",&buff);  printf("\n#q_min %.5E\n", buff); 
  scanf("%lf",&buff); /* printf("\n#q_max %.5E\n", buff); */
  skipline(1);
  q=(double *)malloc(sizeof(double)*(nt+1));
  
    
  for(i=0;i<=nt;i++) scanf("%lf",q+i);
  printf("\n#Q_grid\n"); 
  for(i=0,j=1;i<=nt;i++,j++) 
  {  printf(" %.5E",q[i]); 
     if(j==10) {printf("\n"); j=0;}
  }   
  printf("\n");  

  writeAlpha(stdout,nf,ordr,lambda,NfMx,Mc,Mb,Mt,nt+1,q);

/*
  nf6=nf; L6=lambda;
  if(nf6==6) {nf5=5; L5=findLambda(5,ordr, alpha(6, ordr,L6, 175.) ,175.);}
        else {nf5=nf;L5=lambda;}
  if(nf5==5) {nf4=4; L4=findLambda(4,ordr, alpha(5, ordr,L5, 4.5 ) ,4.5);}
        else {nf4=nf;L4=lambda;}
  if(nf4==4) {nf3=3; L3=findLambda(3,ordr, alpha(4, ordr,L4, 1.4 ) ,1.4);}
        else {nf3=nf;L3=lambda;}

  printf("\n#Alpha\n");
  for(i=0,j=1;i<=nt;i++,j++) 
  {  double al;
     double Q=q[i];
           if(Q<1.4)  al=alpha(nf3, ordr, L3, Q);
     else  if(Q<4.5)  al=alpha(nf4, ordr, L4, Q);
     else  if(Q<175.) al=alpha(nf5, ordr, L5, Q);
     else             al=alpha(nf6, ordr, L6, Q);
     printf(" %.5E",al);
     if(j==10) {printf("\n"); j=0;}
  }   
  printf("\n");  
*/

  skipline(2);
  scanf("%lf",&buff);  printf("\n#x_min %.5E\n", buff);  
  skipline(1);
  printf("\n#X_grid\n");
  for(i=0,j=1;i<=nx;i++,j++) 
  {  scanf("%lf",&buff);     
     printf(" %.5E ",buff); 
     if(j==10) {printf("\n"); j=0;}
  }     
  printf("\n");

  skipline(2);
  
  printf("\n#Interpolation CTEQ%d  %f ",Ncteq,lambda); 
  for(k=0;k<3+NfMx;k++) 
  {  printf("\n#%d-parton\n",k+1);
     for(j=0;j<nt+1;j++) 
     { for(i=0;i<nx+1;i++) {scanf("%lf",&buff); printf(" %.5E",buff);} 
       printf("\n");
     }  
  }
  return 0;
 errExit:
  fprintf(stderr,"%d", ftell(stdin));
       
}
