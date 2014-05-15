/*
  Author Alexander Pukhov.
*/

#include "include/pdt.h"

/* STEQ INTERPOLATION CODES {*/

/*extern double qql[],xxl[],f3[],f7[],*cc3,*ccb,qqlb[],fb[];*/

static double  polint(double x, int n,  double *xa, double *ya)
{  double z[10];
   int i,m;

   for(i=0;i<n;i++) z[i]=ya[i];

   for(m=1;m<n;m++) for(i=0;i<n-m;i++)
   z[i]=(z[i]*(xa[i+m]-x) - z[i+1]*(xa[i]-x))/(xa[i+m]-xa[i]);
   return z[0];
}    


static int  leftX(int dim, double * xa, double x)
{  int k1,k2,k3;

   if(x<xa[0]) return 0;
   if(x>xa[dim-1]) return dim-3;

   k1=0; 
   k2=dim-1;
      
   while(k2-k1>1)
   { k3=(k1+k2)/2;
     if(xa[k3]>x)k2=k3; else k1=k3;
   } 

   k3=k1;
   if(k3<0) k3=0;
   if(k3>dim-2) k3=dim-2; 
   return k3;
}


static double int_cteq4(double x, double q, pdtStr * W)
{ 
  int px = leftX(W->nx, W->x_grid, x);
  double logQ=log(q);
  int pq=leftX(W->nq, W->q_grid, logQ);
  double tmp[3];
  int i;
  
  if(pq >W->nq-3) pq=W->nq-3;
  if(px >W->nx-3) px=W->nx-3;
  
  for(i=0;i<3;i++)tmp[i]=polint(x,3, W->x_grid +px,W->strfun +W->nx*(pq+i)+px);
  return polint(logQ,3,W->q_grid+pq,tmp);
}


static double qSplineCteq6(double x,double *xa, double *ya)
{
  double s12=xa[0]-xa[1], s13=xa[0]-xa[2], s23=xa[1]-xa[2], s24=xa[1]-xa[3], 
         s34=xa[2]-xa[3];
  double sy2=x-xa[1], sy3=x-xa[2];

  double c1=s13/s23, c2=s12/s23, c3=s34/s23, c4=s24/s23;

  double cxx=sy2*sy3/( s12*s34 - (s12+s13)*(s24+s34));
  double c5=(s34*sy2-(s24+s34)*sy3)*cxx/s12, c6=((s12+s13)*sy2-s12*sy3)*cxx/s34;
  double res =

        ( c5*(ya[0] -  ya[1]*c1 +ya[2]*c2 )
         +c6*(ya[3]+ya[1]*c3 - ya[2]*c4) 
         +ya[1]*sy3 - ya[2]*sy2
        )/s23;
  return res;
}

static double int_cteq6(double x, double q, pdtStr * W)
{
  double x3=pow(x,0.3);
  double loglogQ=log(log(q/W->Q0_cteq));
  double * q_grid=W->q_grid_aux;
  double * x_grid=W->x_grid_aux;
  int pq = leftX(W->nq, q_grid, loglogQ)-1;
  int px = leftX(W->nx, W->x_grid_aux, x3)-1;
  int i;
  int qExt=0;
  int xExt=0;
  double tmp[4];
  
  if(pq<0)  { pq=0; qExt=-1;} else if(pq > W->nq-4) { pq=W->nq-4; qExt=1;}
  if(px<=0) { px=0; xExt=-1;} else if(px > W->nx-4) { px=W->nx-4; xExt=1;}

       if(xExt==0) for(i=0;i<4;i++) 
            tmp[i]=qSplineCteq6(x3, x_grid+px, W->strfun+W->nx*(pq+i)+px);
  else if(xExt>0) for(i=0;i<4;i++) 
            tmp[i]=polint(x3,4, x_grid+px, W->strfun+W->nx*(pq+i)+px); 
  else  
  {   double ftmp[4];
      int j;
      ftmp[0]=0;
      for(i=0;i<4;i++)
      {  for(j=1;j<4;j++) ftmp[j]= pow(x_grid[j],2/0.3)*W->strfun[W->nx*(pq+i)+j];
         tmp[i]= polint(x3,4,x_grid,ftmp)/(x*x);
      }         
  }
             
  if(qExt) return polint(loglogQ,4,q_grid+pq,tmp); 
      else return qSplineCteq6(loglogQ,q_grid+pq,tmp);
}

/*} MRST2001 INTERPOLATION CODES {*/
static int locx_(double *xx,int nx,double x)
{
   int jl,ju;
 
   if(x <= xx[0])    return 1;
   if(x >= xx[nx-1])  return  nx - 1;

   for(ju=nx,jl=0;ju - jl > 1;)
   {
      int jm = (ju + jl)/2;
      if(x >= xx[jm-1]) jl = jm; else ju = jm;
   }
   return jl;
} 

static double polderiv_(double x1,double x2,double x3,double y1,double y2,double y3)
{
 return (x3*x3*(y1-y2)-x2*2*(x3*(y1-y2)+x1*(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))
	/((x1-x2)*(x1-x3)*(x2-x3));
}

static double * jeppe1_(int nx,int my,double*xx,double*yy,double*ff)
{  
   int iwt[16*16]
   = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0,
       2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
       0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
       0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1,
       0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1,
      -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,
       9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2,
      -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2,
       2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0,
      -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1,
       4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1 
     };
    int j, k, l, m, n;

    double * cc=(double*)malloc(16*nx*my*sizeof(double));
    double * ff1=(double*)malloc(nx*my*sizeof(double));
    double * ff2=(double*)malloc(nx*my*sizeof(double));
    double * ff12=(double*)malloc(nx*my*sizeof(double));

    for(m=0;m<my;++m)
    {
       ff1[m*nx]=(ff[1+m*nx]-ff[m*nx])/(xx[1]-xx[0]);
       ff1[nx-1+m*nx]=(ff[nx-1+m*nx]-ff[nx-2+m*nx])/(xx[nx-1]-xx[nx-2]);

       for(n=1;n<nx-1;++n)
          ff1[n+m*nx]=polderiv_(xx[n-1],     xx[n],     xx[n+1],
		                ff[n-1+m*nx],ff[n+m*nx],ff[n+1+m*nx]);
    }

    for(n=0;n<nx;++n)
    {
	ff2[n]=(ff[n+nx]-ff[n])/(yy[1]-yy[0]);
	ff2[n+(my-1)*nx]=(ff[n+(my-1)*nx]-ff[n+(my-2)*nx])/(yy[my-1]-yy[my-2]);
	
	for(m=1;m<my-1;++m)
           ff2[n+m*nx]=polderiv_(yy[m-1],       yy[m],     yy[m+1],
	                         ff[n+(m-1)*nx],ff[n+m*nx],ff[n+(m+1)*nx]);
    }

    for(m=0;m<my;++m)
    {
	ff12[m*nx]=(ff2[1+m*nx]-ff2[m*nx])/(xx[1]-xx[0]);
	ff12[nx-1+m*nx]=(ff2[nx-1+m*nx]-ff2[nx-2+m*nx])/(xx[nx-1]-xx[nx-2]);
	for(n=1;n<nx-1;++n)
	   ff12[n+m*nx]=polderiv_(xx[n-1],      xx[n],      xx[n+1],
		                 ff2[n-1+m*nx],ff2[n+m*nx],ff2[n+1+m*nx]);
    }


    for(n=0;n<nx-1;++n)  for(m=0;m<my-1;++m)
    {   int s[4]; int p;
        
        double z[16],cl[16], d1=xx[n+1]-xx[n], d2=yy[m+1]-yy[m];

        s[0]=0;
        s[1]=1;
        s[2]=1+nx;
        s[3]=nx; 


	for(p=n+m*nx, k=0; k<4; ++k)
	{
           z[k]=     ff[p+s[k]];
           z[k+4]=  ff1[p+s[k]]*d1;
           z[k+8]=  ff2[p+s[k]]*d2;
           z[k+12]=ff12[p+s[k]]*d1*d2;
	}

	for(l=0;l<16;++l)
	{  double xxd=0;
           for(k=0;k<16;++k)xxd+=iwt[k+l*16]*z[k];
           cl[l]=xxd;
	}
	l=0;
        p=n+m*nx;	
	for(k=0;k<4;++k) for(j=0;j<4;++j)  cc[p+(k+4*j)*my*nx]=cl[l++];
    }
    free(ff1); free(ff2); free(ff12);
    return cc; 

}

static double jeppe2_(double x,double y,int nx, int my,
                   double *xx,double *yy,double *cc)
{
    int n = locx_(xx, nx, x);
    int m = locx_(yy, my, y);
    double t = (x-xx[n-1])/(xx[n]-xx[n-1]);
    double u = (y-yy[m-1])/(yy[m]-yy[m-1]);
    double z=0;
    int mn=nx*my,l;


    cc+=(n-1)+(m-1)*nx;

    for(l=3;l>=0;l--) 
        z=t*z+((cc[(l+12)*mn]*u+cc[(l+8)*mn])*u+cc[(l+4)*mn])*u+cc[l*mn];
    
    return z;
} 

static double int_mrst2001(double x, double q, pdtStr * W)
{
int i=W->qt0;
return  jeppe2_(log(x),log(q*q),W->nx, W->nq-i,
               /*  xxl  */  W->x_grid_aux    ,
               /*  qqlb */  W->q_grid_aux+i  ,
               /*  ccb  */  W->aux        )/x;
}

/*} END OF INTEPOLATION CODE */


double interFunc(double x, double q, pdtStr * W)
{ 
  if(W->q_grid)
  {  if(q<W->q_threshold) return 0; 
     if(q>W->q_max) W->nLargeQ++; else if(q<W->q_min) W->nSmallQ++;  
  }
  if(x<W->x_min)  W->nSmallX++;
  
  return (W->interpolation)(x,q,W); 
}


double interAlpha(double q, pdtStr * W )
{ 
  double logQ;
  int pq;
  if(!W->alpha) return -1.;
  logQ=log(q);
  pq=leftX(W->nq,W->q_grid,logQ); 
  if(pq>W->nq-3) pq = W->nq-3;
  return polint(logQ, 3, W->q_grid+pq, W->alpha+pq);
}


void freePdtData( pdtStr * data)
{
  if(data->x_grid)     {free(data->x_grid);      data->x_grid=NULL;     }
  if(data->q_grid)     {free(data->q_grid);      data->q_grid=NULL;     }
  if(data->alpha)      {free(data->alpha);       data->alpha =NULL;     }
  if(data->strfun)     {free(data->strfun);      data->strfun=NULL;     }
  if(data->strfun_aux) {free(data->strfun_aux);  data->strfun_aux=NULL; }
  if(data->q_grid_aux) {free(data->q_grid_aux);  data->q_grid_aux=NULL; }
  if(data->x_grid_aux) {free(data->x_grid_aux);  data->x_grid_aux=NULL; }
  if(data->aux)        {free(data->aux);         data->aux=NULL;        }
}


int getPdtData(char * file, int n_parton, pdtStr * data )
{ char pattern[20];
  char buff[100]; 
  char c;
  int  nx=0,nq=1;  
  FILE *f=fopen(file,"r"); 
  int errNo=0;
  
  if(!f) return -1;
  
  data->nq=0;
  data->nx=0;
  data->x_grid=NULL;
  data->q_grid=NULL;
  data->x_grid_aux=NULL;
  data->q_grid_aux=NULL;
  data->aux=NULL;
  data->alpha =NULL;
  data->strfun=NULL;
  data->strfun_aux=NULL;
  data->interpolation=&int_cteq4;
  data->q_threshold=0;
  data->mass=1;
  data->qt0=0;
 
  data->Q0_cteq=0.22;
  data->pow0=0;
  data->pow1=0;
  data->x_min=0;

  data->nSmallX=0;
  data->nSmallQ=0;
  data->nLargeQ=0;
  data->nLargeX=0;
  
  sprintf(pattern,"%d-parton",n_parton);

  while(1==fscanf(f,"%c",&c))   if(c=='#')
  { double qq; int i;
 
    fscanf(f,"%s",buff);
    if(!strcmp(buff,"Mass"))
    {  if(1!=fscanf(f,"%lf",&data->mass)) goto errexit;}  
    else if(!strcmp(buff,"Q_grid"))
    {  long fpos=ftell(f);
       if(data->q_grid || data->strfun) goto errexit;
       nq=0;
       while(fscanf(f,"%lf",&qq)) nq++;
       if(nq<3)  goto errexit; 
       data->nq=nq;
       data->q_grid=malloc(nq*sizeof(double));
       fseek(f,fpos,SEEK_SET);
       for(i=0;i<nq;i++) 
       { if(fscanf(f,"%lf", data->q_grid+i)!=1) goto errexit;
         if(i){ if(data->q_grid[i-1]>=data->q_grid[i]) goto errexit;}
         else if(data->q_grid[0]<=0)  goto errexit;
       }
    } 
    else if(!strcmp(buff,"Interpolation"))
    {                    
       fscanf(f,"%s",buff);
       { if(!strcmp(buff,"CTEQ4"))     data->interpolation=&int_cteq4;
         else if(!strcmp(buff,"CTEQ6"))
         {  data->interpolation= &int_cteq6;
            if(fscanf(f,"%lf",&(data->Q0_cteq))!=1)goto errexit;
         }
         else if(!strcmp(buff,"MRST2001")) data->interpolation=&int_mrst2001;
         else goto errexit; 
       }
    }     
    else if(!strcmp(buff,"X_grid"))
    {  long fpos=ftell(f);
       if(data->x_grid)  goto errexit; 
       while(fscanf(f,"%lf",&qq)) nx++; 
       data->nx=nx;
       if(nx<3)  goto errexit; 
       data->x_grid=malloc(nx*sizeof(double));
       fseek(f,fpos,SEEK_SET);
       for(i=0;i<nx;i++) fscanf(f,"%lf", data->x_grid+i);
       for(i=1;i<nx;i++) if(data->x_grid[i-1]>=data->x_grid[i])
                goto errexit;       
    }
    else if(!strcmp(buff,"Alpha"))    
    {  if(!data->q_grid)  goto errexit; 
       data->alpha=malloc(nq*sizeof(double));
       for(i=0;i<nq;i++)  
       if(fscanf(f,"%lf", data->alpha+i)!=1)  goto errexit; 
       if(fscanf(f,"%lf", &qq)==1)  goto errexit; 
    } else if(!strcmp(buff,pattern))
    {  int nn=nq*nx;

       data->strfun=malloc(nn*sizeof(double));
       for(i=0;i<nn;i++) 
       if(fscanf(f,"%lf", data->strfun+i)!=1) goto errexit; 
       if( fscanf(f,"%lf", &qq)==1)  goto errexit;
       
       break;
    } else if(!strcmp(buff,"q_threshold"))
         { if(fscanf(f,"%lf", &(data->q_threshold))!=1)  goto errexit; } 
      else if(!strcmp(buff,"x_min"))
         { if(fscanf(f,"%lf", &(data->x_min))!=1)  goto errexit; } 
  }
  
  if(data->strfun==NULL) {errNo=-2; goto errexit;}  

  if(data->x_min==0.) data->x_min=data->x_grid[0];
      
  if(data->q_grid)
  { int i;
    data->q_min=data->q_grid[0];
    data->q_max=data->q_grid[data->nq-1];
    for(i=0;i<nq;i++) 
    { if(data->q_grid[i]<data->q_threshold) data->qt0=i+1;
      data->q_grid[i]=log(data->q_grid[i]);
    } 
  }
       if(data->interpolation == &int_cteq4) 
  {if(!data->q_grid) {errNo=-3;goto errexit;}}
  else if(data->interpolation == &int_cteq6) 
  { int i; 
    if(!data->q_grid) { errNo=-3; goto errexit;}
    data->x_grid_aux=(double*)malloc(sizeof(double)*nx);
    if(data->x_grid[0]==0) data->x_grid_aux[0]=0; else data->x_grid_aux[0]=pow( data->x_grid[0], 0.3);
    for(i=1;i<nx;i++) data->x_grid_aux[i]=pow( data->x_grid[i], 0.3);
    data->q_grid_aux=(double*)malloc(sizeof(double)*nq);
    for(i=0;i<nq;i++) data->q_grid_aux[i]=log( data->q_grid[i]-log(data->Q0_cteq) ); 
  }
  else if(data->interpolation == &int_mrst2001)
  { int i;     

    double *f_aux=(double*)malloc(sizeof(double)*nx*nq);
    for(i=0;i<nx;i++)
    { double x=data->x_grid[i];
      int j;
      for(j=0;j<nq;j++)  f_aux[i+nx*j]=x*data->strfun[i+nx*j];
    }  
    

    if(!data->q_grid) {errNo=-3; goto errexit;}
    data->q_grid_aux=(double*)malloc(sizeof(double)*nq);
    for(i=0;i<nq;i++) data->q_grid_aux[i]=2*data->q_grid[i];
    data->x_grid_aux=(double*)malloc(sizeof(double)*nx);
    for(i=0;i<nx;i++) data->x_grid_aux[i]=log(data->x_grid[i]);
    if(data->qt0) data->qt0--;
    data->q_grid_aux[data->qt0]=2*log(data->q_threshold);
    i=data->qt0;
    data->aux=jeppe1_(data->nx,  data->nq-i,
                     /*  xxl */  data->x_grid_aux  ,
                     /* qqlb */  data->q_grid_aux+i , 
                     /* fb   */  f_aux  +i*data->nx);
    free(f_aux);
  }

  fclose(f); 
  return 0;
  errexit:
  { if(errNo==0) errNo=ftell(f); 
    fclose(f);   
    freePdtData(data);
    return errNo;
  } 
}


void delPdtList(pdtList * list)
{ 
  while(list)
  { pdtList * next=list->next;;
    free(list->file);
    free(list->name);
    free(list->partons);
    free(list->items);
    free(list);
    list=next; 
  }
}


long  makePdtList(char * file,  pdtList ** list)
{ char s[100];
  char dName[100];
  FILE * f=fopen(file,"r");  
  int partons_[100], positions_[100];
  int N,L;

  if(!f) return 0;
  while(fscanf(f,"%s",s) ==1 && strcmp(s,"#distribution")!=0) ;
  for(L=0;!feof(f);)
  { long mother;
    int bcount=0,pos;
    char ch;
    L=0;
    if(2!=fscanf(f," \"%[^\"]%*c %ld => ",dName,&mother)) break; 
    for(pos=1,ch=0;ch!='#';)
    {  if(1==fscanf(f," %d ",&N))
       {
         partons_[L]=N;   positions_[L]=pos; L++;       
         if(bcount==0) pos++; 
       }
       else 
       {  if(1!=fscanf(f,"%c",&ch)||!strchr("(#)",ch)
          ||(bcount &&strchr("(#",ch)) || (!bcount && ch==')') )goto exi;
          if(ch=='(') bcount=1; 
          else if(ch==')') { bcount=0; pos++;}
       }      
    } 
    {  pdtList * new=malloc(sizeof(pdtList));
       new->name=malloc(strlen(dName)+1);
       strcpy(new->name,dName);
       new->file=malloc(strlen(file)+1);
       strcpy(new->file,file);
       new->partons=malloc(sizeof(int)*(L+1));
       memcpy(new->partons,partons_,sizeof(int)*L);
       new->partons[L]=0;

       new->items=malloc(sizeof(int)*(L+1));
       memcpy(new->items,positions_,sizeof(int)*L);
       new->items[L]=0;
       new->beamP=mother;
       new->next=*list;
       *list=new;   
    }
    if(1!=fscanf(f,"%s",s)) break;
    if(strcmp(s,"distribution")==0) continue; else { fclose(f); return 0;} 
  }
exi:
  N=ftell(f);
  fclose(f); return N; 
}

int checkPartons( int * pNum, pdtList * L)
{ int *p;
  for(;*pNum; pNum++)
  { 
    if(*pNum==81 || *pNum==83) 
    {  for(p=L->partons ;*p;p++) if(*p==3) break; if(*p==0) return 0;
       for(p=L->partons ;*p;p++) if(*p==1) break; if(*p==0) return 0;
    } else if(*pNum==-81 || *pNum==-83)     
    {  for(p=L->partons ;*p;p++) if(*p==-3) break; if(*p==0) return 0;
       for(p=L->partons ;*p;p++) if(*p==-1) break; if(*p==0) return 0;
    }  else  
    {
      for(p=L->partons ;*p;p++) if(*pNum==*p) break;
      if(*p==0) return 0;
    }  
  }
  return 1;
}

