
#include"micromegas_aux.h"
#include"micromegas.h"

ModelPrtclsStr *OddPrtcls=NULL;
int Nodd=0;

int createTableOddPrtcls(void)
{
   int i;
   if(OddPrtcls) return 0;

   for(i=0,Nodd=0;i<nModelParticles;i++)
   { int o1,o2;
     o1=(ModelPrtcls[i].name[0]=='~');
     o2=(ModelPrtcls[i].aname[0]=='~');
     if(o1&&o2) Nodd++;
     if( (o1&&!o2) || (o2&&!o1) ) return 1;
   }    

   OddPrtcls=( ModelPrtclsStr*)malloc(Nodd*sizeof(ModelPrtclsStr));
   Nodd=0;
   for(i=0,Nodd=0;i<nModelParticles;i++) if(ModelPrtcls[i].name[0]=='~')
         OddPrtcls[Nodd++]=ModelPrtcls[i];
   return 0;
}

void printMasses(FILE * f,int sort)
{
if(f==NULL) return;

fprintf(f,"\nMasses of odd sector Particles:\n");
{ int i,col;
  int *n=malloc(sizeof(int)*Nodd);
  int *nn=malloc(sizeof(int)*Nodd);
  double * mass=malloc(sizeof(double)*Nodd);
  int pow;

  for(pow=0,i=0; i<Nodd; i++) 
  { double  v=findValW(OddPrtcls[i].mass); 
    mass[i]=fabs(v); 
    n[pow]=i;
    nn[pow]=i;
    pow++;
  }
  if(pow>1)
  for(i=0;i<pow-1;)
  {  int i1=i+1;
     if( mass[n[i]] > mass[n[i1]])
     { int k=n[i]; n[i]=n[i1]; n[i1]=k;
       if(i==0) i++; else i--;
     } else i++;
  }

  for(col=0,i=0;i<pow;i++)
  { int k;
    if(sort)k=n[i];else k=nn[i];
  
    fprintf(f,"%-4.4s : %-6.6s= %7.1f ", OddPrtcls[k].name,
           OddPrtcls[k].mass,mass[k]);
                        
    col++;
    if(f)
    { if(col==1 || col==2) fprintf(f,"|| ");
      if(col==3) { col=0; fprintf(f,"\n");}
    }
  }
  fprintf(f,"\n");
  free(n); free(nn); free(mass);

}
}


char * nextOdd(int N,double * Mass)
{
  int i,nn, *n;
  double * mass;
  if(N>Nodd) return NULL;
  n=malloc(sizeof(int)*Nodd);
  mass=malloc(sizeof(double)*Nodd);
  
  for(i=0; i<Nodd; i++) 
  { 
    mass[i]=fabs(findValW(OddPrtcls[i].mass)); 
    n[i]=i;
  }
  for(i=0;i<Nodd-1;)
  {  int i1=i+1;
     if( mass[n[i]]  >  mass[n[i1]] ) 
     { int k=n[i]; n[i]=n[i1]; n[i1]=k;
       if(i==0) i++; else i--;
     } else i++;
  }

  nn=n[N];
  if(Mass) *Mass= mass[nn];  
  free(n); free(mass); 
  
  return OddPrtcls[nn].name;
}
