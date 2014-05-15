#include<stdlib.h>
#include<math.h>
#include"phys_val.h"
#include"interface.h"
#include"histogram.h"
#include"../../include/num_out.h"
#include"num1.h"

#include "event2pyth.h"

hepeup_str E_;
heprup_str R_;
int nout_int;


double sum;


static void wrongParam(int N)
{ 
   if(N) fprintf(stderr,"Wrong parameter %d\n",N);
   else  fprintf(stderr,"Wrong number of parameters \n");  
   fprintf(stderr,  
    "Parameters:\n"
    " 1- name of variable,\n"
    " 2- minimum limit,\n"
    " 3- maximum limit,\n"
    " 4- number of bins(<=300).\n"
    "File with events must be passed to input. For example:\n"
    "   ../bin/lhe2tab \"M(6,-6)\" 1 100 200 < events.lhe >tab.txt\n"
    "So, PDG codes  should be used instead of particle names in \n"
    " standard CalcHEP functions.\n" );
}


static long nEvents=0;
static double totCS=0;

static int skipHeadLine(FILE* flow) {fscanf(flow,"%*[^\n]\n"); return -1;}
static int readCS(FILE* flow) {fscanf(flow,"%lf",&totCS); return 0;}
static int readNEvents(FILE* flow) {fscanf(flow,"%ld",&nEvents); return 0;}



int  main(int argc,char** argv)
{ 
  char varName[NAMELEN];
  int i,k,err;
    
  double minX, maxX;
  int nbin;
  char  key[4];
  physValRec * plist;
  
  double *hist, *dhist;
  double weight,coef;
  long nPoints;  
  double pvect[4*MAXNP];
  double Etot,pmiss[4]={0,0,0,0};
  double qmiss[4]={0,0,0,0};

  int outPos[50];
  char buff[200];

  pinf_int=pinf_ext;
  
 
  if(argc != 5 ) { wrongParam(0); return 1;} 
  
    
       
  if(sscanf(argv[2],"%lf",&minX)!=1){ wrongParam(2); return 1;}
  if(sscanf(argv[3],"%lf",&maxX)!=1 || minX>=maxX){ wrongParam(3); return 1;}
  if(sscanf(argv[4],"%d",&nbin)!=1  || nbin<=0)
    { wrongParam(4); return 1;}



  hist=(double*)malloc(nbin*sizeof(double));
  dhist=(double*)malloc(nbin*sizeof(double));
  for(i=0;i<nbin;i++){hist[i]=0; dhist[i]=0;}

  nPoints=0;
  
  openeventfile_(NULL, -1);
  err=readeventheader_();
     
  nEvents=0;
  while(readevent_()==0)
  { nEvents++;

    for(i=0;i<E_.NUP;i++) if(E_.ISTUP[i]==2)
    { double M=sqrt(fabs(E_.PUP[i][3]*E_.PUP[i][3]-E_.PUP[i][0]*E_.PUP[i][0]-E_.PUP[i][1]*E_.PUP[i][1]-E_.PUP[i][2]*E_.PUP[i][2])); 
        for(k=0;k<4;k++)
        { double dif=E_.PUP[i][k];
          int j;
          for(j=i+1;j<E_.NUP;j++) if(E_.MOTHUP[j][0]==i+1) dif-= E_.PUP[j][k];
          dif=fabs(dif)/M;
          if(dif> qmiss[k]) qmiss[k]=dif;
        }   
    }

    
    for(i=0, nin_int=0; i<E_.NUP; i++)if(E_.ISTUP[i]==-1)
    { double *pp=pvect+4*nin_int;
      for(k=1;k<3;k++) pp[k]=0;
      pp[0]=E_.PUP[i][3];
      pp[3]=E_.PUP[i][2];
      p_masses_[nin_int]=E_.PUP[i][4];
      p_codes_[nin_int]=E_.IDUP[i];
      nin_int++;
    }
      
    for(i=0, nout_int=0; i<E_.NUP; i++)if(E_.ISTUP[i]==1)outPos[nout_int++]=i;
    for(i=0;i!=nout_int-1;)
    { if(E_.IDUP[outPos[i]]>E_.IDUP[outPos[i+1]])
      { int id=outPos[i];
        outPos[i]=outPos[i+1];
        outPos[i+1]=id;
        if(i==0)i++; else i--;
       } else i++;   
    }

    for(i=0; i<nout_int; i++)
    { double* pp=pvect+4*(i+nin_int); 
      pp[0]=E_.PUP[outPos[i]][3];
      for(k=0;k<3;k++) pp[k+1]=E_.PUP[outPos[i]][k];
      p_masses_[nin_int+i]=E_.PUP[outPos[i]][4];
      p_codes_[nin_int+i]=E_.IDUP[outPos[i]];
    }
    
    for(i=0;i<nout_int+nin_int;i++) sprintf(p_names_[i],"%d",p_codes_[i]);

    Etot=pvect[0];
    for(i=1;i<nin_int;i++) Etot+=pvect[4*i];
    for(k=0;k<4;k++)
    { double dif=pvect[k];
      for(i=1;i<nin_int;i++) dif+=pvect[k+i*4];
      for(i=nin_int;i<nin_int+nout_int;i++) dif-=pvect[k+i*4];
      dif=fabs(dif)/Etot;
      if(dif>pmiss[k])  pmiss[k]=dif;
    }     

    
    if(checkPhysValN(argv[1], key, &plist))
    { double z0[100];
      int n0=0;
      physValRec * plist0=plist;
      for(;plist;plist=plist->next) z0[n0++]=calcPhysVal(key[0],plist->pstr,pvect);
      cleanPVlist( plist0);
      switch(key[1])
      { case '^':
         for(k=1;k<n0;k++) if(z0[0]<z0[k]) z0[0]=z0[k];
         n0=1;
         break;
        case '_':
         for(k=1;k<n0;k++) if(z0[0]>z0[k]) z0[0]=z0[k];
         n0=1;
         break;
      }    
      for(k=0;k<n0;k++)
      {
        i=nbin*(z0[k]- minX)/(maxX-minX); 
        if(i>=0 && i<nbin) {hist[i]+=1; nPoints++;}
      }
    }
  }
  
  coef= R_.XSECUP[0]* nbin/(maxX - minX)/nEvents;
  
  for(i=0;i<nbin;i++)
  { dhist[i]=sqrt(hist[i])*coef;
    hist[i]*=coef;
  }


  { char  xname[200], yname[200], xunits[100];
    fprintf(stdout,"#title ?");
/*    for(i=0;i<nin_int+nout_int;i++) 
    { if(i==nin_int) fprintf(stdout," ->"); else if(i) fprintf(stdout,",");
      fprintf(stdout," %s",pinf_int(1,i+1,NULL,NULL));
    } 
*/    
    xUnit(key[0], xunits);
    if(nin_int==2) sprintf(yname,"Diff. cross section [pb/%s]",xunits);
    else        sprintf(yname,"Diff. width [GeV/%s]",xunits);
    fprintf(stdout,"\n");
      

   fprintf(stdout,"#type 1  %%1d-histogram\n");
   fprintf(stdout,"#xName %s\n",argv[1]);
   fprintf(stdout,"#xMin %E\n",minX);
   fprintf(stdout,"#xMax %E\n",maxX);
   fprintf(stdout,"#xDim %d\n",nbin);
   fprintf(stdout,"#yName %s\n",yname);
   fprintf(stdout,"#lost_momenta_max/Etot %.1E %.1E %.1E %.1E\n",pmiss[0],pmiss[1],pmiss[2],pmiss[3]); 
   fprintf(stdout,"#lost_internal_momenta_max/M %.1E %.1E %.1E %.1E\n",qmiss[3],qmiss[0],qmiss[1],qmiss[2]);
  }
  for(i=0;i<nbin;i++) fprintf(stdout,"%-12E  %-12E\n",hist[i],dhist[i]);
  return 0;
}
