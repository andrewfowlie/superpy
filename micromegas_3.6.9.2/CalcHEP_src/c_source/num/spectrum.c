#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
 
#include "VandP.h"
#include "crt_util.h"
#include "file_scr.h"
#include "dynamic_cs.h"
#include "runVegas.h"
#include "vp.h"
#include "spectrum.h"


static void  writeSLHA(void)
{  int i;
   FILE *f;
   char fName[100];
   
   sprintf(fName,"decaySLHA_%d.txt",nSess);
          
   f=fopen(fName,"w");
   
   fprintf(f,"BLOCK ModelParameters\n");
   for(i=0;i<nModelVars;i++)
   fprintf(f," %3d  %16E # %s\n",i+1, (double)varValues[i], varNames[i]);    
   fprintf(f,"#\n");

   for(i=0;i<nModelParticles;i++)
   {  
    fprintf(f,"BLOCK QNUMBERS %d  # %s\n", ModelPrtcls[i].NPDG, ModelPrtcls[i].name);   
    fprintf(f," 1  %d # 3*el.charge\n 2  %d # 2*spin+1\n 3  %d # color dim\n 4  %d # 0={ self-conjugated}\n#\n",
         ModelPrtcls[i].q3, 
         ModelPrtcls[i].spin2+1, 
         ModelPrtcls[i].cdim, 
         strcmp(ModelPrtcls[i].name,ModelPrtcls[i].aname)? 1:0);
   }

   fprintf(f,"BLOCK MASS\n");   
   for(i=0;i<nModelParticles;i++) 
   { char *name=ModelPrtcls[i].name;
     fprintf(f," %d  %E # %s\n",ModelPrtcls[i].NPDG,pMass(name),name);
   }
   fprintf(f,"#\n");

   for(i=0;i<nModelParticles;i++)
   {  txtList all=NULL;
      double mass,width;
      char *name;
      
      if( strcmp(ModelPrtcls[i].mass,"0")==0) continue;
      if( strcmp(ModelPrtcls[i].width,"0")==0) continue;
            
      name=ModelPrtcls[i].name;
      mass=pMass(name);
      if(!mass) continue;
      width=pWidth(name,&all);
      fprintf(f,"DECAY %d  %E # %s\n",ModelPrtcls[i].NPDG,width,name);
      for(;all;all=all->next)
      { int dim; 
        char pn[20], buff[100], *chB,*chE;
        strcpy(buff,all->txt);
        sscanf(buff,"%s", pn);
        chB=strstr(buff,"->");
        chB+=2;
        for(dim=1; ;dim++)
        { 
           chE=strchr(chB,',');
           if(chE)chB=chE+1;else break;
        }
        fprintf(f," %s %d  ",pn,dim);
         
        strcpy(buff,all->txt);        
        chB=strstr(buff,"->");
        chB+=2;
        for(;;)
        { 
           chE=strchr(chB,',');
           if(chE)chE[0]=0;
           sscanf(chB,"%s",pn);
           fprintf(f," %d", qNumbers(pn,NULL,NULL,NULL));
           if(chE)chB=chE+1;else break;
        }
        
        
        chB=strstr(all->txt,"->");   
        fprintf(f,"  # %s \n",chB+2);
      }
      fprintf(f,"#\n");          
   }
   fclose(f);
   { char buff[100];
     sprintf(buff,"See results in file '%s'", fName);
     messanykey(16,5,buff); 
   }                  
}


void show_spectrum(int X, int Y)
{ int i;
  char *menuP=malloc(2+24*(1+nModelParticles));
  int mode=1;  
  menuP[0]=24;
  menuP[1]=0;
  
  sprintf(menuP+1," All Particles -> SLHA  ");
  
  for(i=0;i<nModelParticles;i++)
  { char *mass=ModelPrtcls[i].mass;
    char *name=ModelPrtcls[i].name;
    if(!strcmp(mass,"0")) sprintf(menuP+strlen(menuP)," %-6.6s      Zero       ",name);
    else sprintf(menuP+strlen(menuP)," %-6.6s  %14.4E ",name,pMass(name));
  }
  
  while(mode)
  {  menu1(X,Y,"",menuP,"n_qnumbers",NULL, &mode);
     if(mode==1) writeSLHA(); else
     if(mode)
     { FILE*f=fopen("width.tmp","w");
       txtList LL=NULL;
       char *mass=ModelPrtcls[mode-2].mass;
        
       fprintf(f, "Patricle %s(%s),  PDG = %d,  Mass= ", 
       ModelPrtcls[mode-2].name, ModelPrtcls[mode-2].aname,
       ModelPrtcls[mode-2].NPDG);
       if(strcmp(mass,"0")==0)  fprintf(f, "Zero\n"); else
       { 
         double width;
         fprintf(f,"%.3E ", pMass(ModelPrtcls[mode-2].name));       
         width=pWidth(ModelPrtcls[mode-2].name,&LL); 
         fprintf(f," Width=%.2E\n",width);
       }
       fprintf(f,"Quantum numbers: ");
       { int spin=ModelPrtcls[mode-2].spin2;
         int q3=ModelPrtcls[mode-2].q3; 
         fprintf(f," spin=");
         if(spin&1)  fprintf(f,"%d/2, ",spin); else fprintf(f,"%d, ",spin/2);
         fprintf(f," charge(el.)="); 
         if(q3!=3*(q3/3)) fprintf(f,"%d/3, ",q3);else fprintf(f,"%d ",q3/3); 
         fprintf(f," color=%d\n",ModelPrtcls[mode-2].cdim);
       }
       if(LL) 
       { txtList ll=LL;
          fprintf(f," Branchings & Decay channels:\n");
          for(;ll;ll=ll->next)
          { char buff[100];
            double br;
            sscanf(ll->txt,"%lf %[^\n]",&br,buff);
            fprintf(f," %.2E     %s\n",br,buff);
          }   
       }  
  
       fclose(f);
       f=fopen("width.tmp","r");
       showtext(2,Y,78,  maxRow()-1,"Particle information",f); 
       fclose(f);                   
       unlink("width.tmp");       
     }
  }   
  free(menuP);
  return;
}
