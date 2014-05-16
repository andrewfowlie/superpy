/*
 Copyright (C) 2009, Alexander Pukhov and Neil Christensen
*/
#include "interface.h"
#include<math.h>
#include"comp.h"
#include"parser.h"
#include"rd_num.h"
#include"crt_util.h"
#include"phys_val.h"
#include"read_func.h"
#include"sortarr.h"
#include "VandP.h"

int nComps=0,nCompParts[60];
table compTab={"*** Table ***"," Composites  ",
                   "  Name  |> Comma separated list of particles                           <|",NULL,0}; 
char compName[60][4], compParts[60][60][4];

int wrtcomp_(FILE *nchan){ fprintf(nchan,"\n"); writetable0(&compTab,nchan); return 0;}
int rdrcomp_(FILE *nchan){ 
  int rtrn=0;
  fscanf(nchan,"\n");  
  rtrn=readtable0(&compTab,nchan); 
  if(rtrn) return rtrn; 
  return fillCompositeArray(); 
}


int fillCompositeArray(void)
{
  linelist ln;
  int lineNum=0;
  int i,j,l,pos;
  char compStr[STRSIZ], partsStr[STRSIZ];
  char fieldName[50];
  char *chB;
  
  
  nComps=0;
  for(ln=compTab.strings ;ln;ln=ln->next)
  {  
    int k;
    compStr[0]=0;
    partsStr[0]=0;
    lineNum++;    
    sscanf(ln->line,"%[^|]%*c%[^|]",compStr,partsStr);
    trim(compStr); trim(partsStr);
/*============ Composite name ===========*/
    strcpy(fieldName,"Missing Name.");
    if(!strlen(compStr)) goto errorExit;
    nComps++;nCompParts[nComps-1]=0;
    sprintf(compName[nComps-1],"%s",compStr);
	 
/*================== Particle names ==========*/    
    strcpy(fieldName,"Missing list of particles.");
    if(!strlen(partsStr)) goto errorExit;
    chB=partsStr;
    i=0;
    chB--;
    for(j=0; chB ;chB=strchr(chB,','),j++)
  	 { char pname[100];
  	   int new=1,found=0;
    	chB++;
    	sscanf(chB,"%[^,)]",pname);
    	trim(pname);
    	//Check that pname is in model.
    	if(found==0)for(k=0;k<nModelParticles;k++)if(strcmp(pname,ModelPrtcls[k].name)==0)found++;
    	if(found==0)for(k=0;k<nModelParticles;k++)if(strcmp(pname,ModelPrtcls[k].aname)==0)found++;
    	if(found==0){sprintf(fieldName,"Cannot find %s in model.",pname);goto errorExit;}
    	//Check that pname is not a duplicate.
    	for(l=0;l<nCompParts[nComps-1];l++) if(strcmp(compParts[nComps-1][l],pname)==0)new=0;
    	if(new){
	    	sprintf(compParts[nComps-1][i++],"%s",pname);
   	 	nCompParts[nComps-1]++;
   	}
   	else{sprintf(fieldName,"%s is included more than once.",pname);goto errorExit;}
    }
  } 
  //Test
  /*for(i=0;i<nComps;i++){
  	fprintf(stderr,"%s=%s",compName[i],compParts[i][0]);
  	for(j=1;j<nCompParts[i];j++) fprintf(stderr,",%s",compParts[i][j]);
  	fprintf(stderr,"\n");
  }*/
  return 0;
  errorExit:
    sprintf(errorText," Error in Composite table line %d.\n%s",lineNum,fieldName);
    messanykey(2,10,errorText);
    return 1;  
}



