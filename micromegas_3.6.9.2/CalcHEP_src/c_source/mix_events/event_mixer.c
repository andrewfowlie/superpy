#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<unistd.h>
#include<time.h>
#include"SLHAplus.h"
#include"../../include/version.h"


#include "event2pyth.h"

static  char* pdg2name(int pdg)
{ static char buff[40];
  switch(pdg)
  { 
    case  2212: return "p"; 
    case -2112: return "anti-p";
    case   11: return "e-"; 
    case  -11: return "e+";
    case   13: return "mu-"; 
    case  -13: return "mu+";
    case   22: return "gamma";
    default: sprintf(buff,"PDG(%d)",pdg);
             return buff;
  }
}
             
static char * getField(char * fname, char *field)
{
  static char buff[200];
  FILE* f;
  buff[0]=0;

  f=fopen(fname,"r");
  if(!f) return buff;
  
  for(;;)
  { int r=fscanf(f,"%[^\n]\n",buff);
    if(r==0) { fscanf(f,"\n"); continue;}
    if(r!=1) { buff[0]=0; break;}
    if(strcmp(field, buff) < 0) { int k=strlen(field);
                                  for(;buff[k]==' ';k++) continue;
                                  fclose(f);
                                  return buff+k;
                                }
  }
  fclose(f);
  return buff;  
}

void  upinit_(void)
{ int II;

  R_.NPRUP= 1;
  R_.LPRUP[0]=1;
  for(II=0;II<2;II++) { R_.PDFSUP[II]=-1; R_.PDFGUP[II]=-1;}
  E_.AQEDUP=-1;
  E_.AQCDUP=-1;
  R_.IDWTUP = 3;
  R_.XERRUP[0] = 0.;
  R_.XMAXUP[0] = 1.;
}

int main(int argc,char ** argv)
{
  int N,NEV,MAXEVENTS,II,J,K,err;
  double cs;
  FILE *lun1=NULL;
  long posNevents, posSize;
  int SLHA=0;
  double Lumi;

  if(argc<4){ printf("%s needs arguments: Luminosity[1/fb]  nEvents event_directories ...\n",argv[0]); exit(1); }

  if(sscanf(argv[1],"%lf",&Lumi)!=1) 
            { printf("Can not read first argument. It should be luminosity in [1/fb] units.\n"); exit(1); }

  srand48(time(NULL));
  if(argv[2][0]=='?') NEV=-1; else
  { if(sscanf(argv[2],"%d",&NEV)!=1 || NEV<0)
    {printf("Can not read second argument. Number of events is expected.\n"); exit(1);}  
  } 
  if(access("decaySLHA.txt", R_OK)==0){ readslha_(); SLHA=1;}
  for(II=3;II<argc;II++) scandir_(argv[II],strlen(argv[II])); 
  if( Lumi <= 0 &&  NEV <=0 ) { writeInfo(); exit(0);}

  
  eventstat_(&cs, &MAXEVENTS); 
  if(cs==0) 
  {  printf("There are no events \n");
     exit(1);
  }   

  printf("%.3E  -total cross section[pb]\n", cs);
  printf("%d    -maximum number of events\n", MAXEVENTS);
  if(Lumi >0 && (NEV<=0 || NEV> 1000*Lumi*cs)) NEV=1000*Lumi*cs;

  lun1=fopen("event_mixer.lhe","w");
                                                                     
  upinit_();
      
  if(lun1)
  {
    fprintf(lun1,"<LesHouchesEvents version=\"1.0\">\n");
    fprintf(lun1,"<!--\n");
    fprintf(lun1,"File generated with CalcHEP-PYTHIA interface\n");
    fprintf(lun1,"-->\n");
    fprintf(lun1,"<header>\n");


  fprintf(lun1,"<hepml>\n");
  fprintf(lun1,"<samples xmlns=\"http://mcdb.cern.ch/hepml/0.2/\"\n");
  fprintf(lun1,"    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  fprintf(lun1,"    xsi:schemaLocation=\"http://mcdb.cern.ch/hepml/0.2/ http://mcdb.cern.ch/hepml/0.2/hepml.xsd\">\n");
  fprintf(lun1,"    <files>\n");
  fprintf(lun1,"        <file>\n");
  fprintf(lun1,"            <eventsNumber> ");
  posNevents=ftell(lun1);
  fprintf(lun1,"              </eventsNumber>\n");

  fprintf(lun1,"            <crossSection unit=\"pb\">%f</crossSection>\n",cs);
  fprintf(lun1,"            <fileSize> ");
  posSize=ftell(lun1);
  fprintf(lun1,"                 </fileSize>\n");
  
//  fprintf(lun1,"            <checksum type=\"md5\">d7722af9e2886fa465e3f7b786c3e6e6</checksum>\n");
  fprintf(lun1,"            <comments></comments>\n");
  
  fprintf(lun1,"            <location>\n");
  fprintf(lun1,"              <path/>\n");
  fprintf(lun1,"            </location>\n");
  
  fprintf(lun1,"        </file>\n");
  fprintf(lun1,"    </files>\n");
  
  fprintf(lun1,"    <description>\n");
  
  fprintf(lun1,"        <title> </title>\n");
  fprintf(lun1,"        <abstract></abstract>\n");
  fprintf(lun1,"        <authorComments></authorComments>\n");

  fprintf(lun1,"	<experimentGroup>\n");
  fprintf(lun1,"	    <experiment></experiment>\n");
  fprintf(lun1,"	    <group></group>\n");
  fprintf(lun1,"	    <responsiblePerson></responsiblePerson>\n");
  fprintf(lun1,"	    <description></description>\n");
  fprintf(lun1,"	</experimentGroup>\n");

  fprintf(lun1,"        <generator>\n");
  fprintf(lun1,"            <name>CalcHEP</name>\n");
  fprintf(lun1,"            <version> %s </version>\n",VERSION);
  fprintf(lun1,"            <homepage>http://theory.sinp.msu.ru/~pukhov/calchep.html</homepage>\n");
  fprintf(lun1,"            <description>\n");
  fprintf(lun1," CalcHEP - a package for calculation of Feynman diagrams\n"
               " integration over multi-particle phase space,\n"
               " generation of events and application of decay processes\n");       
  fprintf(lun1,"            </description>\n");
  fprintf(lun1,"        </generator>\n");

  fprintf(lun1,"        <model>\n");
//  fprintf(lun1,"             <name>  %s  </name>\n",getField("run_details.txt","Model:"));
  fprintf(lun1,"             <name></name>\n");
  fprintf(lun1,"             <description> </description>\n",getField("run_details.txt","Model:"));

  fprintf(lun1,"            <parameters>\n");
  if(SLHA) 
  { int i,j;
    char name[50];
    double complex val;
    for(i=1; allBlocks(i,0,name,NULL,NULL,NULL);i++) if(strcmp(name,"MODELPARAMETERS")==0)
    { for(j=1;allBlocks(i,j,name,NULL,NULL,&val);j++)
      { fprintf(lun1,"                <parameter>\n");
        fprintf(lun1,"                    <name>%s</name>\n",slhaComment);
        fprintf(lun1,"                    <value>%f</value>\n",creal(val) );
        fprintf(lun1,"			      <notation>\n");
        fprintf(lun1,"				      <plain></plain>\n");
        fprintf(lun1,"				      <Latex></Latex>\n");
        fprintf(lun1,"			      </notation>\n");
        fprintf(lun1,"			      <description></description>\n");
        fprintf(lun1,"                </parameter>\n");
      }
    } 
  } 
  fprintf(lun1,"            </parameters>\n");

  fprintf(lun1,"        </model>\n");

  fprintf(lun1,"        <process>\n");
  fprintf(lun1,"            <beam1>\n");
  fprintf(lun1,"                <particle KFcode=\"%d\">%s</particle>\n",R_.IDBMUP[0],pdg2name(R_.IDBMUP[0]));
  fprintf(lun1,"                <energy unit=\"GeV\">%f</energy>\n",R_.EBMUP[0]);
//  fprintf(lun1,"                <pdf name=\"CTEQ\" version=\"6m\" PDFLIBset=\"57\" PDFLIBgroup=\"4\" LHAPDFset=\"0\" LHAPDFmember=\"0\" LHAPDFfile=\"\" />\n");
  fprintf(lun1,"                <pdf name= \"%s\"></pdf>\n", getField("run_details.txt","pdf1 :"));  
//  fprintf(lun1,"                <QCDCoupling>\n");
//  fprintf(lun1,"                    <Lambda unit=\"GeV\">0.226200</Lambda>\n");
//  fprintf(lun1,"                    <NFlavours>6</NFlavours>\n");
//  fprintf(lun1,"                    <NLoopsAlphaS>2</NLoopsAlphaS>\n");
//  fprintf(lun1,"                </QCDCoupling>\n");

  fprintf(lun1,"            </beam1>\n");
 
  fprintf(lun1,"            <beam2>\n");
  fprintf(lun1,"                <particle KFcode=\"%d\">%s</particle>\n",R_.IDBMUP[1],pdg2name(R_.IDBMUP[0]));
  fprintf(lun1,"                <energy unit=\"GeV\">%f</energy>\n",R_.EBMUP[1]);
//  fprintf(lun1,"                <pdf name=\"CTEQ\" version=\"6m\" PDFLIBset=\"57\" PDFLIBgroup=\"4\" LHAPDFset=\"0\" LHAPDFmember=\"0\" LHAPDFfile=\"\" />\n");
  fprintf(lun1,"                <pdf name= \"%s\"></pdf>\n", getField("run_details.txt","pdf2 :"));  
  
//  fprintf(lun1,"                <QCDCoupling>\n");
//  fprintf(lun1,"                    <Lambda unit=\"GeV\">0.226200</Lambda>\n");
//  fprintf(lun1,"                    <NFlavours>6</NFlavours>\n");
//  fprintf(lun1,"                    <NLoopsAlphaS>2</NLoopsAlphaS>\n");
//  fprintf(lun1,"                </QCDCoupling>\n");
  fprintf(lun1,"            </beam2>\n");

//All Final state tags are required!!!
  fprintf(lun1,"	    <finalState>\n");
{   char*out=getField("run_details.txt","Process   :");
    out=strstr(out,"->");
    if(out)out+=2; else out="";
  fprintf(lun1,"		<state>%s\n",out);
  fprintf(lun1,"                </state>\n");
  fprintf(lun1,"		<notation>\n");
  fprintf(lun1,"		    <plain>%s\n",out);
}
   fprintf(lun1,"                   </plain>\n");
  fprintf(lun1,"		    <Latex></Latex>\n");
  fprintf(lun1,"		</notation>\n");
  fprintf(lun1,"	    </finalState>\n");

  fprintf(lun1,"            <crossSection unit=\"pb\">%f</crossSection>\n",cs);
 
  fprintf(lun1,"            <subprocesses>\n");
  printProcInfo(lun1);

//  fprintf(lun1,"                  <FactorisationScale>\n");
//  fprintf(lun1,"                      <plain></plain>\n");
//  fprintf(lun1,"                      <Latex></Latex>\n");
//  fprintf(lun1,"                  </FactorisationScale>\n");  

  fprintf(lun1,"            </subprocesses>\n");
  fprintf(lun1,"        </process>\n");
  
  fprintf(lun1,"        <cuts>\n");
  fprintf(lun1,"            <cutSet cutset_id=\"1\">\n");
  fprintf(lun1,"                <cut>\n");
  fprintf(lun1,"                    <object>\n");
  fprintf(lun1,"                        <name> </name>\n");
  fprintf(lun1,"                        <notation>\n");
  fprintf(lun1,"                            <plain> </plain>\n");
  fprintf(lun1,"                            <Latex></Latex>\n");
  fprintf(lun1,"                        </notation>\n");
  fprintf(lun1,"                    </object>\n");
  fprintf(lun1,"                    <minValue></minValue>\n");
  fprintf(lun1,"                    <maxValue></maxValue>\n");
  fprintf(lun1,"                    <logic></logic>\n");
  fprintf(lun1,"                </cut>\n");
  fprintf(lun1,"            </cutSet>\n");
  fprintf(lun1,"        </cuts>\n");
  fprintf(lun1,"        <authors>\n");
  fprintf(lun1,"       	   <author>\n");                                                                                                                            
  fprintf(lun1,"       	      <firstName>CalcHEP</firstName>\n");                                                                                                    
  fprintf(lun1,"              <lastName> </lastName>\n");                                                                                                     
  fprintf(lun1,"              <email>calchep[at]goolegroups.com</email>  \n");                                                                                      
  fprintf(lun1,"              <experiment></experiment>\n");                                                                                                    
  fprintf(lun1,"              <group></group>\n");                                                                                                    
  fprintf(lun1,"       	      <organization></organization>\n");                                                                                             
  fprintf(lun1,"   	   </author>\n");                                                                                                                            
  fprintf(lun1,"         </authors>\n"); 
//  fprintf(lun1,"        <relatedPapers></relatedPapers>\n");
  fprintf(lun1,"    </description>\n");
  fprintf(lun1,"</samples>\n");
  fprintf(lun1,"</hepml>\n");


    fprintf(lun1,"<slha>\n");
    if(SLHA) 
    { int i,j,k,pIn,pOut[20],len;
      int SM[16]={1,2,3,4,5,6,11,12,13,14,15,16,21,22,23,24};
      double width,br;  
      int pdg, eQ3, spinDim, cDim,neutral;

       for(i=1;allQnumbers(i, &pdg,&eQ3,&spinDim,&cDim,&neutral);i++)
       { for(k=0;k<16;k++) if(pdg==SM[k]) break;
         if(k==16) fprintf(lun1,"BLOCK QNUMBERS %d # %s\n"
                      " 1 %d\n 2 %d\n 3 %d\n 4 %d\n\n"   
         , pdg,slhaComment,eQ3,spinDim,cDim,neutral);
       }

       { char blkName[100];   
         int key[100];
         double complex val;
         int first=1;
         for(i=1;allBlocks(i, 0 , blkName, NULL ,NULL, NULL);i++) if(strcmp(blkName,"MASS")==0)
         { for(j=1; allBlocks(i,j, blkName, &len,key, &val);j++) if(len==1)
           { for(k=0;k<16;k++) if(key[0]==SM[k]) break;
             if(k==16)
             {  if(first) { fprintf(lun1,"BLOCK MASS\n"); first=0;}
                fprintf(lun1," %d    %E # \n", key[0], creal(val));
             }
           }              
           break;
         }
         if(!first) fprintf(lun1,"\n");
      }    
      for(i=1; allDecays(i,0,&pIn, &len,pOut,&width,&br);i++)
      { for(k=0;k<16;k++) if(pIn==SM[k]) break;
        if(k==16)
        {  fprintf(lun1,"DECAY %d  %E # %s \n",pIn,width,slhaComment);
           for(j=1; allDecays(i,j,&pIn, &len,pOut,&width,&br);j++)
           { fprintf(lun1," %E  %d ",br, len);
             for(k=0;k<len;k++) fprintf(lun1," %d ",pOut[k]);
             fprintf(lun1," # %s \n",slhaComment);
           }
        }
      }
    }
        
    fprintf(lun1,"</slha>\n");

    if(access("run_details.txt", R_OK)==0) 
    {  FILE * f=fopen("run_details.txt","r");
       int ch;
       for(;;)
       { ch=fgetc(f);
         if(ch== EOF) break;
         fputc(ch,lun1);
       }
      fclose(f);
    }
    
    fprintf(lun1,"</header>\n");	
    fprintf(lun1,"<init>\n");
    fprintf(lun1," %5d %5d %18.11E %18.11E %5d %5d %5d %5d %5d %5d\n",
      R_.IDBMUP[0],R_.IDBMUP[1],R_.EBMUP[0], R_.EBMUP[1], 
      R_.PDFGUP[0],R_.PDFGUP[1],R_.PDFSUP[0],R_.PDFSUP[1],
      R_.IDWTUP,R_.NPRUP);

    fprintf(lun1," %18.11E %18.11E %18.11E %3d\n",
             cs,R_.XERRUP[0],R_.XMAXUP[0], R_.LPRUP[0]); 
    fprintf(lun1,"</init>\n");
  }
 
  for(N=1;N<=NEV;N++)
  {
    while(upevnt_());
    if(E_.NUP==0){printf("Only %d events are generated\n", N-1); break;}
    
    if(lun1)
    {
      fprintf(lun1,"<event>\n");
      
      fprintf(lun1,"%2d %4d %15.7E %15.7E %15.7E %15.7E\n",
          E_.NUP,E_.IDPRUP,E_.XWGTUP,E_.SCALUP,E_.AQEDUP,E_.AQCDUP);
     
      for(II=0;II<E_.NUP;II++)
      { double s;
        fprintf(lun1," %8d %4d",  E_.IDUP[II],E_.ISTUP[II]);
        for(J=0;J<2;J++) fprintf(lun1," %4d", E_.MOTHUP[II][J]);
        for(J=0;J<2;J++) fprintf(lun1," %4d", E_.ICOLUP[II][J]);
        for(J=0;J<4;J++)  fprintf(lun1," %18.11E",E_.PUP[II][J]);
        if(E_.ISTUP[II]==2)
        { 
          s= E_.PUP[II][3]*E_.PUP[II][3] 
          - E_.PUP[II][0]*E_.PUP[II][0] -E_.PUP[II][1]*E_.PUP[II][1]-E_.PUP[II][2]*E_.PUP[II][2];
          if(s<0) s=-sqrt(-s); else s=sqrt(s); 
        } else s=E_.PUP[II][4]; 
         fprintf(lun1," %18.11E",s);
         fprintf(lun1," %11.4E %4.1f\n",E_.VTIMUP[II],E_.SPINUP[II]);
      }
	
      fprintf(lun1,"</event>\n");
    }   
  }       

  if(lun1)
  { long size;
    fprintf(lun1,"</LesHouchesEvents>\n");
    size=ftell(lun1);
    fseek(lun1, posNevents,SEEK_SET);
    fprintf(lun1,"%d",N-1);
    fseek(lun1,posSize,SEEK_SET);
    fprintf(lun1,"%ld", size);
    fclose(lun1);   
  }  
  closeevents_();
  return 0;
}
