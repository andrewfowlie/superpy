#include <ctype.h>
#include <stdio.h>
#include<string.h>
#include<math.h>

void fName2c(char*f_name,char*c_name,int len)
{ int i; for(i=len-1;i>=0 &&f_name[i]==' ';i--);
  c_name[i+1]=0;
  for(;i>=0;i--) c_name[i]=f_name[i];
}

void cName2f(char*c_name,char*f_name,int len)
{ int i; for(i=0;i<len &&c_name[i];i++) f_name[i]=c_name[i];
         for(   ;i<len            ;i++) f_name[i]=' ';
}



static struct{ char name[10];double val;} nameInfo[500];
static int powInfo=0;


void clearinfo_(void){ powInfo=0;}

static double findValW(char*name)
{ 
  int i;
  for(i=0;i<powInfo;i++) if(strcmp(nameInfo[i].name,name)==0) return nameInfo[i].val;
  printf(" variable '%s' is not found\n",name);
  return 0; 
}

static void assignValW(char*name, double val)
{ int i;
  for(i=0;i<powInfo;i++) if(strcmp(nameInfo[i].name,name)==0) 
  { nameInfo[i].val=val;  return;} 
  if(powInfo==500) printf("No space in nameInfo\n"); 
  else
  { strcpy(nameInfo[powInfo].name,name);
    nameInfo[powInfo].val=val;
    powInfo++;
  }
}

void assignvalw_(char*f_name, double *val,int len)
{ 
  char buff[20];
  fName2c(f_name,buff,len);
  assignValW(buff, *val); 
}

double findvalw_(char*f_name,int len) 
{
  char c_name[20];
  fName2c(f_name,c_name,len);
  return findValW(c_name);
}

static int nfscanf(FILE*f,int *n1, int*n2, double* val )
{  
  long pos=ftell(f);
  for(;;)
  { char c;
    int r;

    
    fscanf(f,"%c",&c);
    
    if(c=='#') 
    { 
      do{ r=fscanf(f,"%c",&c);} while( r==1 && c!='\n');   
      if( r!=1 ) return r;
      continue;
    }
    else if(c==' ') 
    {
       if(val)
       {
         if(n1&&n2) r=fscanf(f,"%d %d %lf%*[^\n]",n1,n2,val);
         else if(n1)r=fscanf(f,"%d %lf%*[^\n]",n1,val);
         else r=fscanf(f,"%lf%*[^\n]",val);
       } else    r=fscanf(f,"%d%*[^\n]",n1);  
       fscanf(f,"%c",&c); 
       return r;
    }else
    { fseek(f,pos,SEEK_SET); 
      return 0;    
    }
  }  
}

int readlesh_(char *f_name,int len)
{
  FILE *f;
  char buff[100], name[20];
  int n1,n2,i,err=0;
  double val;
  char fname[100];
  char Zf[3][3]={"Zb","Zt","Zl"};
  int MG1ok=0, MG2ok=0, AmOK=0;

  fName2c(f_name,fname,len);
  f=fopen(fname,"r");
  if(f==NULL) return -2;
  for(;;) 
  { if(fscanf(f,"%s", buff)== EOF) break; 
    if(buff[0]=='#') { fscanf(f,"%*[^\n]\n"); continue;}
    for(i=0;buff[i];i++) buff[i]=toupper(buff[i]);
    if(strcmp(buff,"BLOCK")==0)
    { char rest[200];
      char *c;
      fscanf(f,"%s",buff);
      if(fscanf(f,"%[^\n]",rest))
      { c=strchr(rest,'#');
        if(c) c[0]=0;
        c=strchr(rest,'=');
        if(c && 1==sscanf(c+1,"%lf",&val))assignValW("QSUSY",val);
      } 
      fscanf(f,"%*c");       
      
      for(i=0;buff[i];i++) buff[i]=toupper(buff[i]);
      if(strcmp(buff,"MODSEL")==0)
      {
         for(;nfscanf(f,&n1,NULL,&val)==2;) 
            if(n1==1) assignValW("model",val);
      } 
      else if(strcmp(buff,"SMINPUTS")==0)
      {  
         for(;nfscanf(f,&n1,NULL,&val)==2;)
         {  
           switch(n1)
           {
             case         1   : assignValW("alfEMZ", 1/val); break;
             case         2   : assignValW("GF",val); break; 
             case         3   : assignValW("alfSMZ", val); break;  
             case         4   : assignValW("MZ",  val); break; 
             case         5   : assignValW("MbMb",val); break;  
             case         6   : assignValW("Mtp",val); break;  
             case         7   : assignValW("Ml",val); break;  
           } 
         }
      }   
      else if(strcmp(buff,"MINPAR")==0)
      { int model=findValW("model")+0.1; 
        while(2==nfscanf(f,&n1,NULL,&val))
        switch(model)
        { case 0: if(n1==3) assignValW("tb",val);    break;
          case 1: switch(n1)
          { case 1: assignValW("M0",val);
                    assignValW("gMHu",val);  assignValW("gMHd",val); 
                    assignValW("gMl2",val);  assignValW("gMl3",val);  
                    assignValW("gMr2",val);  assignValW("gMr3",val);
                    assignValW("gMq2",val);  assignValW("gMq3",val);  
                    assignValW("gMu2",val);  assignValW("gMd2",val);  
                    assignValW("gMu3",val);  assignValW("gMd3",val);
            break; 
            case 2: assignValW("Mhlf",val);  
                    assignValW("gMG1",val); assignValW("gMG2",val); assignValW("gMG3",val);
            break;
            
            case 3: assignValW("tb",val);    break;
            case 4: assignValW("sgn",val);   break;
            case 5: assignValW("A0", val);
                    assignValW("gAl",val);   assignValW("gAt",val);  assignValW("gAb",val);
            break;
          } break; 
          case 2: switch(n1)
          {  case 1: assignValW("Lambda",val);    break;
             case 2: assignValW("Mmess",val);    break;
             case 3: assignValW("tb",val);    break;
             case 4: assignValW("sgn",val);    break;
             case 5: assignValW("N5",val);
                     assignValW("N5_1",val);
                     assignValW("N5_2",val);
                     assignValW("N5_3",val); break;
             case 6:         
                     assignValW("Cgrav",val); break;
          } break;  
          case 3: switch(n1)
          { case 1: assignValW("M0",val); break; 
            case 2: assignValW("M32",val);  break;
            case 3: assignValW("tb",val);    break;
            case 4: assignValW("sgn",val);   break;
          } break;   
        }  
      }
      else if(strcmp(buff,"EXTPAR")==0) while(2==nfscanf(f,&n1,NULL,&val))    
          switch(n1)
          {  
           case    1   : assignValW("gMG1"  , val); break;  
           case    2   : assignValW("gMG2"  , val); break;  
           case    3   : assignValW("gMG3"  , val); break;
           case   11   : assignValW("gAt"   , val); break;  
           case   12   : assignValW("gAb"   , val); break;  
           case   13   : assignValW("gAl"   , val); break;  
           case   21   : if(val>0) assignValW("gMHd",sqrt(val));
                         else      assignValW("gMHd",-sqrt(-val));
                                                   break;
           case   22   : if(val>0)assignValW("gMHu", sqrt(val));  
                         else assignValW("gMHu", -sqrt(-val));break;
           case   23   : assignValW("mu",    val); break;
           case   26   : assignValW("MH3",   val); break;
           case   31   : assignValW("gMl1"  , val); break;  
           case   32   : assignValW("gMl2"  , val); break;  
           case   33   : assignValW("gMl3"  , val); break;  
           case   34   : assignValW("gMr1"  , val); break;  
           case   35   : assignValW("gMr2"  , val); break;  
           case   36   : assignValW("gMr3"  , val); break;  
           case   41   : assignValW("gMq1"  , val); break;  
           case   42   : assignValW("gMq2"  , val); break;  
           case   43   : assignValW("gMq3"  , val); break;  
           case   44   : assignValW("gMu1"  , val); break;  
           case   45   : assignValW("gMu2"  , val); break;  
           case   46   : assignValW("gMu3"  , val); break;  
           case   47   : assignValW("gMd1"  , val); break;  
           case   48   : assignValW("gMd2"  , val); break;  
           case   49   : assignValW("gMd3"  , val); break;
           
           case   51   : assignValW("N5_1"  , val); break;
           case   52   : assignValW("N5_2"  , val); break;
           case   53   : assignValW("N5_3"  , val); break;  
          }
    }          
  }

  fclose(f);
  return err;
}



int writelesh_(int * err, char*CODE,char *f_name,int len1,int len2)
{
  FILE *f;
  char name[20];
  int n1,n2;
  double Q;
  char fname[100];
  char codeTxt[100];
  char codeN[100];
  char codeV[100];
  
  fName2c(f_name,fname,len2);
  f=fopen(fname,"w");
  if(f==NULL) return -2;

  fName2c(CODE,codeTxt,len1);

  fprintf(f,"Block SPINFO\n");
  sscanf(codeTxt,"%s %[^\n]",codeN,codeV);
  fprintf(f,"   1  %s\n",codeN);
  fprintf(f,"   2  %s\n",codeV);

  if(*err<0) fprintf(f,"   3  %d\n",-(*err)  );
  if(*err>0) {fprintf(f,"   4  %d\n", (*err)  ); fclose(f); return 1;}

  Q=findValW("QSUSY");
  fprintf(f,"Block MASS   # Mass spectrum\n");
  fprintf(f,"#PDG code      mass           particle\n");
  fprintf(f,"       24   %16.8E    # MW\n",              80.423);
  fprintf(f,"       25   %16.8E    # h0\n",              findValW("Mh")); 
  fprintf(f,"       35   %16.8E    # H0\n",              findValW("MHH" )); 
  fprintf(f,"       36   %16.8E    # A0\n",              findValW("MH3" )); 
  fprintf(f,"       37   %16.8E    # H+\n",              findValW("MHc" )); 
  fprintf(f,"  1000001   %16.8E    # ~d_L\n",            findValW("MSdL")); 
  fprintf(f,"  1000002   %16.8E    # ~u_L\n",            findValW("MSuL")); 
  fprintf(f,"  1000003   %16.8E    # ~s_L\n",            findValW("MSsL")); 
  fprintf(f,"  1000004   %16.8E    # ~c_L\n",            findValW("MScL")); 
  fprintf(f,"  1000005   %16.8E    # ~b_1\n",            findValW("MSb1")); 
  fprintf(f,"  1000006   %16.8E    # ~t_1\n",            findValW("MSt1")); 
  fprintf(f,"  1000011   %16.8E    # ~e_L\n",            findValW("MSeL")); 
  fprintf(f,"  1000012   %16.8E    # ~nue_L\n",          findValW("MSne")); 
  fprintf(f,"  1000013   %16.8E    # ~mu_L\n",           findValW("MSmL")); 
  fprintf(f,"  1000014   %16.8E    # ~numu_L\n",         findValW("MSnm")); 
  fprintf(f,"  1000015   %16.8E    # ~stau_1\n",         findValW("MSl1")); 
  fprintf(f,"  1000016   %16.8E    # ~nu_tau_L\n",       findValW("MSnl")); 
  fprintf(f,"  1000021   %16.8E    # ~g\n",              findValW("MSG")); 
  fprintf(f,"  1000022   %16.8E    # ~neutralino(1)\n",  findValW("MNE1")); 
  fprintf(f,"  1000023   %16.8E    # ~neutralino(2)\n",  findValW("MNE2")); 
  fprintf(f,"  1000024   %16.8E    # ~chargino(1)\n",    findValW("MC1")); 
  fprintf(f,"  1000025   %16.8E    # ~neutralino(3)\n",  findValW("MNE3")); 
  fprintf(f,"  1000035   %16.8E    # ~neutralino(4)\n",  findValW("MNE4")); 
  fprintf(f,"  1000037   %16.8E    # ~chargino(2)\n",    findValW("MC2")); 
  fprintf(f,"  2000001   %16.8E    # ~d_R\n",            findValW("MSdR")); 
  fprintf(f,"  2000002   %16.8E    # ~u_R\n",            findValW("MSuR")); 
  fprintf(f,"  2000003   %16.8E    # ~s_R\n",            findValW("MSsR")); 
  fprintf(f,"  2000004   %16.8E    # ~c_R\n",            findValW("MScR")); 
  fprintf(f,"  2000005   %16.8E    # ~b_2\n",            findValW("MSb2")); 
  fprintf(f,"  2000006   %16.8E    # ~t_2\n",            findValW("MSt2")); 
  fprintf(f,"  2000011   %16.8E    # ~e_R\n",            findValW("MSeR")); 
  fprintf(f,"  2000013   %16.8E    # ~mu_R\n",           findValW("MSmR")); 
  fprintf(f,"  2000015   %16.8E    # ~stau_2\n",         findValW("MSl2")); 

  fprintf(f,"# Higgs mixing\n");
  fprintf(f,"Block ALPHA\n");
  fprintf(f,"     %16.8E\n",findValW("alpha"));

  fprintf(f,"Block GAUGE Q= %16.8E\n",Q);
  fprintf(f,"  1  %16.8E  # U(1)  coupling  \n", findValW("gY") );
  fprintf(f,"  2  %16.8E  # SU(2) coupling\n", findValW("g2") );
  fprintf(f,"  3  %16.8E  # SU(3) coupling\n", findValW("g3") );

  fprintf(f,"Block YU Q= %16.8E\n",Q);
  fprintf(f,"  3   3  %16.8E  # Yt\n", findValW("Yt") );
  fprintf(f,"Block YD Q= %16.8E\n",Q);
  fprintf(f,"  3   3  %16.8E  # Yb\n", findValW("Yb") );
  fprintf(f,"Block YE Q= %16.8E\n",Q);
  fprintf(f,"  3   3  %16.8E  # Yl\n", findValW("Yl") );
    
  fprintf(f,"Block HMIX Q= %16.8E\n",Q);
  fprintf(f,"  1  %16.8E  # mu(Q)MSSM DRbar\n",findValW("mu"));
  fprintf(f,"  2  %16.8E  # tan beta(Q)MSSM DRbar\n",findValW("tb_Q"));
  fprintf(f,"  3  %16.8E  # vev\n",findValW("vev"));
  fprintf(f,"  4  %16.8E  # mA^2\n",findValW("mA_2"));
  
  fprintf(f,"Block MSOFT Q= %16.8E\n",Q);
    fprintf(f,"    1   %16.8E # %s\n", findValW("MG1"  ),"MG1"  );  
    fprintf(f,"    2   %16.8E # %s\n", findValW("MG2"  ),"MG2"  );  
    fprintf(f,"    3   %16.8E # %s\n", findValW("MG3"  ),"MG3"  );  
    fprintf(f,"   21   %16.8E # %s\n", findValW("mH1_2"),"mH1_2");
    fprintf(f,"   22   %16.8E # %s\n", findValW("mH2_2"),"mH2_2");
    fprintf(f,"   31   %16.8E # %s\n", findValW("Ml1"  ),"Ml1"  );  
    fprintf(f,"   32   %16.8E # %s\n", findValW("Ml2"  ),"Ml2"  );  
    fprintf(f,"   33   %16.8E # %s\n", findValW("Ml3"  ),"Ml3"  );  
    fprintf(f,"   34   %16.8E # %s\n", findValW("Mr1"  ),"Mr1"  );  
    fprintf(f,"   35   %16.8E # %s\n", findValW("Mr2"  ),"Mr2"  );  
    fprintf(f,"   36   %16.8E # %s\n", findValW("Mr3"  ),"Mr3"  );  
    fprintf(f,"   41   %16.8E # %s\n", findValW("Mq1"  ),"Mq1"  );  
    fprintf(f,"   42   %16.8E # %s\n", findValW("Mq2"  ),"Mq2"  );  
    fprintf(f,"   43   %16.8E # %s\n", findValW("Mq3"  ),"Mq3"  );  
    fprintf(f,"   44   %16.8E # %s\n", findValW("Mu1"  ),"Mu1"  );  
    fprintf(f,"   45   %16.8E # %s\n", findValW("Mu2"  ),"Mu2"  );  
    fprintf(f,"   46   %16.8E # %s\n", findValW("Mu3"  ),"Mu3"  );  
    fprintf(f,"   47   %16.8E # %s\n", findValW("Md1"  ),"Md1"  );  
    fprintf(f,"   48   %16.8E # %s\n", findValW("Md2"  ),"Md2"  );  
    fprintf(f,"   49   %16.8E # %s\n", findValW("Md3"  ),"Md3"  );  
    
  fprintf(f,"Block STOPMIX\n");
  for(n1=1;n1<=2;n1++)for(n2=1;n2<=2;n2++)
  { sprintf(name,"Zt%d%d",n1,n2);
    fprintf(f,"  %d %d %16.8E # %s\n",n1,n2,findValW(name),name);
  }
  fprintf(f,"Block SBOTMIX\n"); 
  for(n1=1;n1<=2;n1++)for(n2=1;n2<=2;n2++)
  { sprintf(name,"Zb%d%d",n1,n2);
    fprintf(f,"  %d %d %16.8E # %s\n",n1,n2,findValW(name),name);
  }
  fprintf(f,"Block STAUMIX\n");
  for(n1=1;n1<=2;n1++)for(n2=1;n2<=2;n2++)
  { sprintf(name,"Zl%d%d",n1,n2);
    fprintf(f,"  %d %d %16.8E # %s\n",n1,n2,findValW(name),name);
  }
  fprintf(f,"Block NMIX\n");  
    for(n1=1;n1<=4;n1++)for(n2=1;n2<=4;n2++)
  { sprintf(name,"Zn%d%d",n1,n2);
    fprintf(f,"  %d %d %16.8E # %s\n",n1,n2,findValW(name),name);
  }
  fprintf(f,"Block UMIX\n");
  for(n1=1;n1<=2;n1++)for(n2=1;n2<=2;n2++)
  { sprintf(name,"Zu%d%d",n1,n2);
    fprintf(f,"  %d %d %16.8E # %s\n",n1,n2,findValW(name),name);
  }

  fprintf(f,"Block VMIX\n");
  for(n1=1;n1<=2;n1++)for(n2=1;n2<=2;n2++)
  { sprintf(name,"Zv%d%d",n1,n2);
    fprintf(f,"  %d %d %16.8E # %s\n",n1,n2,findValW(name),name);
  }

  fprintf(f,"Block AU Q= %16.8E\n",Q);
  fprintf(f,"  3 3 %16.8E # At\n",findValW("At"));
  fprintf(f,"Block AD Q= %16.8E\n",Q); 
  fprintf(f,"  3 3 %16.8E # Ab\n",findValW("Ab"));
  fprintf(f,"Block AE Q= %16.8E\n",Q);
  fprintf(f,"  3 3 %16.8E # Al\n",findValW("Al"));
/*  fprintf(f,"  2 2 %16.8E # Am\n",findValW("Am")); */
  fprintf(f,"Block SMINPUTS\n"); 
  fprintf(f,"  1  %16.8E # 1/alfEMZ\n",1/ findValW("alfEMZ")); 
  fprintf(f,"  2  %16.8E # G_Fermi\n",1.16637E-5);
  fprintf(f,"  3  %16.8E # alfSMZ\n", findValW("alfSMZ")); 
  fprintf(f,"  5  %16.8E # MbMb\n", findValW("MbMb"));
  fprintf(f,"  6  %16.8E # Mtp\n", findValW("Mtp"));
  fprintf(f,"  7  %16.8E # Mtau\n",findValW("Ml"));
  fprintf(f,"Block MODSEL\n");
  fprintf(f,"  1    0    # General MSSM\n");
  fprintf(f,"Block MINPAR\n"); 
  fprintf(f,"  3  %16.8E # tb\n",findValW("tb")); 
  fclose(f);
  return 0;
}
