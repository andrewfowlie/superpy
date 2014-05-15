#include"syst.h"
#include"histogram.h"
#include"crt_util.h"
#include"edittab.h"
#include"syst.h"
#include"read_func.h"
#include"rd_num.h"
#include"interface.h"
#include"plot.h"
#include"subproc.h" 
#include"phys_val.h"


table histTab={"*** Table ***","Distributions",
                   "Parameter_1|> Min_1  <|> Max_1  <|"
                   "Parameter_2|> Min_2  <|> Max_2  <|",NULL,0};
                   
typedef  struct  histRec 
{ struct histRec * next;
  linelist  mother;
  long  nPoints;
  char key[2][4];               
  physValRec* pList[2];
  char title[2][50];
  double hMin[2],hMax[2];
  double f[900];
  double ff[900];   
} histRec;

static histRec * histPtr=NULL;

int clearHists(void)
{ histRec * hists=histPtr;
  int i;
  if(histPtr==NULL) return 0;
  while(hists)
  { 
    for(i=0;i<900;i++){hists->f[i]=0;hists->ff[i]=0;}
    hists->nPoints=0; 
    hists=hists->next;
  }
  return 1;
}


void  fillHists(double w,double*V)
{ histRec * hists=histPtr;
  int i0,i1,k0,k1,n0,n1,i,k;
  double z0[100],z1[100]; 
  physValRec*p;
    
  for(;hists;hists=hists->next)
  { 
    hists->nPoints++;
    p=hists->pList[0];
    if(!w || !p) continue;
    for(n0=0 ;p;p=p->next,n0++)  z0[n0]=calcPhysVal(hists->key[0][0],p->pstr,V);
    switch(hists->key[0][1])
    { case '^':
        for(k=1;k<n0;k++) if(z0[0]<z0[k]) z0[0]=z0[k];
        n0=1;
        break;
      case '_':
        for(k=1;k<n0;k++) if(z0[0]>z0[k]) z0[0]=z0[k];
        n0=1;
        break;
    }  
     
    if(hists->key[1][0]=='0')      
    { for(k=0;k<n0;k++)
      { 
        i=300*(z0[k] - hists->hMin[0])/(hists->hMax[0] - hists->hMin[0]);
        if(i<0 || i>=300) continue;
        hists->f[i]+=w;
        hists->ff[i]+=w*w;
      }
      continue;
    }

    p=hists->pList[1];
    if(!p) continue;
    for(n1=0 ;p;p=p->next,n1++)  z1[n1]=calcPhysVal(hists->key[1][0],p->pstr,V);
    switch(hists->key[1][1])
    { case '^':
        for(k=1;k<n1;k++) if(z1[0]<z1[k]) z1[0]=z1[k];
        n1=1;
        break;
      case '_':
        for(k=1;k<n1;k++) if(z1[0]>z1[k]) z1[0]=z1[k];
        n1=1;
        break;
    }  
    
    for(k0=0;k0<n0;k0++)
    { i0=30*(z0[k0]-hists->hMin[0])/(hists->hMax[0]-hists->hMin[0]);
      if(i0<0 || i0>=30) continue;
      for(k1=0;k1<n1;k1++)
      { i1=30*(z1[k1]-hists->hMin[1])/(hists->hMax[1]-hists->hMin[1]);
        if(i1<0 || i1>=30) continue;
        hists->f[30*i0+i1]+=w;
        hists->ff[30*i0+i1]+=w*w;
      }
    }
  }  
}

static int approxEq(double x1,double x2)
{ return fabs(x1-x2) <= 1.E-5*(fabs(x1)+fabs(x2)); }

int correctHistList(void)
{ 
   linelist ln;
   histRec * hptr;

   int lineNum;
   int i;
   double rMin[2]={0.,0.}, rMax[2]={0.,0.};
   char  histStr[2][STRSIZ], minStr[2][STRSIZ], maxStr[2][STRSIZ];
   char fieldName[50];    
   errorText[0]=0;
   for( hptr = histPtr ;hptr;hptr=hptr->next) hptr->mother=NULL;
   
   for(ln=histTab.strings,lineNum=1; ln; ln=ln->next,lineNum++)
   {
      sscanf(ln->line,"%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^|]%*c%[^\n]",
              histStr[0],minStr[0],maxStr[0],histStr[1],minStr[1],maxStr[1]);            

      for(i=0;i<2;i++)
      {
        trim(minStr[i]);trim(histStr[i]);trim(maxStr[i]);
/*============ Parameter ===========*/
        if(strlen(histStr[i])==0) 
        { if(i==0)
          { sprintf(fieldName,"Wrong field 'Parameter #%d'",i+1);
            goto errorExit;
          } else { rMin[i]=0;rMax[i]=0;} 
          continue; 
        }
/*================ MIN bound ============*/
        if(calcExpression(minStr[i],rd_num,&rMin[i] )) 
        { 
           sprintf(fieldName,"Wrong field 'Min. #%d bound'",i+1);
           goto errorExit;
        }
/*
        if(1!=sscanf(minStr[i],"%lf%c",rMin+i,&ch))  
        { sprintf(fieldName,"Wrong field 'Min. #%d bound'",i+1);
          goto errorExit;
        }
*/
/*================== MAX bound ==========*/
        if(calcExpression(maxStr[i],rd_num,&rMax[i] )) 
        {       
          sprintf(fieldName,"Wrong field 'Max. #%d bound'",i+1);
          goto errorExit;
        }  

/*
        if(1!=sscanf(maxStr[i],"%lf%c",rMax+i,&ch))  
        { sprintf(fieldName,"Wrong field 'Max. #%d bound'",i+1);
          goto errorExit;
        }
*/
      }
      
      for(hptr = histPtr; hptr; hptr=hptr->next)
       if(hptr->mother==NULL
         &&strcmp(hptr->title[0],histStr[0])==0
         && approxEq(hptr->hMin[0],rMin[0])
         && approxEq(hptr->hMax[0],rMax[0])
         &&strcmp(hptr->title[1],histStr[1])==0
         && approxEq(hptr->hMin[1],rMin[1])
         && approxEq(hptr->hMax[1],rMax[1]) 
         ) { int ok=0;
             cleanPVlist(hptr->pList[0]);
             cleanPVlist(hptr->pList[1]); 
             ok=checkPhysValN(histStr[0],hptr->key[0],hptr->pList); 
             if(strlen(histStr[1]))ok=ok&checkPhysValN(histStr[1],hptr->key[1],hptr->pList+1);
             else { strcpy(hptr->key[1],"0"); hptr->pList[1]=NULL;}
             if(ok) hptr->mother=ln; else hptr=NULL;  
             break; 
            } 
      if(hptr==NULL)
      { 
        physValRec *pL[2];
        char key[2][4];

        if(!checkPhysValN(histStr[0],key[0],pL ))
        { sprintf(fieldName,"Wrong field 'Parameter #1'");
          goto errorExit;
        }
        if(strlen(histStr[1]))
        { 
          if(!checkPhysValN(histStr[1],key[1],pL+1 )) 
          {  sprintf(fieldName,"Wrong field 'Parameter #2'");
             goto errorExit;
          }
        }  
        else {strcpy(key[1],"0") ;pL[1]=NULL;}  
        hptr=malloc(sizeof(histRec));
        hptr->next=histPtr;
        hptr->mother=ln;
        histPtr=hptr;
        for(i=0;i<2;i++)
        { strcpy(hptr->title[i], histStr[i]);
          strcpy(hptr->key[i],key[i]);
          hptr->pList[i]=pL[i]; 
          hptr->hMin[i]=rMin[i]; 
          hptr->hMax[i]=rMax[i];
        } 
        for(i=0;i<900;i++) {hptr->f[i]=0; hptr->ff[i]=0;}
        hptr->nPoints=0; 
      }  
   }      
   hptr = histPtr; 
   histPtr=NULL;  
   while(hptr) 
   { histRec *  hptr_=hptr;
     hptr=hptr->next;
     if(!hptr_->mother) 
     { 
       cleanPVlist(hptr_->pList[0]);
       cleanPVlist(hptr_->pList[1]);     
       free(hptr_);  
     } else {hptr_->next=histPtr; histPtr=hptr_;}
   }
   return 0;
 errorExit:
   { char buff[300];
     sprintf(buff,"Histograms:  Error in  line %d .\n%s. %s",lineNum,fieldName,errorText); 
     messanykey(2,10,buff);
   }  
   return 1;
}

static void writeDistributions(FILE*iprt)
{
  histRec * hists=histPtr;
  int i;
  linelist rec;
  
  for(rec=histTab.strings;rec;rec=rec->next)
  {
     for(hists=histPtr;hists;hists=hists->next) if(hists->mother == rec) 
     {  int Nprint;
        fprintf(iprt, " nPoints=%ld ",hists->nPoints);
        for(i=0;i<2;i++) 
        {  fprintf(iprt, " title=\"%s\"",hists->title[i]);
           fprintf(iprt," hMin=%-12E hMax=%-12E",hists->hMin[i],hists->hMax[i]); 
        }
        if(strcmp(hists->key[1],"0")==0) Nprint=300; else Nprint=900;
        fprintf(iprt,"\n  f: ");   
        for(i=0;i<Nprint;i++) fprintf(iprt," %-12E", hists->f[i]); 
     
        fprintf(iprt,"\n ff: ");
        for(i=0;i<Nprint;i++) fprintf(iprt," %-12E",hists->ff[i]);
     }
     fprintf(iprt,"\n");
  } 
}


static void readDistributions(FILE*iprt)
{
  int i;
  linelist rec;
  for(rec=histTab.strings;rec;rec=rec->next)
  {  int Nread; 
     histRec * hist=malloc(sizeof(histRec));

     hist->next=histPtr;
     hist->mother=rec;
     histPtr=hist;

     fscanf(iprt, " nPoints=%ld",&(hist->nPoints));

     for(i=0;i<2;i++)
     { int nR;
       nR=fscanf(iprt," title=\"%[^\"]",hist->title[i]);

       if(nR==0) strcpy(hist->title[i],"");
       fscanf(iprt,"\" ");       
       if(strlen(hist->title[i]))checkPhysValN(hist->title[i],hist->key[i],hist->pList+i);
       else{ strcpy(hist->key[i],"0");  hist->pList[i]=NULL;} 
       nR=fscanf(iprt,"hMin=%lf hMax=%lf\n",hist->hMin+i,hist->hMax+i); 
     }
    
     if(strcmp(hist->key[1],"0")==0) Nread=300; else Nread=900;     
     fscanf(iprt,"  f: ");
     for(i=0;i<Nread;i++) fscanf(iprt," %lf",hist->f+i); 
     fscanf(iprt," ff: ");
     for(i=0;i<Nread;i++) fscanf(iprt," %lf",hist->ff+i);
     for(;i<900;i++){ hist->f[i]=0; hist->ff[i]=0;}
  } 
  
}

int wrt_hist(FILE *nchan)
{  fprintf(nchan,"\n"); 
   writetable0(&histTab,nchan); 
   writeDistributions(nchan);
   return 0;
}

int rdr_hist(FILE *nchan)
{  fscanf(nchan,"\n"); 
   readtable0(&histTab,nchan);   
   readDistributions(nchan);
   return 0;
}


int wrt_hist2(FILE *fi, char * comment)
{ 
  int i,nameL1,nameL2; 
  linelist ll;

  char minRec1[200],minRec2[200],maxRec1[200],maxRec2[200],nameRec1[200],nameRec2[200];
  if(comment) fprintf(fi,"%s\n",comment); else  fprintf(fi,"\n");

  if(!(fi)) return -1;
  fprintf(fi,"%s\n",histTab.mdlName);
  fprintf(fi,"%s\n",histTab.headln); 

  sscanf(histTab.format,"%[^|]|%*[^|]|%*[^|]|%[^|]",nameRec1,nameRec2 );
 
  nameL1=strlen(nameRec1);
  nameL2=strlen(nameRec2);

  fprintf(fi,"%s|>    Min_1    <|>    Max_1    <|%s|>    Min_1    <|>     Max_1   <|\n",nameRec1,nameRec2); 

  for(ll=histTab.strings;ll;ll=ll->next) 
  { double V;
    char buff[20];

   sscanf(ll->line,"%[^|]|%[^|]|%[^|]|%[^|]|%[^|]|%[^|]",
               nameRec1,minRec1,maxRec1,nameRec2,minRec2,maxRec2);
   fprintf(fi,"%s|", nameRec1);
   calcExpression(minRec1,rd_num,&V);
   sprintf(buff,"%E",V);
   fprintf(fi,"%15.15s|", buff);
   calcExpression(maxRec1,rd_num,&V);
   sprintf(buff,"%E",V);
   fprintf(fi,"%15.15s|", buff);
     
   fprintf(fi,"%s|", nameRec2);
   trim(nameRec2);
   if(strlen(nameRec2)==0) fprintf(fi,"               |\n"); else
   {
      calcExpression(minRec2,rd_num,&V);
      sprintf(buff,"%E",V);
      fprintf(fi,"%15.15s|", buff);
      calcExpression(maxRec2,rd_num,&V);
      sprintf(buff,"%E",V);
      fprintf(fi,"%15.15s\n", buff);   
   }
  }
  
  for(i=0;i<nameL1+nameL1+60;i++) fprintf(fi,"=");
  fprintf(fi,"\n");

/*
   writetable0(&histTab,nchan); 
*/
   writeDistributions(fi);
   return 0;
}


int rdr_hist2(FILE *nchan, char **comment)
{  
    long pos1,pos2;
    pos1=ftell(nchan);
    fscanf(nchan,"%*[^\n]\n");
    pos2=ftell(nchan);
    
    if(comment)
    {  *comment=malloc(pos2-pos1);
       fseek(nchan,pos1,SEEK_SET);
       fscanf(nchan,"%[^\n]\n",*comment);
       pos2=strlen(*comment);
       if(pos2 && (*comment)[pos2-1]=='\n') (*comment)[pos2-1]=0;
    }
    else fscanf(nchan,"%*[^\n]\n"); 

   if(readtable0(&histTab,nchan)) return 1;
   
   readDistributions(nchan);
   return 0;
}


static int strcmp2(char*c1,char*c2)
{ 
  for(;;) if(*c1==' ') c1++;
  else    if(*c2==' ') c2++;
  else    if(*c1==*c2){if(*c1) { c1++; c2++;} else return 0;}
  else    return *c1-*c2; 
}


int add_hist(FILE *f, char **procname)
{  
  table histTab2=histTab;
  histRec * histPtr2=histPtr;
  histRec *r,*r2;
  char*procname2;
  int i,err=0; 
  linelist qLl;
  histRec*histPtr_;
  
  histPtr=NULL;
  histTab.strings=NULL;
  
  if(procname) err=rdr_hist2(f,&procname2); else err=rdr_hist2(f,NULL); 
  if(err) return err;

  if(!histPtr2){ if(procname) *procname=procname2; return 0; }
  if(procname)
  {
     *procname=realloc(*procname,strlen(*procname)+strlen(procname2)+5);
     strcat(*procname,";");
     strcat(*procname,procname2);
     free(procname2);
  }   
  if(!histPtr ){ histTab=histTab2; histPtr=histPtr2;  return 0;}

  if(strcmp(histTab.format,histTab2.format))
  { int l[6],l2[6];
    char *ch,*ch2;
    for(i=0,ch=histTab.format,ch2=histTab2.format; i<6; i++,ch++,ch2++)
    { ch =strchr(ch, '|');  l[i]=ch -histTab.format; 
      ch2=strchr(ch2,'|'); l2[i]=ch2-histTab2.format;
    } 
    
    for(i=5;i;i--) {l[i]-=l[i-1]+1; l2[i]-=l2[i-1]+1;}
    { int blind_=blind;
      char*inkeyString_=inkeyString;
      char buff1[20]="",buff2[20];
      blind=1;     
      for(i=0;i<6;i++,strcat(buff1,"\\09")) 
      { 
        if(l[i]<l2[i]) 
        { sprintf(buff2,"%s\\13%d{}",buff1,l2[i]);
          inkeyString=buff2;
          edittable(1,4,&histTab,1,"n_distrib",0); 
        }else if(l[i]>l2[i])
        { sprintf(buff2,"%s\\13%d{}",buff1,l[i]);
          inkeyString=buff2;
          edittable(1,4,&histTab2,1,"n_distrib",0); 
        }  
      }
      blind=blind_;
      inkeyString=inkeyString_;
      correctHistList();
      { 
        table histTab3=histTab;
        histRec * histPtr3=histPtr;
        histTab=histTab2;
        histPtr=histPtr2;
        correctHistList();
        histTab2=histTab;
        histPtr2=histPtr;
        histTab=histTab3;
        histPtr=histPtr3;    
      } 
    }  
  }
  for(r2=histPtr2;r2;r2=r2->next) for(r=histPtr;r;r=r->next)
  if(r->mother && strcmp2(r->mother->line,r2->mother->line)==0) 
  { double n1=r->nPoints;
    double n2=r2->nPoints;
    double n=n1+n2;
    int Ntot;
    if(n1)
    { if(strcmp(r->key[1],"0")==0) Ntot=300; else Ntot=900;
      if(n2) for(i=0;i<Ntot;i++) 
      {   
        r2->ff[i]=(n/n1)*(n/n1)*(r->ff[i] - r->f[i]*r->f[i]/n1  )
                 +(n/n2)*(n/n2)*(r2->ff[i]- r2->f[i]*r2->f[i]/n2);
        r2->f[i]=(n/n1)*r->f[i]+(n/n2)*r2->f[i];
        r2->ff[i]+=(r2->f[i])*(r2->f[i])/n;
      } 
      else for(i=0;i<Ntot;i++) 
      {   
        r2->f[i]=  r->f[i];
        r2->ff[i]=r2->ff[i];
      }
    }
    r2->nPoints=n;
    if(r->mother->next) r->mother->next->pred = r->mother->pred; 
    if(r->mother->pred) r->mother->pred->next = r->mother->next;
    else                histTab.strings=r->mother->next; 
    free(r->mother);
    r->mother=NULL;    
    break;
  }

  correctHistList();
  
  qLl=histTab.strings;
  if(qLl)
  { 
    for( ;qLl->next; qLl=qLl->next) continue;
    qLl->next=histTab2.strings;
    histTab2.strings=histTab.strings;  
    for(histPtr_=histPtr; histPtr_->next; histPtr_=histPtr_->next) continue;
    histPtr_->next=histPtr2;
    histPtr2=histPtr;
  }  
 
  histTab=histTab2;
  histPtr=histPtr2;
 
  return err;
}


static int nBinMenu(int X, int Y)
{                   

static int kBinMenu=3;
char   strmen[] =
   "\015"
   " 300         "
   " 150         "
   " 100         "
   "  75         "
   "  60         "
   "  50         "
   "  30         " 
   "  25         "
   "  20         "
   "  15         "
   "  12         "
   "  10         "
   "  6          "
   "  5          "
   "  4          "
   "  3          "
   "  2          ";

   void * pscr=NULL;
   
   int n;
   if(!kBinMenu) kBinMenu=3;
           
   menu1(X,Y,"number of bins",strmen,"",&pscr,&kBinMenu);
   if (kBinMenu)
   {
     sscanf(strmen+1+strmen[0]*(kBinMenu-1),"%d",&n);
     put_text(&pscr);
     return n;
   }
   return 0;
}

static int  nBinMenu2(int X, int Y, int*nb1, int*nb2)
{                   

static int nb1_=3, nb2_=3;
char   strmen[] =
   "\011"
   " 30      "
   " 15      "
   " 10      "
   "  6      "
   "  5      "
   "  3      "
   "  2      "
   "  1      ";

   void * pscr1=NULL;
   void * pscr2=NULL;

   if(!nb1_) nb1_=3;
   if(!nb2_) nb2_=3;
 
   menu1(X,Y,"N bin1",strmen,"",&pscr1,&nb1_);
   if (nb1_)
   {
     sscanf(strmen+1+strmen[0]*(nb1_-1),"%d",nb1);
     menu1(X+13,Y,"N bin2",strmen,"",&pscr2,&nb2_);
     put_text(&pscr1);
     if (nb2_) 
     { sscanf(strmen+1+strmen[0]*(nb2_-1),"%d",nb2); 
       put_text(&pscr2);
       return 1;
     }
   }
   return 0;
}

void xUnit(char key, char * units)
{
   switch(key)
   {
   case 'A': sprintf(units,"Deg");  break;
   case 'C': case 'J': case 'P': case 'Y': 
   case 'N': sprintf(units," ");     break; 
   
   case 'T': case 'E': case 'M': 
   case 'Z': sprintf(units,"GeV");  break;                   
   case 'S': sprintf(units,"GeV^2");break;
   default:  sprintf(units,"?");    break;  
   } 
} 


void showHist(int X, int Y,char *title)
{
   char  histStr1[STRSIZ],histStr2[STRSIZ];
   linelist ln=histTab.strings;
   char * menutxt;
   void * pscr=NULL;
   int mode =0;
   int npos=0;
   int width,width1,width2;
   int i,j;
   
   while(ln)
   {  npos++;
      ln=ln->next;     
   }
   if(!npos) return; 

   for(ln=histTab.strings,width1=0,width2=0; ln; ln=ln->next,npos++)
   { sscanf(ln->line,"%[^|]|%*[^|]|%*[^|]|%[^|]",histStr1,histStr2);
     trim(histStr1);
     trim(histStr2);
     { int l1=strlen(histStr1), l2=strlen(histStr2);
       if(width1<l1) width1=l1;
       if(width2<l2) width2=l2;
     }  
   }

   width=width1+width2+3;
   if(width<12) {width1=12-width2-3; width=width1+width2+3;}
   menutxt=malloc(2+width*npos);
   menutxt[0]=width;
   menutxt[1]=0;
   
   for(ln=histTab.strings; ln; ln=ln->next)
   {  
      sscanf(ln->line,"%[^|]|%*[^|]|%*[^|]|%[^|]",histStr1,histStr2);
      trim(histStr1);
      trim(histStr2);
      if(width2) sprintf(menutxt+strlen(menutxt)," %-*.*s| %-*.*s",
                   width1,width1,histStr1,width2,width2,histStr2);
      else       sprintf(menutxt+strlen(menutxt)," %-*.*s  ",
                   width1,width1,histStr1);
   }
   
   for(;;)
   {  
      menu1(X,Y,"Distributions",menutxt,"",&pscr,&mode);

      switch(mode)
      {
      case 0: free(menutxt);return;
      default: 
      {  histRec * hist=histPtr;
         int nBin1,nBin2;
         ln=histTab.strings;
         for(npos=1;npos<mode;npos++) ln=ln->next; 
         for( ;hist && hist->mother!= ln;hist=hist->next){;}
         if(hist)
         {  
            char xname[80],yname[80],units[80];                                    
                                                                                   
            if( hist->nPoints == 0) messanykey(10,10,"Distibution is empty");     
            else
            if(strcmp(hist->key[1],"0")==0)
            while(nBin1=nBinMenu(X,Y+4))                                                                   
            {  double f[300],df[300],coeff;                                        
               int i;
               coeff=nBin1/(hist->nPoints*(hist->hMax[0] - hist->hMin[0]));               

               for(i=0;i<nBin1;i++)                                                 
               {  int k;                                                           
                  f[i]=0;                                                          
                  df[i]=0;                                                         
                  for(k=0;k<300/nBin1;k++)                                          
	          {  f[i] += coeff*hist->f[i*300/nBin1+k];            
                    df[i] += coeff*coeff*hist->ff[i*300/nBin1+k];                   
                  }                                                                
                  df[i]=sqrt(fabs(df[i] - f[i]*f[i]/hist->nPoints));                      
               }                                                                                  
               if(nin_int==2) strcpy(yname,"Diff. cross section [pb");                 
               else       strcpy(yname,"Diff. width [GeV");                        
               xUnit(hist->key[0][0],units); 
               strcpy(xname,hist->title[0]);                          
               if(units[0]) { strcat(yname,"/");strcat(yname,units);}              
               strcat(yname,"]");                                                  
                                                                                   
               plot_1(hist->hMin[0],hist->hMax[0],nBin1,f,df,title,xname,yname);       
            } else
            while(nBinMenu2(X,Y+4,&nBin1,&nBin2))                                                                   
            {  double f[900],df[900],coeff;                                        
               int i,j;                     
               coeff=nBin1*nBin2/(hist->nPoints*(hist->hMax[0]-hist->hMin[0])
                                               *(hist->hMax[1]-hist->hMin[1]));               

               for(i=0;i<nBin1;i++)  for(j=0;j<nBin2;j++)
               {  int k,l,pos;
                  pos=i*nBin2+j;                                                           
                  f[pos]=0;                                                          
                  df[pos]=0;                                                         
                  for(k=0;k<30/nBin1;k++) for(l=0;l<30/nBin2;l++)                                         
	          {  f[pos] += coeff*hist->f[30*(i*30/nBin1+k)+j*30/nBin2+l];
                    df[pos] += coeff*coeff*hist->ff[30*(i*30/nBin1+k)+j*30/nBin2+l];                   
                  }
                  df[pos]=sqrt(df[pos] - f[pos]*f[pos]/hist->nPoints);                      
               } 
                                                                                                                                
               strcpy(xname,hist->title[0]);
               strcpy(yname,hist->title[1]);
               xUnit(hist->key[0][0],units);
               if(units[0])sprintf(xname+strlen(xname),"[%s]",units);
               xUnit(hist->key[1][0],units);
               if(units[0])sprintf(yname+strlen(yname),"[%s]",units);
                                                                                    
               plot_2(hist->hMin[0],hist->hMax[0],nBin1,
                      hist->hMin[1],hist->hMax[1],nBin2,
                       f,df,title,xname,yname);       
            }                                                                                  
         }
      }       
      }
   }
}

void editHist(void) {do  edittable(1,4,&histTab,1,"n_distrib",0); while(correctHistList());} 
