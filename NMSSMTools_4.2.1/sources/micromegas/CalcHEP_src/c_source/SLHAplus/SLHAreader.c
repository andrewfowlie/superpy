#include"SLHAplus.h"

#define StrLn 200
#define BlckLn 50
#define KeyMLn  12
#define DecLen  10

/*
   Different types of  information  are stored in  structures presented below. 
   Any item of information are stored in one  structure. The structures 
   allocated in computer memory and are linked via field 'next'. 
   The types of structures are described below.
*/ 


extern int  FError;

char * slhaComment="";

typedef   struct blockRec
{  struct blockRec* next;
   long double complex val;
   int nkey;
   int  keys[KeyMLn];
   char * txt;
   char * body;
} blockRec;


typedef struct blockStr
{ struct blockStr*  next;
  char name[BlckLn]; 
  double scale;
  char * txt;
  blockRec* dataList;
} blockStr;


typedef struct qNumberStr
{
  struct qNumberStr*  next;
  int pdg;
  int eQ3;
  int cDim;
  int spinDim;
  int anti; 
  char * txt;
} qNumberStr;

typedef   struct decayRec
{  struct decayRec* next;
   double Br;
   int nkey;
   int  pNum[DecLen];
   char * txt;
} decayRec;

typedef struct decayStr
{ struct decayStr*  next;
  int pNum; 
  double pWidth;
  decayRec* dataList;
  char * txt;
} decayStr;


static blockStr* blockList=NULL;
static decayStr* decayList=NULL;
static qNumberStr*qNumberList=NULL;

/*static*/ char*Warnings=NULL;
static int nWarnings=0;
static char creator[StrLn],version[StrLn];


static void cleanBlockRec(blockRec*List)
{
  while(List){ blockRec*list=List; List=List->next; free(list->txt); free(list->body); free(list);}  
} 

static void cleanDecayRec(decayRec*List)
{
  while(List){ decayRec*list=List; List=List->next; free(list->txt);   free(list);}  
} 

/*  release  memory allocated for SLHA data */
void cleanSLHAdata(void)
{ 

  while(blockList)
  { blockStr*block=blockList; blockList=blockList->next; 
    cleanBlockRec(block->dataList);
    free(block->txt);
    free(block); 
  }
  
  while(decayList)
  { decayStr*decay=decayList; decayList=decayList->next; 
    cleanDecayRec(decay->dataList);
    free(decay->txt);
    free(decay); 
  }
  
  while(qNumberList)
  {
    qNumberStr*qnb=qNumberList; qNumberList=qNumberList->next;
    free(qnb->txt);
    free(qnb);
  }
  
  if(Warnings) {free(Warnings); Warnings=NULL;}
  nWarnings=0;
  creator[0]=0;
  version[0]=0;
}

static void addWarning(char * txt)
{
  if(Warnings) Warnings=(char*)realloc(Warnings,strlen(Warnings)+strlen(txt)+4);
  else      {  Warnings=(char*)malloc(strlen(txt)+2); Warnings[0]=0;}
  sscanf(txt,"%[^\n]", Warnings+strlen(Warnings));
  strcat(Warnings,"\n");
   nWarnings++;
}            

static int nLine=0;


static FILE * _f_;
static char * _end_;
static int (*getLnPtr)(int, char*)=NULL;

/* procedure for reading new record of SLHA file 
   for the file open in C
*/

static int getLnC(int N, char * buff)
{ int len;
  char*s=fgets(buff,N-1,_f_);
  if(s==NULL) return -1;
  len=strlen(buff);   
  if(len>0 && buff[len-1]=='\n')
  {  
      buff[len-1]=0;
      if(len>1 && buff[len-2]==13) buff[len-2]=' ';
  }
                              
  return 0;
}

/* skipping of comment lines */

static int readLine(int size, char*buff)
{ int i;
  for(;;)
  {
    int err= getLnPtr(size-1,buff);  
    if(err==-1) return -1;

    nLine++;
    if(err==-2)
    {
       printf("Line %d is too long\n",nLine);  
       cleanSLHAdata();
       return nLine;
    }
    
    for(i=0;buff[i]&&buff[i]==' ';i++);
    if(buff[i]!=0 && buff[i]!='#') 
    { if(_end_ && buff+i == strstr(buff+i,_end_) )  return -1;

/*
      if(i)
      { int k;
        for(k=i;buff[k];k++) buff[k-i]=buff[k];
        buff[k-i]=0; 
      } 
*/           
      return 0;
    }    
  }
}

/* checking that new obtained line is started from keyword BLOCK or DECAY*/

static int findNewBlock(char* buff)
{   const char* block="BLOCK ";
    const char* decay="DECAY ";
    int err,i;
    for(;;)
    {  char * ch;
       err=readLine(StrLn,buff);
       if(err) return 0;
       
       for(ch=buff;*ch && *ch==' '; ch++);
       for(i=0;i<6 && toupper(ch[i])==block[i];i++);
       if(i==6) return 1;
       for(i=0;i<6 && toupper(ch[i])==decay[i];i++);
       if(i==6) return 1;                 
    }
}



static char * getComment(char * comm)
{ char * tret;
  if(comm) 
  { char * end;
    comm++; while(comm[0]==' ') comm++;
    end=comm+strlen(comm)-1;
     while(end>comm && end[0]==' ') end--;
    end[1]=0;            
  }else comm="";
  
  tret=malloc(strlen(comm)+1);
  strcpy(tret,comm);
  return tret;
}
  
/*
  mode= 1*m1+2*m2+4*m4+8*m8+16*m16

  m1  0 overwrite all;  1 keep old data
  m2  0 ignore mistake  1: stop in case of mistake in input file.
  m4  0 read DECAY      1: don't read   Decay 
  m8  0 read BLOCK      1: don't read   Blocks
  m16 0 read QNUMBERS   1: don't read   QNUMBERS 
*/

/* It is a main routine  for SLHA file reading. 
   Its second argument  present a function which 
   reads from a file. It could be function which reads 
   a line from file open in C, or a function which 
   read form a file open in Fortran.
   This functions itself is not presented for user call.
   But it is called from  slhaRead and slhaReadStream .
*/



int slhaBasicReader( int mode, int (*getLnPar)(int, char*), int *anydate,char * end)
{
  char buff[StrLn],name[StrLn],rest[StrLn];
  int n,err,m1,m4,m8,m2,m16;
  double scale;
  char wTxt[100];

  _end_=end;  
  *anydate=0;
  m1=mode&1;
  m2=mode&2;
  m4=mode&4;
  m8=mode&8;
  m16=mode&16;
  
  nLine=0;
  FError=0;
  getLnPtr=getLnPar;
  
  if(m1==0)cleanSLHAdata();
    
  err=readLine(StrLn,buff);
  if(err){FError=1;  if(err==-1 )  return -3 ; else  return err;}
  
  for(;;) 
  { char*block="BLOCK ";
    int L;
    char * bComm, *c;
    
    bComm=strchr(buff,'#'); if(bComm) bComm[0]=0;
 
    for(L=0;L<6 && buff[L]&& toupper(buff[L])==block[L] ;L++) continue;
    if(L==6)
    { int i;
      if( sscanf(buff+6," %s ",name)!=1)
      { 
         sprintf(wTxt,"SLHAreader: Line %d : block name is absent",nLine);
         addWarning(wTxt);
         if(m2) 
         { cleanSLHAdata();
           printf("%s\n",wTxt);
           addWarning(wTxt);
           FError=1;
           return nLine;
         }
         if(!findNewBlock(buff)) return -1;
         continue;   
      }  

      if(strlen(name)>=BlckLn)
      {
         sprintf(wTxt,"SLHAreader: Line %d: Too long name of BLOCK",nLine);
         addWarning(wTxt);
         if(m2) 
         {  cleanSLHAdata();
           printf("%s\n",wTxt);
           addWarning(wTxt);
            FError=1;
            return nLine;
         }
         if(!findNewBlock(buff)) return -1;
         continue; 
      }
      for(i=0;name[i];i++) name[i]=toupper(name[i]);
/* SPINFO & DCINFO */
      if(strcmp(name,"SPINFO")==0 ||strcmp(name,"DCINFO")==0 )
      { 
       *anydate=1;
       for(;;)
       { err=readLine(StrLn,buff); if(err){ return err;}
         if(sscanf(buff,"%d",&n)!=1)  break; 
         if(n==1) sscanf(buff,"%*d %[^\n]",creator);
         else if(n==2) sscanf(buff,"%*d %[^\n]",version);
         else if(n==4)
         {  cleanSLHAdata();  
            addWarning(buff);
            FError=1;                      
            return -2;
         }
         else if(n==3) addWarning(buff);
       }
      }
/* QNUMBERS */      
      else if( strcmp(name,"QNUMBERS")==0) 
      {  *anydate=1;
         if(m16) { if(!findNewBlock(buff)) return -1; else continue;} 
         else  
         {
           int val;
           qNumberStr * newQ;
           
           if( sscanf(buff+6,"%*s %d", &val)!=1) 
           { 
              sprintf(wTxt,"SLHAreader: line %d No pdg code for Qnumbers",nLine);
              addWarning(wTxt);
              if(m2)
              {  cleanSLHAdata();
                 printf("%s\n",wTxt);
                 addWarning(wTxt);
                 FError=1;
                 return nLine;
              }
              if(!findNewBlock(buff)) { return -1;}
              continue;                                                    
           }
           newQ =(qNumberStr*)malloc(sizeof(qNumberStr));
           newQ->next=qNumberList;
           qNumberList=newQ;
           newQ->pdg=val;
           newQ->txt=getComment(bComm);
           newQ->spinDim=newQ->eQ3=newQ->cDim=newQ->anti=-88888888;
           
           for(;;)
           { err=readLine(StrLn,buff); if(err) return err;
             if(sscanf(buff,"%d %d",&n,&val)!=2)  break;
             switch(n)
             {
               case 1: newQ->eQ3=val;      break;
               case 2: newQ->spinDim=val;  break;
               case 3: newQ->cDim=val;     break;
               case 4: newQ->anti=val;  break;
               default:
               {  
                  sprintf(wTxt,"SLHAreader: line %d: unexpected key for QNUMBERS",nLine);
                  addWarning(wTxt);  
               }
             }
           } 
         }   
      }
      else   
/* NORMAL BLOCK */      
      { blockStr*newBlock;
        if( sscanf(buff+6,"%*s %s", rest)==EOF) scale=-1;  
        else if(sscanf(buff+6,"%*s %*s %lf %s", &scale, rest)!=1)
        { 
           sprintf(wTxt,"SLHAreader: line %d: Unexpected BLOCK specification",nLine);
           addWarning(wTxt);
           if(m2)
           {  cleanSLHAdata();
              printf("%s\n",wTxt);
              addWarning(wTxt);
              FError=1;
              return nLine;
           }
           if(!findNewBlock(buff)) { return -1;}
           continue;       
        }

        *anydate=1;         
        if(m8) for(;;)
        {
          err=readLine(StrLn,buff); if(err) return err;
          if(!isdigit(buff[0])) break;
        }
        else
        {
          newBlock=(blockStr*)malloc(sizeof(blockStr));
          newBlock->next=blockList;
          strcpy(newBlock->name,name);         
          newBlock->dataList=NULL;
          newBlock->scale=scale;
          blockList=newBlock;
          newBlock->txt=getComment(bComm);
          
          for(;;)
          { int err=0,nkey=0,k,keys[KeyMLn];
            long double re=0,im=0;
            blockRec*dr;
            
            err=readLine(StrLn,buff); if(err) return err;

            if(buff[0]!=' ')    break;
            
            bComm=strchr(buff,'#'); if(bComm) bComm[0]=0;
            
            for(c=buff+strlen(buff)-1  ;c>buff && c[0]==' ';c--);
            
            c[1]=0;

            dr=(blockRec*)malloc(sizeof(blockRec));
            dr->next=blockList->dataList;
            blockList->dataList=dr;
            dr->val=0;
            dr->nkey=-1;
            dr->txt=getComment(bComm);
            dr->body=malloc(strlen(buff)+1);
            strcpy(dr->body,buff);
            
            if(c[0]==')') 
            {
              for(;c>buff && c[0]!='(';c--);
              if(2!=sscanf(c,"( %Lf , %Lf)", &re,&im))  err=1; 
            } else
            { char rest[StrLn];
              for(;c>buff && c[0]!=' ';c--);
              if(1!=sscanf(c,"%Lf%s",&re,rest))  err=1;
            }
            if(err) continue;
#ifdef OLD            
            { sprintf(wTxt,"SLHAreader: line %d: Unexpected last token " ,nLine);
              addWarning(wTxt);
              if(m2) { cleanSLHAdata(); printf("%s\n",wTxt); FError=1; return nLine; }
              else continue;                                                                       
            } 
#endif
            else 
            c[0]=0;
            
            for(c=strtok(buff," ");c && nkey<KeyMLn;c=strtok(NULL," "),nkey++)
            {   
              if(1!=sscanf(c,"%d%s",keys+nkey,rest)) { err=1; break;}
            }
          
            if(err) continue;
#ifdef OLD            
            {  
              sprintf(wTxt,"SLHAreader: line %d: Unexpected %d token" ,nLine,nkey+1);
              addWarning(wTxt);
              if(m2) { cleanSLHAdata(); printf("%s\n",wTxt); FError=1; return nLine;}
              else continue; 
            }
#endif
            dr->val=re+I*im;
            dr->nkey=nkey;
            for(k=0;k<nkey;k++) dr->keys[k]=keys[k];          
          }
        } 
      }
    }else 
/* DECAY */    
    { char  *decay="DECAY ";
      int pNum;
      double pWidth;
      decayStr*newDecay;
      
      for(L=0;L<6 && buff[L]&& toupper(buff[L])==decay[L] ;L++) continue;
      if(L!=6) 
      { sprintf(wTxt,"SLHAreader line %d: Unexpected first word ",nLine);
        addWarning(wTxt); 
        if(m2) 
        { cleanSLHAdata();
          printf("%s\n",wTxt);
          addWarning(wTxt);
          FError=1; 
          return nLine;
        } 
        if(!findNewBlock(buff)) {return -1;}
        continue;  
      }
      
      if( sscanf(buff+6,"%d %lf %s",&pNum ,&pWidth, rest)!=2)
      {  sprintf(wTxt,"SLHAreader line %d : Unexpected DECAY specification",nLine);
         addWarning(wTxt);
         if(m2) 
         { cleanSLHAdata();
           printf("%s\n",wTxt);
           addWarning(wTxt);
           FError=1;
           return nLine;
         }
         if(!findNewBlock(buff)) {return -1;}
         continue;             
      }

      *anydate=1;           
      if(m4)
      { 
        for(;;)
        { double x;
          err=readLine(StrLn,buff); if(err) return err;
          if(sscanf(buff,"%lf",&x)!=1) break;
        }  
      } else 
      {
        newDecay=(decayStr*)malloc(sizeof(decayStr));
        newDecay->next=decayList;
        newDecay->pNum=pNum;         
        newDecay->dataList=NULL;
        newDecay->pWidth=pWidth;
        newDecay->txt=getComment(bComm);
        decayList=newDecay;
        for(;;)
        { int i;
          int rdn;
          decayRec*dr;
          char end[StrLn];
          
          err=readLine(StrLn,buff); if(err) return err;
          if(isalpha(buff[0])) break;
          bComm=strchr(buff,'#'); if(bComm) bComm[0]=0;
          
          dr=(decayRec*)malloc(sizeof(decayRec));
          dr->next=newDecay->dataList;
          dr->nkey=0;
        
          c=strtok(buff," "); if(c) rdn=sscanf(c,"%lf%s",&(dr->Br),end); else rdn=0;
          if(rdn==1) 
          { c=strtok(NULL," "); 
            if(c) rdn=sscanf(c,"%d%s",&(dr->nkey),end); else rdn=0;
          }
          
          if(rdn==1 && dr->nkey<=DecLen) 
          for(i=0,c=strtok(NULL," ");i<dr->nkey && c&& rdn==1;i++,c=strtok(NULL," ")) 
          { 
             rdn=sscanf(c,"%d%s",dr->pNum+i,end);
          } else i=0;
          if( i!=dr->nkey || c || rdn!=1 )
          { 
            sprintf(wTxt,"SLHAreader line %d :Wrong decay record",nLine);
            addWarning(wTxt);
            free(dr);
            if(m2)
            {
              cleanSLHAdata();
              printf("%s\n",wTxt);
              addWarning(wTxt);
              FError=1; 
              return nLine;
            }   
          } else {dr->txt=getComment(bComm);  newDecay->dataList=dr;}
        }
      }    
    }     
  }
  return 0;
}         

/* see manual arXiv:1008.0181 */

int slhaRead(char *fname, int mode)
{ int err, anydate=0;
  _f_=fopen(fname,"r"); 
  if(_f_==NULL) { FError=1; return -1;}
  err=slhaBasicReader(mode,getLnC,&anydate,NULL);
  fclose(_f_);
  
  if((err==0 || err==-1) && anydate==0) {FError=1; return -3;}
  if(err==-1) err=0; 
  if(err) FError=1;
  
  return err;
}

/* see manual arXiv:1008.0181 */
int slhaReadStream(FILE*f,  int mode, char * end )
{  int err, anydate=0;
  _f_=f;

  err=slhaBasicReader(mode,getLnC,&anydate,end);
  
  if((err==0 || err==-1) && anydate==0) {FError=1; return -3;}
  if(err==-1) err=0; 
  if(err) FError=1;
  return err;
}

static long double complex* slhaValAddress(char * Block, int nKey, int *keys)
{
  char BLOCK[BlckLn];
  blockStr* blck=blockList;
  blockRec * dr;
  int i;
                            
  for(i=0; Block[i];i++)  BLOCK[i]=toupper(Block[i]);
  BLOCK[i]=0;
  
  while(blck && strcmp(BLOCK,blck->name)) blck=blck->next;
  if(!blck) return NULL;
  dr=blck->dataList;
  for(;dr;dr=dr->next)
  { if(dr->nkey!=nKey) continue;
    for(i=0;i<nKey;i++) if(keys[i]!=dr->keys[i]) break;
    if(i==nKey) return &dr->val; else continue;
  }
  return NULL;  
}

/* Below there are functions which print information stored in 
    memory by slhaRead . 
*/

static long double  complex cslhaVal_(char * Block, double Q, int nKey, va_list ap)
{ 
  int keys[KeyMLn];
  int i;
  char BLOCK[BlckLn];
  blockStr* blck=blockList;
  blockRec * dr;

  double  scale[3];
  long double complex val[3];
  int found[3]={0,0,0};

  if(strlen(Block)>19) {FError=1; return 0;}
      
  for(i=0;i<nKey;i++)keys[i]=va_arg(ap, int);

  for(i=0; Block[i];i++)  BLOCK[i]=toupper(Block[i]);
  BLOCK[i]=0;

  for(blck=blockList;  blck; blck=blck->next)if(strcmp(BLOCK,blck->name)==0)
  {  
     int pos;
     if(blck->scale < 0)    { if(found[0]) continue; else pos=0;}
     else if(blck->scale<Q) { if(found[1] && scale[1]>blck->scale) continue; else pos=1;}
     else                   { if(found[2] && scale[2]<blck->scale) continue; else pos=2;}

     dr=blck->dataList;
     for(;dr;dr=dr->next)
     { if(dr->nkey!=nKey) continue;
       for(i=0;i<nKey;i++) if(keys[i]!=dr->keys[i]) break;
       if(i==nKey) 
       {
          found[pos]=1;
          scale[pos]=blck->scale;
          val[pos]=dr->val; 
          break;
       }
     }  
  }
  
  if(found[0]==0 && found[1]==0 && found[2]==0) 
  { printf(" Block '%s', key={",BLOCK);
    for(i=0;i<nKey;i++) printf(" %d",keys[i]);
    printf("} - is absent\n"); 
    FError=1;
    return 0;
  }
  if(found[1]==0 && found[2]==0) return val[0];
  if(found[1]==0) return val[2];
  if(found[2]==0) return val[1];
  return  (val[1]*log(scale[2]/Q)+val[2]*log(Q/scale[1]))/log(scale[2]/scale[1]);  
}


double  slhaValFormat(char * Block, double Q, char * format)
{ 
  int keys[KeyMLn];
  int i;
  char BLOCK[BlckLn];
  blockStr* blck=blockList;
  blockRec * dr;

  double  scale[3];
  long double complex val[3];
  int found[3]={0,0,0};

  if(strlen(Block)>=BlckLn) {FError=1; return 0;}
      
  for(i=0; Block[i];i++)  BLOCK[i]=toupper(Block[i]);
  BLOCK[i]=0;

  for(blck=blockList;  blck; blck=blck->next)
  { 
  
  if(strcmp(BLOCK,blck->name)==0)
  {  char format_[BlckLn +4];
     int pos;
     if(blck->scale < 0)    { if(found[0]) continue; else pos=0;}
     else if(blck->scale<Q) { if(found[1] && scale[1]>blck->scale) continue; else pos=1;}
     else                   { if(found[2] && scale[2]<blck->scale) continue; else pos=2;}

     dr=blck->dataList;
     sprintf(format_," %s %%s",format);
     for(;dr;dr=dr->next)
     { int err;
       double v;
       char body_[StrLn+3];
       char buff[StrLn];
       sprintf(body_,"%s #",dr->body);
       err=sscanf(body_,format_,&v,buff);
       if(err==2 && strcmp(buff,"#")==0) 
       {  found[pos]=1; scale[pos]=blck->scale; val[pos]=v; slhaComment=dr->txt;  break; }
     }  
  }
  }
  if(found[0]==0 && found[1]==0 && found[2]==0) 
  { printf(" Block '%s', data of Format\"%s\" is absent\n", BLOCK,format); 
    FError=1;
    return 0;
  }
  if(found[1]==0 && found[2]==0) return val[0];
  if(found[1]==0) return val[2];
  if(found[2]==0) return val[1];
  return  (val[1]*log(scale[2]/Q)+val[2]*log(Q/scale[1]))/log(scale[2]/scale[1]);  
}

int   slhaSTRFormat(char * Block, char * format, char *txt)
{ 
  int keys[KeyMLn];
  int i;
  char BLOCK[BlckLn];
  blockStr* blck=blockList;
  blockRec * dr;

  if(strlen(Block)>=BlckLn) {FError=1; return 0;}
      
  for(i=0; Block[i];i++)  BLOCK[i]=toupper(Block[i]);
  BLOCK[i]=0;

  for(blck=blockList;  blck; blck=blck->next)
  { 
  
  if(strcmp(BLOCK,blck->name)==0)
  {  char format_[BlckLn +4];

     dr=blck->dataList;
     sprintf(format_," %s %%s",format);
     for(;dr;dr=dr->next)
     { int err;
       double v;
       char body_[StrLn+3];
       char buff[StrLn];
       sprintf(body_,"%s #",dr->body);
       err=sscanf(body_,format_,txt,buff);
//printf("body=|%s| format=|%s|,err=%d\n",body_,format_,err);
       if(err==2 && strcmp(buff,"#")==0) {return 0; }
     }  
  }
  }
  return 1;
}



double complex cslhaVal(char * Block, double Q, int nKey,...)
{
  double complex R;
  va_list ap;
  va_start(ap,nKey);
  R=cslhaVal_(Block,Q,nKey, ap);
  va_end(ap); 
  return R;
}

double slhaVal(char * Block, double Q, int nKey, ...)
{   
  double complex R;
  va_list ap;
  va_start(ap,nKey);
  R=cslhaVal_(Block, Q, nKey, ap);
  va_end(ap); 
  return creal(R);
}


int slhaValExists(char * Block, int nKey, ...)
{ 
  va_list ap; 
  int keys[KeyMLn];
  long double complex* address;
  int i;

  if(nKey>KeyMLn) return 0;  
  va_start(ap,nKey);
  for(i=0;i<nKey;i++)keys[i]=va_arg(ap, int);
  va_end(ap);
  
  address=slhaValAddress(Block,nKey, keys);
  if(address) return 1; else return 0;
}

int slhaWarnings(FILE*f)
{
  if(f&&Warnings) fprintf(f,Warnings);
  return nWarnings;
}


int slhaDecayExists(int pNum)
{ 
   decayStr* decay=decayList;
   for(;decay;decay=decay->next) if( decay->pNum==pNum)
   { decayRec*dr=decay->dataList;
     int n;
     for(n=0; dr; dr=dr->next, n++) continue;
     return n;
   }
   return -1;  
}

double slhaWidth(int pNum)
{  
   decayStr* decay=decayList;
   for(;decay;decay=decay->next) if( decay->pNum==pNum) return decay->pWidth;
   printf("Error: width for particle %d is unknown\n",pNum);
   FError=1;
   return 0;
}

double slhaBranch(int pNum,int N, int * nCh)
{
   decayStr* decay=decayList;
   for(;decay;decay=decay->next) if( decay->pNum==pNum)
   { decayRec*dr=decay->dataList;
     int i;
     if(N<=0){ nCh[0]=0; return 0;}
     for(N--;N&&dr;N--,dr=dr->next) continue;
     if(dr==NULL) { nCh[0]=0;return 0;}
     for(i=0;i<dr->nkey;i++) nCh[i]=dr->pNum[i];
     nCh[i]=0;
     return dr->Br;
   }
   FError=1;
   return 0;
}

double slhaBr(int pNum, int len, ...)
{
  int i,k, list[DecLen], buf[DecLen];
  decayStr*decay=decayList;

  if(len<2 || len > DecLen-1) return 0;
  va_list ap;
  va_start(ap,len);
  for(i=0;i<len;i++) list[i]=va_arg(ap, int);
  va_end(ap);

  for(k=0;k<len-1; )
  { if(list[k]>list[k+1])
    { int mem=list[k+1]; list[k+1]=list[k]; list[k]=mem;
      k--;
      if(k<0) k=1;
    } else k++;   
  }  
 
  for(;decay;decay=decay->next) if( decay->pNum==pNum)
  { decayRec*dr=decay->dataList;
    for(;dr;dr=dr->next)
    { if(dr->nkey != len) continue;
      for(i=0;i<dr->nkey;i++) buf[i]=dr->pNum[i];
      for(k=0;k<len-1; )
      { if(buf[k]>buf[k+1])
        { int mem=buf[k+1]; buf[k+1]=buf[k]; buf[k]=mem;
          k--;
          if(k<0) k=1;
        } else k++;   
      }  
      for(k=0;k<len;k++) if(buf[k]!=list[k]) break;
      if(k==len) return  dr->Br;
    }                           
  }                             
  return 0;                     
}


int findQnumbers(int pdg,int*eQ3,int*spinDim,int*cDim,int*anti)
{ qNumberStr*qNmb=qNumberList;
  for(; qNmb && abs(qNmb->pdg) != abs(pdg); qNmb=qNmb->next);
  if(!qNmb) return 0;

  *spinDim=qNmb->spinDim; 
  *eQ3=qNmb->eQ3; 
  *cDim=qNmb->cDim; 
  *anti=qNmb->anti;
    
  if(pdg!=qNmb->pdg) 
  {  if(*anti==0) return 0;
     (*eQ3)*=-1; if((*cDim)!=8) (*cDim)*=-1;
     return -1;
  }
  return 1;
}

int allQnumbers(int i, int *pdg,int*eQ3,int*spinDim,int*cDim,int*anti)
{ qNumberStr*qNmb=qNumberList;
  int i1;
  if(i<=0) return 0;
  for(i1=0,qNmb=qNumberList; qNmb; qNmb=qNmb->next) i1++;
  if(i1<i) return 0;  
  for(i1-=i,qNmb=qNumberList;i1; i1--) qNmb=qNmb->next;
  if(pdg) *pdg=qNmb->pdg;
  if(spinDim) *spinDim=qNmb->spinDim; 
  if(eQ3)  *eQ3=qNmb->eQ3; 
  if(cDim) *cDim=qNmb->cDim; 
  if(anti) *anti=qNmb->anti;
  slhaComment=qNmb->txt;
  return 1;
}


int allBlocks(int i, int j, char * name, int *keyLen, int * key, double complex*val)
{
  blockStr* bList=blockList;
  blockRec* bRec;
  int i1,j1,k;

  if(i<=0 || j<0) return 0;  
  for(i1=0, bList=blockList; bList; bList=bList->next) i1++;
  if(i>i1) return 0; 
  for(i1-=i, bList=blockList; i1; i1--) bList=bList->next;
  
  if(name)strcpy(name,bList->name);
  if(j==0){ if(keyLen)*keyLen=0; if(val)*val=bList->scale;  slhaComment=bList->txt;  return 1;}
   
  for(j1=0,bRec=bList->dataList; bRec; bRec=bRec->next) j1++;
  if(j>j1) return 0;
  for(j1-=j,bRec=bList->dataList; j1 ;j1--) bRec=bRec->next;
  if(keyLen)*keyLen=bRec->nkey;
  if(key)for(k=0;k<bRec->nkey;k++) key[k]=bRec->keys[k];
  if(val)*val=bRec->val;  
  slhaComment=bRec->txt;   
  return 1;  
}


int allDecays(int i, int j, int * pdg, int * decayLen, int * decay, double*width, double *br)
{
  decayStr* dList=decayList;
  decayRec* dRec;
  int i1,j1,k;
  
  if(i<=0 || j<0) return 0;
  for(i1=0,dList=decayList; dList; dList=dList->next) i1++;
  if(i>i1) return 0;
  for(i1-=i,dList=decayList; i1 ; i1--)dList=dList->next;
  if(pdg) *pdg=dList->pNum;
  if(width) *width=dList->pWidth;

  if(j==0){ if(decayLen) *decayLen=0; slhaComment=dList->txt; return 1;}
  
  for(j1=0,dRec=dList->dataList; dRec ;dRec=dRec->next)j1++;
  if(j>j1) return 0;
  for(j1-=j, dRec=dList->dataList; j1 ;j1--) dRec=dRec->next;
  if(decayLen) *decayLen=dRec->nkey;
  if(decay) for(k=0;k<dRec->nkey;k++) decay[k]=dRec->pNum[k];
  if(br) *br=dRec->Br; 
  slhaComment=dRec->txt;
  return 1;
}


int slhaWrite(char *fname)
{ 
  blockStr* block;
  blockRec* rec;
  decayStr* decay;
  decayRec* dr;
  qNumberStr*qList;
  int i,k,l,err;
  
  FILE*f=fopen(fname,"w");
  if(!f) return 1;
  if(!blockList  && ! decayList) return 2;
  
  if(blockList)
  {
    fprintf(f,"BLOCK SPINFO # General Information\n");
    fprintf(f," 1      %s\n",creator);
    fprintf(f," 2      %s\n",version);
  } else 
  {
    fprintf(f,"BLOCK DCINFO # General Information\n");
    fprintf(f," 1      %s\n",creator);
    fprintf(f," 2      %s\n",version);
  } 
  
   
  if(nWarnings) 
  { char buff[StrLn]; 
    char *  c;
    for(c=Warnings;;)
    {  
       sscanf(c,"%[^\n]",buff); 
       fprintf(f," 3 %s \n",buff);
       c=strchr(c,'\n');
       if(!c || c[1]==0) break;
       c++;
    } 
  } 
#define NEW
#ifdef NEW    
{
  blockStr* bList=blockList;
  blockRec* bRec;
  int i1,j1,k,l,i;


  for(i1=0, bList=blockList; bList; bList=bList->next) i1++;
  for(k=0;k<i1;k++) 
  { int i;
    bList=blockList;
    for(i=0;i<k;i++) bList=bList->next;
    fprintf(f,"BLOCK %s ",bList->name);
    if(bList->scale >0 ) fprintf(f,"Q= %E ",bList->scale);
    if(bList->txt) fprintf(f,"# %s\n",bList->txt); else fprintf(f,"\n");    

    for(j1=0,bRec=bList->dataList; bRec; bRec=bRec->next) j1++;
    for(l=0;l<j1;l++)
    { bRec=bList->dataList;
      for(i=0;i<l;i++) bRec=bRec->next;
      fprintf(f," %s",bRec->body);
      if(bRec->txt) fprintf(f," # %s\n",bRec->txt); else fprintf(f,"\n");
    }
  }  
}
#else
  for(k=1;;k++)
  {
    double complex val;
    double scale;
    char name[BlckLn];
    int len, key[KeyMLn];

    if(0==allBlocks(k,0,name,&len, key, &val)) break;
    
    fprintf(f,"BLOCK %s " ,name);
    scale=creal(val);
    if(scale>0) fprintf(f,"Q= %E ",scale);
    fprintf(f," # %s\n",slhaComment);

    for(l=1;;l++)
    {
       if(0==allBlocks(k,l,name,&len, key, &val))break;
       for(i=0;i<len; i++) fprintf(f," %3d",key[i]);
       if(cimagl(val)==0)   fprintf(f,"  %16LE ",creall(val));
                     else   fprintf(f,"  (%16LE,%16LE) ",creall(rec->val), cimagl(rec->val)); 
       fprintf(f," # %s\n",slhaComment);
    }            
  }

#endif 


  for(k=1;;k++)
  { int pdg,eQ3,spinDim,cDim,anti;
    if(0==allQnumbers(k,&pdg,&eQ3,&spinDim,&cDim,&anti))break;
    fprintf(f,"BLOCK QNUMBERS  %d # %s\n" , pdg,slhaComment);
    fprintf(f," 1  %d # 3*(electric charge)\n" , eQ3);
    fprintf(f," 2  %d # 2*spin+1 \n" , spinDim);
    fprintf(f," 3  %d # dimention of color \n" ,  cDim);
    fprintf(f," 4  %d # 0 for self conjugated\n" , anti);
  }

  for(k=1;;k++)
  { double width,br;
    int pdg,len, out[DecLen];
    if(0==allDecays(k,0, &pdg, &len, out, &width, &br))break;
    fprintf(f,"DECAY %d %E # %s\n" ,pdg,width,slhaComment); 
    for(l=1;;l++)
    {  if(0==allDecays(k,l, &pdg, &len, out, &width, &br))break;
       fprintf(f," %16E  %2d ", br,len);
       for(i=0;i<len; i++) fprintf(f," %10d",out[i]);
              fprintf(f," # %s\n",slhaComment);
    }            
  }  
  fclose(f);
  return 0;
}
