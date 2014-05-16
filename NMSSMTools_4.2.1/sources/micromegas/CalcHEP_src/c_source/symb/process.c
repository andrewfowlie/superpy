/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include "chep_crt.h"
#include "s_files.h"
#include "screen.h"
#include "file_scr.h"
#include "read_mdl.h"
#include "process.h"

whohow     liminsp, LimQ;
whohow     limout;

int ZWmax=1000, ZWmin=-1000;

int  nin, nout, n_x;   /* Number of X-particles */
shortstr  processch="";
char limpch[STRSIZ]="", deloutch[STRSIZ]="";

hadron hadrons[MAXINOUT];

#define ycons 19
#define errtxt  "This particle is absent in the model"


void  nilprtcl(whohow p_list)
{int  i; 
   for (i = 0; i < whohowMAX; i++) {p_list[i].who = 0; p_list[i].how = 0;}  
} 

int polarized(int p, int Prtcl)
{
  char name[20];
  char* Pos;
  if(nin<2 || p>2) return 0;
  sprintf(name,"%s%%",prtclbase1[Prtcl].name);
  Pos=strstr(hadrons[p-1].contents,name);
  if(Pos==NULL) { return 0;}
  if(Pos==hadrons[p-1].contents)  { return 1;}
  if(Pos[-1]==' ' || Pos[-1]==','){ return 1;}
}  

static void  addlim(whohow p_list,int j,int k, int anti)
{int i; 

   if (anti && prtclbase[j-1].anti < j) j = prtclbase[j-1].anti; 

   for (i = 0; i < whohowMAX-1; i++) 
   { 
      if (p_list[i].who == j)
      { if(  ( k>0 && p_list[i].how>0 && k<p_list[i].how) 
           ||( k<0 && p_list[i].how<0 && k>p_list[i].how)  
          ) { p_list[i].how=k; return; }
      }
      if(p_list[i].who == 0) 
      {  p_list[i].who = j; 
         p_list[i].how = k; 
         p_list[i+1].who = 0; 
         return;
      } 
   } 
} 


static char ** stritems(char * format, char * s)
{
  int i, space=1, item=0, l=strlen(s);
  char ** res=NULL;

  for(i=0;i<l;i++) 
  if ( !strchr(format,s[i])) {if(space) { space=0;  item++; }}
  else space=1;
 
  res=malloc((item+1)*sizeof(char *));
  res[item]=NULL;
  
  item=0;
  space=1;
  for(i=0;i<l;i++) 
  if( !strchr(format,s[i]))
  {if (space) { space=0; res[item++]=s+i;}}
  else space=1;  

  return res;
}  

static void  prtcllist(int  key)
{
 char         fullname[STRSIZ];
 char         hlp[60];
 char         p1[60], p2[60];
 int         i, j, pnum;
 linelist     ln;
 int  tabMax,tabSz;

 static int    nTot,nFirst;

	tabMax=ycons -7;
	if (key==0)
	{
		scrcolor(FGmain,BGmain);
		for (i = 2; i <= 24; i++)
		{  goto_xy(1,i);
			clr_eol();
		}
		goto_xy(14,3);
		scrcolor(Blue,BGmain);
		print("List of particles (antiparticles)");
		nTot=0;
		nFirst=1;
	}
	else
	{
		if (nTot <= 3 *tabMax )   return;
		switch (key)
		{
		  case KB_DOWN : nFirst+=3;         break;
		  case KB_UP   : nFirst-=3;         break;
		  case KB_PAGED: nFirst +=3*tabMax; break;
		  case KB_PAGEU: nFirst -=3*tabMax; break;
		}
		if (nFirst <1) nFirst=1;
		if (nTot-nFirst+3<3*tabMax )  nFirst=1+3*((nTot+2)/3) -3*tabMax;
		clrbox(1,4,79,5+tabMax);
	}
	goto_xy(3,5); scrcolor(FGmain,BGmain);
	for(i=0,ln=prtcls_tab.strings;  ln; ln=ln->next )
	{  sscanf(ln->line,"%[^|]%*c%[^|]%*c%[^|]%*c%*[^|]%*c%*[^|]%*c%*[^|]%*c%*[^|]%*c%*[^|]%*c%[^|]",
			 fullname,p1,p2,hlp);
		trim(p1);
		locateinbase(p1,&pnum);
		trim(hlp);
		if (prtclbase[pnum-1].top != NULL && strcmp(hlp,"*") != 0)
		{
			i++;
			if (i>=nFirst && (i-nFirst)/3 <tabMax )
			{
				print("%s",p1);
				if (strcmp(p1,p2) == 0) print("     "); else print("(%s)",p2);
				trim(fullname);
				print("- %s",fullname);
				j = i % 3;
				if (j == 0)	goto_xy(3,where_y() + 1);
						else	goto_xy(3 + 26 * j,where_y());
			}
		}
		
	}
	nTot=i;
	tabSz=MIN((nTot+2)/3,tabMax);
	chepbox(1,4,79,5+tabSz);

	if (nFirst >1 ) { goto_xy(72,4); print("PgUp");  }

	if (nFirst+3*tabSz <= nTot    ) { goto_xy(72,5+tabMax); print("PgDn");  }

	scrcolor(FGmain,BGmain);
}




static char * errTxt=NULL;
static char err_Txt[40];

static int input(int y0, char*hlp, char*directive, char*text, int cur, int lim)
{
   if(errTxt)
   {
      goto_xy(1,y0+1); 
      scrcolor(FGmain,BGmain);
      clr_eol();
      print("%s","Error: ");
      scrcolor(Red,BGmain);
      print("%s",errTxt);
      scrcolor(FGmain,BGmain);
      if(blind) { printf("%s\n",errTxt); sortie(110);}
      be_be();
      errTxt=NULL;
   }
   for(;;) 
   { int rc;
     goto_xy(1,y0);  
     scrcolor(FGmain,BGmain);
     print("%s",directive);
     rc=str_redact(text,cur, lim);
     switch(rc)
     { case  KB_UP:
       case  KB_DOWN: 
       case  KB_PAGED:
       case  KB_PAGEU:  prtcllist(rc); continue;
       case  KB_F1:     show_help(hlp); continue; 
     }
     return rc;                            
   }
}

static int enter_h(int * y,char* name,int num,int scat)
{ int      i,m,j=0;
  int      redres;
  char     hadrch[STRSIZ];  

  char ** items;
  char * errpos=NULL;
  locateinbase(name,&j);
  if(j)
  { 
    if(pseudop(j)) return -1;
    strcpy(hadrons[num].name,name);
    strcpy(hadrons[num].contents,name); 
    hadrons[num].parton[0] = j; 
    hadrons[num].pow = 1;
    if(!scat && num==0 && strcmp(prtclbase1[j].massidnt,"0")==0) 
    { errTxt="Decay of massless particle.";
      return -1;
    }   
    if(name[strlen(name)-1]=='%')
    { 
      if( !scat || nout||strcmp(prtclbase1[j].massidnt,"0") || 
          (prtclbase1[j].spin!=1 && prtclbase1[j].spin!=2) || 
          strchr("LR",prtclbase1[j].hlp)
        ) 
       { errTxt="This particle can not be polarized";   
         return -1;
       } else  hadrons[num].polarized[0]=1; 
    } else hadrons[num].polarized[0]=0;
    return 0;
  }

  for(i=0;i<num;i++) 
  if(!strcmp(hadrons[i].name,name))
  {  strcpy(hadrons[num].name,name);
     strcpy(hadrons[num].contents,hadrons[i].contents);
     hadrons[num].pow=hadrons[i].pow;
     for(j=0;j<hadrons[num].pow;j++)
     { hadrons[num].parton[j]=hadrons[i].parton[j];  
       if(nout) hadrons[num].polarized[j]=0;
       else hadrons[num].polarized[j]=hadrons[i].polarized[j];
     }  
     return 0;
  }

  hadrch[0]=0; 

  for(i=num;i<MAXINOUT;i++) if(!strcmp(hadrons[i].name,name))
  {  
     strcpy(hadrch,hadrons[i].contents);
     break;
  } 

  if(*y>=maxRow()-1) { goto_xy(1,*y); clr_eol();} else (*y)++;  

  do
  {
     m=errpos? errpos-hadrch+1: 0;
  
     do 
     {  char direction[100];
        sprintf(direction,"composite '%s'  consists of: ",name);
        redres=input(*y, "s_ent_2", direction,  hadrch, m , STRSIZ-1);
        if(redres==KB_ESC) return 1;               
     }  while (redres!=KB_ENTER && redres!=KB_ESC);

     if (redres == KB_ESC || strcmp(hadrch,"") == 0)  return 1;
  
      
     items=stritems(" ,",hadrch);
     for(m=0,hadrons[num].pow=0; items[m]; m++) 
     { char  name[100];
       if(hadrons[num].pow>=100) {errTxt="too many partons";break;}   
       sscanf(items[m],"%[^ ,]",name);
       locateinbase(name,&j);  
       if (j==0 || pseudop(j)) 
       { errTxt= "This particle is absent in the model"; break;}

       if(!scat && num==0 && strcmp(prtclbase1[j].massidnt,"0")==0)
       { errTxt="Decay of massless particle.";
         break;
       }
       if(name[strlen(name)-1]=='%')
       { 
         if( !scat || nout||strcmp(prtclbase1[j].massidnt,"0") || 
              (prtclbase1[j].spin!=1 && prtclbase1[j].spin!=2) ||
              strchr("LR",prtclbase1[j].hlp)
           ) 
         { errTxt="This particle can not be polarized";   
              break;
         } else  hadrons[num].polarized[hadrons[num].pow]=1; 
       } else hadrons[num].polarized[hadrons[num].pow]=0; 
       hadrons[num].parton[hadrons[num].pow++]=j;
     }

     errpos=items[m]; 
     if(!errpos)
     {  for(i=0;i<hadrons[num].pow;i++)
        for(j=i+1;j<hadrons[num].pow;j++)
        if(hadrons[num].parton[i]==hadrons[num].parton[j])
        {  errpos=items[j];
           errTxt="duplicate parton";
        }
     }
     strcpy(hadrons[num].name,name);
     strcpy(hadrons[num].contents,hadrch);
     free(items);
  }
  while(errpos);   
  return 0;
}


static int restrict_p(int y0, int anti, char * mess, char * hlp, 
char *  inputstr, whohow  liminsp)
{ 
   int m, j, k;
   int ntot,ntotQ,forQ;

   char ** items;   
   char *errpos=NULL;

   int r; 

do{
    for(r=1;r!=KB_ENTER;) 
    { 
       r=input(y0, hlp,  mess,inputstr, errpos? 1+errpos-inputstr:1,STRSIZ-1);
       if(r==KB_ESC) return 1;
    }
    trim(inputstr);
    nilprtcl(liminsp);
    nilprtcl(LimQ);
    ntot=0;
    ntotQ=0;
    items=stritems(",",inputstr);
    for(m=0 ; items[m]; m++)
    {  char frgm[100];
       char *n;
       forQ=0;
       sscanf(items[m],"%[^,]",frgm);
             
      if(n=strstr(frgm,"!="))       
      { if(!anti) {errTxt="wrong restriction"; break;}
        if(sscanf(n+2,"%d",&k)!=1) {errTxt="wrong number"; break;}
        else  { n[0]=0; k++; forQ=1;}
      }
      else if(n=strstr(frgm,"<"))
      { if(!anti) {errTxt="wrong restriction"; break;}
        if(sscanf(n+1,"%d",&k)!=1) {errTxt="wrong number"; break;}
        else {n[0]=0; k=-k-1;}
      }
      else if(n=strstr(frgm,">"))
      {
        if(sscanf(n+1,"%d",&k)!=1) {errTxt="wrong number"; break;}
        else {n[0]=0;k++;}
      }
      else k=1;
      trim(frgm);

      if(strcmp(frgm,"%Z+W")==0)
      {  
         if(k>0) ZWmax=k-1;
         if(k<0) ZWmin=-k-1;
         continue;
      }
                                                 

      locateinbase(frgm,&j);
      if ((j == 0) || (!anti&&pseudop(j)))
      {  errTxt="wrong limit statement";
         break;
      }
      if(forQ) 
      { ntotQ++; 
        if(ntotQ>=whohowMAX) { errTxt="Too many items"; break;}
        addlim(LimQ,j,k,anti);
      }else  
      { ntot++;
        if(ntot>=whohowMAX) { errTxt="Too many items"; break;}
        addlim(liminsp,j,k,anti);
      }  
    }
    errpos=items[m];
    free(items);
  }while(errpos);
   return 0;
}


int enter(void)
{   
   int  i, y0;
   int scat;
   int  redres;

   char ** items=NULL;
   char * errpos=NULL;
   char * arrpos=NULL;

   int curh=0;

   ZWmax= 1000;
   ZWmin=-1000;
     
   prtcllist(0);
   scrcolor(Red,BGmain);
   y0 = ycons;
   errTxt=NULL;
label_1:
   if(y0<maxRow()) y0++;  
   for(;y0>ycons;y0--){goto_xy(1,y0); clr_eol();} 

   if(arrpos) strncpy(arrpos,"->",2);

   redres=input(y0,"s_ent_1", "Enter  process: ",processch,errpos?errpos-processch+1:1,SSTRLEN-1); 
   switch (redres)
   {
	case KB_F1:   /*  Help  */
	       show_help("s_ent_1");
	       goto label_1;
	case KB_ESC: 
	       clrbox(1,2,maxCol(),24);
               return 1;
   }   /*  Case  */

   
   arrpos=strstr(processch,"->");
   if(arrpos)  strncpy(arrpos,", ",2); 
   else { errTxt="'->' is absent "; goto label_1; }
 
   if(items) free(items);
   items=stritems(" ,",processch);
   if( items[0] && items[1] &&items[1] <arrpos ) scat=1; else scat=0;
    
   for(i=0,nin=0,nout=0,n_x=0,curh=0; items[i]; i++)   
   {  char name[100];

      sscanf(items[i],"%[^, ]", name); 
      if(strlen(name)>7) {errTxt="Too long name"; break;}
      if(strlen(name)==0){errTxt="Empty item"; break;}
      if(items[i]>arrpos && strlen(name) == 3 && name[1] == '*' &&
             (name[2] == 'X' || name[2] == 'x')  &&
             isdigit(name[0])) {n_x += name[0]-'0'; nout +=name[0]-'0';}   
      else
      {
        if(curh==MAXINOUT-1){errTxt="Too many particles"; break;}
        if(items[i]>arrpos) nout++; else nin++; 
        if(nin>2){errTxt="Too many incoming particles"; break;}                  
        if(enter_h(&y0,name,curh,scat)) break;
 
        curh++;
      }
   }
   
   errpos=items[i]; 
   free(items); items=NULL;
   if(errpos && !errTxt) 
   { sprintf(err_Txt,"wrong definition of %d-th particle",i+1);
      errTxt=err_Txt;
   }   
   if(errTxt) goto label_1;    
   if(nout+nin<3)
   {  errTxt="The total number of particles should  exceed 2";
      goto label_1;
   }
   if(nin<1)
   {  errTxt="Incoming  particle(s) is not defined";
      goto label_1; 
   }

   if(errpos||nin+nout>MAXINOUT) 
   {if(!errTxt) errTxt=errtxt; goto label_1; }   
    
   y0++;
   if(y0>=maxRow()-1) {y0=maxRow()-1; goto_xy(1,y0); clr_eol();}

   if(restrict_p(y0,1, "Exclude diagrams with ","s_ent_4",limpch,liminsp)
     ) goto label_1;

   for( i=0;LimQ[i].who;i++)  addlim(liminsp,LimQ[i].who, LimQ[i].how,0);

   y0++;
   if(y0>=maxRow()-1) {y0=maxRow()-1; goto_xy(1,y0); clr_eol();} 
   if(n_x && restrict_p(y0,0,"Exclude X-particles ","s_ent_5",deloutch,limout)
     ) goto label_1;   
          
   for(i=nin+nout;i<MAXINOUT;i++) hadrons[i].name[0]=0;
  
   strncpy(arrpos,"->",2);

   return 0;
}
