/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include<string.h>
#include "syst.h"
#include "syst2.h"
#include "physics.h"
#include "procvar.h"
#include "reader_c.h"

FILE * ext_h;

int checkNaN;

static void * bact5(char ch,void * mm1,void * mm2)
{  int dbl;
   char   *m1, *m2, *ans;
   char   r_n,p_m;

   m1 = (char *) mm1;
   m2 = (char *) mm2;
   if (ch=='+' || ch=='-') p_m='P'; else p_m='M';
   
   if( (ch=='+'||ch=='-') && m1[1]=='N'&& m2[1]=='N') r_n='N'; else r_n='R'; 
   if( ch=='/' && m1[1]=='N' && m2[1]=='N') dbl=1; else dbl=0;
   
   if (m1[0] == 'M' || ch =='+'|| ch =='-')    lShift(m1,3); else { lShift(m1,2);m1[0]='(';strcat(m1,")"); }
   if ((m2[0] == 'M' || ch =='+')&& ch != '/') lShift(m2,3); else { lShift(m2,2);m2[0]='(';strcat(m2,")"); }

   ans= m_alloc(strlen(m1)+strlen(m2)+30);

   switch (ch)
   {
      case '+': 
              if (m1[0] == '-')        sprintf(ans,"%c%c|%s%s",p_m,r_n,m2,m1);
	      else  if (m2[0] == '-')  sprintf(ans,"%c%c|%s%s",p_m,r_n,m1,m2);
              else                     sprintf(ans,"%c%c|%s+%s",p_m,r_n,m1,m2);
      break;

      case '-': 
	      if (m2[0] == '-')  sprintf(ans,"%c%c|%s+%s",p_m,r_n,m1,m2+1);
              else               sprintf(ans,"%c%c|%s-%s",p_m,r_n,m1,m2);
      break;
      

      case '*':
              if (m2[0] != '-')        sprintf(ans,"%c%c|%s*%s",p_m,r_n,m1,m2);
              else if (m1[0] != '-')   sprintf(ans,"%c%c|%s*%s",p_m,r_n,m2,m1);
              else                     sprintf(ans,"%c%c|%s*%s",p_m,r_n,m1+1,m2+1);
      break;

      case '/': if(dbl) sprintf(ans,"%c%c|%s/(double)(%s)",p_m,r_n,m1,m2); 
		else if(m2[0] != '-')   sprintf(ans,"%c%c|%s/%s",p_m,r_n,m1,m2);
		else
                {  if (m1[0] == '-') sprintf(ans,"%c%c|%s/%s",p_m,r_n,m1+1,m2+1);
                   else    sprintf(ans,"%c%c|-%s/%s",p_m,r_n,m1,m2+1);
                }
                checkNaN=1;
      break;

      case '^':
                 sprintf(ans,"%c%c|pow(%s,%s)",p_m,r_n,m1,m2);
                 checkNaN=1;
      break;
      default:  checkNaN=1;                 
   }   /* Case */
	return (void *) ans;
}


static void * uact5(char* ch,void * mm)
{ char  *m, *ans;
  m = (char *) mm;
  ans=m_alloc(strlen(m)+30);

  if (strcmp(ch,"-") == 0)
  {
     if(m[0] == 'M')
     {  if (m[3] == '-')  sprintf(ans,"M%c|%s",m[1],m+4);
	else              sprintf(ans,"M%c|-%s",m[1],m+3);
     } else sprintf(ans,"M%c|-(%s)",m[1],m+3);
  } 
  return (void *) ans;
}

void * act_c(char * name,int n, void ** args)
{ int l,i;
  char * ans;
  char tp='R';
  checkNaN=1;
  if(!isalpha(name[0]))
  { 
    if(n==1) return uact5(name,args[0]);
    if(n==2) return bact5(name[0],args[0],args[1]);
  } else if(strcmp(name,"one")==0)
  { ans=m_alloc(6); sprintf(ans,"MN|1"); return ans; }
  
  if(ext_h) 
  { 
    char _name_[100];
    
    sprintf(_name_," %s ",name);
    
    checkNaN=1;
    if(strstr(EXTFunc,_name_)==NULL  )
    { fprintf(ext_h, " extern double %s(",name);
      if(n==0) fprintf(ext_h,"void);\n");
      else for(i=1;i<=n;i++)
      { if(((char*)args[i-1])[1]=='S') fprintf(ext_h,"char*" );
          else fprintf(ext_h,"double");
          if(i==n)fprintf(ext_h,");\n"); else fprintf(ext_h,",");
      }
    }                                                                                     
  } 
  
  if(strcmp(name,"aWidth")==0) name="aWidth_ext";
  l=n+10+strlen(name);
  for(i=0;i<n;i++) l+=strlen((char*)args[i]); 
  ans=m_alloc(l);
  if(!strcmp(name,"if") && n==3)
  {
    if( ((char*)args[0])[1]=='N' && ((char*)args[1])[1]=='N') tp='N';
     
    sprintf(ans,"M%c|(%s>0 ? %s : %s)",tp, (char*)args[0]+3, (char*)args[1]+3,
               (char*)args[2]+3);
  }             
  else     
  { sprintf(ans,"M%c|%s(",tp,name);
    for(i=0;i<n;i++) 
    { 
      strcat(ans,(char*)args[i]+3);
      strcat(ans,",");
    }
    if(n) ans[strlen(ans)-1]=')'; else strcat(ans,")");
  }
  return ans;
}

void *  rd_c(char* s)
{  char      *p;
   int        l;
   
   if(s[0]=='"'){ p = m_alloc(strlen(s)+10); sprintf(p,"MS|%s",s);} else
   {
     p = m_alloc(40);
     if (isdigit(s[0])) 
     {   sprintf(p,"MN|%s",s); 
         if( strchr(s,'.')||strchr(s,'E')||strchr(s,'e')) p[1]='R';   
     }              
     else
     {
       for(l=1;l<=nmodelvar;l++)
       {
         if(!strcmp(s,modelvars[l].varname))
         {
            sprintf(p,"MR|%s",vararr[l].alias);
            return (void *) p;
         }
       }
     }  
   }
   return (void *) p;
}
