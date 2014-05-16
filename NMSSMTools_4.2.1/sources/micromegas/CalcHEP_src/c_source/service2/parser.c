/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#define ILEN 200
#define MAXFUN 200
#include "parser.h"
#include<stdio.h>
#include <ctype.h>
#include <string.h> 
#include <stdlib.h>

int  rderrcode,   rderrpos;

static char       * source;
static rdelement  rd0;
static int        count;
static operation  act0;
static char       keyChar[]=".^/*+-),|%";
static void **    toDel;
static int        nDel, delFlag;


static void addToDel(void * res)
{ toDel[nDel]=res;
  nDel++;
  toDel=realloc(toDel, (nDel+1)*sizeof(void*));
}
            
static void * act(char * operation,int n, void ** args)
{ void * res=act0(operation,n,args);
  if(!res && !rderrcode) rderrcode=cannotevaluate;
  if(rderrcode) return NULL;
  if(delFlag) addToDel(res);
  return res;
}
  
static void * rd(char * name)
{ void * res=rd0(name);
  if(!res){if(!rderrcode) rderrcode=cannotread;rderrpos=count-strlen(name)+1;}
  if(rderrcode) { return NULL;} 
  if(delFlag) addToDel(res);  
  return res;
}

static void*readSum(void);

static void  skip(void) { while(source[count] == ' ') count++;}


static void * readSmplTerm(void)
{
   int      m1;
   char tmp[ILEN+1];
   int len;
   void * s_res=NULL;
   int n=0;                                                                      
   void * tmp_arg[MAXFUN];                                                       
   char ch;
   
   skip();
   if(source[count] == '(') 
   { 
      count++;
      if((s_res=readSum()))
      {if(source[count] != ')')  rderrcode = braketexpected; else count++;}
      return s_res;
   }
   ch=source[count]; 
   if(!isalnum(ch)&& ch!='_'  && ch!='"') {rderrcode=unexpectedcharacter; return NULL;}

   m1 = count; 
   if(isdigit(source[count]))
   {
      do count++; while (isdigit(source[count]));
      if(source[count]=='.') do count++; while (isdigit(source[count]));
      if(source[count]=='E'||source[count]=='e')
      { count++;
        if(source[count]=='+'||source[count]=='-') count++;
        if(!isdigit(source[count])) 
        { rderrcode=unexpectedcharacter; 
          return NULL;
        }
        else do count++; while (isdigit(source[count]));
      }
   }                 
   else if(source[count]=='"') 
   {count++; while(count-m1<=ILEN &&  (source[count]!='"' || source[count-1]=='\\') ) count++; count++;} 
   else {count++; while(isalnum(source[count])||source[count]=='_') count++;}

   len=count-m1;
   if(len>ILEN) {rderrcode=toolongidentifier; return NULL;} 
   else         sprintf(tmp,"%*.*s",len,len,source+m1); 
   
   if(tmp[0]=='0' && !strchr(tmp,'.') && len>1 )
   { rderrcode=unexpectedcharacter; return NULL;}
   
   if(isdigit(tmp[0])|| tmp[0]=='"' ||source[count] != '(') return rd(tmp);

   do                                                                            
   { count++;  skip(); if(source[count]==')')break;                                                                   
     if(n==MAXFUN) {rderrcode=toomanyagruments; break;}                          
     if(!(tmp_arg[n]=readSum()))  break;                                                      
     n++;                                                                        
   } while (source[count]==',');                                                 
   if(!rderrcode) 
   {if(source[count] != ')') rderrcode = braketexpected; else count++;}  

   if(!rderrcode && !(s_res=act(tmp,n,tmp_arg)))  rderrpos = m1+1;
   return s_res;
}   /*  ReadSmplTerm  */ 


static void * readSmplTermPlus(void)
{ void * s_res; 
  s_res=readSmplTerm();
  if(s_res) 
  { skip();
    if(!strchr(keyChar,source[count])) 
    {rderrcode = operationexpected; return NULL;}
  } 
  return s_res;                              
}

static void * readMonom(int level)
{  void * m_res;
   char com[2]="?";
   skip();
   m_res= level ? readMonom(level-1):readSmplTermPlus(); 
   
   if(m_res)
   {
      while(!strchr(keyChar+1+level,source[count]))
      {  int  sgn_pos;
         void*  tmp_arg[2]; 

         sgn_pos = ++count;  
         tmp_arg[1]= level ? readMonom(level-1):readSmplTermPlus();         
         if(!tmp_arg[1])  return NULL;
         tmp_arg[0]=m_res;
         com[0]=keyChar[level];
         m_res=act(com,2,tmp_arg);
         if(!m_res) rderrpos = sgn_pos; 
      }
   } 
   return m_res; 
}

static void *  readSum()
{  char   sgn;
   void*  tmp_arg[2];
   int    sgn_pos;
   void*  s_res;

   skip();
   sgn = source[count];
   sgn_pos = count+1;
   if (sgn == '-' || sgn == '+') count++; else sgn = '+';
   if(!(s_res=readMonom(3)))  return NULL;
   
   if(sgn=='-' && !(s_res=act("-",1,&s_res))){rderrpos=sgn_pos;return NULL;}
   tmp_arg[0]=s_res;

   while (!strchr(",|%)",source[count]) )
   {
      sgn = source[count];
      sgn_pos = ++count;
      if(!(s_res=readMonom(3))) return NULL;
      if (sgn=='-' && !(s_res=act("-",1,&s_res))){rderrpos=sgn_pos+1; return NULL;}
      tmp_arg[1]=s_res;   
      s_res = act("+",2,tmp_arg);
      if (!s_res) { rderrpos = sgn_pos; return NULL; }
      tmp_arg[0]=s_res;
   }
   return s_res;
}  

void* readExpression(char* source1,rdelement rd1,operation act1,clear del1)
{  void*  rslt;
   source=source1; 
   act0=act1,
   rd0=rd1;
   delFlag=(del1!=NULL);

   if(delFlag) {nDel=0; toDel=malloc(sizeof(void*));}   
   count=0;
   rderrpos=0;
   rderrcode=0;
   rslt= readSum();
   if(rderrcode==0 && ( source[count]==',' || source[count]==')') ) 
   {rslt=NULL; rderrcode=unexpectedcharacter;}
    
   if(!rslt && !rderrpos) rderrpos=count+1;

   if(delFlag)
   {
   
     while(--nDel>=0) if(toDel[nDel]!=rslt) del1(toDel[nDel]);
      free(toDel); 
   }
   return  rslt;
}

char * errmesstxt(int errcode)
{
   switch (rderrcode)
   {
    case braketexpected:      return "Bracket expected";          
    case unexpectedcharacter: return "Unexpected character";      
    case operationexpected:   return "Operation expected";        
    case toolongidentifier:   return "Too long identifier/number";
    case toomanyagruments:    return "Too many arguments";        
    case cannotread:          return "Undefined identifier";      
    case cannotevaluate:      return "Can not evaluate function"; 
    case unknownidentifier:   return "Unknown identifier";        
    case unexpectedoperation: return "Unexpected operator/function"; 
    case unknownfunction:     return "Unknown function";          
    case wrongnumberofagr:    return "Wrong number of arguments"; 
    case typemismatch:        return "Type mismatch";             
    case naninoperation:      return "NAN in operation";          
    case toolargenumber:      return "Too large or float number"; 
    case indexuncompatibility:return "Index uncompatibility";     
    case rangecheckerror:     return "Range check error";         
         default:             return "Unknown error code";        
   }

}
