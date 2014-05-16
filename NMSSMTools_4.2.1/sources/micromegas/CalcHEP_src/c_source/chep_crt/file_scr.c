/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "crt.h"
#include "crt_util.h"
#include "help.h"

#include "file_scr.h"


char  errorText[256] = "";

static char buff[90]=
"                                                                                ";

static int fgets1(char *  c,int  n, FILE * fp)
{  int l;
   if(fgets(c,n,fp)==NULL) return 0;	
   l= strlen(c) - 1;
   if (c[l] == '\n')  { c[l] = '\0'; return 1;}  
   else 
   { fscanf(fp,"%*[^\n]");
     fscanf(fp,"%*c");
     return -1;
   }  	
}



static int readnewline(FILE* fileptr , int  lwidth , linelist * lastln)
{
  linelist newline;
  int i,l;
  if(*lastln && (*lastln)->next) return 0;
  newline=(linelist)m_alloc(sizeof(linerec)-STRSIZ+lwidth+2);
  if (newline == NULL)
  { if (*lastln) strcpy((*lastln)->line," -- FILE TOO LARGE --");
     return 0;
  }
  else
  {
     if (!fgets(newline->line,lwidth+2,fileptr) )
     {  free(newline);
	return 0;
     }
     newline->next=NULL;
     newline->pred=*lastln;
     if(newline->line[strlen(newline->line)-1]=='\n')
	           newline->line[strlen(newline->line)-1]=0;  
     if (strlen(newline->line)==lwidth+1)
     {  ungetc(newline->line[lwidth],fileptr);
	newline->line[lwidth]=0;
     }
     
     l=strlen(newline->line);
     for (i=0;i<l;i++) 
     if (newline->line[i]<32)   
     {
        sprintf(newline->line,"%.*s",lwidth -1,"Unprinted symbols!");
        i=l;
        fseek(fileptr,0,SEEK_END);
     }
     if(*lastln) (*lastln)->next=newline; else *lastln=newline;
  }
  return 1;
}

static void dellinelist(linelist  ln)
{
  linelist lndel;
  while (ln !=NULL)
  {
	 lndel=ln;
	 ln=ln->next;
	 free(lndel);
  }

}

void cleartab(table * tab)
{
  dellinelist(tab->strings);
  tab->strings=NULL;
}


void showtext(int x1, int y1, int x2, int y2, char* headline,
 FILE * fileptr )
{
   int   key;
   void *  pscr;
   int   i,l;
   linelist  currentline=NULL,lastline;
   int width=x2-x1-1;
   int fc=fColor,bc=bColor;	

   if (fileptr==NULL || !readnewline(fileptr,width,&currentline))
   {
      messanykey(10,10," File is empty");
      return;
   }
	       
   if(y2<=y1+1) 
      for (y2=y1+2;y2<maxRow()&&readnewline(fileptr,width,&currentline);y2++)
         currentline=currentline->next;	      
            
   get_text(x1,y1,x2,y2,&pscr);

   scrcolor(Black,White);
   chepbox(x1,y1,x2,y2);
   scrcolor(Red,White);	
   goto_xy(x1+1,y1); print("*");
   l=strlen(headline);
   if ( l < (x2-x1) )
   {
      goto_xy(x1+(x2-x1-l)/2,y1);
      scrcolor(White,Black);
      print(headline);
      scrcolor(Black,White);
   }

   key = KB_PAGEU;
   do
   {
      switch (key)
      {
	case KB_PAGEU:
	   for (i=y1+1;i<y2;i++)
	     if (currentline->pred !=NULL)  currentline=currentline->pred;
           break;

	case KB_UP:
	   for(i=1;i<(y2-y1)/2;i++)
              if (currentline->pred !=NULL)  currentline=currentline->pred;
	   break;

	case KB_PAGED: 
	   for (i=y1+1;i<y2;i++)
	   {
	      readnewline(fileptr,width,&currentline);
	      if (currentline->next !=NULL)  currentline=currentline->next;
	   }
	   break;

	case KB_DOWN:
	   for(i=1;i<(y2-y1)/2;i++)
           {
	      readnewline(fileptr ,width,&currentline);
	      if (currentline->next !=NULL)  currentline=currentline->next;
	   }
			break;

      }  /*  Case  */

/*    display */

      lastline=currentline;
      goto_xy(x2-4,y1); 
      if( lastline->pred == NULL) for (i=0;i<4;i++) print("%c",boxFrame[1]);
              else                print("PgUp");
      for (i =y1+1; i <y2; i++)
      {
	   goto_xy(x1+1,i);
	   readnewline(fileptr,width,&lastline);
	   if ( lastline !=NULL ) { print("%s",lastline->line);
	   lastline=lastline->next;}
	   print(buff+where_x()+80-x2);
      }
      goto_xy(x2-4,y2);
      if (lastline == NULL) for (i=0;i<4;i++) print("%c",boxFrame[5]);
              else                            print("PgDn");   
      key = inkey();
/* mouse filter */
       if ( (key == KB_MOUSE)&&(mouse_info.but1==2)
         &&(mouse_info.col >= x1)&&(mouse_info.col <= x2)  )
       { 
         if (x2 - mouse_info.col <4 ) 
         { if(mouse_info.row == y1) key=KB_PAGEU; else            
           if(mouse_info.row == y2) key=KB_PAGED;
         }   
         if ( (mouse_info.row == y1)&&(mouse_info.col - x1 <3)) key=KB_ESC;
         if ((mouse_info.row > y1)&&(mouse_info.row < y2))                        
         { if (mouse_info.row < (y1+y2)/2)  key=KB_UP; else key=KB_DOWN;}
       }                   
/* end of filter */                                        
   }  while ((key != KB_ESC)&&(key !=KB_BACKSP));
     
   scrcolor(fc,bc);
   put_text(&pscr);
   while (currentline->pred != NULL )  currentline=currentline->pred;
   dellinelist(currentline);

}


int  readtable0(table*  tab,FILE * txt)
{  
  linelist lastln,newline;                                                     
  int len,i;
  linerec zeroline;                                                            
                                                                          
  if(txt==NULL) return 1;                                                      
  cleartab(tab);                                                               
  if(!fgets1(tab->mdlName,STRSIZ,txt)) return 1;                                             
  if(!fgets1(tab->headln,70,txt))return 1;
  if(fgets1(tab->format,STRSIZ,txt)<=0)return -2;                                              
  len=strlen(tab->format);
  while( len && tab->format[len] !='|')  len--;  
  if(len==0) return 1;                                    
  tab->format[len+1]=0;                                                        
                                                                          
   zeroline.next=NULL;                                                         
   lastln= &zeroline;                                                          
                                                                          
   newline=(linelist)m_alloc(sizeof(linerec)-STRSIZ+len+2);                    
   while (fgets1(newline->line ,len+1,txt) && newline->line[0]!='=')           
   {                                                                           
       newline->line[len+1]=0;                                                 
       for (i=strlen(newline->line);i<=len;i++) newline->line[i]=' ';          
       lastln->next=newline;                                                   
       newline->next=NULL;                                                     
       newline->pred=lastln;                                                   
       lastln=newline;                                                         
       newline=(linelist)m_alloc(sizeof(linerec)-STRSIZ+len+2);                
   }                                                                           
   free( newline);                                                             
   tab->strings=zeroline.next;                                                 
   if (tab->strings != NULL)tab->strings->pred=NULL;                           
   return 0;
}

int  readtable(table*  tab, char* fname)
{
  FILE * txt;
  int err; 
  txt=fopen(fname,"r");
  if(txt==NULL) return -1;
  err=readtable0(tab,txt);
  fclose(txt);   
  return err;
}



void  writetable0(table*  tab, FILE * fi)
{ int i;

  linelist ll=tab->strings;
  if(!(fi)) return;
  f_printf(fi,"%s\n",tab->mdlName);
  f_printf(fi,"%s\n",tab->headln);
  f_printf(fi,"%s\n",tab->format);
  while (ll) 
  { i=strlen(ll->line)-1; 
    while(i>=0 && ll->line[i]==' ')i--;      
    f_printf(fi,"%*.*s\n",i+1,i+1,ll->line); 
    ll=ll->next;
  }
  for(i=0;i<strlen(tab->format);i++) f_printf(fi,"=");
  f_printf(fi,"\n");
}

void  writetable(table*  tab,char*  fname)
{
  FILE * fi;
  fi=fopen(fname,"w");
  writetable0(tab, fi);
  fclose(fi);
}
