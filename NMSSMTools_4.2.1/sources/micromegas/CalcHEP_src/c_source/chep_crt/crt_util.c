/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "crt_util.h"
#include "help.h"

int      yesnokey(void)
{
 int         x0, y0;
 char         ch;
 int      key;
 int      res;

   y0 = where_y();
   x0 = where_x();
   print(" Y/N ?");

label_1: be_be();
   key = inkey();
   if (key < 0) goto label_1;
   ch = (char)key;
   switch (ch)
   {
      case 'Y':   res = 1;  break;
      case 'y':   res = 1;  break;
      case 'N':   res = 0;  break;
      case 'n':   res = 0;  break;
      default:    goto label_1;
   }
   goto_xy(x0,y0);
   print("      "); 
   return res;
} 

static int maxline(int* lines,char* txt)
{ int      tmp=0, i=0, k; 
  
    *lines = 0; 
    for (k = 1; k <= strlen(txt); k++) 
      if (txt[k-1] == '\n') 
      { ++(*lines); 
        if ((k - i) > tmp) tmp = k - i;
        i = k; 
      }
    if ((k - i) > tmp) tmp = k - i; 
    return tmp; 
} 


static int  message(int x1,int y1,char* txt1)
{  char     txt[STRSIZ];
   int      mess, marg, i=0, x2, y2;
   void *      dump;
   int         xold, yold;
   char * c;   

     strcpy(txt,txt1);
     c=txt;
     xold = where_x();
     yold = where_y();  
     x2 = x1 + maxline(&y2,txt) /*+4*/;
     y2 = y1 + y2 + 1;
     get_text(x1,y1,x2,y2,&dump); 
     chepbox(x1,y1,x2,y2);     
     be_be();
     c=strtok(txt,"\n");
     
     while(c)
     {
        ++(i);
        marg = (x2 - x1 - strlen(c)) / 2;
        if (y1+i == y2) 
        {  goto_xy(x1+1+marg ,y1 + i);
           print("%s",c);
        }else
        {  char buff[STRSIZ];
           memset(buff,' ',x2-x1 -1);
           memcpy(buff+marg,c,strlen(c));
           buff[x2-x1-1]=0;
           goto_xy(x1+1,y1+i);
           print(buff);
        }
        c=strtok(NULL,"\n");   
     }
ret_:
     if (blind) mess='Y' ; else  do{mess = inkey();} while (mess == KB_SIZE);
/* mouse filter */
     if (mess == KB_MOUSE)
     {
        if (  (mouse_info.but1 !=2) ||
              (mouse_info.row < y1) || (mouse_info.row > y2) ||
              (mouse_info.col < x1) || (mouse_info.col > x2)   ) goto ret_;
        if (mouse_info.row == y2 )
        {  if (mouse_info.col < (x1+x2)/2 /* -1*/ ) mess='Y';
           if (mouse_info.col > (x1+x2)/2 /*+1 */) mess='N';
        }       
     }   
/* end of filter */         
     put_text(&dump); 
     goto_xy(xold,yold);     
     return mess;
}   /* message */

int     mess_y_n(int x1,int y1,char* txtstr)
{
  int      key;
  int      fcolor_tmp=fColor;
  int      bcolor_tmp=bColor;
  char  newtext[STRSIZ];

  scrcolor(White,Red);
  sprintf(newtext,"%s\n( Y / N ?) ",txtstr);       

  for(;;)
  {
    key = message(x1,y1,newtext);
/*    if(blind) printf(" %s   YES\n",newtext); */
    switch (key)
    {
      case 'N':
      case 'n': scrcolor(fcolor_tmp,bcolor_tmp);return 0;
      case 'Y':
      case 'y': scrcolor(fcolor_tmp,bcolor_tmp);return 1;
    }
  }
}

void  messanykey(int x1,int y1,char* txtstr)
{  
 int  fcolor_tmp=fColor;
 int  bcolor_tmp=bColor;
 char  newtext[STRSIZ];     
   scrcolor(White,Blue);
   sprintf(newtext,"%s\nPress any key ",txtstr);       
   message(x1,y1,newtext);
   scrcolor(fcolor_tmp,bcolor_tmp);
}


int informline(long curent, long total)
{  int  xx;
   int res=0;
   static int Y;
   static int X;
   int xm=where_x(), ym=where_y(), fc=fColor,bc=bColor;
  
   if (curent == 0)
   {
      Y=maxRow();
      X=15;
/*      goto_xy(15,Y); scrcolor(Black,Red);
      memset(b,' ',50); b[50] = '\0';
*/
     goto_xy(15,Y); scrcolor(White,Red);
     print(" Calculation in progress. Calculation in progress.");                                                                     
/*            12345678901234567890123456789012345678901234567890 */
/*      print("%s",b); */
      goto_xy(15,Y);
      scrcolor(fc,bc);
      escpressed();
      goto exi;;
   }

   if(X>65 || X<15) X=15;
   
   xx = 15 + (50 * curent) / total;
   if (xx > 65) xx = 65;
   scrcolor(Black,Red);
   goto_xy(X,Y);
   while (where_x() < xx) print("%c",'X');
   if(xx !=X) { X=xx; res=escpressed();}
   if (curent >= total)
   {
      scrcolor(White,Black);
      goto_xy(1,Y); clr_eol();
      goto exi;
   }
   
   exi:
   
   scrcolor(fc,bc); goto_xy(xm,ym);   
   return res;
   
   
}


void chepbox(int x1,int y1,int x2,int y2)
{int  i,y;
 char  b[STRSIZ];

   i = x2 - x1;

   b[0]   = boxFrame[0];             /* ACS_ULCORNER */
   memset(b+1,boxFrame[1] ,i-1);     /* ACS_HLINE */
   b[i]   =  boxFrame[2];            /* ACS_URCORNER */
   b[i+1] = '\0';
   goto_xy(x1,y1);print(b);

   b[0]   = boxFrame[6];
   memset(b+1,boxFrame[5] ,i-1);
   b[i]   =  boxFrame[4];
   b[i+1] = '\0';
   goto_xy(x1,y2);  print(b);

   b[0]   = boxFrame[7];
   memset(b+1,' ',i-1);
   b[i]   =  boxFrame[3];
   b[i+1] = '\0';

   for (y = y1 + 1; y < y2; y++) 
   {  goto_xy(x1,y); print("%c",boxFrame[7]);
      goto_xy(x2,y); print("%c",boxFrame[3]); 
   }
      

}

void  menu0(int col,int row,char* label, char* menstr ,
	  void (**funcKey)(int) ,char** funcKeyMess, void ** hscr, int* kk)
{  int    i, j, k, col1, npage,lastrow;
   long	  lastpage;
   int    ink;
   int    ncol;
   void * pscr;
   int  fkPos[11];
   int  height;
   char fmt[20];
   int  lastLine;

/* colors */
   int label_fg  =Yellow;
   int label_bg  =Blue;
   int help_fg1  =White;
   int help_fg2  =Black;
   int help_bg   =DarkGray;
   int box_fg    =White;
   int box_bg    =DarkGray;
   int star_fg   =Red;
   int page_fg   =Black;
   int actFunc_fg=Black;
   int actFunc_bg=White;
   int pasFunc_fg=White;
   int pasFunc_bg=DarkGray;

/* save colors */
   int      fcolor_tmp=fColor;
   int      bcolor_tmp=bColor;
   void *hscr_=NULL;
   
   if(hscr==NULL) hscr=&hscr_;
   lastLine=maxRow();
   if (funcKey == NULL) for (i=0;i<11;i++) fkPos[i]=0; else
   {  int xx;
      scrcolor(FGmain,BGmain);
      goto_xy(1,lastLine); clr_eol();
      xx=0;
      for (j=0;j<10;j++) { if(funcKey[j] && funcKeyMess[j]) 
                                          xx=xx+4+strlen(funcKeyMess[j]);}
      xx= (80 - xx )/2 ;
      goto_xy(xx,lastLine);
      for (i=0;i<10;i++)
      { fkPos[i]=where_x();
        if (funcKey[i] && funcKeyMess[i])
        { scrcolor(help_fg1,help_bg); print(" F%i-",i+1);
          scrcolor(help_fg2,help_bg); print(funcKeyMess[i]);
        }
      }
      fkPos[10]=where_x();
   }

   clearTypeAhead();

   if (*kk < 0) *kk = -(*kk);
      ncol=menstr[0];
      sprintf(fmt,"%%%d.%ds",ncol,ncol);
      height=strlen(menstr)/ncol;
      if(height==0) { *kk=0; return; }
      if (row+height+1 >lastLine-2) height=lastLine-3-row;
      lastpage = 1+    (strlen(menstr)/ncol -1)/height ;
   if(label[0] ==0 || row == 1) 
   { if (*hscr == NULL)  get_text(col,row,col+ncol+1,row+2,hscr);} 
   else
   {  char label_[STRSIZ];
      int shft,sz;
      if (*hscr == NULL) get_text(col,row-1,col+ncol+1,row+2,hscr); 
      for(i=0;i<ncol+2;i++) label_[i]=' ';
      label_[ncol+2]=0;  
      sz=strlen(label);
      if(sz >ncol+2) {shft=0;sz=ncol+2;} else shft=(ncol+2 -sz)/2;  
      memcpy(label_+shft,label,sz); 
      scrcolor(label_fg,label_bg);
      goto_xy(col,row-1);
      print(label_);
   }

   get_text(col,row + 3,col + ncol + 1,row + height + 1,&pscr);

   if (*kk <= 0  || *kk > lastpage * height   )
   {  npage = 1;
      k = 1;
   }
   else {
      k = ((*kk) - 1) % height + 1;
      npage = ((*kk) - 1) / height + 1;
   }
   col1 = col + 1;

label_1:
      scrcolor(box_fg,box_bg);
      chepbox(col,row,col + ncol + 1,row + height + 1);
      scrcolor(star_fg,box_bg);
      goto_xy(col+1,row); print("<");
/*      goto_xy(col+1,row + height+1);  print("?"); */
      scrcolor(page_fg,box_bg);
      if (npage > 1)
      {
         goto_xy(col + ncol - 2,row);
         print("PgUp");
      }
      if (npage < lastpage)
      {
         goto_xy(col + ncol - 2,row + height + 1);
         print("PgDn");
      }

      if(npage<lastpage) lastrow=height;
		  else   lastrow = (strlen(menstr)/ncol)%height;

   lastrow=0;
   scrcolor(pasFunc_fg,pasFunc_bg);
   for (j = 1; j <= height; j++)
   {  int shift;
      goto_xy(col + 1,row + j);
      shift=1+(j-1 + (npage-1)*height)*ncol;
      if(shift<strlen(menstr)) {print(fmt,menstr+shift );lastrow++;}
		 else           print(fmt," ");

   }

   scrcolor(actFunc_fg,actFunc_bg);
   if (k > lastrow) k = lastrow;
   goto_xy(col + 1,row + k);
   if (lastrow) print(fmt,menstr+1+(k-1+(npage-1)*height)*ncol);
   while (1)
   {  int jump=1,mousePos;

      scrcolor(pasFunc_fg,pasFunc_bg); 

      ink = inkey();
/* mouse filter */
      if ((ink==KB_MOUSE)&&(mouse_info.but1 == 2))
      {
         if (mouse_info.row == lastLine )
         for(i=0; i<10;i++)
         if ((mouse_info.col > fkPos[i]) && (mouse_info.col < fkPos[i+1]))
         {  if (i==9)ink='0'; else ink='1'+i;}

         if ( (mouse_info.col >= col ) && (mouse_info.col <=col+ncol+1) )
         {  mousePos = mouse_info.row - row;

            if (col+ncol+1-mouse_info.col <4)
            {
               if (mousePos==0)        ink=KB_PAGEU;
               if (mousePos==height+1) ink=KB_PAGED;
            }

            if ((mousePos == 0 ) && ( mouse_info.col-col <4)) ink=KB_ESC;

            if ((mousePos < 0)&&(mousePos >= height))
            {
               if (mousePos > k)  {ink=KB_DOWN; jump=mousePos - k;}
               if (mousePos < k)  {ink=KB_UP;   jump=k - mousePos;}
               if (mousePos==k )   ink=KB_ENTER;
            }
         }
      }
/* end of filter */
      if (lastrow == 0) goto label_3;
label_4:
      switch (ink)
      {
        case KB_MOUSE:
        if (mouse_info.but1 != 2) break;
        if (mouse_info.row == lastLine )
        for(i=0; i<10;i++)
          if ((mouse_info.col > fkPos[i]) && (mouse_info.col < fkPos[i+1]))
          {  if (i==9)ink='0'; else ink='1'+i;
             goto label_4;
          }

        if ( (mouse_info.col < col ) || (mouse_info.col >col+ncol+1) ||
             (mouse_info.row < row ) || (mouse_info.row >row+height+1) ) break;

           mousePos = mouse_info.row - row;
           if ((mousePos == 0 ) && ( mouse_info.col-col <4)) ink=KB_ESC;
           if ((mousePos != 0)&&(mousePos != height+1))
           {
              if  (mousePos > k)  { ink=KB_DOWN; jump=mousePos - k;}
              if (mousePos < k      ) { ink=KB_UP;   jump=k - mousePos;}
           }
           if (col+ncol+1-mouse_info.col <4)
           {
              if (mousePos==0)        ink=KB_PAGEU;
              if (mousePos==height+1) ink=KB_PAGED;
           }
           if (mousePos==k       ) ink=KB_ENTER;
           if (ink!=KB_MOUSE) goto label_4;

          break;

		  case  KB_DOWN: 
           if(k==lastrow)
           { 
              if(npage < lastpage)
              {  k=1;
                 npage++;
                 goto label_1;
              } else { be_be(); break; }           
           }else {ink= KB_RIGHT; goto label_4;}

		  case  KB_UP: 
           if(k==1)
           {
              if (npage > 1)
	      {
                 k=height;
		 npage--;
		 goto label_1;
              }
              else{ be_be(); break; }           
           }
           else {ink= KB_LEFT; goto label_4;}

		  case KB_LEFT: 
            goto_xy(col1,row + k);
	    print(fmt,menstr+1+(k-1+(npage-1)*height)*ncol);
            k = k - jump;
            if (k == 0) k = lastrow;
            scrcolor(actFunc_fg,actFunc_bg);
            goto_xy(col1,row + k);
	    print(fmt,menstr+1+(k-1+(npage-1)*height)*ncol);
         break;

                  case KB_RIGHT:
	    goto_xy(col1,row + k);
	    print(fmt,menstr+1+(k-1+(npage-1)*height)*ncol);
            k = k + jump ;
            if(k > lastrow) k = 1;
            scrcolor(actFunc_fg,actFunc_bg);   
            goto_xy(col1,row + k);
	    print(fmt,menstr+1+(k-1+(npage-1)*height)*ncol);
         break;

                  case KB_ENTER:
            scrcolor(box_fg,box_bg);
            chepbox(col,row,col+ncol+1,row+2);
            put_text(&pscr);
            goto_xy(col1,row + 1);
            scrcolor(actFunc_fg,actFunc_bg);
	    print(fmt,menstr+1+(k-1+(npage-1)*height)*ncol);
            *kk = (npage - 1) * height + k;
            goto_xy(1,lastLine);scrcolor(FGmain,BGmain); clr_eol();
            if(hscr==&hscr_) put_text(hscr); 
            refresh_scr();
            escpressed();
            return;

			case  KB_BACKSP:
			case  KB_ESC:
            goto label_3;

			case KB_PAGEU:
            if (npage > 1)
            {
               npage--;
               goto label_1;
            }
            else
               be_be();
         break;

			case KB_PAGED:
            if (npage < lastpage)
            {
               npage++;
               goto label_1;
            }
            else
               be_be();
	    break;

	case  '1':  case'2':	   case  '3':	case  '4':	case  '5':
	case '6':   case '7':	case  '8':	case  '9':	case  '0':
	case KB_F1: case KB_F2: case KB_F3: case KB_F4: case KB_F5:
	case KB_F6: case KB_F7: case KB_F8: case KB_F9: case KB_F10:
	{  int fk;
           if (funcKey==NULL) break;
	   if ( ink>='0' && ink <='9') { fk=ink-'0';if (fk==0) fk=10;}
	      else fk=ink-KB_F1+1;
	   if ((funcKey[fk-1]) != NULL)
           {  void * saveHlp;
	      get_text(1,lastLine,maxCol(),lastLine,&saveHlp);
	      scrcolor(box_fg,box_bg);
	      goto_xy(col+1,row); print("%c",boxFrame[1]);
	      goto_xy(1,lastLine);scrcolor(FGmain,BGmain);clr_eol();
	      (*funcKey[fk-1])((npage-1)*height+k);	 
	      put_text(&saveHlp);
	      scrcolor(star_fg,box_bg);
	      goto_xy(col+1,row); print("<");
	   }
           break;
        }
        case 6:  /* ^F */
        case 'f':
        case 'F':
        {  char name[32]="";
           int key=correctStr(5,23,"Enter name(Esc for exit):",name,30,1);
           if(key)
           { char * x=strstr(menstr+1,name);
             if(x==NULL) 
             {
               if(blind) sortie(121); else messanykey(10,10, "Not detected"); 
               break;
             }
             k= ((x-menstr)-1)/menstr[0];
             npage = k/ height + 1;
             k = k% height + 1;
             goto label_1;
           }
        } break; 
     }
  }

label_3:
   put_text(hscr);
   put_text(&pscr);
   scrcolor(fcolor_tmp,bcolor_tmp);
   goto_xy(1,lastLine); clr_eol();
   *kk = 0;
}


static char * f_hlp=NULL;

static void f1_key(int x)
{
  int len;
  char ss[STRSIZ];
  if (f_hlp == NULL)  goto exi;
  len=strlen(f_hlp)-1;
  strcpy(ss,f_hlp);
  if (ss[len] =='*') sprintf(ss+len,"%d",x);
  if( show_help(ss) ) return;
            
  exi:  messanykey(10,10,"Help is absent");return;
}


static void f2_key(int x) {  show_help("man");  }

char *f3_mess[8]      ={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
void (*f3_key[8])(int)={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
            

void menu1(int col, int row,char* label,char* menustr,char* help,
                                                      void* hscr,int * kk)
{
void   (*funcKey[11])(int);
char * funcKeyMess[11];

int i;
static char f1_mess[10]="Help";
static char f2_mess[10]="Man";


char * help_tmp;

help_tmp=f_hlp;
f_hlp=help;

funcKey[0]=  f1_key; funcKeyMess[0]= f1_mess;
funcKey[1]=  f2_key; funcKeyMess[1]= f2_mess;

for(i=0;i<8;i++){ funcKey[2+i]=f3_key[i]; funcKeyMess[2+i]= f3_mess[i];}

menu0(col,row,label,menustr,funcKey,funcKeyMess,hscr,kk);

f_hlp=help_tmp;

}


int  str_redact(char* txt,int npos, int maxLen)
{  int    i, x, x0, y0;
   int key;
   int first;
   char    txttmp[STRSIZ];
   int t_color=Yellow;   
   int saveFcolor=fColor, saveBcolor=bColor;
   clearTypeAhead();

   strcpy(txttmp,txt);
   if(strlen(txt)>maxLen) txt[maxLen]=0;
   first = 1;
   
   scrcolor(t_color,Black);
   y0 = where_y();
   x0 = where_x();
   print("%-*s",maxLen,txt);
      if (npos<1) npos=1;
   if (npos>strlen(txt)) npos=strlen(txt)+1;
   x = x0 + npos - 1;
   scrcolor(Black,White); 
   goto_xy(x,y0);
   if(npos>strlen(txt)) print("%c",' ');else print("%c",txt[npos-1]);
   scrcolor(t_color,Black);
   goto_xy(x,y0);

   for(;;)
   {
      key = inkey();

      switch (key)
      {	
         case KB_MOUSE:
         if ( mouse_info.row == where_y() && mouse_info.but1 == 2 ) 
         {  npos=mouse_info.col-x0+1; }
         break;
         	
         case 9:      /*  Tab  */
            do
            {
               npos++;
               if (npos > strlen(txt)) npos = 0;
            }  while (!(npos==0 || txt[npos-1] == ','));
            x = x0 + i;
         break;

         case KB_END:
            npos=strlen(txt)+1;   /* End    */
         break;

         case KB_HOME:
            npos=1;   /*  Home  */
         break;


         case KB_RIGHT:
            npos++;
         break;

         case KB_LEFT:
            npos--;
         break;


         case KB_BACKSP:
	   if (npos > 1) npos--; else break;

         case KB_DC:
	   if (npos<=strlen(txt)) for(i=npos-1; txt[i]; i++) txt[i]=txt[i+1];
         break;


         case KB_ESC:
            strcpy(txt,txttmp);
            npos=1;

         case KB_ENTER:
         case KB_F1:
	 case KB_PAGEU:
	 case KB_PAGED:
	 case KB_UP:
	 case KB_DOWN:
            goto label_1; 

         default:
	    if (key > 31 && key <128 && npos<maxLen)                                         
            {                                                               
              if (first)                                                    
              {                                                             
                 sprintf(txt,"%c",key);                                     
                 npos = 2;                                                  
                 x = x0 + npos - 1;                                         
              }                                                             
              else                                                          
              {                                                             
                 if(strlen(txt)>=maxLen-1) txt[maxLen-2]=0; 
                 for (i=MIN(strlen(txt)+1,maxLen-1);i>= npos;i--) 
                                            txt[i] = txt[i-1]; 
                 txt[npos-1] = key; 
                 npos++;
              }                                                             

            }                                                               
      }

      goto_xy(x0,y0); print("%-*s",maxLen,txt);

      if(npos<=0) npos=1;
      if(npos>strlen(txt)) npos=strlen(txt)+1;  
      x=x0+npos-1;
      scrcolor(Black,White); 
      goto_xy(x,y0);        
      if(npos>strlen(txt)) print("%c",' ');else print("%c",txt[npos-1]);
      scrcolor(t_color,Black);
      goto_xy(x,y0);      
   
      first = 0;
   }

label_1: 
   goto_xy(x0,y0);
   if(key!=KB_ESC) print("%s",txt);
   scrcolor(saveFcolor,saveBcolor);
   return key;
}


static int  correctVar(char* txt,int x,int y, int  type, void * var,int clear)
{  int         npos,key,err;
   int xx;
   void *     pscr; 
   char * buff;
   int *  I;
   long * L;
   double * D;
   int maxLen;

   get_text(x,y,maxCol(),y,&pscr);   
   scrcolor(White,Black);
   goto_xy(x,y);
   print(txt);
   xx=where_x();
   

   if(type<0) {buff=var; maxLen=-type;} else 
   { buff=malloc(100); 
     switch (type)
     {
        case 'L':  L=var;  sprintf(buff,"%ld",*L); break;
        case 'D':  D=var;  sprintf(buff,"%lg",*D); break;
        case 'I':  I=var;  sprintf(buff,"%d" ,*I); break;
     }
     maxLen=MAX(15,strlen(buff));
   }    
   npos = 1;      

   do
   { 
      goto_xy(xx,y);
      key  = str_redact(buff,npos,maxLen);

      if (key == KB_ESC)
      {  
         put_text(&pscr);
         if(type>=0) free(buff);
         return 0;
      }
      err=0;
      if (key == KB_ENTER)
      {
/*         trim(buff); */
         if(type<0) err=1; else
         {
            switch (type)
            {   
              case 'L':   err=sscanf(buff,"%ld",L); break;
              case 'D':   err=sscanf(buff,"%lf",D); break;
              case 'I':   err=sscanf(buff,"%d" ,I); break;
            }
         }                                     
         if (err<=0)   messanykey(10,10," incorrect number"); 
      }
      
   }  while (err<=0);
   
   if(clear)put_text(&pscr); else free(pscr); 
   if(type>=0) free(buff);   
   return 1;   
} 


int improveStr(char * str, char * mark, char * format, ...)
{ va_list args;
  int i;
  char  buff[STRSIZ];
  char * sub;
          
  sub=strstr(str,mark);
  if(sub == NULL)  return -2;
  va_start(args, format);
  /*r=*/ vsprintf(buff,format,args);
  va_end(args);
  memcpy(sub,buff,strlen(buff));
  for(i=strlen(buff);i<strlen(mark);i++) sub[i]=' ';
  return 0;
}




int  correctDouble (int x,int y,char* txt,double * var,int clear)
{ return correctVar (txt,x,y,'D',var,clear); }

int  correctLong (int x,int y,char* txt,long * var,int clear)
{ return correctVar (txt,x,y,'L',var,clear); }

int  correctInt (int x,int y,char* txt,int * var,int clear)
{  return correctVar (txt,x,y,'I',var,clear); }

int  correctStr (int x,int y,char* txt,char * var, int maxLen, int clear)
{  return correctVar (txt,x,y,-maxLen,var,clear); }



