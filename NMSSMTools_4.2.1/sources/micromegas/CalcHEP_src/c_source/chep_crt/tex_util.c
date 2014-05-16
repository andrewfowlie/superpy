/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "syst.h"
#include "crt.h"
#include "crt_util.h"
#include "help.h"
#include "tex_util.h"



void texPicture(int x1,int y1,int x2,int y2,int hsizePt,int  vsizePt)
{
  texxshift=-x1;
  texxscale=(double)hsizePt/(x2-x1);
  texymax1=y2;
  texyscale=(double)vsizePt/(y2-y1);
}



void texStart(char * f_name, char* upstr,char * charSize)
{
   out_tex=fopen(f_name,"w");

   texflag=1;

   f_printf(out_tex,"%% %s\n",upstr);
   f_printf(out_tex,"%%\\documentstyle[axodraw]{article}\n");
   f_printf(out_tex,"%%\\begin{document}\n");
   
   f_printf(out_tex,"{\n");

   f_printf(out_tex,"\\unitlength=1.0 pt\n");
   f_printf(out_tex,"\\SetScale{1.0}\n");
   f_printf(out_tex,"\\SetWidth{0.7}      %% line    size control\n");

   if (setTexCharSize(charSize)) f_printf(out_tex,"\\%s",charSize);
   else f_printf(out_tex,"\\scriptsize");
   f_printf(out_tex,"    %%  letter  size control\n");
}

void texFinish(void)
{
   f_printf(out_tex,"}\n%%\\end{document}\n");
   fclose(out_tex);
   texflag=0;
}




int texmenu(int * pictureX,int * pictureY,char * letterSize)
{

   char sizestr[]="\016"
        " tiny         "
        " scriptsize   "
        " footnotesize "
        " small        "
        " normalsize   "
        " large        "
        " Large        ";

   char  menustr[100];
   void* pscr=NULL;
   void* pscr2=NULL;
   int k=1;
   int l=1;

   for(;;)
   {

       sprintf(menustr,"\026"
       " picture(%4d,%4d)   "
       " Letter = %-12.12s"
       " Write to file        ",*pictureX,*pictureY,letterSize);

          menu1( 33,11,"LaTex menu",menustr,"h_tex1",&pscr,&k);
          switch(k)
          { case 0: return 0;
            case 1:  correctInt (15,20,"Picture X-size=",pictureX,1);
                     correctInt (15,20,"Picture Y-size=",pictureY,1);
                     break;
            case 2:  menu0( 12,12,"",sizestr,NULL,NULL,&pscr2,&l);
                     if(l!=0)
                     {  put_text(&pscr2);
                        sprintf(letterSize,"%12.12s",sizestr+2+(l-1)*14);
                      }
                      break;
            case 3:  put_text(&pscr);
                     return 1;
          }
   }

}
