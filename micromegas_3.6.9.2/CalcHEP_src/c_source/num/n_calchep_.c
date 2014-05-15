/*
 Copyright (C) 1997,2006  Alexander Pukhov
*/

#include "chep_crt.h"
#include "files.h" 
#include "rw_sess.h"
#include "mc_menu.h" 
#include "param.h"
#include "interface.h"

#include"n_calchep_.h"

#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>
              
#include "../../include/num_out.h"
#include "dynamic_cs.h"



void n_comphep(void)
{
  clr_scr(FGmain,BGmain); 

  while(checkParam()) 
   if(mess_y_n(15,15, "Quit the session?"))
     {w_sess__(NULL); return;} else change_parameter(54,7,0);
  do
  { int err=monte_carlo_menu();
    switch(err)
    { case 1:printf("Energy is too small!\n"); sortie(123);
      case 2:printf("Can not evaluate cuts limlts\n"); sortie(124);
      case 3:printf("Can not evaluate regularization paremeters"); sortie(125);
    }
  }
  while(!mess_y_n(15,15,"Quit session?"));
 
  w_sess__(NULL);
}
