/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include"err_code.h"
#include "crt_util.h"
int err_code=0;

void  errormessage(void)
{
   switch(err_code)
   {
      case 1:
         messanykey(10,10,"Lost of precision in SQME ");
         break;
      case 2:
         messanykey(10,10,"Zero denominator");
      break;
      
      case 3: 
         messanykey(10,10,"Too many points for integration.");
      break;
       
      case 4:
         messanykey(10,10,"Energy too small");
      break;

      case 5:
         messanykey(10,10,"Can not evaluate constraints");
      break;
      
      case 10:
         messanykey(10,10,"User Break");
      break;

      default:
         messanykey(10,10," Error ?? ");
   }
}
