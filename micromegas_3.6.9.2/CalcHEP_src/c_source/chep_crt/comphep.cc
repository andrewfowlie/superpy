
#include"chep_crt.h"

void * pscr=NULL;
int k=1;
int main(int argc,char** argv)
{

   start1("QQ=QQ","", "");


menu0(10,20,"my first menu","\24"
                            "menu function 1     "
                            "menu function 2     "
                            "menu function 3     "
                            "menu function 4     "
                            "menu function 5     ",NULL,NULL,&pscr,&k);
                           
  
 inkey();



 finish();
 return 0;
}
