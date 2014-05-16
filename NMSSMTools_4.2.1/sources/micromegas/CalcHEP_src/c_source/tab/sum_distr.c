#include "histogram.h"

 
int main(int np, char ** par)
{
  int n;
  char *process=NULL;
  FILE*f;


  if(np<3) 
  { printf(" This  routine is intended to sum the  distributions produced\n" 
           " in  calcHEP numerical sessions (files dist_#).\n"
           " The names of files must be submitted as parameters\n"
           " Resulting distribution is directed to stdout(screen)\n");
    return 1;
  }

  for(n=1;n<np;n++)
  { 
    f=fopen(par[n],"r"); 
    if(!f) return 2;
    if(add_hist(f,&process)) {fclose(f);return 3;}
    fclose(f);
  } 
  wrt_hist2(stdout,process);

  return 0;
}
