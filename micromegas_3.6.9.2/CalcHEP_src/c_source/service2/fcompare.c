#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <unistd.h>
                     
#include "fcompare.h"

#define SIZE 4096

int fcompare(char *f1, char*f2)
{
  int  d1,d2,r1,r2,d,err;
  char buf1[SIZE], buf2[SIZE];

 
  d1=open(f1,O_RDONLY); if(d1<0 ) return 1;
  d2=open(f2,O_RDONLY); if(d2<0 ) {close(d1); return 2;}
  
  for(err=0;;)
  {
    r1=read(d1,buf1, SIZE);
    r2=read(d2,buf2, SIZE); 
    if(r1!=r2 ) {err=3;break;} 
    if(r1==0)   break;
    
    d=memcmp(buf1, buf2, r1);

    if(d) {err=3; break;}
    if(r1<SIZE) break;
  }  
  close(d1); close(d2); return err;
}
