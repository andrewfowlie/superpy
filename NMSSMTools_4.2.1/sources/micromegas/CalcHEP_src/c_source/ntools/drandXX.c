/*
Restored version of standard drand48 and seed48 procedures 
*/

#include<stdlib.h>
#include<stdio.h>
#include"drandXX.h" 

static double float48=1./(((double) 0x10000)*((double) 0x10000)*((double) 0x10000));

static unsigned long Xlong=0x1234ABCD;
static unsigned long Xshort=0x330E;
#define Along   0x5DEEC
#define Ashort  0xE66D
#define Along16 0xDEEC0000
#define  Clong  0
#define  Cshort 0xB   
#define FIRST16 0xFFFF


double drandXX(void)
{
  unsigned  long  bot=Xshort*Ashort; 
  unsigned  long  top= Clong+ (bot>>16);

  bot=(bot&FIRST16)+Cshort;

  Xlong=(top + (bot>>16)+Along*Xshort+Ashort*Xlong+((Along16*Xlong)))&0xFFFFFFFF;
  Xshort=bot&FIRST16;
  return ( (double)Xshort + Xlong*(double)0x10000 )*float48; 
}

char * seedXX(char * init)
{
  unsigned long Xlong_,Xshort_;
  static char cbuff[13];
  
  sprintf(cbuff,"%08X%04X",Xlong,Xshort);
  
  if(init)
  { if(2==sscanf(init,"%8lX%4lX",&Xlong_,&Xshort_))
     { Xlong=Xlong_;
       Xshort=Xshort_;
     } else return NULL;
  }     
  return cbuff;
}
/*
#include<limits.h>
int main(void)
{ int i;
  int N=1000000;
  double sum=0;
   
  
 unsigned short seed16v[3];
 unsigned short * buff; 
 int z[3];

seed16v[0]=Xshort;
seed16v[1]=Xlong&0xFFFF;
seed16v[2]=Xlong>>16; 

printf(" %u %u %u\n",seed16v[0],seed16v[1],seed16v[2]);


 buff=seed48(seed16v);

printf(" %u %u %u\n",seed16v[0],seed16v[1],seed16v[2]);
 
 
 buff=seed48(seed16v);
 for(i=0;i<3;i++) z[i]=buff[i];
 printf(" %d %d %d\n",z[0],z[1],z[2]);
 
   
printf("%E %E\n", erand48(seed16v),drandXX());

printf("sizeoflong=%d\n",sizeof(long));
  for(sum=0,i=0;i<N;i++)   if( fabs(erand48(seed16v )-drandXX())> 1.E-10) printf("N=%d\n",i); 
printf("check OK\n");
}
*/
