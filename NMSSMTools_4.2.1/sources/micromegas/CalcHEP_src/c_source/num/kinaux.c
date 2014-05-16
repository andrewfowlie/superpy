/*
 Copyright (C) 1997, Alexander Pukhov 
*/
#include "interface.h"
#include "syst.h"
#include"kinaux.h"
#include"4_vector.h"
#include"sortarr.h"

kinmtc_ kinmtc_1[10];
/*
static void lvprn(char * mess, char * lv)
{ int i;
  printf("%s=(", mess);
  for(i=0;lv[i];i++) printf("%d",lv[i]);
  printf(")");
}
*/


void lvmirr_(char *lv)
{
  char buff[PLISTLEN];
  int i1,i2=0;

  strcpy(buff,lv);
  for(i1=1;i1<=nin_int+nout_int;i1++) { if( !strchr(buff,i1)) lv[i2++]=i1;} 
  lv[i2]=0; 
} 

void coninv_(char *lv)
{
    if (lv[0] == 0) return;
    SORTARR(lv,strlen(lv));

    if( 2*strlen(lv) < nin_int +nout_int)  return;
    if( (2*strlen(lv) == nin_int +nout_int) && lv[0]==1) return;
    
    lvmirr_(lv);
} 

int eqvect_(char *lv1, char *lv2)
{
    char lvbuf1[10], lvbuf2[10];
    strcpy(lvbuf1,lv1); coninv_(lvbuf1);
    strcpy(lvbuf2,lv2); coninv_(lvbuf2);
    return !strcmp(lvbuf1,lvbuf2);
} 
  
int spole_(char *lv)
{
   int i1=0,i2=0,n,i=0;
   
   while(n=lv[i++]){ if(n>nin_int) i1++; else i2++;}  
   if(i2==nin_int) i2=0;
   return  ! (i1 && i2);
} 


static int lvdiff_(char *lv1, char *lv2, char *lvres)
{
    int i1,i,nc=0,r1=0;
   
    i1=strlen(lv1);
    nc=0;
    for(i=0;i<i1;i++) if(strchr(lv2,lv1[i])) nc++; else lvres[r1++]=lv1[i];
    
    if(nc==strlen(lv2))
    {
      lvres[r1]=0;
      return 0;
    } else if(nc==0)
    { 
       sprintf(lvres,"%s%s",lv1,lv2); 
       lvmirr_(lvres);
       return 0;
    } else return 1;
    
} 


static int isknown(char *lv,int level)
{
    char j, k, n, len;
    int numbx=nin_int +nout_int +1;

    if (eqvect_(lv, kinmtc_1[level-1].lvin))  return 0;
    
    len = strlen(lv);

    for(n=0; n<len; n++)
    {  int icmp=1;
       for (k=0; k<2; k++) for (j=0; j<level-1; j++) 
       {
          if (kinmtc_1[j].lvout[k][1] == 0 && 
              kinmtc_1[j].lvout[k][0] == lv[n])  icmp = 0;
       }
	
       for (k = 1; k <= nin_int; ++k) if (k == lv[n] || lv[n] == numbx) icmp=0;
       if (icmp) return 0;
    }
    return 1;
}


void sngpos_(char *lv, int * ndec, int * nclust, char *lvaux)
{
    int i, k,ln[2];
    int ierr[2];
    char lvaux_cp[2][10];

    if (lv[0] == 1 && lv[1] == 2 && lv[2] == 0 && nin_int == 2) 
    {
	lvaux[0] = 0;
	*ndec = 0;
	*nclust = 0;
	return; 
    }
    
    for (i = 0; i < nout_int - 1; ++i)
    {  
      for (k = 0; k < 2; k++) 
      {
         ierr[k]=lvdiff_(lv,kinmtc_1[i].lvout[k],lvaux_cp[k]);
         if (!ierr[k]) ierr[k] = !isknown(lvaux_cp[k],i+1); 
         if (!ierr[k]) ln[k]=strlen(kinmtc_1[i].lvout[k])+strlen(lvaux_cp[k]);
      } 
      if (! (ierr[0]&&ierr[1]) )  
      {
         *ndec = i+1;
         if((!ierr[0]) && (!ierr[1])) k=(ln[0]<ln[1])? 1:2;
         else k=ierr[0]?2:1;
         *nclust = k;
         strcpy( lvaux,lvaux_cp[k-1]);
         return;
      }
    }
    printf("error in sngpos_\n"); sortie(51);
}
