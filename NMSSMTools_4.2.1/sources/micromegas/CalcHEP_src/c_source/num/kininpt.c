/*
 Copyright (C) 1997, Dmitry Kovalenko
*/

#include "interface.h"
#include "syst.h"
#include "crt_util.h"
#include "kinaux.h"
#include"kininpt.h"

static void wrt_(char *mes, char *lvv)
{
   int i;

   scrcolor(Blue,BGmain); print(mes);
   scrcolor(FGmain,BGmain); for(i=0; lvv[i];i++ )  print("%d",lvv[i]);
} 

static int fill_(char *lvar)
{
  int j;
  int nn, ncr;
 
  char rdar[STRSIZ];

  rdar[0]=0;
  if(str_redact(rdar,1,10)==KB_ESC) return 0 ;

  lvar[0] = 0;
  ncr = 0;
    
  for (j = 0; j <strlen(rdar); ++j) 
  {
	nn =  rdar[j]   - '0';
	if (nn > nin_int && nn <= nin_int+nout_int && !strchr(lvar,nn)) 
	{
	    lvar[ncr] = nn;
	    lvar[++ncr]=0;
	}
  }
  if(ncr == 0)
  {   messanykey(10,10,"ERROR: 1st cluster is empty!");
        return 0;
  }
  return 1;
} 

static void uskin_(void)
{
    int i,l1,l2,nc,i1;
    char icmp[PLISTLEN];
    int i_moth[10],k_moth[10];

    i = 0;
L1:

       
/* * Fill the vector of the in-particles for this decay */
    if (i == 0) 
    {
	kinmtc_1[0].lvin[0] = 1;
	if (nin_int == 2)  kinmtc_1[0].lvin[1] = 2;
	kinmtc_1[0].lvin[nin_int]=0;
    } else  strcpy(kinmtc_1[i].lvin,kinmtc_1[i_moth[i]].lvout[k_moth[i]]);
L5:
    goto_xy(3,11+i); print("%40.40s","");
    goto_xy(3,11+i); wrt_("in=  ", kinmtc_1[i].lvin);

    if (i > 0 && kinmtc_1[i].lvin[2] == 0) 
    {
	kinmtc_1[i].lvout[0][0] = kinmtc_1[i].lvin[0];
	kinmtc_1[i].lvout[1][0] = kinmtc_1[i].lvin[1];
	kinmtc_1[i].lvout[0][1] = 0;
	kinmtc_1[i].lvout[1][1] = 0;
    }
    else if (i == 0 && nout_int == 2) 
    {
	kinmtc_1[0].lvout[0][0] = nin_int + 1;
	kinmtc_1[0].lvout[1][0] = nin_int + 2;
        kinmtc_1[0].lvout[0][1] = 0;
        kinmtc_1[0].lvout[1][1] = 0;
    }
    else
    {  int j;
       goto_xy(13,11+i); scrcolor(Blue,BGmain); print("-> out1= ");
       if (!fill_(kinmtc_1[i].lvout[0]))  
       {
          if (i>0) i--;
          goto_xy(1,where_y()); clr_eol();
          goto L1;
       }
           
       nc=0;
       i1=0;
       strcpy(icmp,kinmtc_1[i].lvin);
       if (i == 0) lvmirr_(icmp);


       for(j=0;j<strlen(icmp);j++)
       if(strchr(kinmtc_1[i].lvout[0],icmp[j])) nc++;else 
          kinmtc_1[i].lvout[1][i1++]=icmp[j];

       kinmtc_1[i].lvout[1][i1]=0;  
    
       if(nc != strlen( kinmtc_1[i].lvout[0]))     
       {
          messanykey(10,10, "ERROR: particle(s) have to be from out state");
          goto L5;
       }

       if (!strlen(kinmtc_1[i].lvout[1])) 
       {
          messanykey(10,20,"ERROR: 2-nd cluster is empty");
	  goto L5;
       }
    }
    goto_xy( 13,11+i);  wrt_("-> out1= ", kinmtc_1[i].lvout[0]);
    goto_xy( 25,11+i);  wrt_("out2= ", kinmtc_1[i].lvout[1]);
    
    l1=strlen(kinmtc_1[i].lvout[0]);
    l2=strlen(kinmtc_1[i].lvout[1]);

    if (l1>1) {i_moth[i+1] =i; k_moth[i+1]=0;}
    if (l2>1) {i_moth[i+l1]=i; k_moth[i+l1]=1;}
   
    if (i < nout_int - 2)   { ++i; goto L1; }
    return;
} /* uskin_ */

int entkin_(void)
{
    int i__1,i;
    void * pscr=NULL;
    get_text(1,5,52,20,&pscr);
    clrbox(1,5,52,20);    

    goto_xy(3,6); scrcolor(Red,BGmain);
     print("========= Current kinematical scheme =========\n");
    i__1 = nout_int-1;
    scrcolor(FGmain,BGmain);
    for (i = 0; i < i__1; ++i)
    {   goto_xy(3,7+i); print("%40.40s","");
        goto_xy(3,7+i); wrt_("in= ",kinmtc_1[i].lvin);
        goto_xy(13,7+i);wrt_( "-> out1= ", kinmtc_1[i].lvout[0]);
        goto_xy(25,7+i);wrt_("out2= ", kinmtc_1[i].lvout[1] );
    }
    goto_xy(3,6+i__1+1); scrcolor(Red,BGmain);
    print("==============================================\n");

    if( mess_y_n(10,13,"Input new kinematics?"))
    {  uskin_();
       put_text(&pscr);
       return 2;
    } else
    {  put_text(&pscr);
       return 0;
    }  
} /* entkin_ */

void stdkin_(void)
{
    int i__1, i__2;
    int i, j;

    i__1 = nout_int - 1;
    
    for (i = 0; i < i__1; ++i) 
    {
	kinmtc_1[i].lvout[0][0] = i + nin_int+1;
	kinmtc_1[i].lvout[0][1] = 0;

	i__2 = nout_int - i - 1;
	for (j = 0; j < i__2; ++j) kinmtc_1[i].lvout[1][j] = i + nin_int + j + 2;
	kinmtc_1[i].lvout[1][i__2]=0;
    }

    kinmtc_1[0].lvin[0] = 1;
    if (nin_int == 2) kinmtc_1[0].lvin[1] = 2;
    kinmtc_1[0].lvin[nin_int]=0;
    
    for (i = 1; i < i__1; ++i) strcpy(kinmtc_1[i].lvin,kinmtc_1[i-1].lvout[1]);
} /* stdkin_ */


int wrtkin_(FILE *nchan)
{
    int i__1;
    int i, k, c;

    int ndec=nout_int-1;
    
   fprintf(nchan,"\n");
   i__1 = ndec;
   for (i = 0; i < i__1; ++i) 
   {
     for (k=0;c=kinmtc_1[i].lvin[k];k++) fprintf(nchan,"%d",c);
     fprintf(nchan," -> ");
     for (k = 0; c=kinmtc_1[i].lvout[0][k]; ++k)  fprintf(nchan,"%d",c);
     fprintf(nchan," , ");
     for (k = 0; c=kinmtc_1[i].lvout[1][k] ; ++k) fprintf(nchan,"%d",c);
     fprintf(nchan,"\n");
    }	
    return 0;
}

int rdrkin_(FILE *nchan)
{
    int i;
    char strin[10],strout1[10],strout2[10];

    for (i = 0; i <  nout_int-1; ++i) 
    { int l,k,c;
      fscanf(nchan,"%s -> %s , %s", strin,strout1,strout2); 

      for(k=0,l=0;c=strin[k];k++)  if(c!=' ')kinmtc_1[i].lvin[l++]=c-'0';
      kinmtc_1[i].lvin[l]=0; 

      for(k=0,l=0;c=strout1[k];k++)if(c!=' ')kinmtc_1[i].lvout[0][l++]=c-'0';
      kinmtc_1[i].lvout[0][l]=0; 
      
      for(k=0,l=0;c=strout2[k];k++)if(c!=' ')kinmtc_1[i].lvout[1][l++]=c-'0';
      kinmtc_1[i].lvout[1][l]=0; 

    }

    return 0;
}
