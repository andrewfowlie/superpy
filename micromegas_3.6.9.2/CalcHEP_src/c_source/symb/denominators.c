/*
Copyright (C) 1997, Alexander Pukhov
*/
#include"s_files.h"
#include"out_serv.h"
#include"denominators.h"
#include"saveres.h"
#include"prepdiag.h"
#include"pvars.h"
#include"process.h"
#include"c_out.h"


static int  nincount(int v,int l)
{int  i, vv, ll, summ;

   if (IN_PRTCL  & vcs.vertlist[v-1][l-1].prop) return 1;
   if (OUT_PRTCL & vcs.vertlist[v-1][l-1].prop) return 0;
   summ = 0;
   vv = vcs.vertlist[v-1][l-1].nextvert.vno;
   ll = vcs.vertlist[v-1][l-1].nextvert.edno;
   for (i = 1; i <= vcs.valence[vv-1]; i++)
      if (i != ll)
         summ += nincount(vv,i);
   return summ;
}

static int  checkOut( int p, int v, int l)
{ int i,vv,ll,s;

  if( vcs.vertlist[v][l].prop&IN_PRTCL) return -1;
  if( vcs.vertlist[v][l].prop&OUT_PRTCL)
  { int pt=vcs.vertlist[v][l].partcl; 
    if( pt==p) return 1; 
    else 
    {  long N=prtclbase[pt-1].N;
       if(N==22 || N==21) return 0;
       else return 2;
    } 
  }
  
  vv=vcs.vertlist[v][l].nextvert.vno-1;
  ll=vcs.vertlist[v][l].nextvert.edno-1;                                      


  for(s=0,i=0;i< vcs.valence[vv];i++) if(i!=ll)
  { int r= checkOut(p,vv,i);                                       
    if(r<0) return -1;
    if(r>1) return  2;
    s+=r;
  } 
  if(s) return s;
  
}    

static int rtypepropag(int v, int l)
{
  int p,r;
  int vv,ll;
  
  p=vcs.vertlist[v-1][l-1].partcl;

  r=checkOut(p,v-1,l-1);
  if(r==-1) 
  {   vv=vcs.vertlist[v-1][l-1].nextvert.vno-1; 
      ll=vcs.vertlist[v-1][l-1].nextvert.edno-1;
      p=vcs.vertlist[vv][ll].partcl;
      r=checkOut(p,vv,ll);
  }
  if(r!=1) r=0;
  return r;
}



int  ttypepropag(int v,int l)
{  int r;
   if (nin == 2) 
   { int r=nincount(v,l);
     if(r==1) return 1;
   }   
   r=rtypepropag(v, l);
   if(r==1) return 1;
   return 0;
}

      
      
static int stype(char * momStr)
{ 
  int c,i; 
  if(nin==1) return 1;
  for(i=0,c=0; momStr[i];i++) if (momStr[i]<=2) c++; 
  if(c&1) return 0; else return 1;
}       


/* momdep from prepdiag.h must be calculated before */
void  calcdenominators(vcsect vcs )
{ int  v, l, k;
  char buff[MAXINOUT+1];



   denrno = 0;
   for (v = 1; v <= vcs.sizet; v++)
   for (l = 1; l <= vcs.valence[v-1]; l++)
   {  edgeinvert *ln = &vcs.vertlist[v-1][l-1];
      if(!(ln->moment<0||((IN_PRTCL|OUT_PRTCL)&ln->prop)||pseudop(ln->partcl)))
      {
         for( k=1;k<=momdep[ln->moment-1][0];k++) buff[k-1]=momdep[ln->moment-1][k];
         buff[k-1]=0;

         k=-1; while(buff[++k]) if(buff[k]<0) buff[k]=-buff[k];
         if( (2*strlen(buff) > nin+nout) ||
             (2*strlen(buff)==nin+nout && !strchr(buff,1))
           )
         { int ll=0;
           char buff2[MAXINOUT+1];
           for(k=1;k<=nin+nout;k++){ if(!strchr(buff,k)) buff2[ll++]=k;}
           buff2[ll]=0;
           strcpy(buff,buff2);
         }
         k=0;
         while(buff[k])
         {  if(!k) k++;
            if(buff[k]<buff[k-1])
            { int c=buff[k];
              buff[k]=buff[k-1];
              buff[k-1]=c;
              k--;
            } else k++;
         }

         strcpy(denom[denrno].momStr,buff);
         if(v <= vcs.sizel) denom[denrno].power = 1; else denom[denrno].power = -1;

         denom[denrno].mass=modelVarPos(prtclbase[ln->partcl-1].massidnt);
         if(ttypepropag(v,l)&&!tWidths) denom[denrno].width = 0; else 
         denom[denrno].width=modelVarPos(prtclbase[ln->partcl-1].imassidnt);

         for (k = 0; k < denrno; k++)
         if ( !strcmp(denom[denrno].momStr,denom[k].momStr) &&
               denom[denrno].mass  ==  denom[k].mass &&
               denom[denrno].width ==  denom[k].width )
         {   denom[k].power=2;    goto label_1;}
         denrno++;
label_1:;
      }
   }
}


void  denominatorStatistic(int nsub, 
   int * n_swidth, int *n_twidth, int * n_0width, denlist * allDenominators, 
   FILE * fd)
{ 
   int i;
   catrec    cr;
   denlist    den_, den_tmp;
   deninforec   dendescript={0};
    
   (*n_swidth)  = 0;
   (*n_twidth)  = 0;
   (*n_0width) = 0;

   den_ =NULL;
   
   fseek(catalog,0,SEEK_SET);
   while (FREAD1(cr,catalog))
   {
      
      if (cr.nsub_ == nsub)
      { 
         whichArchive(cr.nFile,'r');
         dendescript.cr_pos = ftell(catalog) - sizeof(cr);

         fseek(archiv,cr.denompos,SEEK_SET);
         readDenominators();
         dendescript.tot_den=denrno; 

         for (i = 0; i < dendescript.tot_den; i++)
         {  
            dendescript.denarr[i].power=denom[i].power;
            dendescript.denarr[i].width=denom[i].width;
            den_tmp = den_;  
            while (den_tmp != NULL &&
              (  strcmp(denom[i].momStr,den_tmp->momStr)
              ||  denom[i].mass!=den_tmp->mass 
              ||  denom[i].width!=den_tmp->width ) ) den_tmp=den_tmp->next;
            if(den_tmp == NULL)
            {  
               den_tmp = (denlist)getmem_((unsigned)sizeof(denlistrec));
               den_tmp->next = den_;
               strcpy(den_tmp->momStr,denom[i].momStr);
               den_tmp->mass=denom[i].mass;
               den_tmp->width=denom[i].width;
               den_tmp->stype= stype(denom[i].momStr);
               den_ = den_tmp;
               if(denom[i].width) 
               { if(den_tmp->stype) den_tmp->order_num= ++(*n_swidth);
                         else       den_tmp->order_num= ++(*n_twidth);
               }  else              den_tmp->order_num= ++(*n_0width);
            }
            dendescript.denarr[i].order_num=den_tmp->order_num;
            dendescript.denarr[i].stype=den_tmp->stype;
         }      
         if(fd) FWRITE1(dendescript,fd);
      }  /* if CR.nsub_ =nsub */
   }
   
/*   if(ArchNum) fclose(archiv);  */
   whichArchive(0,0);
  *allDenominators=den_;
}
