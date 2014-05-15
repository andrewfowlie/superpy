/*
 Copyright (C) 1997, Alexander Pukhov 
*/

#include <unistd.h>
#include <sys/file.h>
#include <sys/utsname.h>
#include <signal.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include"syst.h"



char *  trim(char* p)
{  int n1=0, n2, k=-1;
   n2=(int)strlen(p)-1;
   if(n2<0) return p;
   while(!isgraph(p[n1]) && n1 <= n2) n1++;
   while(!isgraph(p[n2]) && n1 <  n2) n2--;
   while(++k < n2-n1+1)  p[k] = p[k+n1];
   p[k] = '\0';
   return p;
}


void  (*diskerror)(void) = NULL;

int f_printf(FILE *fp,char * format, ... )
{ va_list args;
  char dump[STRSIZ];
  int  r;
  va_start(args, format);

        vsprintf(dump,format,args);
	va_end(args);
	r = fputs(dump,fp);
	if (r == EOF)
	{ if(diskerror)(*diskerror)(); else sortie(65); }
	return r;
}


size_t f_write(void* ptr,  size_t size,  size_t n,  FILE * fp)
{ size_t nn;
   if ((size == 0)||(n == 0)) return 0;
   nn= fwrite(ptr,size,n,fp);
   if (nn != n) {if(diskerror) (*diskerror)() ; else sortie(65);}
   return nn;
}

void  unLockFile(int fd) 
{ 
if(fd<=0) return;

#ifdef LINK

return;
#endif 
lseek(fd,SEEK_SET,0);  


#ifdef FCNTL
{ 
  struct flock myLock;
  myLock.l_type= F_UNLCK;       /*  F_RDLCK ||  F_WRLCK || F_UNLCK */
  myLock.l_whence=SEEK_SET;
  myLock.l_start=0;
  myLock.l_len=10;
  ret=fcntl(fd, F_SETLK, &myLock);
}
#endif

#ifdef FLOCK
  ret = flock(fd,LOCK_UN);
#endif
 
#ifdef LOCKF
  ret = lockf(fd,F_ULOCK,10);
#endif

lseek(fd,SEEK_SET,0);
write(fd,"                              ",30);
  close(fd);
}

int setLockFile(char * fname)
{ int fd,ret;
  char lockTxt[100]; 

#ifdef LINK
{  int nlk;
   struct stat finfo;
   if(AuxLock[0]==0)
   {
     strcpy(AuxLock,'auxLock_');
     uname(AuxLock+8);
     sprintf(AuxLock+strlen(AuxLock),"_%d",getpid());
   }

   if(access(AuxLock),F_OK)
   { 
      FILE *f fopen(AuxLock,"w");
      fprintf(f,"%s",AuxLock);
      fclose(f);
   }

   stat(AuxLock,&finfo);
   nlk=finfo.st_nlink;
   if(!link(AuxLock, fname)) return 1;
   stat(AuxLock,&finfo); 
   if(finfo.st_nlink==nlk+1) return 1;
   return 0;
}
#endif



#if defined (FLOCK) || defined (LOCKF) || defined (FCNTL) 
  fd=open(fname,O_WRONLY|O_CREAT,0666);
  if(fd<0) return -1;
  lseek(fd,SEEK_SET,0); 
#else 
  return -1;
#endif


#ifdef FCNTL
{ 
  struct flock myLock;
  myLock.l_type= F_WRLCK;       /*  F_RDLCK ||  F_WRLCK || F_UNLCK */
  myLock.l_whence=SEEK_SET;
  myLock.l_start=0;
  myLock.l_len=0;
  ret=fcntl(fd, F_SETLK, &myLock);
}
#endif

#ifdef FLOCK  
  ret=flock(fd,LOCK_EX|LOCK_NB);   
#endif

#ifdef LOCKF
  ret = lockf(fd, F_TLOCK, 10);
#endif

  
  if(ret) { close(fd); return 0;}
  {
      struct utsname buff;
      uname(&buff);
      sprintf(lockTxt," process %d on %s\n",getpid(),buff.nodename);
  }                       
  write(fd,lockTxt,strlen(lockTxt));
  lseek(fd,SEEK_SET,0);
  return fd;
}


static char * lockFileName=NULL;
static int   lockFileFD=0;

int writeLockFile(char * fname)
{   
  if(lockFileFD) return lockFileFD;
  lockFileFD=setLockFile(fname);
  if(lockFileFD>0)
  { lockFileName=malloc(strlen(fname)+1);
    strcpy(lockFileName,fname);
  }
  return lockFileFD;
}

void sortie(int code) 
{ 
  if(lockFileFD>0)
  {  unLockFile(lockFileFD);
     unlink(lockFileName);
  }   
  exit(code);
}

