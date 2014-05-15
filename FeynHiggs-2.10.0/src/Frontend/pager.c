/*
	pager.c
		starts the pager to view FeynHiggs
		Frontend output
		this file is part of FeynHiggs
		last modified 20 Jul 11 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>


#if NOUNDERSCORE
#define pageron_ pageron
#define pageroff_ pageroff
#endif


static int pagerpid = -1;
static int oldout = -1, olderr = -1;


void pageron_()
{
  if( isatty(1) ) {
    int handle[2];

    const char *pager = getenv("PAGER");
    if( pager == NULL ) pager = "less";
    else if( *pager == 0 ) return;

    if( oldout == -1 ) {
      oldout = dup(1);
      olderr = dup(2);
    }
    pipe(handle);

    pagerpid = fork();
    switch( pagerpid ) {
    case -1:
      return;

    case 0:
      dup2(handle[0], 0);
      close(handle[0]);
      close(handle[1]);
#if VT100
      putenv("LESS=-R");
#endif
      execlp(pager, pager, NULL);
      exit(-1);

    default:
      dup2(handle[1], 1);
      dup2(handle[1], 2);
      close(handle[0]);
      close(handle[1]);
      break;
    }
  }
}


void pageroff_()
{
  if( pagerpid != -1 ) {
    int status;
    close(1);
    close(2);
    wait(&status);
    dup2(oldout, 1);
    dup2(olderr, 2);
    pagerpid = -1;
  }
}

