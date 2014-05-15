/*
	logfile.c
		I/O redirection for logfiles
		this file is part of FormCalc
		last modified 6 May 04 th
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

#ifdef HAVE_UNDERSCORE
#define openlog openlog_
#define closelog closelog_
#define flush flush_
#endif


static int truestdout = -1;


int openlog(char *dir, const int *serial, const int len)
{
  int logfile;
  struct stat filestat;
  char filename[512];
  char *last = memccpy(filename, dir, ' ', len);
  *--last = 0;

  if( stat(filename, &filestat) == 0 ) {
    if( !S_ISDIR(filestat.st_mode) ) {
      fprintf(stderr, "%s already exists\n", filename);
      exit(1);
    }
  }
  else if( mkdir(filename, 0755) == -1 ) {
    fprintf(stderr, "Cannot create directory %s\n", dir);
    exit(1);
  }

  sprintf(last, "/%06d", *serial);
  if( stat(filename, &filestat) == 0 && (filestat.st_mode & 0100) == 0 ) {
    printf("%s already complete\n", filename);
    return 1;
  }

  logfile = creat(filename, 0744);
  if( logfile == -1 ) {
    fprintf(stderr, "Cannot create %s\n", filename);
    exit(1);
  }

  { int fortranstdout = 6; flush(&fortranstdout); }
  puts(filename);
  fflush(stdout);

  if( truestdout == -1 ) truestdout = dup(1);
  dup2(logfile, 1);
  close(logfile);

  return 0;
}


void closelog()
{
  { int fortranstdout = 6; flush(&fortranstdout); }
  fchmod(1, 0644);
  dup2(truestdout, 1);
}

