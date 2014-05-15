/*
	table.c
		reorders FeynHiggs output for plotting purposes
		this file is part of FeynHiggs
		last modified 29 Feb 12 th

	Syntax: FeynHiggs file.in flags | table var1 var2... > file.out


Note: Variable names including a `-' are assumed to be decay channels
      or couplings.  Unlike other names, these can be valid items even
      though they may not be present in a particular FeynHiggs run if
      the channel is not open.  Secondly, the decay/coupling information
      may have a varying number of fields.

      To make the column arrangement more predictable, table *always*
      fills in three fields (0 if not present) when such a variable is
      requested.
*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#if VT100
#define BOLD "\e[1m\e[34m"
#define ERROR "\e[1m\e[31m"
#define RESET "\e[0m"
#else
#define BOLD
#define ERROR
#define RESET
#endif

#define MAXIN 2048
#define MAXOUT 20480


typedef struct {
  char *data;
  int ncur, nmin, ntot;
} value;


static void writeout(int argc, value *val)
{
  char tmp[MAXOUT], *t = tmp;
  int c, p;

  for( c = 1; c < argc; ++c ) {
    if( val[c].nmin == 0 ) val[c].nmin = val[c].ncur;
    if( val[c].ncur ) {
      if( val[c].ntot == 0 ) val[c].ntot = val[c].ncur;
      t += sprintf(t, "%s", val[c].data);
    }
    for( p = val[c].ncur; p < val[c].nmin; ++p ) {
      *t++ = '\t';
      *t++ = '0';
    }
    val[c].ncur = 0;
  }
  *t++ = '\n';
  *t = 0;
  fputs(tmp, stdout);
}


int main(int argc, char **argv)
{
  char out[MAXOUT], *d = out;
  value val[argc];
  int c = 0, term = 0, col;

  memset(val, 0, sizeof val);
	/* always three output fields for decays and couplings */
  for( c = 1; c < argc; ++c )
    if( strchr(argv[c], '-') ) val[c].nmin = 3;

  while( !feof(stdin) ) {
    char in[MAXIN], *s = in;
    *in = 0;
    fgets(in, sizeof in, stdin);

    if( strstr(in, "error") || strstr(in, "warning") ) {
      fputs(in, stderr);
      continue;
    }

    if( strstr(s, "END OF OUTPUT") ) {
      writeout(argc, val);
      d = out;
      continue;
    }

    if( *s == '%' ) ++s;
    if( *s != '|' ) continue;
    ++s;

    if( term != '\\' ) {
      char var[128];
      int next = -1;
      sscanf(s, "%s =%n", var, &next);
      if( next == -1 ) continue;
      s += next;
      for( c = 1; c < argc; ++c )
        if( strcmp(var, argv[c]) == 0 ) break;
      if( c == argc ) continue;
      val[c].data = d++;
      val[c].ncur = 0;
    }

    term = 0;
    do {
      int before, behind = -1;
      d[-1] = '\t';
      sscanf(s, " %n%s%n", &before, d, &behind);
      if( behind == -1 || (term = *d) == '\\' ) {
        d[-1] = 0;
        break;
      }
      ++val[c].ncur;
      d += behind - before + 1;
      s += behind;
    } while( *s != '\n' );
  }

  fputs(BOLD "Column Layout:\n" RESET, stderr);

  col = 0;
  for( c = 1; c < argc; ++c )
    switch( val[c].nmin ) {
      int i;
    case 0:
      fprintf(stderr, ERROR "%-26s not found" RESET "\n", argv[c]);
      break;
    case 1:
      fprintf(stderr, "%-26s Column %d\n", argv[c], ++col);
      break;
    default:
      d = out;
      for( i = 0; i < val[c].nmin; ++i )
        d += sprintf(d, " %d", ++col);
      if( val[c].ntot == 0 )
        sprintf(d, "\t" ERROR "Warning: all zero" RESET);
      fprintf(stderr, "%-26s Columns%s\n", argv[c], out);
      break;
    }

  return 0;
}

