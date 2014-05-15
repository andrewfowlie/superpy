/*
	unescape.c
		Removes escaped characters from a string
		and splits the string at the & character.
		Useful for CGI scripts.
		Last modified 9 May 05 th

Usage:

1. unescape string1 string2 ...
   Works on string1 string2 ...

2. unescape
   Reads the input strings from stdin.

In bash do "eval `unescape $QUERY_STRING`", this should
work with both, METHOD=PUT and METHOD=GET.

Any characters the shell might hiccup on are replaced
by "_" on the l.h.s. of every item, and the r.h.s. is
wrapped in single quotes, e.g.

#abc=def	becomes	_abc='def'
7file=passwd	becomes	_7file='passwd'
xxx		becomes	xxx='1'

*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>


static inline int digit(const char c)
{
  return (unsigned)(c - '0') <= 9;
}

static inline int letter(const char c)
{
  return (unsigned)((c & 0xdf) - 'A') <= 'Z';
}


static void print( const char *s )
{
  char new[2*strlen(s)], *d = new, *lf = new, c;
  int quote = 0;

  while( (c = *s++) ) {
    if( c == '&' || c == '\n' ) {
      if( quote & 2 ) *d++ = '\'';
      else if( quote ) {
        *d++ = '=';
        *d++ = '1';
      }
      quote = 0;
      c = '\n';
      lf = d + 1;
    }
    else {
      quote |= 1;
      if( c == '+' ) c = ' ';
      else if( c == '%' ) {
        char hex[4];
        hex[0] = *s++;
        hex[1] = *s++;
        hex[2] = 0;
        c = strtol(hex, NULL, 16);
      }
      if( (quote & 2) == 0 ) {
        if( c == '=' ) {
          *d++ = c;
          c = '\'';
          quote = 2;
        }
        else if( digit(c) ) {
          if( d == lf ) *d++ = '_';
        }
        else if( !letter(c) ) c = '_';
      }
    }
    *d++ = c;
  }

  if( quote & 2 ) *d++ = '\'';
  else if( quote ) {
    *d++ = '=';
    *d++ = '1';
  }
  *d = 0;

  puts(new);
}


int main( int argc, char **argv )
{
  if( argc > 1 )
    while( --argc ) print(*++argv);
  else
    while( !feof(stdin) ) {
      char line[4096];
      *line = 0;
      fgets(line, sizeof(line), stdin);
      print(line);
    }

  return 0;
}

