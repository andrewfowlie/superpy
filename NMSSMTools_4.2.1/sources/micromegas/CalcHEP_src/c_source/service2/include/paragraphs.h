#ifndef __PARAGRAPHS__
#define __PARAGRAPHS__

#include<stdio.h>

typedef struct rw_paragraph {char * keyword;
                             int(*rw_command)(FILE *) ;
                            } rw_paragraph;  

extern void  readParagraphs(FILE * f, int n, rw_paragraph *  rd_array);
extern void writeParagraphs(FILE * f, int n, rw_paragraph * wrt_array);

#endif
