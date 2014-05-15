
#include<string.h>
#include"paragraphs.h"

void readParagraphs(FILE * f, int n, rw_paragraph * rd_array)
{
   int i;
   char word[100];

   while(EOF!=fscanf(f,"%*[^#]"))
   { fscanf(f,"#%s",word);
     for(i=0;i<n;i++) if(strcmp(word,rd_array[i].keyword)==0) 
     {  int rd_code;
        if(!(rd_array[i].rw_command)) break;
	fscanf(f,"%*c");
        rd_code=rd_array[i].rw_command(f);
        if(rd_code<0) return;
        if(rd_code) fprintf(stderr,"Error in reading '%s'\n",word); 
        break;
     }
     if(i==n) 
     {  if(strcmp(word,"END")==0) return; else
        fprintf(stderr,"Unknown keyword '%s' The item is passed by.\n ",word);
     }
   }
}


void writeParagraphs(FILE * f, int n, rw_paragraph * wrt_array)
{
   int i;

   for(i=0;i<n;i++)
   { if(wrt_array[i].keyword&&strlen(wrt_array[i].keyword) 
       && wrt_array[i].rw_command )
     {
      fprintf(f,"\n#%s ",wrt_array[i].keyword);
      wrt_array[i].rw_command(f);
     } 
   }    
    fprintf(f,"\n#END\n"); 
}

