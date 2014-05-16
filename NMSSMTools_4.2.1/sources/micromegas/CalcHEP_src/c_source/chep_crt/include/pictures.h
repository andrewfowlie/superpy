#ifndef __PICTURES__ 
#define __PICTURES__

extern void  pictures(
  int  n_diag_tot,                  /* total number of pictures */ 
  int sizeX,int sizeY,              /* size of a picture */
  void (*pict)(int n,int x,int y), /* prorgam, which draw picture number 'n'
                                   starting from left-top corner (x,y) */  
  char* commandStr,                /* Auxilary commands, separeted by comma */                                      
  int (*command)(int n, char key), /* program to perform the auxilary command,
                                   'n'- is a current picture number,
                                   'key' is a capital letter from the command
                                    name*/    
char * help                        /* help file */   
                         );
#endif
