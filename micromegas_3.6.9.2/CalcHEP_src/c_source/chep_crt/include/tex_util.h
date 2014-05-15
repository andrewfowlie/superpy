#ifndef __TEX_UTIL__
#define __TEX_UTIL__
extern  int texmenu(int * pictureX,int * pictureY,char * letterSize);
extern  void texPicture(int x1,int y1,int x2,int y2,int vsizePt,int hsizePt);
extern  void texStart(char * fname, char* upstr,char * charScale);
extern  void texFinish(void);

#endif
