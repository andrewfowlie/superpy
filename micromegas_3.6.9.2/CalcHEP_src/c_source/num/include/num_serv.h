#ifndef __NUM_SERV_
#define __NUM_SERV_
  typedef double  (*r_func)(void);
  extern void  paramdependence(r_func ff, char * procname,char * resultname);
#endif

