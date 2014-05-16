#ifndef __EVENT_TO_PYTH__
#define __EVENT_TO_PYTH__

#define MAXNUP 500
typedef struct  
{    int     NUP,IDPRUP;
     double  XWGTUP,SCALUP,AQEDUP,AQCDUP;
     int     IDUP[MAXNUP],ISTUP[MAXNUP],MOTHUP[MAXNUP][2],ICOLUP[MAXNUP][2];
     double  PUP[MAXNUP][5], VTIMUP[MAXNUP],SPINUP[MAXNUP];
} hepeup_str;


#define MAXPUP 100
typedef struct 
{  int IDBMUP[2];
   double  EBMUP[2];
   int PDFGUP[2],PDFSUP[2];
   int  IDWTUP,NPRUP;
   double XSECUP[MAXPUP],XERRUP[MAXPUP],XMAXUP[MAXPUP];
   int LPRUP[MAXPUP];
} heprup_str;

#define R_ heprup_
#define E_ hepeup_


extern heprup_str R_;
extern hepeup_str E_;

extern int readslha_(void);
extern int scandir_( char * dirname, int len);
extern void upinit_(void);
extern int upevnt_(void);
extern void closeevents_(void);
extern void eventstat_(double * cs, int * nevents);
extern void writeInfo(void);

extern  int  openeventfile_(char *fname, int len);
extern void  closeeventfile_(void);
extern  int  readeventheader_(void);
extern  int  readevent_(void);

extern void  printProcInfo(FILE*f);
#endif
