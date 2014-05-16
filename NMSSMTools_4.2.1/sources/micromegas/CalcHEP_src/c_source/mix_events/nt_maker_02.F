#define NEW  
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
    
      character Line*10000
      Integer N1,N2,K

      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      INTEGER I,NEVENT

#ifdef NEW
      INTEGER ERR
      INTEGER  OpenEventFile,ReadEventHeader,readEvent
      EXTERNAL OpenEventFile,ReadEventHeader,readEvent
#endif      
ccccccccccccccccccccccccccccccccccccc 
      CHARACTER*40 NTNAME
      PARAMETER(NWPAWC = 500000)
      COMMON/QUEST/IQUEST(100)
      COMMON/PAWC/PAW(NWPAWC)
      INTEGER IQUEST
      INTEGER*4 NPRT,IKF(100),ISTA(100),IPAR(100)
      REAL*4    CS,XPUP(100),YPUP(100),ZPUP(100),EPUP(100)
      COMMON/NTVAR/ NPRT,IKF,IPAR,ISTA,CS,XPUP,YPUP,ZPUP,EPUP 
      character*20 NAME

      SAVE
                 
      NUMBER=0
      NQUANT=50000
      CALL HLIMIT(NWPAWC)    
      
  1   NUMBER=NUMBER+1 !number of ntuples for a givent QUANT
c      WRITE (NAME, '(I10)') NUMBER
c      name=adjustl(name)

      IF(NUMBER.lt.10) THEN
        WRITE (NAME,fmt='(A12,I1,A3)')  'event_mixer_',NUMBER,'.nt'
      ELSEIF (NUMBER.lt.100) THEN
        WRITE (NAME,fmt='(A12,I2,A3)')  'event_mixer_',NUMBER,'.nt'
      ELSEIF (NUMBER.lt.1000) THEN
        WRITE (NAME,fmt='(A12,I3,A3)')  'event_mixer_',NUMBER,'.nt'
      ELSEIF (NUMBER.lt.10000) THEN
         WRITE (NAME,fmt='(A12,I4,A3)')  'event_mixer_',NUMBER,'.nt'
      ELSE
         WRITE (NAME,*)  'event_mixer_',NUMBER,'.nt'
      ENDIF
      
      CALL HROPEN(10,'calcnt',NAME, 'N',1024,ISTAT)

c      CALL HROPEN(10,'calcnt',
c     &            'event_mixer_'//NAME(2:LEN_TRIM(NAME))//'.nt',
c     &            'N',1024,ISTAT)

c      print *, LEN_TRIM(NAME), name,NUMBER
c      print *,'**',trim(name),'***',name
      CALL HBNT(10,'calcnt hbnt',' ')
      CALL HBNAME(10,'data',NPRT,'NPRT[0,100],
     >IKF(NPRT),IPAR(NPRT),ISTA(NPRT),CS,XPUP(NPRT),
     >YPUP(NPRT),ZPUP(NPRT),EPUP(NPRT)')
ccccccccccccccccccccccccccccccccccccc 
      
      IF(NUMBER.gt.1) goto 20
#ifdef NEW 
      if(iargc().ne.1) then 
        write(*,*)'Name of LH event file should be passed as parameter'
        stop
      endif
      call getarg(1,Line)
      Err=OpenEventFile(Line)
      if(Err.ne.0) then 
         write(*,*) 'Can not open file ', Line
         stop
      endif
      Err=ReadEventHeader()
      if(Err.ne.0) then 
        write(*,*) 'Can not filed HEPRUP block in input file'
        stop
      endif
#else       
      OPEN(99,FILE='event_mixer.out',STATUS='OLD')
c      OPEN(89,FILE='paw_input.out'  ,STATUS='UNKNOWN')
     
 10   read(99,'(2000A)',err=400,END=400) Line
      IF(Line(1:1).eq.'#') goto 10
      
      read(Line(1:144),'(2I6,2D18.10,6I6,3D18.10,I6)') 
     &IDBMUP,EBMUP,PDFGUP,PDFSUP,
     &IDWTUP,NPRUP,XSECUP(1),XERRUP(1),XMAXUP(1),
     &LPRUP(1)     
#endif      
     
 20   NEVENT=0


600   Continue
#ifdef NEW 
      if(readEvent().ne.0) goto 400
#else 
      read(99,'(2000A)',err=400,END=400) Line
#endif
      NEVENT=NEVENT+1
      
      if(NEVENT.gt.NQUANT)  THEN
        call hrout(10,icycle,' ')    
        call hrend('calcnt')
	goto 1
      endif
 
      IF(NEVENT.eq.1000*(NEVENT/1000)) print *, NEVENT , NUMBER   
#ifndef NEW
      read(Line(1:84),300,err=400)   NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDU
 300  format(2I6,4D18.10)
      N2=84
      DO K=1,NUP
       N1=N2+1
       N2=N1+15
       read(Line(N1:N2),100) IDUP(K)
      ENDDO
      
      DO K=1,NUP
       N1=N2+1
       N2=N1+5
       read(Line(N1:N2),100) ISTUP(K)
      ENDDO
      DO K=1,NUP
       DO J=1,2
       N1=N2+1
       N2=N1+5
         read(Line(N1:N2),100) MOTHUP(J,K)
       ENDDO
      ENDDO
      DO K=1,NUP
       DO J=1,2
       N1=N2+1
       N2=N1+5
       read(Line(N1:N2),100) ICOLUP(J,K)
       ENDDO
      ENDDO
      DO K=1,NUP
       DO J=1,5
       N1=N2+1
       N2=N1+17
       read(Line(N1:N2),200) PUP(J,K)
       ENDDO
      ENDDO
      DO K=1,NUP
       N1=N2+1
       N2=N1+17
       read(Line(N1:N2),200) VTIMUP(K)
      ENDDO
      DO K=1,NUP-1
       read(Line(N1:N2),200) SPINUP(K)
       N1=N2+1
       N2=N1+17
      ENDDO
       N1=N2+1
       N2=N1+17
       read(Line(N1:N2),200) SPINUP(NUP)
  
 100  FORMAT(I16)
 200  FORMAT(D18.10)
 
#endif 
ccccc fill ntuple cccccccccccccccc
      NPRT=NUP
      CS=XSECUP(1)
      DO I=1,NUP
c---------------------              
	  IKF(I)=IDUP(I)
	  IPAR(I)=MOTHUP(1,I)
	  ISTA(I)=ISTUP(I)
	  XPUP(I)=PUP(1,I)
	  YPUP(I)=PUP(2,I)
	  ZPUP(I)=PUP(3,I)
	  EPUP(I)=PUP(4,I)
c---------------------              

c        IF(ISTUP(I).eq.1) THEN
c	  write(89,'(I10,$)') IDUP(I)
c	  DO K=1,4
c            write(89,'(D14.6,$)') PUP(K,I)
c	  ENDDO
c	ENDIF
      ENDDO
c      write(89,'(D14.6)') XSECUP(1)
    
      CALL HFNT(10)

cccccccccccccccccccccccccccccccccc 
      
      goto 600
 400  NUP=0
      
      call hrout(10,icycle,' ')    
      call hrend('calcnt')
#ifdef NEW
      Call CloseEventFile
#else 
      CLOSE(99)
      CLOSE(89)
#endif

      STOP
      END




