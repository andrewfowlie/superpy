double I1(msq,mq,mne)
{

double dsqn=(msq-mne)*(msq+mne);
      DSQN=MSQ**2-MNE**2
double dqn=(mq-mne)*(mq+mne)
      DQN =MQ**2 -MNE**2
double dsqq=(msq-mq)*(msq+mq)      
      DSQQ=MSQ**2-MQ**2
      
      D3M=  msq**2+mq**2-MNE**2
      D3M1= msq**2-mq**2-MNE**2
      
      R1  =DSQQ/MNE**2
      R2  =DQN /MSQ**2
      R3  =DSQN/MQ**2
      
      DEL =2d0*MNE**2*(MQ**2+MSQ**2)-MNE**4-(MSQ**2-MQ**2)**2
      ADE =abs(DEL)
C$      IF(ADE.lt.1E-20)PRINT *,'DEL=0(!!!)'
      LAM=0d0
      IF(DEL.GE.0)LAM=2d0*ATAN(sqrt(ADE)/D3M)/sqrt(ADE)
      IF(DEL.LT.0)LAM=LOG((D3M+sqrt(ADE))/(D3M-sqrt(ADE)))/sqrt(ADE)

      CMD=0d0
      
      IF    (I.eq.1) THEN
        
	CMD=1d0/DEL*(R2/3d0 - 2d0/3d0*R3-5d0/3d0+
     $              (2d0*msq**2-2d0/3d0*mne**2)*LAM)
c       print *,'CI1:',DEL, R2/3d0,- 2d0/3d0*R3, 5d0/3d0,
c     $                (2d0*msq**2-2d0/3d0*mne**2)*LAM
	
==============================



      REAL*8 FUNCTION CMD(I,MSQ,MQ,MNE)
      IMPLICIT DOUBLE PRECISION (A-H,M,O-Z)
      REAL*8 LAM
      
      common/check/p1,p2

      DSQN=MSQ**2-MNE**2
      DQN =MQ**2 -MNE**2
      DSQQ=MSQ**2-MQ**2
      
      D3M=  msq**2+mq**2-MNE**2
      D3M1= msq**2-mq**2-MNE**2
      
      R1  =DSQQ/MNE**2
      R2  =DQN /MSQ**2
      R3  =DSQN/MQ**2
      
      DEL =2d0*MNE**2*(MQ**2+MSQ**2)-MNE**4-(MSQ**2-MQ**2)**2
      ADE =abs(DEL)
C$      IF(ADE.lt.1E-20)PRINT *,'DEL=0(!!!)'
      LAM=0d0
      IF(DEL.GE.0)LAM=2d0*ATAN(sqrt(ADE)/D3M)/sqrt(ADE)
      IF(DEL.LT.0)LAM=LOG((D3M+sqrt(ADE))/(D3M-sqrt(ADE)))/sqrt(ADE)

      CMD=0d0
      
      IF    (I.eq.1) THEN
        
	CMD=1d0/DEL*(R2/3d0 - 2d0/3d0*R3-5d0/3d0+
     $              (2d0*msq**2-2d0/3d0*mne**2)*LAM)
c       print *,'CI1:',DEL, R2/3d0,- 2d0/3d0*R3, 5d0/3d0,
c     $                (2d0*msq**2-2d0/3d0*mne**2)*LAM
	
      ELSEIF(I.eq.2) THEN
      
        CMD=(LOG(MSQ**2/MQ**2)-D3M1*LAM)/2d0/MNE**4+
     $ (((mq**4-mq**2*msq**2)/mne**2-7./3.*mq**2-2./3.*DSQN)*LAM+
     $ R2/3.+R1+2./3.)/DEL
      
      ELSEIF(I.eq.3) THEN
      
        CMD=-3d0/DEL**2*D3M+LAM/DEL*(-1d0+6d0*mq**2*msq**2/del)
	
      ELSEIF(I.eq.4) THEN
      
        CMD=
     $   ((2d0*log(msq/mq) - D3M1*LAM)/2d0/mne**2
     $   -1d0/msq**2
     $   -mq**2*d3m1/DEL*LAM)/mne**4
     $ +((mq**2/mne**4-(1d0-mq**2/mne**2)**2/msq**2+1d0/2d0/mne**2)
     $   +3d0*mq**2/DEL*
     $   (1d0 +  R1 + (-R1*mq**2-2d0*mq**2-msq**2+mne**2)*LAM))/DEL
     
c      p1=del*lam
           
c      p2=(-mq**2*d3m1)/mne**4

      
      ELSEIF(I.eq.5) THEN
    	
         CMD=1d0/2d0/mne**4*(log(msq**2/mq**2)-D3M1*LAM)-
     $      1d0/DEL*(LAM*(2d0*DSQN+3d0*MQ**2+R1*mq**2)-3d0-R1)
     
      ELSE
      
C$      print *,'wrong "I" IN CMD(I,...) !!!! check I:=',I
     
      ENDIF
      
      RETURN
      END
