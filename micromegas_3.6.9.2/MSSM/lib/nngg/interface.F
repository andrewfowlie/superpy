      real*8 function  vcsnngg(v)
      IMPLICIT NONE
      real*8 v
      real*8 vcs1,res(2),coss

      data coss/0/

      call nnaaini(1,1,1,1,v)
      call differentialcs(res,coss)
      vcs1=res(2)
      call nnaaini(1,1,-1,-1,v)
      call differentialcs(res,coss)
      vcs1=vcs1+res(2)
      vcsnngg=vcs1*v
      end


