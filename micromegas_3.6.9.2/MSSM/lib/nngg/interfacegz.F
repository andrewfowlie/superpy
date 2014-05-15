      real*8 function vcsnngz(v)
      implicit none 
      real*8 v,vcs,res(2),coss
      integer ok,nnazini
      coss=0
      vcsnngz=0
      ok=nnazini(1,1,1,1,v)
      if(ok.eq.0) return
      call gzdifferentialcs(res,coss)
      vcs=res(2)
      ok=nnazini(1,1,-1,-1,v)
      call gzdifferentialcs(res,coss)
      vcs=vcs+res(2)
      vcsnngz=vcs*v 
      end
