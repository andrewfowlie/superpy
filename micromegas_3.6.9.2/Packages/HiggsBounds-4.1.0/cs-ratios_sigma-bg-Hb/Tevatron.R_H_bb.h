!******************************************************
        real*8 function tev_rH_bb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: K_ggH_TEV = 3.6
!* /!\: For the fit log_tev_rH_bb2 it is is important that the 
!*      coefficients are evaluated in double precision.
!*      Otherwise the result would be off by a factor larger than one!!
!*      It seems there are strong cancellations occuring in this
!*      expressions.
!******************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.3%
        implicit none
        real*8 x,log_tev_rH_bb,log_tev_rH_bb1,log_tev_rH_bb2,log_tev_rH_bb3
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1
        real*8 a2,b2,c2,d2,e2,f2,g2,h2,i2,j2,k2
        real*8 a3,b3,c3,d3,e3,f3,g3,h3,i3,j3,k3

        i1 = -40214.5223235366
        j1 = -1015.82748577826
        k1 = 182.869391507808
        a1 = -9.93339956971163
        b1 = 0.13148103336327
        c1 = -0.00183575723834311
        d1 = 1.45447917443553e-05
        e1 = -7.10199642418681e-08
        f1 = 2.08956850813578e-10
        g1 = -3.39191964925567e-13
        h1 = 2.32720334359468e-16

        i2 = -63262438241.0588
        j2 = -63146342615.3088
        k2 = 2604904414.22459
        a2 = -43169576.3385549
        b2 = 395121.031747372
        c2 = -2243.50437140996
        d2 = 8.27973195689632
        e2 = -0.0199851080043984
        f2 = 3.05478972125454e-05
        g2 = -2.69076223528998e-08
        h2 = 1.0427465954974e-11

        i3 = -1068990.36335925
        j3 = -116790.739045306
        k3 = 498392.529304714
        a3 = -3776261.14931063
        b3 = 70868.7892325959
        c3 = -570.267632826798
        d3 = 2.54878006194773
        e3 = -0.0068325177492403
        f3 = 1.09851398575053e-05
        g3 = -9.80790240576833e-09
        h3 = 3.75133090525309e-12
                                                                                                                
        log_tev_rH_bb1=i1/x**3+j1/x**2+k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        log_tev_rH_bb2=i2/x**3+j2/x**2+k2/x+a2+b2*x+c2*x**2+d2*x**3+e2*x**4+f2*x**5+g2*x**6+h2*x**7
        log_tev_rH_bb3=i3/x**3+j3/x**2+k3/x+a3+b3*x+c3*x**2+d3*x**3+e3*x**4+f3*x**5+g3*x**6+h3*x**7

        if(x .lt. 50d0) then
        write(*,*)'function tev_rH_bb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 290d0)) then
        log_tev_rH_bb=log_tev_rH_bb1
        endif
        if((x .gt. 290d0) .and. (x .le. 340d0)) then
        log_tev_rH_bb=log_tev_rH_bb2
        endif
        if((x .gt. 340d0) .and. (x .le. 400d0)) then
        log_tev_rH_bb=log_tev_rH_bb3
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rH_bb might not be a good fit (m_H > 400 GeV)'
        endif
        
        tev_rH_bb=exp(log_tev_rH_bb)
        
        end function

        
