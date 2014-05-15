!******************************************************
        real*8 function tev_rHZ_ddb(x)
!******************************************************
!* PDF used: MRST NNLO 2006
!* K-factors used: drop out
!******************************************************
!c x : Higgs mass
!c fit: valid in range [50:400], deviations from data below 0.4%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x

        i = 70218.5986224673
        j = -5029.65551079772
        k = 149.029291563756
        a = -2.03481900046978
        b = 0.0217025637309646
        c = -0.000132490305931743
        d = 4.82474402412335e-07
        e = -9.96206232060581e-10
        f = 9.22849382455714e-13
        g = 1.20550194299852e-16
        h = -6.04287787982214e-19
!        FIT_CONVERGED = 1
!        FIT_NDF = 60
!        FIT_STDFIT = 0.000224135498729174
!        FIT_WSSR = 3.01420330743452e-06

        if(x .lt. 50d0) then
        write(*,*)'function tev_rHZ_ddb might not be a good fit (m_H < 50 GeV)'
        endif
        if((x .ge. 50d0) .and. (x .le. 400d0)) then
        tev_rHZ_ddb=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        endif
        if(x .gt. 400d0) then
        write(*,*)'function tev_rHZ_ddb might not be a good fit (m_H > 400 GeV)'
        endif     
                                                                                                                     
        end function

