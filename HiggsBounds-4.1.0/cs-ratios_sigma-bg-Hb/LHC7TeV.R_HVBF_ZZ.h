!******************************************************
        real*8 function lhc7_rHVBF_ZZ(x)
!******************************************************
!* Used HAWK with settings of input parameters of arXiv:1101.0593 [hep-ph].
!* (PDF: MSTW 2008 NNLO)
!* Calculated Born CS with only t-channel graphs (*),
!* (a) with only ZZ-fusion graphs and 
!* (b) with all WW/ZZ fusion graphs, and calculated 
!* ratio (a)/(b).
!* (*): s-channel graphs are neglected because they are 
!*      highly suppressed after VBF cuts.
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:1000], deviations from data below 0.2%
!*      slight extrapolation to range [70:1020] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x
        real*8 a0,b0,c0,d0

        a0 = 0.237215317453696
        b0 = 0.000102401893913729
        c0 = -8.00259372973129e-08
        d0 = -1.1874470382624e-11

        i = 5170.13979049644
        j = 169.605166589112
        k = -7.27222629415314
        a = 0.32245038001123
        b = -0.000343674727914338
        c = 9.99626245300363e-07
        d = -7.7551434795945e-10
        e = -1.73547041303278e-12
        f = 4.41364218643659e-15
        g = -3.67002024110552e-18
        h = 1.08578419805009e-21

        lhc7_rHVBF_ZZ=0d0

        if(x .lt. 70d0) then
        write(*,*)'function lhc7_rHVBF_ZZ might not be a good fit (m_H < 70 GeV)'
        endif
        if((x .ge. 70d0) .and. (x .lt. 90d0)) then
        lhc7_rHVBF_ZZ=a0+b0*x+c0*x**2+d0*x**3
        endif
        if((x .ge. 90d0) .and. (x .le. 1020d0)) then
        lhc7_rHVBF_ZZ=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        endif
        if(x .gt. 1020d0) then
        write(*,*)'function lhc7_rHVBF_ZZ might not be a good fit (m_H > 1020 GeV)'
        endif
                                                                                                                                        
        end function
        
