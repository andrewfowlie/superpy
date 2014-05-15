!******************************************************
        real*8 function lhc7_rHZ_bbb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.2%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,x,log_lhc7_rHZ_bbb

        a1 = -3.92833920968973
        b1 = -0.00168025141081915
        c1 = -3.48443687089216e-05
        d1 = 1.52750406535549e-07
        e1 = -1.94872842022944e-10

        lhc7_rHZ_bbb=0d0
        
        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHZ_bbb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc7_rHZ_bbb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc7_rHZ_bbb=exp(log_lhc7_rHZ_bbb)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHZ_bbb might not be a good fit (m_H > 310 GeV)'
        endif
                                                                                                                                        
        end function
        
