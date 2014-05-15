!******************************************************
        real*8 function lhc7_rHZ_ssb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.25%
!*      slight extrapolation to range [80:310] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a1,b1,c1,d1,e1,x,log_lhc7_rHZ_ssb

        a1 = -2.3749351144681
        b1 = -0.000699294911440839
        c1 = -3.56629276999138e-05
        d1 = 1.52410935751157e-07
        e1 = -2.00074252300485e-10
        
        lhc7_rHZ_ssb=0d0
        
        if(x .lt. 80d0) then
        write(*,*)'function lhc7_rHZ_ssb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        log_lhc7_rHZ_ssb=a1+b1*x+c1*x**2+d1*x**3+e1*x**4
        lhc7_rHZ_ssb=exp(log_lhc7_rHZ_ssb)
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc7_rHZ_ssb might not be a good fit (m_H > 310 GeV)'
        endif


                                                                                                                        
                                                                                                                                        
        end function
        
