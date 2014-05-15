!******************************************************
        real*8 function lhc8_rHZ_ddb(x)
!******************************************************
!* Used LO calculation and NNLO QCD and EW K-factors.
!* This reproduces the numbers of arXiv:1101.0593 [hep-ph] 
!* to better than 1%. (PDF: MSTW 2008 NNLO)
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:300], deviations from data below 0.2%
!*      slight extrapolation to range [80:310] is allowed.
!* 24/9/2012, Oliver Brein / Oscar Stål
!******************************
        implicit none
        real*8 a,b,c,d,e,x

        a = 0.37158485156533D0
        b = 2.74861371238629D-05
        c = 2.15319965594886D-06
        d = -1.98315785119313D-09
        e = -8.13856600329891D-11
		
        lhc8_rHZ_ddb=0d0

        if(x .lt. 80d0) then
        write(*,*)'function lhc8_rHZ_ddb might not be a good fit (m_H < 80 GeV)'
        endif
        if((x .ge. 80d0) .and. (x .le. 310d0)) then
        lhc8_rHZ_ddb=a+b*(x-200d0)+c*(x-200d0)**2+d*(x-200d0)**3+e*(x-200d0)**4
        endif
        if(x .gt. 310d0) then
        write(*,*)'function lhc8_rHZ_ddb might not be a good fit (m_H > 310 GeV)'
        endif

                                                                                                                                        
        end function
        
