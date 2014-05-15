!******************************************************
        real*8 function lhc7_rHVBF_WW(x)
!******************************************************
!* Used HAWK with settings of input parameters of arXiv:1101.0593 [hep-ph].
!* (PDF: MSTW 2008 NNLO)
!* Calculated Born CS with only t-channel graphs (*),
!* (a) with only WW-fusion graphs and 
!* (b) with all WW/ZZ fusion graphs, and calculated 
!* ratio (a)/(b).
!* (*): s-channel graphs are neglected because they are 
!*      highly suppressed after VBF cuts.
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:1000], deviations from data below 0.3%
!*      slight extrapolation to range [70:1020] is allowed.
!* 27/4/2011, Oliver Brein (revised: 29/4/2011)
!******************************
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x
        real*8 a0,b0,c0,d0

        a0 = 0.764629881486091
        b0 = -0.000128901014699432
        c0 = 2.248155381437e-07
        d0 = -2.38464077139572e-10

        i = 353411.391484508
        j = -13532.7423394166
        k = 216.248242966815
        a = -1.13447368036044
        b = 0.0100380489651014
        c = -3.45357638530139e-05
        d = 7.69000469151772e-08
        e = -1.10305464638084e-10
        f = 9.79724660472145e-14
        g = -4.88552475279702e-17
        h = 1.04278241770435e-20

        lhc7_rHVBF_WW=0d0
                                                                                                                                                                                
        if(x .lt. 70d0) then
        write(*,*)'function lhc7_rHVBF_WW might not be a good fit (m_H < 70 GeV)'
        endif
        if((x .ge. 70d0) .and. (x .lt. 200d0)) then
        lhc7_rHVBF_WW=a0+b0*x+c0*x**2+d0*x**3
        endif
        if((x .ge. 200d0) .and. (x .le. 1020d0)) then
        lhc7_rHVBF_WW=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
        endif
        if(x .gt. 1020d0) then
        write(*,*)'function lhc7_rHVBF_WW might not be a good fit (m_H > 1020 GeV)'
        endif


                                                                                                                                        
        end function
        
