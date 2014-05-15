!******************************************************
        real*8 function lhc8_rH_bb(x)
!******************************************************
!* ggH process : consensus numbers used (PDF: MSTW 2008 NNLO)
!*		 see: arXiv:1101.0593 [hep-ph] 
!* 		 numbers taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt7TeV
!* bbH process : numbers in range [70:1020] generated with bbh@nnlo v. 1.3 (PDF: MSTW 2008 NNLO) 
!* 		 with input parameters according to arXiv:1101.0593 [hep-ph]
!******************************************************
!* x : Higgs mass in GeV
!* fit: strictly valid in range [90:1000], deviations from data below 0.9%
!*      slight extrapolation to range [70:1020] is allowed.
!* 5/9/2012, Oliver Brein (implemented by OS)
!******************************
        implicit none
        real*8 x,log_lhc8_rH_bb
        real*8 a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1
        real*8 a2,b2,d2,f2,h2,i2,k2
        real*8 a3,b3,d3,f3,h3,i3,k3
        real*8 a4,b4,c4,d4,e4,f4,g4,h4,i4,j4,k4

        i1 = 50289348.6325164D0
        j1 = -5457532.5725116D0
        k1 = 229808.827286398D0
        a1 = -5224.3873291114D0
        b1 = 72.5564132006601D0
        c1 = -0.6525299930268D0
        d1 = 0.00387236419984D0
        e1 = -1.5051035421818D-5
        f1 = 3.68196043426451D-8
        g1 = -5.1376601793327D-11
        h1 = 3.11529787208313D-14
        
        
        i2 = 12182.6157148462D0
        k2 = -8767.3919801225D0
        a2 = 106.303172943418D0
        b2 = -0.4182926536687D0
        d2 = 2.3418996986372D-06
        f2 = -1.0432598252495D-11
        h2 = 1.9711083344475D-17
        
        i3 = -3281.81931165533D0
        k3 = 7488.3694050317D0
        a3 = -51.390894803456D0
        b3 = 0.0766433197729768D0
        d3 = -9.62237115264618D-08
        f3 = 8.81655982140704D-14
        h3 = -2.81620641485924D-20
        
        i4 = -9807081985.35512D0
        j4 = 166792119.488553D0
        k4 = -1229001.04292627D0
        a4 = 5222.71458350981D0
        b4 = -14.3131068007305D0
        c4 = 0.0263893978685678D0
        d4 = -3.32580710408263D-05
        e4 = 2.83322845999699D-08
        f4 = -1.56343097264634D-11
        g4 = 5.05161631215867D-15
        h4 = -7.26365833097573D-19



        log_lhc8_rH_bb=0d0

        if(x .lt. 70d0) then
        write(*,*)'function lhc8_rH_bb might not be a good fit (m_H < 70 GeV)'
        endif
        if((x .ge. 70d0) .and. (x .le. 250d0)) then
        log_lhc8_rH_bb=i1/x**3+j1/x**2+k1/x+a1+b1*x+c1*x**2+d1*x**3+e1*x**4+f1*x**5+g1*x**6+h1*x**7
        endif
        if((x .gt. 250d0) .and. (x .le. 340d0)) then
        log_lhc8_rH_bb=i2/x**3+k2/x+a2+b2*x+d2*x**3+f2*x**5+h2*x**7

        endif
        if((x .gt. 340d0) .and. (x .le. 400d0)) then
        log_lhc8_rH_bb=i3/x**3+k3/x+a3+b3*x+d3*x**3+f3*x**5+h3*x**7
        endif
        if((x .gt. 400d0) .and. (x .le. 1020d0)) then
        log_lhc8_rH_bb=i4/x**3+j4/x**2+k4/x+a4+b4*x+c4*x**2+d4*x**3+e4*x**4+f4*x**5+g4*x**6+h4*x**7
        endif
        if(x .gt. 1020d0) then
        write(*,*)'function lhc8_rH_bb might not be a good fit (m_H > 1020 GeV)'
        endif
        
        lhc8_rH_bb=exp(log_lhc8_rH_bb)
        
        end function
        
