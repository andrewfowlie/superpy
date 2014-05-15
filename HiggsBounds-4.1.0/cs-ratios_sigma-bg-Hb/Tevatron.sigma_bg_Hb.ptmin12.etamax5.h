!******************************************************
        real*8 function tev_cs_bg_Hb_c4_SM(x)
!******************************************************
!* Fit for sigma(pp-bar -> b g -> H b) in pikobarn as a function of Higgs mass in GeV
!* However, sigma(pp-bar -> b g -> H b)+sigma(pp-bar -> b-bar g -> H b-bar) is 
!* needed here. Hence, the factor of 2 below.
!*
!* PDF used: MSTW NNLO 2008
!* K-factors used: none. 
!* However, used: a) scale choice MU_R=MU_F=M_h/4
!*                b) 1-loop running bottom mass instead of pole mass
!* This resembles pretty closely the NLO result of Maltoni et al. (hep-ph/0204093)
!* cuts used: ptmin=12 GeV, etamax=5 (c4)
!* /!\: using ptmin as low as 12 GeV may render the prediction less accurate
!******************************************************
!c x : Higgs mass
!c fit: valid in range [10:400], deviations from data below 0.6%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x,log_tev_cs_bg_Hb_c4

        i = 858.031227860774d0
        j = -258.119990404928d0
        k = 27.1187184643582d0
        a = -0.228082879265122d0
        b = -0.0758457063398605d0
        c = 0.000321567763038373d0
        d = -1.28719909200427d-06
        e = 3.44247069068616d-09
        f = -5.5337665142955d-12
        g = 4.58899173696998d-15
        h = -1.3449829067125d-18

       log_tev_cs_bg_Hb_c4=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
       tev_cs_bg_Hb_c4_SM=2d0*exp(log_tev_cs_bg_Hb_c4)

       end function
        
