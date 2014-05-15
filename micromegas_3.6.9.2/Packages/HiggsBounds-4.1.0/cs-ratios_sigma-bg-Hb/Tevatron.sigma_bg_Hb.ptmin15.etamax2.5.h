!******************************************************
        real*8 function SMCS_tev_bg_Hb_c3(x)
!******************************************************
!* Fit for sigma(pp-bar -> b g -> H b) in pikobarn as a function of Higgs mass in GeV
!* However, sigma(pp-bar -> b g -> H b)+sigma(pp-bar -> b-bar g -> H b-bar) is 
!* needed here. Hence, the factor of 2 below.
!*
!* PDF used: MRST NNLO 2006
!* K-factors used: none. 
!* However, used: a) scale choice MU_R=MU_F=M_h/4
!*                b) 1-loop running bottom mass instead of pole mass
!* This resembles pretty closely the NLO result of Maltoni et al. (hep-ph/0204093)
!* cuts used: ptmin=15 GeV, etamax=2.5 (c3)
!******************************************************
!c x : Higgs mass
!c fit: valid in range [10:400], deviations from data below 0.6%
        implicit none
        real*8 a,b,c,d,e,f,g,h,i,j,k,x,log_tev_cs_bg_Hb_c3

        i = 975.922569925692
        j = -262.219528073721
        k = 23.4358421199538
        a = -0.827082997657741
        b = -0.0695990730309521
        c = 0.000294055015901147
        d = -1.22958624351404e-06
        e = 3.63264089688587e-09
        f = -6.84539813558712e-12
        g = 7.28039954882708e-15
        h = -3.29940845455953e-18
!        FIT_CONVERGED = 1
!        FIT_NDF = 68
!        FIT_STDFIT = 0.00191372574757335
!        FIT_WSSR = 0.000249039544110912

       log_tev_cs_bg_Hb_c3=i/x**3+j/x**2+k/x+a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7
       SMCS_tev_bg_Hb_c3=2d0*exp(log_tev_cs_bg_Hb_c3)

       end function
        