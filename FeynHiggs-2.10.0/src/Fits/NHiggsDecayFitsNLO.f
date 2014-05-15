
	if( m .lt. 110.D0 )
     &    Warning('Extrapolating GammaSM(H0VV('//Digit(h)//',3)) in MHiggs')

	if( m.lt.161.D0 ) then

        GammaSM(H0VV(h,3)) = 
     &  exp(-65.10019525061139D0 + 
     &    m*(0.949915980247105D0 + 
     &       m*(-0.005454555928428384D0 + 0.000011493488849045667D0*m))
     &    )

	else

        GammaSM(H0VV(h,3)) = 
     &  exp((1.58417879062567D11 + 
     &      m*(-1.1891121875847956D10 + 
     &         m*(1.5493674059984714D8 + 
     &            m*(-733118.4403864534D0 + 1177.118512098065D0*m))))/
     &    (-3.549562407540174D11 + 
     &      m*(6.26673912938523D9 + 
     &         m*(-3.3482554942270122D7 + 
     &            m*(35157.54099115199D0 + 102.241073227734D0*m)))))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'GammaSM(H0VV(h,3)) =', GammaSM(H0VV(h,3)) ENDL
#endif


	if( m .lt. 110.D0 )
     &    Warning('Extrapolating GammaSM(H0VV('//Digit(h)//',4)) in MHiggs')

	if( m.lt.161.D0 ) then

        GammaSM(H0VV(h,4)) = 
     &  exp(-112.75895906661255D0 + 
     &    m*(2.209012186197331D0 + 
     &       m*(-0.015885711991485498D0 + 0.00003984296368949889D0*m)))

	else

        GammaSM(H0VV(h,4)) = 
     &  exp((1.21831447488354D9 + 
     &      m*(-2.0732520918997794D7 + 
     &         m*(128474.34303947467D0 + 
     &            m*(-357.8081500249324D0 + 0.40839544958693924D0*m))))
     &     /(-1.4885590393233255D8 + 
     &      m*(964350.3775372116D0 + 
     &         m*(2390.3346143235276D0 + 
     &            m*(-23.298635407525836D0 + 0.04837647108675248D0*m)))
     &      ))

	endif

#ifdef DETAILED_DEBUG
	DPROD 'GammaSM(H0VV(h,4)) =', GammaSM(H0VV(h,4)) ENDL
#endif

