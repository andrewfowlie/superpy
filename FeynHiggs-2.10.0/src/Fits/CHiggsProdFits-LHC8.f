
	if( sqrtm .lt. 14.142135623730951D0 )
     &    Warning('Extrapolating tHm2 in MHiggs')
	if( TBeff .lt. 30.D0 )
     &    Warning('Extrapolating tHm2 in TBeff')

        tHm2 = exp(7.021818834219673D0 - 
     &     sqrtm*(0.24264801383854634D0 + 
     &        0.0024692174421524853D0*sqrtm + 
     &        2.4536969668640866D-9*TBeff) + 
     &     TBeff*(0.08830669851003883D0 - 
     &        0.00047277613771311733D0*TBeff))

#ifdef DETAILED_DEBUG
	DPROD "tHm2 =", tHm2 ENDL
#endif


	if( sqrtm .lt. 14.142135623730951D0 )
     &    Warning('Extrapolating tHm2lo in MHiggs')
	if( TBeff .lt. 30.D0 )
     &    Warning('Extrapolating tHm2lo in TBeff')

        tHm2lo = exp(6.706421238043519D0 - 
     &     sqrtm*(0.24470193778367805D0 + 
     &        0.0025272634073900285D0*sqrtm + 
     &        1.7694439279323158D-9*TBeff) + 
     &     TBeff*(0.08830668570333006D0 - 
     &        0.00047277612263374864D0*TBeff))

#ifdef DETAILED_DEBUG
	DPROD "tHm2lo =", tHm2lo ENDL
#endif


	if( sqrtm .lt. 14.142135623730951D0 )
     &    Warning('Extrapolating tHm2hi in MHiggs')
	if( TBeff .lt. 30.D0 )
     &    Warning('Extrapolating tHm2hi in TBeff')

        tHm2hi = exp(7.281658836714188D0 - 
     &     sqrtm*(0.24156900820289887D0 + 
     &        0.0024717298679832666D0*sqrtm + 
     &        1.3392486950262758D-9*TBeff) + 
     &     TBeff*(0.08830670751479178D0 - 
     &        0.0004727764088130346D0*TBeff))

#ifdef DETAILED_DEBUG
	DPROD "tHm2hi =", tHm2hi ENDL
#endif

