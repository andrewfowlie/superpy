
	if( sqrtm .lt. 10.D0 )
     &    Warning('Extrapolating btagbhSM('//Digit(h)//') in MHiggs')

        btagbhSM(h) = exp(11.656268920612366D0 + 
     &     sqrtm*(-0.6201183772901906D0 + 0.005467860164016601D0*sqrtm)
     &     )

#ifdef DETAILED_DEBUG
	DPROD "btagbhSM(h) =", btagbhSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating tthSM('//Digit(h)//') in MHiggs')

        tthSM(h) = exp(7.668104619282241D0 + 
     &     sqrtm*(0.9233583568187562D0 + 
     &        sqrtm*(-0.1657563305726811D0 + 
     &           sqrtm*(0.007941294236898499D0 - 
     &              0.00012509592650744857D0*sqrtm))))

#ifdef DETAILED_DEBUG
	DPROD "tthSM(h) =", tthSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating qqhSM('//Digit(h)//') in MHiggs')

        qqhSM(h) = exp(9.477144600773386D0 + 
     &     sqrtm*(0.01917853525435689D0 + 
     &        sqrtm*(-0.01586840079453014D0 + 
     &           sqrtm*(0.0005156314311776905D0 - 
     &              5.642390955525177D-6*sqrtm))))

#ifdef DETAILED_DEBUG
	DPROD "qqhSM(h) =", qqhSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating WhSM('//Digit(h)//') in MHiggs')

        WhSM(h) = exp(8.071713542712953D0 + 
     &     sqrtm*(1.4254134929944258D0 + 
     &        sqrtm*(-0.25289140687250683D0 + 
     &           sqrtm*(0.013523476235153417D0 - 
     &              0.00025504685613550714D0*sqrtm))))

#ifdef DETAILED_DEBUG
	DPROD "WhSM(h) =", WhSM(h) ENDL
#endif


	if( sqrtm .lt. 9.486832980505138D0 )
     &    Warning('Extrapolating ZhSM('//Digit(h)//') in MHiggs')

        ZhSM(h) = exp(13.52369429479672D0 + 
     &     sqrtm*(-0.6710557531650944D0 + 0.006119295282304429D0*sqrtm)
     &     )

#ifdef DETAILED_DEBUG
	DPROD "ZhSM(h) =", ZhSM(h) ENDL
#endif


	if( sqrtm .lt. 10.488088481701515D0 )
     &    Warning('Extrapolating StSth('//Digit(h)//') in MHiggs')
	if( mst1 .lt. 150.D0 )
     &    Warning('Extrapolating StSth('//Digit(h)//') in MStop1')

        StSth(h) = exp((7.393265441716565D7 + 
     &       mst1*(-334902.238690035D0 - 1152.3708860300499D0*mst1 - 
     &          12194.806996326963D0*sqrtm) + 
     &       sqrtm*(-8.326881834682285D6 + 
     &          126542.25115702617D0*sqrtm))/
     &     (3.9518151195296748D6 + 
     &       mst1*(55658.75817701751D0 + 34.18404384496346D0*mst1 - 
     &          167.12163223629233D0*sqrtm) + 
     &       sqrtm*(361192.31847366353D0 - 16169.109899792384D0*sqrtm))
     &     )

#ifdef DETAILED_DEBUG
	DPROD "StSth(h) =", StSth(h) ENDL
#endif

