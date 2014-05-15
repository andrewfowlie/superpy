		PROGRAM AllAnalyses
		 use S95tables_type1  
       	 use S95tables 
		 use theory_BRfunctions, only : setup_BRSM,deallocate_BRSM
 		 use channels, only : setup_channels,check_channels
         use store_pathname
			 
		   character(len=100),allocatable :: filename(:)
		  character(LEN=pathname_length+150) :: fullfilename 
		 integer analysis, n, k, x, col
		 integer n_datafiles
		 character(LEN=80) c
		 character(LEN=80) :: datafile(500)
		 
		 type(table1),allocatable :: S95_tnew(:)

		 call setup_BRSM 

		 call setup_S95tables      
		PRINT *, ""
		PRINT *, " LISTING ALL HIGGSBOUNDS ANALYSES"
		PRINT *, ""
		PRINT *, " Tables of Type I"
		PRINT *, ""
		PRINT *, "  #        ID  expt        E       lumi  part  label                                " &
		// "            xmin   xmax     step  SM"
		PRINT *, " -----------------------------------------------------------------------------------" &
		// "------------------------------------"
		 DO I=lbound(S95_t1,dim=1), ubound(S95_t1,dim=1)
           WRITE(*,'(I5,I10,A,A,A,A,F6.2,A,F8.3,I4,A,A,3F8.2,I3)') I, S95_t1(I)%id," ", adjustl(S95_t1(I)%expt),  "  ",&
           "   ", &
           S95_t1(I)%energy,"   ",  &
           S95_t1(I)%lumi,S95_t1(I)%particle_x,"   ", ADJUSTL(S95_t1(I)%label), &
           S95_t1(I)%xmin,S95_t1(I)%xmax,S95_t1(I)%sep,S95_t1(I)%SMlike
		 ENDDO
		PRINT *, " -----------------------------------------------------------------------------------" &
		// "------------------------------------"


		PRINT *, ""
		PRINT *, " Tables of Type II"
		PRINT *, " ------------------------------------------------------------------"
	     DO I=lbound(S95_t2,dim=1), ubound(S95_t2,dim=1)
		 	WRITE(*,'(I5,I10,A,A,A,A)') I,  S95_t2(I)%id," ", adjustl(S95_t2(I)%expt), "   ", ADJUSTL(S95_t2(I)%label)
		 ENDDO
		 PRINT *, " ------------------------------------------------------------------"

		PRINT *, ""
		
	
		 END		 
		
