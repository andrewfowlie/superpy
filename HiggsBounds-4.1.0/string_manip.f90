module string_manip
  implicit none 

contains

 !********************************************************
 function strtolcase(str,n)
 !********************************************************
    implicit none
    integer :: strcapA = 65    ! from WRITE(*,*) "A:",ICHAR("A")
    integer :: strcapZ = 90    ! etc
    !integer :: strlowa = 97
    !integer :: strlowz = 122
    integer :: casediff = 32   ! = 97-65  = 122-90
    character(len=n) :: strtolcase
    character(len=n) :: str
    integer :: i,n,ich
    !-----------------------

    strtolcase = str

    do i=1,n
       ich = ichar(strtolcase(i:i))
       if( ich>=strcapA .and. ich<=strcapZ ) strtolcase(i:i) = char(ich+casediff)
    enddo

 end function strtolcase

 !********************************************************
 function remove_from_character_onwards(string_in,n,charact)
 !********************************************************
  implicit none
  character(len=n),intent(in)::string_in
  character(len=n) ::remove_from_character_onwards
  integer,intent(in) :: n
  integer :: position_of_charact
  character(len=1) :: charact

  position_of_charact=index(string_in, charact, .false.)!.false. means returns the position of the first substring from the left
  !write(*,*)'start_of_comment',start_of_comment

  if(position_of_charact.eq.0)then
    remove_from_character_onwards=string_in
  else
    remove_from_character_onwards=string_in(:position_of_charact-1)
  endif
 end function remove_from_character_onwards

 !********************************************************
 function remove_from_character_backwards(string_in,n,charact)
 !********************************************************
  implicit none
  character(len=n),intent(in)::string_in
  character(len=n) ::remove_from_character_backwards
  integer,intent(in) :: n
  integer :: position_of_charact
  character(len=1) :: charact

  position_of_charact=index(string_in, charact, .false.)!.false. means returns the position of the first substring from the right
  !write(*,*)'position_of_charact',position_of_charact

  if(position_of_charact.eq.0)then
    remove_from_character_backwards=string_in
  else
    remove_from_character_backwards=string_in(position_of_charact+1:)
  endif
 end function remove_from_character_backwards


 !********************************************************
 function strip_off_comment(string_in,n)
 !********************************************************
  implicit none
  character(len=n),intent(in)::string_in
  character(len=n) ::strip_off_comment
  integer,intent(in) :: n

  strip_off_comment=remove_from_character_onwards(string_in,n,'#')

 end function strip_off_comment


 !********************************************************
 subroutine split_into_col(str,len_str,columns,len_columns)
 !reads the first n columns where n=ubound(columns,dim=1)
 !if there are m columns and m<n, columns(m+1)...columns(n) will be empty 
 !********************************************************
  implicit none
  character(len=len_str),intent(in) :: str
  character(len=len_str) :: temp_line,bit_of_line
  character(len=len_columns) :: columns(:)
  integer,intent(in) :: len_str,len_columns
  integer :: n,i,position_of_first_space

  columns=''
  n=size(columns,dim=1) 
 
  temp_line=str

  do i=1,n
 
    temp_line=adjustl(temp_line)

    position_of_first_space=index(temp_line, ' ', .false.)

    if(position_of_first_space.gt.1)then
      bit_of_line= temp_line(:position_of_first_space-1)

      !if(len(trim(bit_of_line)).gt.len_columns)then
      ! stop'string columns is not big enough'
      !endif

      columns(i)=bit_of_line
  
      temp_line = temp_line(position_of_first_space:)
    endif
  enddo

  !if(trim(temp_line).ne.'')then
  !  stop'there is still something left in line'
  !endif

 end subroutine split_into_col
 !********************************************************
 subroutine saferead_int(str,integ)
 !********************************************************
  implicit none

  character(len=*),intent(in) :: str
  integer :: integ,integ_initial
  integer :: ios

  integ_initial=integ
  read(str,*,iostat=ios)integ
  if(ios.ne.0)integ=integ_initial
  
 end subroutine saferead_int

 !********************************************************
 subroutine saferead_dble(str,dbleprec)
 !********************************************************
  implicit none

  character(len=*),intent(in) :: str
  double precision :: dbleprec,dbleprec_initial
  integer :: ios

  dbleprec_initial=dbleprec
  read(str,*,iostat=ios)dbleprec
  if(ios.ne.0)dbleprec=dbleprec_initial

 end subroutine saferead_dble
 !********************************************************
end module string_manip
