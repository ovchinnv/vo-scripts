module output
 use size, only : me
 
 int, save :: fout=-1
 bool, save :: output_initialized=.false.
!
 int, private, parameter :: minerrorlev=0, minwarnlev=0
 int, private :: l
 character*9, private :: stat
   
 contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine error(whoami,msg,level)
  implicit none
  int :: level
  character*(*) :: whoami, msg
!  
  if (.not.output_initialized) call output_init()
  if (level.lt.minerrorlev) then
   stat=' FATAL'; l=6
  else
   stat=' NONFATAL'; l=9
  endif
!
  if (me.le.0) then
   write(fout,'(5A)') stat(1:l),' ERROR (',whoami,'): ', msg
  endif
!
  if (level.lt.minerrorlev) stop
!
 end subroutine error
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine warning(whoami,msg,level)
  implicit none
  int :: level
  character*(*) :: whoami, msg
!  
  if (.not.output_initialized) call output_init()
  if (level.lt.minerrorlev) then
   stat=' FATAL'; l=6
  else
   stat=' NONFATAL'; l=9
  endif
!
  if (me.le.0) then
   write(fout,'(5A)') stat(1:l),' WARNING (',whoami,'): ', msg
  endif
!
  if (level.lt.minwarnlev) stop
!
 end subroutine warning
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 subroutine message(whoami,msg)
  implicit none
  character*(*) :: whoami, msg
!  
  if (.not.output_initialized) call output_init()
!
  if (me.le.0) then
   write(fout,'(4A)') ' MESSAGE (',whoami,'): ', msg
  endif
!
 end subroutine message
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
 subroutine output_init(filename)
 implicit none
 character*(*), optional :: filename
 character*400 :: fname
 int :: flen
 
 call output_done()
!
 if (present(filename)) then
  fname=adjustl(filename)
  flen=len_trim(fname)
  fout=987
  open(unit=fout, file=fname(1:flen), form='FORMATTED', status='UNKNOWN')
 else
  fout=6
 endif
! 
 output_initialized=.true. 
 end subroutine output_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
 subroutine output_done()
 implicit none
 if (output_initialized) then
  if (fout.ne.0.and.fout.ne.6.and.fout.ne.5) close(fout)  
  output_initialized=.false.
  fout=-1
 endif
 end subroutine output_done
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module output
