
 implicit none
 integer :: i
 character*26, parameter  :: lower='abcdefghijklmnopqrstuvwxyz'
 character*26, parameter  :: upper='ABCDEFGHIJKLNMOPQRSTUVWXYZ'
 character*10, parameter  :: digits='0123456789'
 character, parameter :: lower2(26)=(/ (lower(i:i),i=1,26)/)
 character, parameter :: upper2(26)=(/ (upper(i:i),i=1,26)/)
 character, parameter :: digits2(10)=(/ (digits(i:i),i=1,10)/)
 character, parameter :: tilde='~'
 character, parameter :: decimal='.'
 character, parameter :: underscore='_'
 character, parameter :: hyphen='-'
 character, parameter :: slash='/'
 character, parameter :: space(2) = (/' ', '\t'/)
 character, parameter :: comment(2) = (/'*','!'/)
 character, parameter :: equals(1) = (/'='/)
 character*200, parameter :: allowed=upper//lower//digits//decimal//underscore//hyphen//slash//tilde
!
 integer :: allowed_len, lower_len, upper_len, digits_len
 character :: allowed2(200),c
 logical :: qtag, qeq, qval, qerror
!
 integer :: ioerr
 integer :: l=0, ltag=0, lval=0
 character*400 :: cmdline
 character*200 :: tag, val

 allowed_len=len_trim(allowed)
 do i=1,allowed_len
  allowed2(i)=allowed(i:i)
 enddo 


 open(1, file='0readme',form='formatted', IOSTAT=ioerr)
!
 qerror=.false.
 do while (.true.)
  read(1,'(A)',IOSTAT=ioerr) cmdline
  if (ioerr.eq.0) then
!   write(0,*) cmdline

! `read' from the line
   cmdline=adjustl(cmdline)
   l=len_trim(cmdline)
!   write(0,*) cmdline(1:l)
! add comment character to know when to stop below
   if (l.lt.200) l=l+1
   cmdline(l:l)='*'
   if (any(comment.eq.cmdline(1:1))) l=1 ! skip lines that are comments

   i=0
   qtag=.true. ; ltag=0 ! each line is required to begin with a tag
   qval=.false.; lval=0
   qeq=.false.
!
   do while (l.gt.1)
     i=i+1
     
     c=cmdline(i:i)

!     write(0,*) tag(1:ltag), val(1:lval),i

     if (any(comment.eq.c.or.i.eq.l)) then ! end of command line
      if ((qtag.and.i.gt.1).or.qeq) then
       write(0,*) 'ERROR: UNEXPECTED END OF LINE'
       qerror=.true.; exit 
      elseif (qval) then
       if (lval.eq.0) then
        write(0,*) 'ERROR: MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'"; SKIPPING LINE.'
        qerror=.true.; exit 
       else ! store tag and val
        write(0,*) tag(1:ltag),'==>',val(1:lval)
        ltag=0; lval=0
        qerror=.true.; exit
       endif
      else
       exit
      endif ! qtag
!
     elseif (any(space.eq.c)) then ! completed a tag or value
!
      if (qtag) then
       qtag=.false.
       qeq=.true.
      elseif (qeq) then
!      nothing
      elseif (qval) then
       qval=.false. 
       qtag=.true.
       if (lval.eq.0) then
        write(0,*) 'ERROR: MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'"; SKIPPING LINE.'
        qerror=.true.; exit 
       else ! store tag and val
        write(0,*) tag(1:ltag),'==>',val(1:lval)
        ltag=0; lval=0
        cmdline=adjustl(cmdline(i:l)) ! remove tag/val pair from string
        l=len_trim(cmdline)
        i=0
       endif ! lval
       ltag=0
      endif ! qtag
!
     elseif (any(equals.eq.c)) then ! completed a tag or value
!
      if (qtag) then
       qtag=.false.
       if (ltag.eq.0) then
        write(0,*) 'ERROR: MISSING PARAMETER NAME. SKIPPING LINE '
        qerror=.true.; exit 
       endif
!
       qval=.true. ; lval=0
       cmdline=adjustl(cmdline(i+1:l)) ! remove 'tag=' and leading spaces pair from string
       l=len_trim(cmdline)
       i=0
!
       qeq=.false.
      elseif (qeq) then
       qeq=.false.
       qval=.true. ; lval=0
       cmdline=adjustl(cmdline(i+1:l)) ! remove 'tag=' and leading spaces pair from string
       l=len_trim(cmdline)
       i=0
      elseif (qval) then
       if (lval.eq.0) then
        write(0,*) 'ERROR: MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'". SKIPPING LINE.'
        qerror=.true.; exit 
       else ! store tag and val
        write(0,*) tag(1:ltag),'==>',val(1:lval)
        ltag=0; lval=0
        cmdline=adjustl(cmdline(i:l)) ! remove tag/val pair from string
        l=len_trim(cmdline)
        i=0
       endif ! lval
       qtag=.true. ; ltag=0
       qval=.false.
      endif ! qtag       
!
     elseif (any(allowed2(1:allowed_len).eq.c)) then ! check that the characters are allowed
!
      if (qtag) then
       if (all(lower2(1:26).ne.c).and.all(upper2(1:26).ne.c)) then
        if (ltag.eq.0) then 
         write(0,*) 'ERROR: PARAMETER NAMES MUST START WITH A LETTER. SKIPPING LINE.'
         qerror=.true.; exit
        elseif (all(digits2(1:10).ne.c).and.underscore.ne.c) then
         write(0,*) 'ERROR: ILLEGAL CHARACTER IN PARAMETER NAME. SKIPPING LINE.'
         qerror=.true.; exit
        endif
       endif
       ltag=ltag+1; tag(ltag:ltag)=c
      elseif (qeq) then
       write(0,*) 'ERROR: MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'". SKIPPING LINE.'
       qerror=.true.; exit 
      elseif (qval) then
       lval=lval+1; val(lval:lval)=c
      else
       write(0,*) 'ERROR: UNKNOWN ERROR. SKIPPING LINE'
       qerror=.true.; exit 
      endif
!      
     else
       write(0,*) 'ERROR: UNRECOGNIZED CHARACTER "'//c//'". SKIPPING LINE.'
       exit 
     endif ! character loop
    enddo ! while l.gt.0
   else ! end of file
    exit
   endif
  enddo ! over all lines in the file
!
  if (qerror) then 
   write(0,*) 'ERROR(S) FOUND IN INPUT. ABORT.'
   stop
  endif
 end
 
 