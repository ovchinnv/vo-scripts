      module parser
      private
!     read, parse & store input file; 
!
!    define a derived type to store simulation parameters
      type params
       character*200, dimension(:), pointer :: tag, val
       int, dimension(:), pointer :: tlen, vlen
       int :: length ! length of the vector
       int :: last ! index of last element
       bool :: initialized=.false. ! has the vector been initialized
      end type params
!
      int, parameter, private :: expand_incr=200
!
      bool, save :: parser_initialized=.false. ! set to true after parse_file is called successfully; private
      type (params) :: parameters
!
      public atoi     	! convert string to intger 
      public atod	! convert string to double
      public atol	! convert string to bool
      public getval	! return tag value
      public parse_file ! read input file and store all parameters
      public list_params ! list parameters
!
      contains
!******************************************** implement data routines*********************
       subroutine params_init( v )
       implicit none
       type (params) :: v
c       if (associated(v%...)) deallocate(v%...) ! testing unassigned pointer is an error!
       allocate(v%tag(expand_incr),v%val(expand_incr),v%tlen(expand_incr),v%vlen(expand_incr))
c
       v%tag=''; v%value=''; v%tlen=0; v%vlen=0
       v%length=expand_incr
       v%last=0
       v%initialized=.true.
       end subroutine params_init
!ccccc
       subroutine params_done( v )
       implicit none
       type (params) :: v
       if (associated(v%tag)) deallocate(v%tag)
       if (associated(v%val)) deallocate(v%val)
       if (associated(v%tlen)) deallocate(v%tlen)
       if (associated(v%vlen)) deallocate(v%vlen)
       v%length=0
       v%last=0
       v%initialized=.false.
       end subroutine params_done
!
       subroutine params_expand( v )
       implicit none
       type (params) :: v
       int :: newlength
       character*200, dimension(:), allocatable, target :: ntag, nval
       int, dimension(:), allocatable, target :: ntlen, nvlen
!
       if (.not.v%initialized) then 
        call params_init(v) 
       else
!    assume length is valid
        newlength=v%length+expand_incr ! temporary storage space
        allocate(ntag(newlength),nvalue(newlength),ntlen(newlength),nvlen(newlength)) ! copy old data
        ntag(1:v%length)=v%tag; nval(1:v%length)=v%val; ntlen(1:v%length)=v%tlen; nvlen(1:v%length)=v%vlen ! deallocate old array
        deallocate(v%tag,v%val,v%tlen,v%vlen) ! point to new data
        v%tag=>ntag; v%val=>nval; v%tlen=>ntlen; v%vlen=>nvlen
        v%length=newlength
       endif 
       end subroutine params_expand
!ccccc
       function params_add(v,newtag,newval,ltag,lval) ! add a new element to the list (not necessarily unique) 
!                                       and return its index
       implicit none
       type (params) :: v
       int :: params_add, ltag, lval
       character*(*) :: newtag, newval
!
       if (.not.v%initialized) call params_init(v) 
!    add element to the list
       if (v%last.eq.v%length) call params_expand(v)       
       j=v%last+1
       v%tag(j)=newtag(1:ltag); v%tlen(j)=ltag
       v%val(j)=newval(1:lval); v%vlen(j)=lval
       v%last=j
       params_add=j
       end function params_add
!ccccc
       function params_uadd(v,newtag,newval,ltag,lval) ! add a UNIQUE new element to the list and return its index
!                                                         if the element already exists, overwrite and warn       
       implicit none
       type (params) :: v
       int :: j, params_uadd
!
       if (.not.v%initialized) call params_init(v) 
       do j=1,v%last
        if (all(v%tag(j).eq.newtag) then
!       found element
         params_uadd=j
         call warning('PARAMS_UADD','Parameter '//newtag(1:ltag)//' already present and has the value '//v%val(1:v%vlen(j))//'. Will overwrite.',0)
         v%val(j)=newval(1:lval); v%vlen(j)=lval
         return
        endif
       enddo
!    add element to the list: use regular routine
       params_uadd=params_add(v,newtag,newval,ltag,lval)
!
       end function params_uadd
!ccccc
       function params_getval( v,atag ) ! returns v%i(j) if j is valid
       implicit none
       type (params) :: v
       character*200 :: params_getval, atag
       int :: j
!
       params_getval=''
       if (.not.v%initialized) then 
!
       elseif (j.gt.v%last.or.j.le.0) then 
!
       else
        do j=1,v%last
         if (all(v%tag(j).eq.atag) then
!       found element
          params_getval=v%val(1:v%vlen(j))
          return
         endif
        enddo
       endif
!
       end function params_getval
!****************************************** end of data routines **************************************
       function getval(atag)
       implicit none
       character*200 :: value, atag
!
       command=''
       if (.not.parser_initialized) then  
        call error('GETVAL','Parameter '//trim(atag)//' Not found. Cannot continue.',-1)
       else
        value=params_getval(parameters,atag)
        if (len(trim(value)).eq.0) call error('GETVAL','Parameter '//trim(atag)//' Not found. Cannot continue.',-1)
       endif
!
       end function getval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine parser_init()
       implicit none
!    initialize params structure
       call params_done(parameters)
       call params_init(parameters)
       parser_initialized=.true.
       end subroutine parser_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine parser_done()
       implicit none
       call params_done(parameters)
       parser_initialized=.false.
       end subroutine parser_done
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       subroutine parse_file(fid)
       use size, only: me, communicator
       use fout
       implicit none
!
#ifdef PARALLEL
       INCLUDE 'mpif.h'
#endif
       int :: fid ! input file handle
!
       int :: i, j
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
       character*10, parameter :: whoami = 'PARSE_FILE'
!
       int :: allowed_len, lower_len, upper_len, digits_len
       character :: allowed2(200),c
       logical :: qtag, qeq, qval, qerror
!
       int :: ioerr, ierr
       int :: l=0, ltag=0, lval=0
       character*400 :: cmdline
       character*200 :: tag, val
!
       allowed_len=len_trim(allowed)
       do i=1,allowed_len
        allowed2(i)=allowed(i:i)
       enddo 
!
!     do work
!
       qerror=.false.
       do while (.true.)
        if (me.le.0) read(1,'(A)',IOSTAT=ioerr) cmdline ! if running in parallel, then only the root node is passed a valid handle
#ifdef PARALLEL
        if (communicator.ne.MPI_COMM_NULL) call MPI_BCAST(ioerr,1,MPI_INT,0,communicator,ierr)
#endif
        if (ioerr.eq.0) then
#ifdef PARALLEL
        if (communicator.ne.MPI_COMM_NULL) call MPI_BCAST(cmdline,len(cmdline),MPI_BYTE,0,communicator,ierr) ! broadcast to all CPUs
#endif
!   write(0,*) cmdline
! read from the line
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
             call warning(whoami,'UNEXPECTED END OF LINE',0)
             qerror=.true.; exit 
            elseif (qval) then
             if (lval.eq.0) then
              call warning(whoami, 'MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'"; SKIPPING LINE.',0)
              qerror=.true.; exit 
             else ! store tag and val
              write(output,'(A)') tag(1:ltag)//' = '//val(1:lval))
              j=params_uadd(parameters,tag,val,ltag,lval)
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
              call warning(whoami, 'MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'"; SKIPPING LINE.',0)
              qerror=.true.; exit 
             else ! store tag and val
              write(output,'(A)') tag(1:ltag)//' = '//val(1:lval))
              j=params_uadd(parameters,tag,val,ltag,lval)
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
              call warning(whoami, 'MISSING PARAMETER NAME. SKIPPING LINE ',0)
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
              call warning(whoami, 'MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'". SKIPPING LINE.',0)
              qerror=.true.; exit 
             else ! store tag and val
              write(output,'(A)') tag(1:ltag)//' = '//val(1:lval))
              j=params_uadd(parameters,tag,val,ltag,lval)
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
               call warning(whoami, 'PARAMETER NAMES MUST START WITH A LETTER. SKIPPING LINE.',0)
               qerror=.true.; exit
              elseif (all(digits2(1:10).ne.c).and.underscore.ne.c) then
               call warning(whoami, 'ILLEGAL CHARACTER IN PARAMETER NAME. SKIPPING LINE.',0)
               qerror=.true.; exit
              endif
             endif
             ltag=ltag+1; tag(ltag:ltag)=c
            elseif (qeq) then
             call warning(whoami, 'MISSING VALUE FOR PARAMETER "'//tag(1:ltag)//'". SKIPPING LINE.',0)
             qerror=.true.; exit 
            elseif (qval) then
             lval=lval+1; val(lval:lval)=c
            else
             call warning(whoami, 'UNKNOWN ERROR. SKIPPING LINE',0)
             qerror=.true.; exit 
            endif
!      
           else
             call warning(whoami, 'UNRECOGNIZED CHARACTER "'//c//'". SKIPPING LINE.',0)
             exit 
           endif ! character loop
         enddo ! while l.gt.1
        else ! end of file
         exit
        endif
       enddo ! over all lines in the file
!
       endif ! me.le.0
!
#ifdef PARALLEL
       call MPI_BCAST(error,MPI_BOOL,1,0,communicator,ioerror)
#endif
       if (qerror) then 
         call error(whoami, 'ERROR(S) FOUND IN INPUT. CANNOT CONTINUE.',-1)
       else
         call message(whoami, 'INPUT FILE READ.')
       endif
!
       end subroutine parse_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine list_params()
! 
       use size, only: me
       use fout
       implicit none
       character*11, parameter :: whoami = 'LIST_PARAMS'
       int :: i
!
       if (.not.parser_initialized) call parser_init()
!
       if (me.le.0) then
        write(output,'(A)') whoami//' : THE FOLLOWING PARAMETERS ARE DEFINED'
        do i=1,parameters%last
         write(output,'(A)') '\t'//parameters%tag(i)(1:parameters%tlen(i))//' == '//parameters%val(i)(1:parameters%vlen(i))
        enddo
       endif
       end subroutine list_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! auxiliary functions (they need not be part of this module)
       function atoi(a)
! NOTE: no overflow check yet
       use fout, only: error
       implicit none
       int :: atoi, i, l, j, k, sgn, base
       character*200 :: a, b
       character*10, parameter  :: digits='0123456789'
       character, parameter :: hyphen='-'
       character, parameter :: digits2(10)=(/ (digits(i:i),i=1,10)/)
       character*4, parameter :: whoami = 'ATOI'
! convert string to integer
       i=0
!
       b=adjustl(a)
       l=len_trim(b)
       if (l.ge.1) then
        if (l(1:1).eq.hyphen) then
         sgn=-1
         b(1:l-1)=b(2:l); l=l-1
        else
         sgn=1
        endif
       endif ! l.ge.1
!       
       base=1
       do j=l,1,-1
        k=maxloc(digits2.eq.b(l))(1)-1
        if (k.lt.0) then 
         call error(whoami, 'ERROR CONVERTING STRING "'//a//'" TO INTEGER.',-1)
         i=-999; sgn=1;
         exit
        else
         i=i+base*k
         base=base*10
        endif
       enddo
       atoi=i*sgn
       end function atoi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atof(a)
       use fout, only: error
       implicit none
       float :: atof,f 
       int :: i, l, j, k, sgn, base
       character*200 :: a, b
       character*10, parameter  :: digits='0123456789'
       character, parameter :: hyphen='-'
       character, parameter :: decimal='.'
       character, parameter :: digits2(10)=(/ (digits(i:i),i=1,10)/)
       character*4, parameter :: whoami = 'ATOF'
       bool :: fraction=.false.
! convert string to floating point number
       f=0
!
       b=adjustl(a)
       l=len_trim(b)
       if (l.ge.1) then
        if (l(1:1).eq.hyphen) then
         sgn=-1
         b(1:l-1)=b(2:l); l=l-1
        else
         sgn=1
        endif
       endif ! l.ge.1
!       
       base=1
       do j=l,1,-1
        if (b(l).eq.decimal) then
         if (fraction) then ! two decimal points are invalid
          call error(whoami, 'ERROR CONVERTING STRING "'//a//'" TO REAL.',-1)
          f=NaN; sgn=1;
          exit
         else
          fraction=.true.
          f=10.0*f/base
          base=1
         endif
        else	 
         k=maxloc(digits2.eq.b(l))(1)-1
         if (k.lt.0) then 
          call error(whoami, 'ERROR CONVERTING STRING "'//a//'" TO REAL.',-1)
          f=NaN; sgn=1;
          exit
         else
          f=f+1.0*base*k
          base=base*10
         endif
        endif
       enddo
       atoi=f*sgn
       end function atof
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function atol(a)
       implicit none
       bool :: atol, l
!
       end function atol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module parser
