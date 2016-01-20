      module size
!     contains basic information about the size of the calculation      
       private ! keep mpif.h hidden here
#ifdef PARALLEL
       include 'mpif.h'
#endif
!
       public size_initialize
       int, save, public ::  nx=0, ny=0, nz=0, nzg=0
       int, save, public :: ix1, ix2, jy1, jy2, kz1, kz2
!       
       int, save, public ::  communicator, ncpu=-1, me=-1
       logical, save, public ::  size_initialized

      contains
       subroutine size_initialize(mx, my, mz, comm)
       implicit none
       int :: mx, my, mz, mzg
       integer*4, optional :: comm
       int :: bug
       character*(8) :: whoami
       data whoami /'SET_SIZE'/
!
       size_initialized=.false.
!
#ifdef PARALLEL
!
       communicator=MPI_COMM_NULL
!
       if (present(comm)) then 
        if (comm.ne.MPI_COMM_NULL) then
         call MPI_SIZE(comm, ncpu, bug)
         call MPI_RANK(comm, me, bug)
         communicator=comm
        endif
       endif
#else
       me=0
       ncpu=1
       communicator=-1
#endif ! PARALLEL
!
       if (mx.gt.0) then
        nx=mx
       else
        call error(0, whoami, ' NX </= ZERO. ABORTING.')
        return
       endif
!  
       if (my.gt.0) then
        if (mod(my-2,ncpu).ne.0) then 
         call error(0, whoami, ' NY-2 NOT DIVISIBLE BY NCPU. ABORTING.')
         return
        else
         ny=my
        endif
       else
        call error(0, whoami, ' NY </= ZERO. ABORTING.')
        return
       endif
!       
       if (mz.gt.0) then
        if (mod(mz-2,ncpu).ne.0) then 
         call error(0, whoami, ' NZ-2 NOT DIVISIBLE BY NCPU. ABORTING.')
         return
        else
         nzg=mz
         nz=(nzg-2)/ncpu+2
        endif
       else
        call error(0, whoami, ' NZ </= ZERO. ABORTING.')
        return
       endif
!
!    indices:
       ix1=1;    jy1=1;     kz1=1;
       ix2=nx-1; jy2=ny-1;  kz2=nz-1;
!
       size_initialized=.true.
!
       end subroutine size_initialize
      end module size
!             