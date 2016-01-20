      PROGRAM eddie_driver
      use eddie
      implicit none      
      integer :: communicator=-1 ! default (serial case)
#ifdef PARALLEL
      include 'mpif.h'
      integer :: ierr
      call mpi_init(ierr)
      communicator=MPI_COMM_WORLD
      call eddy_init(communicator)
#else
      call eddy_init()
#endif
!
      call eddie()
!
#ifdef PARALLEL
      call mpi_finalize(Ierr)
#endif
!
      end
!
