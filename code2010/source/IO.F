      module IO
      use size
      private
!
      public IOVECTOR ! read/write velocity field from restart file
      public IOSCALAR ! read/write once component of a velocity field (or pressure field)
      public IOVSS
      public IORESTART ! read restart files
      public IOVECTOR2D
      
      contains
!     ***************************************************************
!     subroutines for I/O
!     ***************************************************************
       SUBROUTINE IOVECTOR(NAME,UO,VO,WO,DIR,TLEVEL)
!
       IMPLICIT NONE
!
#ifdef PARALLEL
       include 'mpif.h'
#endif
       int ::  dir,jp
       integer*4 :: status(mpi_status_size)
       character*(*) name
       int :: i,j,k,ierr
       float, dimension(nx,ny,nz) :: uo,vo,wo,u,v,w
       float :: tlevel
      
       character*(8) :: whoami
       data whoami /'IOVECTOR'/
!
!
       if (.not.size_initialized) then 
        call error(whoami, ' SIZE NOT INITIALIZED. NOTHING DONE.',0)
        return
       endif
!
!     if dir = -1 read else write
!
       if(dir .eq. -1) then
         write(output,*) me, ' is reading ',name
         open(19,file=name,status='old',form='unformatted')
         read(19) tlevel
         if(ncpu .ne. 1) then
            do k=1,(me)*(nzg-2)/ncpu
               read(19)
               read(19)
               read(19)
            enddo
         endif
         do k=1,kz2+1
            read(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            read(19) ((vo(i,j,k),i=1,nx),j=1,ny)
            read(19) ((wo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
!	    WRITE(100+me,*) 'INITIAL DATA READ from ',name
!	    CLOSE(100+me)
      elseif(dir .eq. 1) then
!     note that all files read the input, but only root writes to output
         if(me.eq.0) then
            open(19,file=name,status='unknown',form='unformatted')
            write(19) tlevel
            do k=1,kz2
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
            enddo
            if(ncpu .eq. 1) then
               k=kz2+1
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
            endif               
#if PARALLEL
            do jp=1,ncpu-1
               call mpi_recv(u(1,1,1),nx*ny*nz,mpi_float,jp,0,
     &              communicator,status,ierr)
               call mpi_recv(v(1,1,1),nx*ny*nz,mpi_float,jp,0,
     &              communicator,status,ierr)
               call mpi_recv(w(1,1,1),nx*ny*nz,mpi_float,jp,0,
     &              communicator,status,ierr)
               do k=2,kz2+(jp+1)/ncpu
                  write(19) ((u(i,j,k),i=1,nx),j=1,ny)
                  write(19) ((v(i,j,k),i=1,nx),j=1,ny)
                  write(19) ((w(i,j,k),i=1,nx),j=1,ny)
               enddo
            enddo
#endif
            close(19)
         else
#if PARALLEL
            call mpi_send(uo(1,1,1),nx*ny*nz,mpi_float,0,0,
     &           communicator,ierr)
            call mpi_send(vo(1,1,1),nx*ny*nz,mpi_float,0,0,
     &           communicator,ierr)
            call mpi_send(wo(1,1,1),nx*ny*nz,mpi_float,0,0,
     &           communicator,ierr)
#endif
         endif
!
       elseif(dir .eq. 2) then
!     note that all files  write to output
         open(19,file=name//CHAR(me+48),
     $        status='unknown',form='unformatted')
!
         do k=kz1,kz2+1
            write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
       endif
!
       end subroutine IOVECTOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       SUBROUTINE IOSCALAR(NAME,UO,DIR,TLEVEL)
!
       IMPLICIT NONE
#ifdef PARALLEL
       include 'mpif.h'
#endif
!
       int :: dir,jp
       integer*4 :: status(mpi_status_size)
       character*(*) name
       int :: i,j,k,ierr
       float ::  uo, u
       float :: tlevel
!
       character*(8) :: whoami
       data whoami /'IOSCALAR'/
!
       if (.not.size_initialized) then 
        call error(me, whoami, ' SIZE NOT INITIALIZED. NOTHING DONE.')
        return
       endif
!
!     IF DIR = -1 READ ELSE WRITE
       if(dir.eq.-1) then
         write(0,*) me, ' is reading ',name
         open(19,file=name,status='old',form='unformatted')
         read(19) tlevel
         if(ncpu .ne. 1) then
            do k=1,(me)*(nzg-2)/ncpu
               read(19)
            enddo
         endif
         do k=1,kz2+1
            read(19) ((uo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
       elseif(dir.eq.1) then
!     note that all files read the input, but only root writes to output
         if(me.eq.0) then
            open(19,file=name,status='unknown',form='unformatted')
               write(19) tlevel
            do k=1,kz2
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            enddo
            if(ncpu .eq. 1) then
               k=kz2+1
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            endif               
#if PARALLEL
            do jp=1,ncpu-1
               call mpi_recv(u(1,1,1),nx*ny*nz,mpi_float,jp,0,
     &              communicator,status,ierr)
               do k=kz1+1,kz2+(jp+1)/ncpu
                  write(19) ((u(i,j,k),i=1,nx),j=1,ny)
               enddo
            enddo
#endif
            close(19)

         else
#if PARALLEL
            call mpi_send(uo(1,1,1),nx*ny*nz,mpi_float,0,0,
     &           communicator,ierr)
#endif
         endif
!
       elseif(dir .eq. 2) then
!     note that all files  write to output
         open(19,file=name//CHAR(me+48),
     &        status='unknown',form='unformatted')
            write(19) tlevel
         do k=kz1,kz2+1
            write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
       endif
!
       end subroutine IOSCALAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
       SUBROUTINE IORESTART(NAME,UO,VO,WO,P,UAB1,VAB1,WAB1,DIR,TLEVEL)
!     read old-style restart file
       IMPLICIT NONE
#if PARALLEL
       include 'mpif.h'
#endif
!
       int :: dir,jp
       integer*4 :: status(mpi_status_size)
       character*(*) name
       int :: i, j, k, n, ierr
       float, dimension(nx,ny,nz) ::  
     &   uo, vo, wo, p, uab1, vab1, wab1, p
       float :: tlevel, dummy1
!
       character*(9) :: whoami
       data whoami /'IORESTART'/
!
       if (.not.size_initialized) then 
        call error(me, whoami, ' SIZE NOT INITIALIZED. NOTHING DONE.')
        return
       endif
!
!     if dir = -1 read else write
!
       if(dir .eq. -1) then
         open(19,file=name,status='unknown',form='unformatted')
         read(19) tlevel
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((uo(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
!
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((vo(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
!
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((wo(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
!
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((p(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
!
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((uab1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
!
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((vab1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
!
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,me)
         read(19) (((wab1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=me,ncpu-2)          
         close(19)
       endif
!
      end subroutine IORESTART
!ccccccccccccccccccccccccc routine to read from a 2D field
       SUBROUTINE IOVECTOR2D(NAME,UO,VO,WO,DIR,TLEVEL)
!
       IMPLICIT NONE
!
       int :: dir,jp
       integer*4 :: status(mpi_status_size)
       character*(*) name
       integer :: i,j,k,ierr
      
       float, dimension(nx,ny,nz) :: uo, vo, wo
       float tlevel
!
       character*(10) :: whoami
       data whoami /'IOVECTOR2D'/
!
       if (.not.size_initialized) then 
        call error(me, whoami, ' SIZE NOT INITIALIZED. NOTHING DONE.')
        return
       endif
!
!
!     if dir = -1 read else write
!
       if(dir .eq. -1) then
         write(0,*) me, ' is reading ',name
         open(19,file=name,status='old',form='unformatted')
         read(19) tlevel
         do k=1,(me)*(nzg-2)/ncpu
           read(19)
           read(19)
           read(19)
         enddo
         do k=1,kz2+1
            read(19) ((uo(i,j,k),i=1,nx),j=1,3)
            read(19) ((vo(i,j,k),i=1,nx),j=1,3)
            read(19) ((wo(i,j,k),i=1,nx),j=1,3)
         enddo
         close(19)
c	    WRITE(100+me,*) 'INITIAL DATA READ from ',name
c	    CLOSE(100+me)
       endif
!      fill other spanwise positions
       do j=3,ny
        uo(:,j,:)=uo(:,2,:)
        vo(:,j,:)=vo(:,2,:)
        wo(:,j,:)=wo(:,2,:)
       enddo
!
       end subroutine IOVECTOR2D
!
      end module IO
!