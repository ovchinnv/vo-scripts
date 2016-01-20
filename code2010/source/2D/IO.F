c
c     ***************************************************************
c
C     subroutine for I/O - read or write a vector quantity
c
c     ***************************************************************
c
      SUBROUTINE IOVECTOR(NAME,UO,VO,WO,NNX,NNY,NNZ,NNZG,DIR,TLEVEL)

      IMPLICIT NONE
      include 'mpif.h'
      include 'headers/common.h'
      include 'headers/dimension.h'

      integer dir,nnx,nny,nnz,nnzg,type,jp,status(mpi_status_size)
      character*(*) name
      integer i,j,k,ierr
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real u (nx,ny,nz),v (nx,ny,nz),w (nx,ny,nz)
      real tlevel
      
      integer nfile

#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

c     if dir = -1 read else write

      if(dir .eq. -1) then
         write(0,*) myrank, ' is reading ',name
         open(19,file=name,status='old',form='unformatted')
         read(19) tlevel
         if(mysize .ne. 1) then
            do k=1,(myrank)*(nzg-2)/mysize
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
	    WRITE(100+myrank,*) 'INITIAL DATA READ from ',name
	    CLOSE(100+myrank)
      elseif(dir .eq. 1) then
c     note that all files read the input, but only root writes to output
         if(myrank.eq.0) then
            open(19,file=name,status='unknown',form='unformatted')
            write(19) tlevel
            do k=1,kz2
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
            enddo
            if(mysize .eq. 1) then
               k=kz2+1
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
               write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
            endif               
            do jp=1,mysize-1
               call mpi_recv(u(1,1,1),nx*ny*nz,type,jp,0,
     &              mpi_comm_hydro,status,ierr)
               call mpi_recv(v(1,1,1),nx*ny*nz,type,jp,0,
     &              mpi_comm_hydro,status,ierr)
               call mpi_recv(w(1,1,1),nx*ny*nz,type,jp,0,
     &              mpi_comm_hydro,status,ierr)
               do k=2,kz2+(jp+1)/mysize
                  write(19) ((u(i,j,k),i=1,nx),j=1,ny)
                  write(19) ((v(i,j,k),i=1,nx),j=1,ny)
                  write(19) ((w(i,j,k),i=1,nx),j=1,ny)
               enddo
            enddo
            close(19)
         else
            call mpi_send(uo(1,1,1),nx*ny*nz,type,0,0,
     &           mpi_comm_hydro,ierr)
            call mpi_send(vo(1,1,1),nx*ny*nz,type,0,0,
     &           mpi_comm_hydro,ierr)
            call mpi_send(wo(1,1,1),nx*ny*nz,type,0,0,
     &           mpi_comm_hydro,ierr)
         endif
      elseif(dir .eq. 2) then
c     note that all files  write to output
         open(19,file=name//CHAR(myrank+48),
     $        status='unknown',form='unformatted')
c
c$$$         write(6,*) myrank,nwrite,name,
c$$$     $        CHAR(myrank+48),name//CHAR(myrank+48)

         if(ifield.eq.2) then
            write(19) tlevel,myrank
         endif
         do k=kz1,kz2+1
            write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
      endif
c
      return
      end
c
c     ***************************************************************
c
C     subroutine for I/O - read or write a scalar quantity
c
c     ***************************************************************
c
      SUBROUTINE IOSCALAR(NAME,UO,NNX,NNY,NNZ,NNZG,DIR,TLEVEL)
c
      IMPLICIT NONE
      include 'mpif.h'
      include 'headers/common.h'
      include 'headers/dimension.h'
c
      integer dir,nnx,nny,nnz,nnzg,type,jp,status(mpi_status_size)
      character*(*) name
      integer i,j,k,ierr
      real uo(nx,ny,nz)
      real u (nx,ny,nz)
      real tlevel
#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif


C     IF DIR = -1 READ ELSE WRITE

      if(dir.eq.-1) then
         write(0,*) myrank, ' is reading ',name
         open(19,file=name,status='old',form='unformatted')
         read(19) tlevel
         if(mysize .ne. 1) then
            do k=1,(myrank)*(nzg-2)/mysize
               read(19)
            enddo
         endif
         do k=1,kz2+1
            read(19) ((uo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
      elseif(dir.eq.1) then
c     note that all files read the input, but only root writes to output
         if(myrank.eq.0) then
            open(19,file=name,status='unknown',form='unformatted')
c            if(ifield.eq.2) then
               write(19) tlevel
c            endif
            do k=1,kz2
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            enddo
            if(mysize .eq. 1) then
               k=kz2+1
               write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            endif               
            do jp=1,mysize-1
               call mpi_recv(u(1,1,1),nx*ny*nz,type,jp,0,
     &              mpi_comm_hydro,status,ierr)
               do k=kz1+1,kz2+(jp+1)/mysize
                  write(19) ((u(i,j,k),i=1,nx),j=1,ny)
               enddo
c               if(jp .eq. mysize-1) then
c                  k=kz2+1
c                  write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
c               endif
            enddo
            close(19)

         else
            call mpi_send(uo(1,1,1),nx*ny*nz,type,0,0,
     &           mpi_comm_hydro,ierr)
         endif
      elseif(dir .eq. 2) then
c     note that all files  write to output
         open(19,file=name//CHAR(myrank+48),
     $        status='unknown',form='unformatted')
         if(ifield.eq.2) then
            write(19) tlevel
         endif
         do k=kz1,kz2+1
            write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
      endif
      return
      end
c
c     ***************************************************************
c
C     subroutine for I/O - write a vector and two scalars
c
c     ***************************************************************
c
      SUBROUTINE IOVSS(NAME,UO,VO,WO,PO,TV,NX,NY,NZ,NZG,DIR,TLEVEL)

      IMPLICIT NONE
      include 'mpif.h'
      include 'headers/common.h'

      integer dir,nx,ny,nz,nzg,type,jp,status(mpi_status_size)
      character*(*) name
      integer i,j,k,ierr
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real po(nx,ny,nz),tv(nx,ny,nz)
      real tlevel
c      
      integer nfile
c
#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
      if(dir .eq. 2) then
c     note that all files  write to output
         open(19,file=name//CHAR(myrank+48),
     $        status='unknown',form='unformatted')
         if(ifield.eq.2) then
            write(19) tlevel,myrank
         endif
         do k=kz1,kz2+1
            write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((po(i,j,k),i=1,nx),j=1,ny)
            write(19) ((tv(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
      endif
c
      return
      end
c
      SUBROUTINE IORESTART(NAME,UO,VO,WO,P,UAB1,VAB1,WAB1,NX,NY,NZ,NZG,DIR,TLEVEL)

      IMPLICIT NONE
      include 'mpif.h'
      include 'headers/common.h'

      integer dir,nx,ny,nz,nzg,type,jp,status(mpi_status_size)
      character*(*) name
      integer i,j,k,n,ierr
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real uab1(nx,ny,nz),vab1(nx,ny,nz),wab1(nx,ny,nz),p(nx,ny,nz)
      real p1(nx,ny,nz)
      real u (nx,ny,nz),v (nx,ny,nz),w (nx,ny,nz)
      real u1 (nx,ny,nz),v1 (nx,ny,nz),w1 (nx,ny,nz)      
      real tlevel,dummy1
c      
      integer nfile

#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

c     if dir = -1 read else write

      if(dir .eq. -1) then
         open(19,file=name,status='unknown',form='unformatted')
         if(ifield.eq.2) then
            read(19) tlevel
         endif
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((uo(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
c
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((vo(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
c
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((wo(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
c
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((p(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
c
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((uab1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
c
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((vab1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
c
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=1,myrank)
         read(19) (((wab1(i,j,k),i=1,nx),j=1,ny),k=1,nz)
         read(19) ((((dummy1,i=1,nx),j=1,ny),k=1,nz),n=myrank,mysize-2)          
         close(19)
c
c
      elseif(dir .eq. 2) then
c     note that all files  write to output
         open(19,file=name//CHAR(myrank+48),
     $        status='unknown',form='unformatted')
c
c$$$         write(6,*) myrank,nwrite,name,
c$$$     $        CHAR(myrank+48),name//CHAR(myrank+48)

         if(ifield.eq.2) then
            write(19) tlevel,myrank
         endif
         do k=kz1,kz2+1
            write(19) ((uo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((vo(i,j,k),i=1,nx),j=1,ny)
            write(19) ((wo(i,j,k),i=1,nx),j=1,ny)
         enddo
         close(19)
      endif
c
      return
      end
c
