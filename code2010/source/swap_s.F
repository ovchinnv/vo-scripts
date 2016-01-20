      subroutine swap_xz(a,at,jmax,kmax,mysize,myrank,MPI_COMM_HYDRO)

      include 'headers/dimension.h'
     
      include "mpif.h"
      parameter(nproc_max=1024)
      integer MPI_COMM_HYDRO
#ifdef SINGLE
      real at(jmax/mysize,*)
      real a(jmax,*)
      real aux((nx-2)*(ny-2)*(nzg-2))
#endif
#ifdef DOUBLE
      real*8 at(jmax/mysize,*)
      real*8 a(jmax,*)
      real*8 aux((nx-2)*(ny-2)*(nzg-2))
#endif
      integer sendcounts(0:nproc_max-1),sdispl(0:nproc_max-1)
      integer reccounts(0:nproc_max-1),rdispl(0:nproc_max-1)
      
      logical flag

      save flag,sendcounts,sdispl,reccounts,rdispl

c      data flag/.true./

c#ifdef TIME_SWAP
c#endif
c     
c     We initialize the strides

c      if (flag) then
c         flag=.false.
      
         do j=0,mysize-1
            sendcounts(j)=jmax*kmax/mysize/mysize
            sdispl(j)=j*jmax*kmax/mysize/mysize
            
            reccounts(j)=jmax*kmax/mysize/mysize
            rdispl(j)=j*jmax*kmax/mysize/mysize
         enddo
c      endif

c   rearrange data before sending
      ind=1
      do jblk=0,mysize-1

       do k=1,kmax/mysize
        do j=1+jblk*jmax/mysize,(jblk+1)*jmax/mysize
      	 aux(ind)=a(j,k)
      	 ind=ind+1
      	enddo
       enddo	 

      enddo
cccccccccccccccccccccccc
c     The bulk of the job
#ifdef DOUBLE         
         call mpi_alltoallv(aux,sendcounts,sdispl,
     &        MPI_DOUBLE_PRECISION,
     &        at,reccounts,rdispl,MPI_DOUBLE_PRECISION,
     &        MPI_COMM_HYDRO,ierr)

#endif
#ifdef SINGLE         
         call mpi_alltoallv(aux,sendcounts,sdispl,
     &        MPI_REAL,
     &        at,reccounts,rdispl,MPI_REAL,
     &        MPI_COMM_HYDRO,ierr)

#endif

      return
      end


      subroutine swap_zx(at,a,jmax,kmax,mysize,myrank,MPI_COMM_HYDRO)
      
      include 'headers/dimension.h'
      include "mpif.h"
      
      integer MPI_COMM_HYDRO
      parameter(nproc_max=1024)
#ifdef SINGLE
      real at(jmax/mysize,*)
      real a(jmax,*)
      real aux((nx-2)*(ny-2)*(nzg-2))
#endif
#ifdef DOUBLE
      real*8 at(jmax/mysize,*)
      real*8 a(jmax,*)
      real*8 aux((nx-2)*(ny-2)*(nzg-2))
#endif
      integer ind

      integer sendcounts(0:nproc_max-1),sdispl(0:nproc_max-1)
      integer reccounts(0:nproc_max-1),rdispl(0:nproc_max-1)
      
      logical flag

      save flag,sendcounts,sdispl,reccounts,rdispl

c      data flag/.true./
      
c     
c     We initialize the strides (for symmetry, now the sends are 
c     really recvs and viceversa)
c      if (flag) then
c         flag=.false.
      
         do j=0,mysize-1
            sendcounts(j)=jmax*kmax/mysize/mysize
            sdispl(j)=j*jmax*kmax/mysize/mysize
            
            reccounts(j)=jmax*kmax/mysize/mysize
            rdispl(j)=j*jmax*kmax/mysize/mysize
         enddo
c      endif

c     The bulk of the job
#ifdef DOUBLE         
         call mpi_alltoallv(
     &        at,reccounts,rdispl,MPI_DOUBLE_PRECISION,
     &        aux,sendcounts,sdispl,
     &        MPI_DOUBLE_PRECISION,
     &        MPI_COMM_HYDRO,ierr)
#endif

#ifdef SINGLE
         call mpi_alltoallv(
     &        at,reccounts,rdispl,MPI_REAL,
     &        aux,sendcounts,sdispl,
     &        MPI_REAL,
     &        MPI_COMM_HYDRO,ierr)

#endif
c   rearrange data
      ind=1
      do jblk=0,mysize-1

       do k=1,kmax/mysize
        do j=1+jblk*jmax/mysize,(jblk+1)*jmax/mysize
      	 a(j,k)=aux(ind)
      	 ind=ind+1
      	enddo
       enddo	 

      enddo
c  
      return
      end
      
c
      subroutine swap_yz(a,at,imax,jmax,kmax,mysize,myrank,MPI_COMM_HYDRO)

      include "mpif.h"
      parameter(nproc_max=1024)
      integer MPI_COMM_HYDRO
#ifdef SINGLE
      real at(imax/mysize,jmax,*)
      real a(imax,jmax,*)
#endif
#ifdef DOUBLE
      real*8 at(imax/mysize,jmax,*)
      real*8 a(imax,jmax,*)
#endif
      integer sendcounts(0:nproc_max-1),sdispl(0:nproc_max-1)
      integer reccounts(0:nproc_max-1),rdispl(0:nproc_max-1)
      
      logical flag

      save flag,sendcounts,sdispl,reccounts,rdispl

         do j=i,mysize-1
            sendcounts(i)=imax/mysize
            sdispl(i)=i*imax/mysize
            
            reccounts(i)=imax/mysize
            rdispl(i)=i*imax/mysize*jmax*kmax/mysize
         enddo

c     The bulk of the job
      do k=1,kmax/mysize
       do j=1,jmax
#ifdef DOUBLE         
         call mpi_alltoallv(a(1,j,k),sendcounts,sdispl,
     &        MPI_DOUBLE_PRECISION,
     &        at(1,j,k),reccounts,rdispl,MPI_DOUBLE_PRECISION,
     &        MPI_COMM_HYDRO,ierr)

#endif
#ifdef SINGLE         
         call mpi_alltoallv(a(1,j,k),sendcounts,sdispl,
     &        MPI_REAL,
     &        at(1,j,k),reccounts,rdispl,MPI_REAL,
     &        MPI_COMM_HYDRO,ierr)

#endif
       enddo
      enddo      

      return
      end


      subroutine swap_zy(at,a,imax,jmax,kmax,mysize,myrank,MPI_COMM_HYDRO)
      
      include "mpif.h"

      integer MPI_COMM_HYDRO
      parameter(nproc_max=1024)
#ifdef SINGLE
      real at(imax/mysize,jmax,*)
      real a(imax,jmax,*)
#endif
#ifdef DOUBLE
      real*8 at(imax/mysize,jmax,*)
      real*8 a(imax,jmax,*)
#endif

      integer sendcounts(0:nproc_max-1),sdispl(0:nproc_max-1)
      integer reccounts(0:nproc_max-1),rdispl(0:nproc_max-1)
      
      logical flag

      save flag,sendcounts,sdispl,reccounts,rdispl

         do i=0,mysize-1
            sendcounts(i)=imax/mysize
            sdispl(i)=i*imax/mysize
            
            reccounts(i)=imax/mysize
            rdispl(i)=i*imax/mysize*jmax*kmax/mysize
         enddo
c      endif

c     The bulk of the job
      do k=1,kmax/mysize
       do j=1,jmax
#ifdef DOUBLE         
         call mpi_alltoallv(
     &        at(1,j,k),reccounts,rdispl,MPI_DOUBLE_PRECISION,
     &         a(1,j,k),sendcounts,sdispl,
     &        MPI_DOUBLE_PRECISION,
     &        MPI_COMM_HYDRO,ierr)

#endif
#ifdef SINGLE
         call mpi_alltoallv(
     &        at(1,j,k),reccounts,rdispl,MPI_REAL,
     &         a(1,j,k),sendcounts,sdispl,
     &        MPI_REAL,
     &        MPI_COMM_HYDRO,ierr)

#endif
       enddo
      enddo      

      return
      end

