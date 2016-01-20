      subroutine pfft_yz(p,dp,dummyt,dummy,dummy2,au,ap,dy,dz,me,comm,datatype)
c     V. Ovchinnikov     8/06/03
c     DO NOT MODIFY!!! No warranty whatsoever
      implicit none
      include 'headers/dimension.h'
      include 'mpif.h'
c     
      integer me,comm,datatype,ierr,status(MPI_STATUS_SIZE)
      real p(nx,ny,nz),dp(nx,ny,nz),
     &     dummy(nx-2,ny-2,nz-2),dummyt(nx-2,(ny-2)/ncpu,nzg-2),
     &     dummy2((ny-2)/ncpu,nx-2,nzg-2),ftemp(ny,nx-2)
c     
      real dy,dz,au(nx),ap(nx)
      real*8 ky(ny/2), kz(nzg-2), kabs(ny/2,nzg-2)
      real*8 pi,fac,dy2,dz2
      integer i,j,k,ii,jj,kk,jbeg,jend
      logical getKabs
c     work arrays for fft routines
      integer ifaxy(13),ifaxz(13)
      real trigsy(3*(ny-2)/2+1),trigsz(2*(nzg-2))
      real work(ny*nx+2*nx*nzg) ! over-dimensioned
c     matrix coefficients and auxiliary arrays    
      real a(nx-2),b(nx-2),c(nx-2),r(nx-2),u(nx-2)
      real oddball(2*(nzg-2),nx-2)
      real zero(2*(nx-2))
c     
      save kabs,getKabs,ifaxy,ifaxz,trigsy,trigsz,a,c,dy2,dz2
      data getKabs /.true./
c
      if (getKabs) then
         if (me.eq.0) write(6,*) 'Initializing poisson solver...'
         getKabs=.false.
         pi=4.0*atan(1.0)
	 dy2=dy*dy
	 dz2=dz*dz
c     fill wavenumber arrays
c     write(6,*) 'generating wavenumber arrays'
c     kx
         ky(1)=0      
         do j=2,ny/2
            ky(j)=ky(j-1)+2.*pi/(ny-2)
         enddo
c     
c     Commented below because the function is real and the transform is symmetric
c     thus only half of the wavenumbers are needed for the first component of the frequency
c     
c     do j=ny/2+2,ny
c     ky(j)=ky(j-1)+2*pi/ny-2
c     enddo
c     
c     kz
         kz(1)=0      
         do k=2,nzg/2-1
            kz(k)=kz(k-1)+2.*pi/(nzg-2)	    
         enddo
         kz(nzg/2)=-pi
         do k=nz/2+1,nzg-2
           kz(k)=kz(k-1)+2.*pi/(nzg-2)
         enddo
         do k=1,nzg-2
            do j=1,ny/2
                  kabs(j,k)=2.*((cos(ky(j))-1)*dy2+(cos(kz(k))-1)*dz2)
            enddo
         enddo
      kabs(1,1)=1. ! to avoid overflow on mode 0
c     also initialize FFT routines
      call fftfax(ny-2,ifaxy,trigsy)
      call cftfax(nzg-2,ifaxz,trigsz)
c     partially set up matrix coefficients
        do i=1,nx-2
	    a(i)=au(i)*ap(i+1)
	    c(i)=au(i+1)*ap(i+1)
	enddo    
c     enforce correct pressure boundary condition
      do k=2,nz-1
         do i=2,nx-1
            p(i,ny,k)=p(i,2,k)
            p(i,1,k)=p(i,ny-1,k)
         enddo
      enddo
      do k=2,nz-1
         do j=2,ny-1
            p(nx,j,k)=p(nx-1,j,k)
            p(1,j,k)=p(2,j,k)
         enddo
      enddo
      if (ncpu.eq.1) then 
        do j=2,ny-1
           do i=2,nx-1
              p(i,j,1)=p(i,j,nz-1)
              p(i,j,nz)=p(i,j,2)
           enddo
        enddo
      else	
        call refreshbc(p,nx*ny,nz-2,ncpu,me,comm)
c       now enforce periodic boundary conditions in z     
        if (me.eq.0) then
          call MPI_SEND(p (1,1,2),nx*ny,datatype,ncpu-1,2,comm,ierr)
          call MPI_RECV(p (1,1,1),nx*ny,datatype,ncpu-1,4,comm,status,ierr)
        elseif (me.eq.ncpu-1) then  	
          call MPI_RECV(p (1,1,nz  ),nx*ny,datatype,0,2,comm,status,ierr)
          call MPI_SEND(p (1,1,nz-1),nx*ny,datatype,0,4,comm,ierr)
        endif
      endif ! ncpu	
c
      endif ! getKabs
c     
c     enforce periodic boundary conditions in y
      do i=1,nx
       do k=1,nz
        dp(i,ny,k)=dp(i,2,k)
	dp(i,1,k)=dp(i,ny-1,k)
       enddo
      enddo 	    
c       
c     now take the transform
      do k=2,nz-1
c      swap x and y direction in a dummy array
       do i=1,nx-2
        do j=1,ny
	 ftemp(j,i)=dp(i+1,j,k)
	enddo
       enddo    	 
c     
       call fft99(ftemp,work,trigsy,ifaxy,1,ny,ny-2,nx-2,-1)
c
c     shift data before swapping & store oddball with 0th mode 
       do i=1,nx-2
        do j=1,ny-2
         dummy(i,j,k-1)=ftemp(j,i)
	enddo
	 dummy(i,2,k-1)=ftemp(ny-1,i)
       enddo
      enddo 	 
c     swap x and z directions
      call swap_xz(dummy,dummyt,(ny-2)*(nx-2),nzg-2,ncpu,me,comm)
c     remove oddball 
      if (me.eq.0) then
        do k=1,nzg-2
	 do i=1,nx-2
	  oddball(2*k-1,i)=dummyt(i,2,k)
	  oddball(2*k,i)=0.
	  dummyt(i,2,k)=0.
	 enddo
	enddo  
      endif	
c     switch x and y dimensions for transform
      do k=1,nzg-2
        do i=1,nx-2
	  do j=1,(ny-2)/ncpu
	    dummy2(j,i,k)=dummyt(i,j,k)
	  enddo
	enddo
      enddo	      
c
c     call second fft in z (all modes except oddball)	
      do j=1,(ny-2)/ncpu,2
        call cfft99(dummy2(j,1,1),work,trigsz,ifaxz,(nx-2)*(ny-2)/ncpu/2,(ny-2)/ncpu/2,nzg-2,nx-2,-1)
      enddo

c     now the oddball
      if (me.eq.0)
     &   call cfft99(oddball(1,1),work,trigsz,ifaxz,1,(nzg-2),nzg-2,nx-2,-1)
c     solve tridiagonal system
c     save 0th mode
       if (me.eq.0) then 
         do i=1,nx-2
	   zero(2*i-1)=dummy2(1,i,1)
	   zero(2*i  )=dummy2(2,i,1)
	 enddo
       endif
c       	   
      fac=1./(nzg-2) ! this factor is gained in cfft99
c
      do k=1,nzg-2
        do j=1,(ny-2)/ncpu,2
c  real part	  
          do i=1,nx-2
            r(i)=fac*dummy2(j,i,k)
	    b(i)=-(a(i)+c(i))+kabs((j+1+(ny-2)*me/ncpu)/2,k)
	  enddo
 	    b(1)=b(1)+a(1) !neumann
	    b(nx-2)=b(nx-2)+c(nx-2) !neumann
c	  call tridiag(a,b,c,r,u,nx-2,nx-2)
	  call tridag(a,b,c,r,u,nx-2,nx-2)
	  do i=1,nx-2
	    dummy2(j,i,k)=u(i)
	  enddo
c     imaginary part	
          do i=1,nx-2
            r(i)=fac*dummy2(j+1,i,k)
	  enddo
c	  call tridiag(a,b,c,r,u,nx-2,nx-2)
	  call tridag(a,b,c,r,u,nx-2,nx-2)
	  do i=1,nx-2
	    dummy2(j+1,i,k)=u(i)
	  enddo
	enddo
c now the oddball
       if (me.eq.0) then 
          do i=1,nx-2
            r(i)=fac*oddball(2*k-1,i)
	    b(i)=-(a(i)+c(i))+kabs(ny/2,k)
	  enddo
	  b(1)=b(1)+a(1) !neumann
	  b(nx-2)=b(nx-2)+c(nx-2) !neumann
	  call tridag(a,b,c,r,u,nx-2,nx-2)
	  do i=1,nx-2
	    oddball(2*k-1,i)=u(i)
	  enddo
c  imaginary part
          do i=1,nx-2
            r(i)=fac*oddball(2*k,i)
	  enddo
	  call tridag(a,b,c,r,u,nx-2,nx-2)
	  do i=1,nx-2
	    oddball(2*k,i)=u(i)
	  enddo
       endif
      enddo
c  0,0 mode
        if (me.eq.0) then 
	 do i=1,nx-2
          r(i)=fac*zero(2*i-1)
          b(i)=-(a(i)+c(i))
        enddo
          b(1)=b(1)+a(1) !neumann
          b(nx-2)=b(nx-2)+c(nx-2) !neumann
c
	  u(nx-2)=0.
          call tridag(a(1),b(1),c(1),r(1),u(1),nx-3,nx-3)
c	  u(2)=(r(1)-b(1)*u(1))/c(1)
c	  do i=3,nx-2
c	   u(i)=(r(i-1)-a(i-1)*u(i-2)-b(i-1)*u(i-1))/c(i-1)
c          enddo
	  do i=1,nx-2
	    dummy2(1,i,1)=u(i)
	  enddo
c     imaginary part
	 do i=1,nx-2
          r(i)=fac*zero(2*i)
        enddo
	  u(nx-2)=0.
c	  u(2)=(r(1)-b(1)*u(1))/c(1)
c	  do i=3,nx-2
c	   u(i)=(r(i-1)-a(i-1)*u(i-2)-b(i-1)*u(i-1))/c(i-1)
c          enddo
          call tridag(a(1),b(1),c(1),r(1),u(1),nx-3,nx-3)
	  do i=1,nx-2
	    dummy2(2,i,1)=u(i)
	  enddo
	 endif

c     fft back      
c     call inv fft in z (all modes except oddball)	
      do j=1,(ny-2)/ncpu,2
        call cfft99(dummy2(j,1,1),work,trigsz,ifaxz,(nx-2)*(ny-2)/ncpu/2,(ny-2)/ncpu/2,nzg-2,nx-2,1)
      enddo	  
c     now the oddball
      if (me.eq.0)
     &   call cfft99(oddball(1,1),work,trigsz,ifaxz,1,(nzg-2),nzg-2,nx-2,1)
c     swap x and y dimensions before calling swap
      do k=1,nzg-2
        do i=1,nx-2
	  do j=1,(ny-2)/ncpu
	    dummyt(i,j,k)=dummy2(j,i,k)
	  enddo
	enddo
      enddo	      
      if (me.eq.0) then
        do k=1,nzg-2
	 do i=1,nx-2
	  dummyt(i,2,k)=(oddball(2*k-1,i))
	 enddo
	enddo  
      endif	
c       
c     swap back z and x directions
      call swap_zx(dummyt,dummy,(ny-2)*(nx-2),nzg-2,ncpu,me,comm)
c
      do k=2,nz-1
       do i=1,nx-2
        do j=1,ny-2
         ftemp(j,i)=dummy(i,j,k-1)
	enddo
	 ftemp(ny-1,i)=dummy(i,2,k-1)
	 ftemp(ny,i)=0.
         ftemp(2,i)=0.
       enddo
       call fft99(ftemp,work,trigsy,ifaxy,1,ny,ny-2,nx-2,1)
       
c     now take the transform
c      swap x and y direction in a dummy array
       do i=1,nx-2
        do j=1,ny
	 dp(i+1,j,k)=ftemp(j,i)
	enddo
       enddo    	 
      enddo

c     apply boundary conditions to dp at ghost points
      do k=2,nz-1
         do i=2,nx-1
            dp(i,ny,k)=dp(i,2,k)
c            p(i,ny,k)=p(i,2,k)
            dp(i,1,k)=dp(i,ny-1,k)
c            p(i,1,k)=p(i,ny-1,k)
         enddo
      enddo
      do k=2,nz-1
         do j=2,ny-1
            dp(nx,j,k)=dp(nx-1,j,k)
c            p(nx,j,k)=p(nx-1,j,k)
            dp(1,j,k)=dp(2,j,k)
c            p(1,j,k)=p(2,j,k)
         enddo
      enddo
      if (ncpu.eq.1) then 
        do j=2,ny-1
           do i=2,nx-1
              dp(i,j,1)=dp(i,j,nz-1)
              dp(i,j,nz)=dp(i,j,2)
           enddo
        enddo
      else	
        call refreshbc(dp,nx*ny,nz-2,ncpu,me,comm)
c        call refreshbc(p,nx*ny,nz-2,ncpu,me,comm)
c        now enforce periodic boundary conditions in z     
        if (me.eq.0) then
          call MPI_SEND(dp(1,1,2),nx*ny,datatype,ncpu-1,1,comm,ierr)
c          call MPI_SEND(p (1,1,2),nx*ny,datatype,ncpu-1,2,comm,ierr)
          call MPI_RECV(dp(1,1,1),nx*ny,datatype,ncpu-1,3,comm,status,ierr)
c          call MPI_RECV(p (1,1,1),nx*ny,datatype,ncpu-1,4,comm,status,ierr)
        else if (me.eq.ncpu-1) then  	
          call MPI_RECV(dp(1,1,nz  ),nx*ny,datatype,0,1,comm,status,ierr)
c          call MPI_RECV(p (1,1,nz  ),nx*ny,datatype,0,2,comm,status,ierr)
          call MPI_SEND(dp(1,1,nz-1),nx*ny,datatype,0,3,comm,ierr)
c          call MPI_SEND(p (1,1,nz-1),nx*ny,datatype,0,4,comm,ierr)
        endif
      endif !ncpu	
c     update pressure
      do i=1,nx
       do j=1,ny      
        do k=1,nz
	 p(i,j,k)=p(i,j,k)+dp(i,j,k)
	enddo
       enddo
      enddo
c       
      return
      end
c     