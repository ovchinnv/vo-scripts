c.......................................................................
c
      subroutine PresDirect(divv,dp,p,xu,yv,zw,au,cw,icycle,
     &                      myrank,mpi_comm_hydro,time)

      implicit none
c
      include 'headers/dimension.h'
c
c
c     ...Parameters...
c
      integer mysize,myrank,n,m,kkk,nwk
      parameter (mysize=ncpu)
c
c     ...Scalar Arguments...
c
      integer  ii,kk
      integer  icycle
c
c     ...Array Arguments...
c
      real     xu(nx), yv(ny), zw(nzg)
      real     au(nx), cw(nzg)
      real     dp(nx,ny,(nzg-2)/mysize+2),
     &     divv(nx,ny,(nzg-2)/mysize+2), 
     &     p(nx,ny,(nzg-2)/mysize+2)
      integer mpi_comm_hydro
      real time(3)
      

      real*8 cpu1,cpu2,cpu3
c
c     ...External Subroutine...
c
      external norm
c
c     ...Intrinsic Functions...
c
      intrinsic int, log10, max, real
c
c     ...Local Scalars...
c
      integer          bcpn, bcps, bcpw, bcpe
      integer          i, iflg, j, jl, k
      integer          l, ierror
      integer          nzl,nyl
      double precision dely
      double precision dx, dz, xc, zc, lx, lz, mean, meane
      integer          ierr, kerr
      double precision eps, sol, err2, errm, diff
      double precision normre, normx, normb
c
c     ...Local Arrays...
c
      double precision ak(ny),am(nx-2),bm(nx-2),cm(nx-2)
      double precision an(nzg-2),bn(nzg-2),cn(nzg-2)
      double precision xr(ny,nx-2), work(ny-1,nx-2)
      double precision dummy(nx-2,ny-2,(nzg-2)/mysize),dummyt(nx-2,(ny-2)/mysize,nzg-2)
      double precision w(40060)
      double precision rhsp(nx-2,nzg-2)
      double precision bb(nx-2,nzg-2)

      double precision exact(nx-2,nzg-2)
c
c     ...Local Arrays for fft...
c
      integer ifax(13)
      double precision trigx(3*(ny-2)/2+1)
c
c     ...Local variables for NAG routine ...
c
      integer          n1, n2, n1m, itmax
      parameter (itmax=500)
      integer          ndir, ixn, iyn, ifail, itcoun, itused
      double precision aparam, conres, conchn
      double precision a(nx-2,nzg-2),b(nx-2,nzg-2),c(nx-2,nzg-2),d(nx-2,nzg-2),e(nx-2,nzg-2)
      double precision q(nx-2,nzg-2)
      double precision wrksp1(nx-2,nzg-2),wrksp2(nx-2,nzg-2),wrksp3(nx-2,nzg-2)
      double precision resids(itmax), chngs(itmax)

      double precision modeOavg
C
C     ...Save statement...
C
      save aparam, ndir, ixn, iyn, conres, conchn
c.......................................................................
c     Executable Statements
c.......................................................................
C.......................................................................
C allocate variables unsaved local variables (deallocated at the end!)
C.......................................................................
      nzl=(nzg-2)/mysize+2
      nyl=(ny-2)/mysize+2
      
      n1=nx-2
      n2=nzg-2
      n1m=n1
c      itmax=100
      if((mod(n2,mysize).gt.0).or.(mod(ny-2,mysize).gt.0)) goto 997
c
cTHIS IS THE F77 VERSION!      
c      allocate(xr(ny,nx-2),work(ny-1,nx-2))
c      allocate(rhsp(nx-2,nzg-2))
c      allocate(bb(nx-2,nzg-2))
c      allocate(exact(nx-2,nzg-2))
c      allocate(q(n1m,n2),wrksp1(n1m,n2),wrksp2(n1m,n2))
c      allocate(wrksp3(n1m,n2),resids(itmax),chngs(itmax))
c      allocate(dummy(nx-2,ny-2,nzl-2))
c      allocate(dummyt(nx-2,nyl-2,nzg-2))
c      allocate(trigx(3*(ny-2)/2+1))
      

c.......................................................................
c     Initialize the damn thing
c.......................................................................
      if (icycle .eq. 1) then
c.......................................................................
c     Initialize direct solver
c.......................................................................
c     Dimension of array for blktri work space and saved arrays
c.......................................................................
         n = nzg-2
         m = nx-2
         k = int(log10(real(n))/log10(real(2)))+1
         l = 2**(k+1)
         nwk =(k-2)*l+k+5+max(2*n,6*m) 
         if (nwk.gt.40060) then 
	   write(6,*) 'increase size of work array "w" to ',nwk
	   stop
	 endif
c         allocate(w(nwk))
c         allocate(ak(ny),am(nx-2),bm(nx-2),cm(nx-2))
c         allocate(an(nzg-2),bn(nzg-2),cn(nzg-2))
c         allocate(a(n1m,n2),b(n1m,n2),c(n1m,n2),d(n1m,n2),e(n1m,n2))
      endif
      do k = 1, nzl-2
         do j = 1, ny-2
            do i = 1, nx-2
               dummy(i,j,k) = 0.0
            end do
         end do
      end do
      call peqina(xu,zw,au,cw,am,bm,cm,an,bn,cn,nx,ny,nzg)
c     
c     Set boundary conditions
c
      bcpn = 1
      bcps = 1
      bcpe = 1
      bcpw = 1
      dely  = yv(3)-yv(2)
      call peqinb(divv,dp,am,bm,cm,an,bn,cn,ak,dely,
     &     bcpn,bcps,bcpw,bcpe,nx,ny,nzg)

c     
      if(icycle.eq.1) then
c     blktri needs to be called once to initialize the working array!
         iflg = 0
         call blktri(iflg,1,nzg-2,an,bn,cn,1,nx-2,am,bm,cm,nx-2,rhsp,
     *        ierror,w)
c.......................................................................
c     Initialization for NAG routine 
c.......................................................................
         do k=1, nzg-2
            do i=1, nx-2
               a(i,k) = an(k)
               b(i,k) = am(i)
               c(i,k) = bm(i) + bn(k)
               d(i,k) = cm(i)
               e(i,k) = cn(k)
            end do
         end do 
c     iteration acceleration factor
         aparam = ((n1-1)**2+(n2-1)**2)/3.0d0 
         if (bcpn*bcps*bcpw*bcpe .eq. 1) then
c     An arbitrary constant can be added to the solution.
            ndir   = 0             
            ixn    = n1/2       ! i at which the solution is to be set to zero.
            iyn    = n2/2       ! k at which the solution is to be set to zero.
         else
c     The system of equations has a unique solution.
            ndir   = 1
            ixn    = 0
            iyn    = 0
         end if
         conres = 1.0e-7        ! convergence criterion
         conchn = 10.*conres    ! convergence criterion
      endif
c      	 
c     Skip fft if in 2D
      if (ny.eq.3) then 
c     fill dummy directly
        do k=1, nzl-2
         do j = 1, ny-2 ! i.e. from 1 to 1
            do i = 1, nx-2
               dummy(i,j,k) = divv(i+1,j+1,k+1)
            end do 
         end do  
        end do
      else 
c.......................................................................
c     Initialize fft
c.......................................................................
      if (icycle.eq.1) call fftfax_4p(ny-2,ifax,trigx)
c.......................................................................
c.......................................................................
c     Prepare for fft; apply fft in the y direction
c.......................................................................
c      cpu1=mpi_wtime()
       do k = 2, nzl-1
         do i = 2, nx-1
            do j = 2, ny-1
               xr(j,i-1) = divv(i,j,k)
            end do
            xr(1 ,i-1) = xr(ny-1,i-1)
            xr(ny,i-1) = xr(2   ,i-1)
         end do
c     
         call fft99_4p(xr,work,trigx,ifax,1,ny,ny-2,nx-2,-1)
c.......................................................................
c     Fast Fourier Transform of the RHS  stored in dummy(i,j,k).
c     Each couple j,j+1 corrsponds to the real (j odd) and imaginary 
c     (j even) part of the mode j/2 (j=1 => mode 0)
c     Because  j=2 corresponds to the imaginary part of the mode 0, 
c     it is identically zero. Similarly, j=ny is zero. Therefore
c     we replace j=2 with the content at j=ny-1 (the nyquist freq.)
c.......................................................................
         do j = 1, ny-2
            do i = 1, nx-2
               dummy(i,j,k-1) = xr(j,i)
            end do 
         end do  
         do i=1,nx-2
c            dummy(i,2,k-1)=xr(ny-1,i)
         enddo
       end do
      endif !2D
c.......................................................................
c     Swap
c.......................................................................
      call swap_xz(dummy,dummyt,(nx-2)*(ny-2),
     &     nzg-2,mysize,myrank,MPI_COMM_HYDRO)
c.......................................................................
c     From now on variables are in the Fourier space
c.......................................................................
c     Invert the system of equations for each l and determine the 
c     unknown Fourier coefficients of the transformed pressure. 
c     The coefficients are ny-2, each determined by the composition 
c     of a couple of the values stored in dummy, so that:
c     Qr(0)          = dummy(1)+I*dummy(2)=dummy(1)
c     Qr(1)          = dummy(3)+I*dummy(4)
c     ...
c     Qr((ny-2)/2-1) = dummy((ny-2)-1)+I*dummy((ny-2)-0)
c     Qr((ny-2)/2)   = dummy(ny-1)+I*0
c     Qr((ny-2)/2+1) = dummy((ny-2)-2)-I*dummy((ny-2)-1)
c     ...
c     Qr((ny-2)-1)   = dummy(3)-I*dummy(4)
c     
c     By definition, the second (j=2) and the last (j=ny) are zero.
c.......................................................................
c     cpu3=mpi_wtime()
c      
      do jl = 1, nyl-2          ! local index

         j=(ny-2)/mysize*myrank+jl ! global index
         
c     
         if (j .eq. 1) then     ! mode zero, special treatment
c.......................................................................
c         First mode
c.......................................................................
c     Assign the right side of equation

            do k = 1, nzg-2
               do i = 1, nx-2
                  rhsp(i,k) = dummyt(i,jl,k)
                  bb(i,k)   = rhsp(i,k) ! temporary
               end do
            end do
      

c.......................................................................
c     Call the direct solver for first attempt solution
c.......................................................................
            iflg = 1 
            call blktri(iflg,1,nzg-2,an,bn,cn,1,nx-2,am,bm,cm,nx-2,rhsp,
     *           ierror,w)
      goto 234 !VO 7/15/04
            if (ierror .ne. 0) write(*,*) j,ierror
c.......................................................................
c     Check it out !
c.......................................................................
c            call norm(am,bm,cm,an,bn,cn,
c     $           rhsp,bb,nx,nzg,normx,normb,normre)
c            write(*,2000) j, normb, normx, normre, normre/(normb+1.e-12)


c.......................................................................
c     Call the iterative solver to improve the solution
c.......................................................................
            do k=1, nzg-2
               do i=1, nx-2
                  q(i,k) = dummyt(i,jl,k)
               end do
            end do
            ifail  = 1
c     iteration counter
            itcoun = 0      
            call d03ebf(n1,n2,n1m,a,b,c,d,e,q,rhsp,aparam,
     *           itmax,itcoun,itused,ndir,ixn,iyn,
     *           conres,conchn,resids,chngs,
     *           wrksp1,wrksp2,wrksp3,ifail)
c            write (*,99999) (i,resids(i),chngs(i),i=1,itused,50)
c            write (*,99999) (i,resids(i),chngs(i),i=itused,itused)
99999       format (' iter, residual, change',2x,i4,10x,1pd11.4,4x,
     $           1pd11.4)
c     Check error flag
            if (ifail.ne.0 .and. ifail.ne.5) then
               write (*,99997) 'Error in D03EBF  ifail =', ifail
            end if
            if (ifail.eq.5) then
               write (0,99997)'0th mode for p: 
     *              no convergence after itmax', 
     *              itmax,'; it could be good enough ...'
            end if
99997       format (1x,a,1x,i4,1x,a)
c.......................................................................
c     Check it out !
c.......................................................................
c            call norm(am,bm,cm,an,bn,cn,
c     $           rhsp,bb,nx,nzg,normx,normb,normre)
c            write(*,2000) j, normb, normx, normre, normre/(normb+1.e-12)


c     call norm(am,bm,cm,an,bn,cn,rhsp,bb,nx,nzg,normx,normb,normre)
c     write(*,2000) j, normb, normx, normre, normre/(normb+1.e-12)
c.......................................................................
c     Dump the solution into dummy
c.......................................................................
 234        continue
c      goto 345
            modeOavg=0.
            do k=1,nzg-2
	     do i=1,nx-2
	      modeOavg=modeOavg+rhsp(i,k)
	     enddo
	    enddo
	    modeOavg=modeOavg/(nx-2)/(nzg-2)  

            do k = 1, nzg-2
               do i = 1, nx-2
                  dummyt(i,jl,k) = rhsp(i,k)-modeOavg
               end do
            end do 
 345   continue	    
c.......................................................................
c     There are still (ny-2)/2-2 couples of (real and imag)
c     Fourier coefficients to determine. 
c.......................................................................
         elseif(j .eq.2 ) then  ! this is the Nyquist freq.!
c.......................................................................
c     Blktri for other modes.
c     Change the diagonal coefficients for j=3,..., ny-1.
c.......................................................................
            l = (ny-2)/2
            do i = 1, nx-2
               bm(i) = bm(i) - ak(l)
            end do
c     Assign the RHS
            do k = 1, nzg-2
               do i = 1, nx-2
                  rhsp(i,k) = dummyt(i,jl,k)
                  bb(i,k)   =  rhsp(i,k)              ! temporary
               end do
            end do
c.......................................................................
c     Call the solver
c.......................................................................
            iflg = 1 
            call blktri(iflg,1,nzg-2,an,bn,cn,1,nx-2,am,bm,cm,nx-2,rhsp,
     *           ierror,w)
            if (ierror .ne. 0) write(*,*) j,ierror
c.......................................................................
c     Check it out !
c.......................................................................
c           call norm(am,bm,cm,an,bn,cn,rhsp,bb,nx,nzg,normx,normb,normre)
c            write(*,2000) j, normb, normx, normre, normre/normb
c.......................................................................
c     Fix the diagonal coefficient
c.......................................................................
            do i = 1, nx-2
               bm(i) = bm(i) + ak(l)
            end do
c.......................................................................
c     Dump the solution into dummy
c.......................................................................
            do k=1,nzg-2
               do i=1,nx-2
                  dummyt(i,jl,k) = rhsp(i,k)
               end do
            end do
         elseif(j .gt. 2) then  ! all the others
c.......................................................................
c     Blktri for other modes.
c     Change the diagonal coefficients for j=3,..., ny-1.
c.......................................................................
            l = (j-1)/2

            do i = 1, nx-2
               bm(i) = bm(i) - ak(l)
            end do
c     Assign the RHS
            do k = 1, nzg-2
               do i = 1, nx-2
                  rhsp(i,k) = dummyt(i,jl,k)
                  bb(i,k)   =  rhsp(i,k)              ! temporary
               end do
            end do
c.......................................................................
c     Call the solver
c.......................................................................
            iflg = 1 
            call blktri(iflg,1,nzg-2,an,bn,cn,1,nx-2,am,bm,cm,nx-2,rhsp,
     *           ierror,w)
            if (ierror .ne. 0) write(*,*) j,ierror
c.......................................................................
c     Check it out !
c.......................................................................
c           call norm(am,bm,cm,an,bn,cn,rhsp,bb,nx,nzg,normx,normb,normre)
c            write(*,2000) j, normb, normx, normre, normre/normb
c.......................................................................
c     Fix the diagonal coefficient
c.......................................................................
            do i = 1, nx-2
               bm(i) = bm(i) + ak(l)
            end do
c.......................................................................
c     Dump the solution into dummy
c.......................................................................
            do k=1,nzg-2
               do i=1,nx-2
                  dummyt(i,jl,k) = rhsp(i,k)
               end do
            end do
            
         end if
      end do
c.......................................................................
c     swap back
c.......................................................................
c      cpu2=mpi_wtime()
       call swap_zx(dummyt,dummy,(nx-2)*(ny-2),
     &     nzg-2,mysize,myrank,MPI_COMM_HYDRO)
c
c     Skip fft if in 2D
      if (ny.eq.3) then 
c     fill dp from dummy directly
        do k=1, nzl-2
         do j = 1, ny-2 ! i.e. from 1 to 1
            do i = 1, nx-2
               dp(i+1,j+1,k+1) = dummy(i,j,k)
            end do 
         end do  
        end do
      else 
c      
c      time(3)=mpi_wtime()-cpu3
c      time(2)=time(2)+(mpi_wtime()-cpu2)
      


c     Transform back into physical space.
c     There are ny coefficients of the Fourier transformed solution. 
c.......................................................................
c      cpu1=mpi_wtime()
       do k = 1, nzl-2
        do i = 1, nx-2
          do j = 1, ny-2
            xr(j,i) = dummy(i,j,k)
          end do
          xr(2,i)=0.
c          xr(ny-1,i)=dummy(i,2,k)
          xr(ny,i)=0.
        end do
        call fft99_4p(xr,work,trigx,ifax,1,ny,ny-2,nx-2,+1)
c.......................................................................
c       Fill of the solution in the interior points
c.......................................................................
        do i = 1, nx-2
          do j = 1, ny-2
            dp(i+1,j+1,k+1) = xr(j+1,i)
          end do
        end do
       end do
c      time(1)=time(1)+(mpi_wtime()-cpu1)
      endif
C.......................................................................
C     Update pressure field
C.......................................................................
      do k = 2, nzl-1
        do j = 2, ny-1
          do i = 2, nx-1
            p(i,j,k) = p(i,j,k) + dp(i,j,k)
          end do
        end do
      end do
c.......................................................................
c     Cyclic b.c. in y direction
c.......................................................................
      do k = 1, nzl
        do i = 1, nx
          p(i,1,k)   = p(i,ny-1,k)
          p(i,ny,k)  = p(i,2,k)
          dp(i,1,k)  = dp(i,ny-1,k)
          dp(i,ny,k) = dp(i,2,k)
        end do
      end do

      call refreshbc(dp,nx*ny,nzl-2,mysize,myrank,MPI_COMM_HYDRO)
      call refreshbc(p,nx*ny,nzl-2,mysize,myrank,MPI_COMM_HYDRO)
c      deallocate(xr,work,rhsp,exact,q)
c      deallocate(wrksp1,wrksp2,wrksp3)
c      deallocate(dummy,dummyt)
c      deallocate(resids,chngs)
      return
c
 997  if(myrank.eq.0) then
         write(6,998) 6,nzg-2,ny-2,mysize
      endif
 998  format('ierror=',I2,';',i3,' or',i3,' Not divisible by', i3)
      call mpi_finalize()
      stop


 1980 format(10f10.4)
 1997 format(1p8e11.3)
 1999 format(10f10.4)
 2000 format(1x,'j=',i2,1x,'nob=',e12.5,1x,'nox=',e12.5,1x,
     *       'nore=',e12.5,1x,'nore/nob=',e12.5)
c
      return
      end
c
c.......................................................................
c.......................................................................
c
      subroutine norm(am, bm, cm, an, bn, cn, x, b, nx, nzg,
     *                normx, normb, normr)
c
c.......................................................................
c
c     ...Scalar Arguments...
c
      integer  nx, nzg
      double precision normx, normb, normr
c
c     ...Array Arguments...
c
      double precision  b(nx-2,nzg-2)
      double precision  x(nx-2,nzg-2)
c
c     ...Intrinsic Functions...
c
      intrinsic sqrt
c
c     ...Local Scalars...
c
      integer i, k
      double precision res
c
c     ...Array Arguments...
c
      double precision am(nx-2), bm(nx-2), cm(nx-2)
      double precision an(nzg-2), bn(nzg-2), cn(nzg-2)

      normx = 0.0
      normb = 0.0
      normr = 0.0
      do k = 1, nzg-2
        do i = 1, nx-2
          res = bm(i)*x(i,k) + bn(k)*x(i,k) 
          if (i.gt.1)    res = res + am(i)*x(i-1,k )
          if (i.lt.nx-2) res = res + cm(i)*x(i+1,k )
          if (k.gt.1)    res = res + an(k)*x(i ,k-1)
          if (k.lt.nzg-2) res = res + cn(k)*x(i ,k+1)
          res = res - b(i,k)
          normx = normx + x(i,k)*x(i,k)
          normb = normb + b(i,k)*b(i,k)
          normr = normr + res*res
        end do
      end do
      normx = sqrt(normx/(nx-2)/(nzg-2))
      normb = sqrt(normb/(nx-2)/(nzg-2))
      normr = sqrt(normr/(nx-2)/(nzg-2))

      return
      end
c
c.......................................................................
C#######################################################################
C     BEGIN TEMPORARY STUFF
C
C     DO K = 1, nzg-2
C       DO I = 1, NX-2
C         DO J = 1, NY
C           XR(J,I) = DUMMY(I,J,K)
C         END DO
C       END DO
C       CALL FFT99_4P(XR,WORK,TRIGX,IFAX,1,NY,NY-2,NX-2,+1)
C       DO I = 1, NX-2
C         DO J = 1, NY-2
C           DP(I+1,J+1,K+1) = XR(J+1,I)
C         END DO
C       END DO
C     END DO
C
C     K = 10 
C     I = 10
C     DO J = 1, NY
C       WRITE(*,*)  J, DIVV(I,J,K), DP(I,J,K) 
C     END DO
C
C     STOP
C
C     END TEMPORARY STUFF
C#######################################################################
C#######################################################################
C         BEGIN TEMPORARY STUFF
C         EPS   = 1.E-8
C         LX = XU(NX-1)-XU(1)
C         LZ = ZW(nzg-1)-ZW(1)
C         DO K = 1, nzg-2
C           ZC = 0.5*(ZW(K+1)+ZW(K))
C           DO I = 1, NX-2
C             XC = 0.5*(XU(I+1)+XU(I))
C             EXACT(I,K) = XC*XC + ZC*ZC 
C           END DO
C         END DO
C         DO K = 1, nzg-2
C           DO I = 1, NX-2
C             RHSP(I,K) = BM(I)*EXACT(I,K) + BN(K)*EXACT(I,K) 
C             IF(I.GT.1)    RHSP(I,K) = RHSP(I,K) + AM(I)*EXACT(I-1,K) 
C             IF(I.LT.NX-2) RHSP(I,K) = RHSP(I,K) + CM(I)*EXACT(I+1,K)
C             IF(k.GT.1)    RHSP(I,K) = RHSP(I,K) + AN(K)*EXACT(I,K-1)
C             IF(K.LT.nzg-2) RHSP(I,K) = RHSP(I,K) + CN(K)*EXACT(I,K+1)
C           END DO
C         END DO
C         DO K=1, nzg-2
C            WRITE(7,1000) K, AN(K),BN(K),CN(K)
C         END DO
C         DO I=1, NX-2
C            WRITE(8,1000) I, AM(I),BM(I),CM(I)
C         END DO
C 1000    FORMAT(1X,I3,3X,3E14.7)
C         CALL BLKTRI(1,1,nzg-2,AN,BN,CN,1,NX-2,AM,BM,CM,NX-2,RHSP,
C    *                IERROR,W)
C         IF (IERROR .NE. 0) WRITE(*,*) 'ERROR IN BLKTRI:',J,IERROR
C         DIFF = RHSP(10,1) - EXACT(10,1)
C         DO K = 1, nzg-2
C           DO I = 1, NX-2
C             RHSP(I,K) = RHSP(I,K) - DIFF
C           END DO
C         END DO
C
C         COMPUTE ERROR 
C
C         SOL  = 0.0
C         ERR2  = 0.0
C         ERRM = -999
C         DO K = 1, nzg-2
C           DO I = 1, NX-2
C             SOL  = SOL + RHSP(I,K)**2
C             ERR2 = ERR2 + (RHSP(I,K)-EXACT(I,K))**2
C             IF (ABS(RHSP(I,K)-EXACT(I,K)) .GT. ERRM) THEN
C               IERR = I
C               KERR = K
C               ERRM = ABS(RHSP(I,K)-EXACT(I,K))
C             END IF
C           END DO
C         END DO
C         ERR2 = SQRT(ERR2/(NX-2)/(nzg-2))
C         SOL  = SQRT(ERR2/(NX-2)/(nzg-2))
C
C         WRITE(*,*)
C         WRITE(*,*) ' L2  error ',err2,' normalized ',err2/(sol+eps)
C         WRITE(*,*) ' max error ',errm,' normalized ',errm/(sol+eps)
C         WRITE(*,*) ' at        ',ierr,kerr
C
C         I = 10
C         DO K = 1, nzg-2
C          WRITE(*,*)  K, RHSP(I,K), EXACT(I,K)
C         END DO
C
C         STOP
C
C         END TEMPORARY STUFF
C#######################################################################
