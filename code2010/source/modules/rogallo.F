       real xref_,x_,dref, L_,L11,stretch,Re_,zlen
       integer nose,kose
       parameter (xref_=50.) !reference location
!       parameter (x_=15.17) !location of the inflow w/r to the reference location
       parameter (dref=0.)  ! disp. thickness at reference location
       parameter (nose=1000) ! number of points for OSE discretization
       parameter (kose=64/2) ! resolution
       parameter (stretch=1.00) ! grid stretching (for ose routine)
       parameter (zlen=30.0) ! side length of noise box
       
       common /ose/ L_, Re_
       

!       parameter (xref_=0.)
!       parameter (x_=32.59)
!       parameter (dref=0.)

!       parameter (x_=-10.) ! for WBL calculation only
       parameter (x_=-50.) 


        real zstart, zend ! locations of inflection points of the noise cutoff functions
        real shape1, shape2
	parameter (zstart=1,zend=45)
 	parameter (shape1=7.0,shape2=7.0) ! sharpness of htan jump at above locations
	
	integer nmodes   ! number of modes in the representation of rogallo-type inflow
	parameter (nmodes=2000)
	


!     OSE/ROGALLO FUNCTIONALITY FLIP CASE
      include 'headers/ose.h'
!      arrays for isotropic inflow     
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      real uose(ny,nx),vose(ny,nx),wose(ny,nx)
      real cutoffu(nx),cutoffw(nx),shape
      integer jcounter, count(nmodes), is(nmodes,kose) 
!     other vars for isotropic inflow
      integer seed,kkmax ! kose defined in ose.h
!      
!cccccccccccccccccccccccccccccccccccccccccccccccccc
!     GRID DATA
      real eta(nose),detacen(nose),detacor(nose)
!     OSE INTERPOLATION DATA
      real a_w(nx),b_w(nx),a_u(nx),b_u(nx)
      integer ps_w(nx),ps_u(nx),nxloc
!cccc BLASIUS VARIABLES      
      complex U_(nose),d2U(nose)
      real    Ur(2*nose),d2Ur(2*nose)
!cccc OSE VARIABLES
      complex phi(nose),phip(nose)
      real phir(2*nose),phipr(2*nose)
!cccc 







! 
      real uosm,vosm,wosm ! for rescaling inflow noise on the fly
! 







      if (jdose.and.rogal) then 
       write(0,*) 'error: jd-ose and rogallo noise are both on!'
       return
      endif
      if (box2.and.(jdose.or.rogal)) then 
       write(0,*) 'error: box2 and jd-ose (or rogallo) noise are both on!'
       return
      endif
!       
!      if (rogal) call rogal_new_setup(jcounter,kx,ky,kz,uh,vh,wh,cutoffu,cutoffw,uscale,xu,xcf,yv)
      if (rogal) call rogal_old_setup(jcounter,kx,ky,kz,count,is,uh,vh,wh,cutoffu,cutoffw,uscale,xu,xcf,yv)
      if (jdose) call jdose_new_setup(jcounter,kx,ky,kz,uh,vh,wh,cutoffu,cutoffw,uscale,eta,detacor,d2U,U_,Ur,xu,xcf,yv,
     &                                a_w,b_w,a_u,b_u,ps_w,ps_u)
      if (jdose.or.rogal) then 
         open(101,file='pscale.dat',form='formatted',status='old')
	   read(101,*) puscale, pvscale, pwscale
	   read(101,*) psample
	 close(101)  
!      target turbulence level
         plevel=u_rms ! u_rms is in the input file
	 pind=1
	 do while (zwg(pind).le.u_rms_loc)
	   pind=pind+1
	 enddo
      endif	

      subroutine grid_ose(y,dycen,dycor,ny,len,stretch,L)
      implicit none
c
      integer i,ny
      real y(ny),dycen(ny),dycor(ny),dy,scale,L,len,leny,stretch
c      write(6,*) 'L=',L,'Re=',Re
      leny=len/L
c      write(6,*) L_,leny,stretch
cccc  generate exponential grid ccccccccccccccccccccccccccccc
      y(1)=0.
      dy=stretch
cccc  
      do i=2,ny
         y(i)=y(i-1)+dy
         dy=dy*stretch
      enddo
      scale=leny/y(ny)
c      write(6,*) scale
      do i=2,ny
         y(i)=y(i)*scale
      enddo
cccc  create metrics
      do i=1,ny-1
         dycen(i)=1./(y(i+1)-y(i))
      enddo
      dycen(ny)=dycen(ny-1)
      do i=2,ny-1
         dycor(i)=1./(0.5*(y(i+1)-y(i-1)))
      enddo
      dycor(1)=dycen(1)
      dycor(ny)=dycen(ny-1)
c     
      open(1,file='grid.dat',form='formatted',status='unknown')
      do i=1,ny
         write(1,'(I3, 3G30.20)') i,y(i),dycen(i),dycor(i)
      enddo
      close(1)
      return
      end
ccccccccccccccccccccccccccccccccccccccccc



      subroutine setinterp_ose(source,target,ns,nt,a,b,pos,L)
      implicit none
      integer ns,nt
      real source(ns), target(nt),a(nt),b(nt),L
      integer pos(nt)
      integer i,ii
c      write(6,*) L*source(1),L*source(2),target(1)
      if (abs(target(nt)-source(ns)).lt.1E-10) target(nt)=source(ns) ! correct for almost equal domains
c      if (L*source(1).gt.target(1).or.L*source(ns).lt.target(nt)) stop 'domain mismatch...aborting'
      ii=1
      do i=1,nt
       do while (.not.(target(i).ge.L*source(ii).and.target(i).le.L*source(ii+1)))

          ii=ii+1
          if (ii.gt.ns-1) then 
	    ii=ii-1
c	  stop "please debug 'setinterp'"
            write(6,*) 'setinterp: warning: target domain not within source domain'
	    write(6,*) target (i),ii,source(ii)*L
	    a(i)=0.
	    b(i)=0.
	    pos(i)=1
	    goto 99
	  endif
       enddo   
       pos(i)=ii
       a(i)=(L*source(ii+1)-target(i))/(source(ii+1)-source(ii))/L
       b(i)=1.-a(i)
 99   enddo     
      return
      end



c     DO NOT MODIFY! 
      subroutine setup_ose(y,dycor,d2U,U,Ur,ny,nu) ! nu is given as ru1 in the main 
      implicit none
      include 'headers/ose.h'
c     
      integer ny
      real V(ny)
      complex U(ny)
      complex d2U(ny)
      real y(ny),dycor(ny)
      real Ur(2*ny),dU(ny),nu
      integer i
c     
ccccccccccccccccccccccccccccccccccccccccc
c     write(6,*) L, Uinf, nu, Re
      call BLASIUS_ose(y,Ur,V,dU,ny,ny,ny,nu,1.,x_,xref_,dref)
c     spread data in Ur to make complex
      do i=ny,1,-1
         U(i)=cmplx(Ur(i),0)
      enddo
c     compute 2nd derivative of U (at cell corner)
      do i=2,ny-1
         d2U(i)=cmplx((dU(i+1)-dU(i))*dycor(i),0)
      enddo
c     debug   write(2,*) (i,y(i),dble(U(i)),dU(i),(U(i+1)-U(i))*dycen(i),char(10),i=1,ny-1)
      return
      end




      subroutine reshuffle(jcounter,uh,vh,wh,ky,kz,seed)
      implicit none
      include 'headers/ose.h'
c      
      integer jcounter,kk,rnd
      real ky(nmodes), kz(nmodes)
      complex uh(nmodes), vh(nmodes), wh(nmodes)
      real ran1
      integer seed
c      
      write(0,*) 'reshuffling phases'
      do kk=1,jcounter
       rnd=2*ran1(seed)
       uh(kk)=rnd*uh(kk)+(1-rnd)*conjg(uh(kk))
       vh(kk)=rnd*vh(kk)+(1-rnd)*conjg(vh(kk))
       wh(kk)=rnd*wh(kk)+(1-rnd)*conjg(wh(kk))
c
c       rnd=2*rnd-1       ! get -1 or 1
c       ky(kk)=rnd*ky(kk)
c       kz(kk)=rnd*kz(kk)
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reshuffle_old(jcounter,uh,vh,wh,count,is,seed)
      implicit none
      include 'headers/ose.h'
c      
      integer jcounter,kk,i,j
      complex uh(nmodes), vh(nmodes), wh(nmodes), temp
      real ran1
      integer seed
      integer count(nmodes), is(nmodes,kose)
c      
      do kk=1,jcounter
        if (count(kk).ne.0) then 
c u
	    i=int(ran1(seed)*count(kk))+1
	    j=is(kk,i)
	    temp=uh(kk)
	    uh(kk)=uh(j)
	    uh(j)=temp	
c v
	    i=int(ran1(seed)*count(kk))+1
	    j=is(kk,i)
	    temp=vh(kk)
	    vh(kk)=vh(j)
	    vh(j)=temp	
c w
	    i=int(ran1(seed)*count(kk))+1
	    j=is(kk,i)
	    temp=wh(kk)
	    wh(kk)=wh(j)
	    wh(j)=temp
        endif
      enddo
      return
      end
       
cccccccccccccccccccccccccccccccccccccc
      subroutine rogal_new_setup(jcounter,kx,ky,kz,uh,vh,wh,cutoffu,cutoffw,uscale,xu,xcf,yv)
c
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/ose.h'
c
      real pi, ftfac
      integer count
      real cutoffu(nx), cutoffw(nx)
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      
      integer jcounter, kkmax
      integer i,j,k
      
      real xu(nx), yv(ny), xcf(nx)
      real uscale
c
      pi=4.*atan(1.0)
      ftfac=(1./(2.*kose))**3  !1./nx/nx/nx  -- ft factor
c     cutoff function
      do i=1,nx
        cutoffu(i)=(1-(tanh(shape1 * (zstart-xu (i)) )+1)/(tanh(shape1)+1))*
     &             (1-(tanh(shape2 * (xu (i)-zend) )+1)/(tanh(shape2)+1))

        cutoffw(i)=(1-(tanh(shape1 * (zstart-xcf(i)) )+1)/(tanh(shape1)+1))*
     &	           (1-(tanh(shape2 * (xcf(i)-zend) )+1)/(tanh(shape2)+1))
     
      enddo 	 
cccccccccc hardwire	 
	 
         write(456,'(2 G25.15)') (xu (i),cutoffu(i),i=1,nx)
      	 write(456,'(2 G25.15)') (xcf(i),cutoffw(i),i=1,nx)
	 close(456)
c	 return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read wavenumbers
         open(101,file=NOISEFILE,form='formatted',status='old') ! NOISEFILE is in the input
         read(101,'(I5)') jcounter,kkmax
	 if (kkmax.ne.kose) stop 'number of modes allocated is different from that specified (jd-ose [rogallo])'
         if (jcounter.gt.nmodes) stop 'increase space for inflow wavenumbers (jd-ose [rogallo])'
         do k=1,jcounter 
          read(101,'(2 I5,9G25.15)') i,count,kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
c          read(101,*) i,count,kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
         enddo
         close(101)
c
	 uscale=ftfac*1. ! disturbance scaling parameter, do _not_ remove ftfac
c
c     set up grid
         if (abs(zlen-yv(ny-1)).gt.1E7) stop 'OSE [Rogallo]: spanwise domain size not correct (see ose.h)'
      return
      end
      
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
      subroutine rogal_old_setup(jcounter,kx,ky,kz,count,is,uh,vh,wh,cutoffu,cutoffw,uscale,xu,xcf,yv)
c
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/ose.h'
c
      real pi, scx
      real cutoffu(nx), cutoffw(nx)
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      
      integer jcounter, kkmax, count(nmodes)
      integer i,j,k,ipos, is(nmodes,kose)
      
      real xu(nx), yv(ny), xcf(nx)
      real uscale
c
      pi=4.*atan(1.0)
      scx=2.*pi/zlen
c     cutoff function
      do i=1,nx
        cutoffu(i)=(1-(tanh(shape1 * (zstart-xu (i)) )+1)/(tanh(shape1)+1))*
     &             (1-(tanh(shape2 * (xu (i)-zend) )+1)/(tanh(shape2)+1))

        cutoffw(i)=(1-(tanh(shape1 * (zstart-xcf(i)) )+1)/(tanh(shape1)+1))*
     &	           (1-(tanh(shape2 * (xcf(i)-zend) )+1)/(tanh(shape2)+1))
     
      enddo 	 
cccccccccc hardwire	 
	 
         write(456,'(2 G25.15)') (xu (i),cutoffu(i),i=1,nx)
      	 write(456,'(2 G25.15)') (xcf(i),cutoffw(i),i=1,nx)
	 close(456)
c	 return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read wavenumbers
         open(101,file=NOISEFILE,form='formatted',status='old') ! NOISEFILE is in the input
         read(101,'(I5)') jcounter,kkmax
	 if (kkmax.ne.kose) stop 'number of modes allocated is different from that specified (jd-ose [rogallo])'
         if (jcounter.gt.nmodes) stop 'increase space for inflow wavenumbers (jd-ose [rogallo])'
         do k=1,jcounter 
          read(101,'(2 I5,9G25.15)') i,count(k),kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
          do ipos=1,count(k)
           read(101,'(I5)') is(k,ipos)
          enddo	  	       
         enddo
         close(101)
c
	 uscale=sqrt(scx)    ! disturbance scaling parameter
c
c     set up grid
         if (abs(zlen-yv(ny-1)).gt.1E7) stop 'OSE [Rogallo]: spanwise domain size not correct (see ose.h)'
      return
      end
      
cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
      subroutine jdose_new_setup(jcounter,kx,ky,kz,uh,vh,wh,cutoffu,cutoffw,uscale,eta,detacor,d2U,U_,Ur,xu,xcf,yv,
     &                           a_w,b_w,a_u,b_u,ps_w,ps_u)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/ose.h'
c
      real pi, ftfac, count
      real cutoffu(nx), cutoffw(nx)
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      
      integer jcounter, kkmax
      integer i,j,k
      
      real uscale
      real xu(nx), xcf(nx), yv(ny)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     GRID DATA
      real eta(nose),detacen(nose),detacor(nose)
c     OSE INTERPOLATION DATA
      real a_w(nx),b_w(nx),a_u(nx),b_u(nx)
      integer ps_w(nx),ps_u(nx)
ccccc BLASIUS VARIABLES      
      complex U_(nose),d2U(nose)
      real    Ur(2*nose)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      pi=4.*atan(1.0)
      ftfac=(1./(2.*kose))**3  !1./nx/nx/nx  -- ft factor
c     cutoff function
      do i=1,nx
        cutoffu(i)=(1-(tanh(shape2 * (xu (i)-zend) )+1)/(tanh(shape2)+1))
        cutoffw(i)=(1-(tanh(shape2 * (xcf(i)-zend) )+1)/(tanh(shape2)+1))
      enddo 	 
cccccccccc hardwire	 
	 
         write(456,'(2 G25.15)') (xu (i),cutoffu(i),i=1,nx)
      	 write(456,'(2 G25.15)') (xcf(i),cutoffw(i),i=1,nx)
	 close(456)
c	 return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read wavenumbers
         open(101,file=NOISEFILE,form='formatted',status='old') ! NOISEFILE is in the input
         read(101,'(I5)') jcounter,kkmax
	 if (kkmax.ne.kose) stop 'number of modes allocated is different from that specified (jd-ose [rogallo])'
         if (jcounter.gt.nmodes) stop 'increase space for inflow wavenumbers (jd-ose [rogallo])'
         do k=1,jcounter 
          read(101,'(2 I5,9G25.15)') i,count,kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
         enddo
         close(101)
c
	 uscale=ftfac*1. ! disturbance scaling parameter, do _not_ remove ftfac
c
         if (-x_.ne.x_ref) stop 'xref mismatch between input file and ose.h' ! x_ obsolescent
	 
         L_=(ru1*(x_-xref_)/uinfty)**0.5 ! in ose.h ! x_ obsolescent
         Re_=uinfty*L_/ru1
c     set up grid
         if (abs(zlen-yv(ny-1)).gt.1E7) stop 'OSE: spanwise domain size not correct (see ose.h)'
c
         call grid_ose(eta,detacen,detacor,nose,zlen,stretch,L_) ! stretch and nose in ose.h
c
         call setinterp_ose(eta,xu(2),nose,nx-1,a_w(2),b_w(2),ps_w(2),L_) 
         call setinterp_ose(eta,xcf(2),nose,nx-1,a_u(2),b_u(2),ps_u(2),L_) ! because the 1st point is ghost
         call setup_ose(eta,detacor,d2U,U_,Ur,nose,ru1)
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rogallo(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type)
      implicit none
      include 'headers/dimension.h'
      include 'headers/common.h'
      include 'headers/ose.h'
      include 'mpif.h'
      
c     scalars for isotropic inflow
      real x_u,x_w,y_v,y_u,z_u,z_w
      complex zh1,zh2,xh_u,xh_w,yh_v,yh_u,v1,v2,v3
c
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      real uose(ny,nx),vose(ny,nx),wose(ny,nx)
      integer jcounter,i,iii,j,jjj,k,kk
      real tlevel
      integer ierr,status(MPI_STATUS_SIZE)
      integer nxloc
      integer type
      
      real yv(ny),ycf(ny),xu(nx),xcf(nx),zwg(nzg)
c 
c     generate inflow with fourier wall-normal basis
      do i=2,nx-1
       do j=2,ny-1
        uose(j,i)=0.    
        vose(j,i)=0.    
        wose(j,i)=0.
       enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      x_u=uinfty*tlevel+zwg(1)!      /uinfty
      x_w=uinfty*tlevel+0.5*(3.*zwg(1)-zwg(2))!   /uinfty

      do kk=1,jcounter
      
         xh_u=exp((0.,-1.)*kx(kk)*x_u)
         xh_w=exp((0.,-1.)*kx(kk)*x_w)
c
c     NOTE:FLIPPED CASE -- distribute load between cpus
         nxloc=(nx-2)/mysize
	 if (mod(nx-2,mysize).ne.0) nxloc=nxloc+1
c	 
         do iii=nxloc*myrank+2,min(nxloc*(myrank+1)+1,nx-1)
	  z_u=xcf(iii)          !z_u (streamwise v. coordinate)
	  z_w=xu(iii)           !z_w (transverse v. coordinate)  
          zh1=exp((0.,-1.)*kz(kk)*z_u)
          zh2=exp((0.,-1.)*kz(kk)*z_w)

             do jjj=2,ny-1
                y_v=yv(jjj)
		y_u=ycf(jjj)
                yh_v=exp((0.,-1.)*ky(kk)*y_v)
                yh_u=exp((0.,-1.)*ky(kk)*y_u)
	        v1=xh_u*yh_u*zh1
		v2=xh_w*yh_v*zh1
		v3=xh_w*yh_u*zh2
                wose(jjj,iii)=wose(jjj,iii)+2.*real(uh(kk)*v1)
                vose(jjj,iii)=vose(jjj,iii)+2.*real(vh(kk)*v2)
                uose(jjj,iii)=uose(jjj,iii)+2.*real(wh(kk)*v3)
             enddo
         enddo
       enddo  
c
         if (myrank.eq.0) then 
           do k=1,mysize-1
	    if (k*nxloc+2.lt.nx) then  
	     call MPI_RECV(uose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,0,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(vose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,1,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(wose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,2,MPI_COMM_HYDRO,status,ierr)
            endif
	   enddo
	 else
            if (myrank*nxloc+2.lt.nx) then  
	     call MPI_SEND(uose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,0,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(vose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,1,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(wose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,2,MPI_COMM_HYDRO,ierr)
            endif
         endif         
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rogallo_old(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type)
      implicit none
      include 'headers/dimension.h'
      include 'headers/common.h'
      include 'headers/ose.h'
      include 'mpif.h'
      
c     scalars for isotropic inflow
      real x_u,x_w,y_v,y_u,z_u,z_w
      complex zh1,zh2,xh_u,xh_w,yh_v,yh_u,v1,v2,v3
c
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      real uose(ny,nx),vose(ny,nx),wose(ny,nx)
      integer jcounter,i,iii,j,jjj,k,kk
      real tlevel
      integer ierr,status(MPI_STATUS_SIZE)
      integer nxloc
      integer type
      
      real yv(ny),ycf(ny),xu(nx),xcf(nx),zwg(nzg)
c 
c     generate inflow with fourier wall-normal basis
      do i=2,nx-1
       do j=2,ny-1
        uose(j,i)=0.    
        vose(j,i)=0.    
        wose(j,i)=0.
       enddo
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      x_u=uinfty*tlevel+zwg(1)!      /uinfty
      x_w=uinfty*tlevel+0.5*(3.*zwg(1)-zwg(2))!   /uinfty

      do kk=1,jcounter
      
         xh_u=exp((0.,1.)*kx(kk)*x_u)
         xh_w=exp((0.,1.)*kx(kk)*x_w)
c
c     NOTE:FLIPPED CASE -- distribute load between cpus
         nxloc=(nx-2)/mysize
	 if (mod(nx-2,mysize).ne.0) nxloc=nxloc+1
c	 
         do iii=nxloc*myrank+2,min(nxloc*(myrank+1)+1,nx-1)
	  z_u=xcf(iii)-50.!(full plate start)        !z_u (streamwise v. coordinate)
	  z_w=xu(iii)-50.                            !z_w (transverse v. coordinate)  
          zh1=exp((0.,1.)*kz(kk)*z_u)
          zh2=exp((0.,1.)*kz(kk)*z_w)

             do jjj=2,ny-1
                y_v=yv(jjj)
		y_u=ycf(jjj)
                yh_v=exp((0.,1.)*ky(kk)*y_v)
                yh_u=exp((0.,1.)*ky(kk)*y_u)
	        v1=xh_u*yh_u*zh1
		v2=xh_w*yh_v*zh1
		v3=xh_w*yh_u*zh2
                wose(jjj,iii)=wose(jjj,iii)+2.*real(uh(kk)*v1)
                vose(jjj,iii)=vose(jjj,iii)+2.*real(vh(kk)*v2)
                uose(jjj,iii)=uose(jjj,iii)+2.*real(wh(kk)*v3)
             enddo
         enddo
       enddo  
c
         if (myrank.eq.0) then 
           do k=1,mysize-1
	    if (k*nxloc+2.lt.nx) then  
	     call MPI_RECV(uose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,0,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(vose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,1,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(wose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,2,MPI_COMM_HYDRO,status,ierr)
            endif
	   enddo
	 else
            if (myrank*nxloc+2.lt.nx) then  
	     call MPI_SEND(uose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,0,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(vose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,1,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(wose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,2,MPI_COMM_HYDRO,ierr)
            endif
         endif         
c
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jdorrsom(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type,eta,detacen,detacor,U_,d2U,
     &                 a_w,b_w,a_u,b_u,ps_w,ps_u)
      implicit none
      include 'headers/dimension.h'
      include 'mpif.h'
      include 'headers/ose.h'
      include 'headers/common.h'
      
c     scalars for isotropic inflow
      real x_u,x_w,y_v,y_u,z_u,z_w
      complex zh1,zh2,xh_u,xh_w,yh_v,yh_u,v1,v2,v3
c
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      real uose(ny,nx),vose(ny,nx),wose(ny,nx)
      integer jcounter,i,iii,j,jjj,k,kk
      real tlevel,zwg(nzg)
      integer ierr,status(MPI_STATUS_SIZE)
      integer nxloc
      integer type
      
      real yv(ny),ycf(ny),xu(nx),xcf(nx)

c     GRID DATA
      real eta(nose),detacen(nose),detacor(nose)
c     OSE INTERPOLATION DATA
      real a_w(nx),b_w(nx),a_u(nx),b_u(nx)
      integer ps_w(nx),ps_u(nx)
ccccc BLASIUS VARIABLES      
      complex U_(nose),d2U(nose)
      complex phi(nose),phip(nose)

c     generate inflow with os wall-normal basis
      do i=2,nx-1
       do j=2,ny-1
        uose(j,i)=0.    
        vose(j,i)=0.    
        wose(j,i)=0.
       enddo
      enddo
c       	  
      x_u=uinfty*tlevel+zwg(1)!      /uinfty
      x_w=uinfty*tlevel+0.5*(3.*zwg(1)-zwg(2))!       /uinfty

      do kk=1,jcounter
c
        xh_u=exp((0.,-1.)*kx(kk)*x_u)
        xh_w=exp((0.,-1.)*kx(kk)*x_w)
c	        
        call ose_(Re_,kx(kk)*L_/uinfty,kz(kk)*L_,U_,d2U,phi,phip,eta,detacen,detacor,nose)
c
c     NOTE:FLIPPED CASE -- distribute load between cpus
        nxloc=(nx-2)/mysize
	if (mod(nx-2,mysize).ne.0) nxloc=nxloc+1
c
        do iii=nxloc*myrank+2,min(nxloc*(myrank+1)+1,nx-1)
	  i=ps_u(iii) ! u and v coords
	  zh1=(a_u(iii)*phip(i)+b_u(iii)*phip(i+1))/L_
	  i=ps_w(iii) ! w coords
	  zh2=(a_w(iii)* phi(i)+b_w(iii)* phi(i+1))
             do jjj=2,ny-1
                y_v=yv(jjj)
		y_u=ycf(jjj)
                yh_v=exp((0.,-1.)*ky(kk)*y_v)
                yh_u=exp((0.,-1.)*ky(kk)*y_u)
	        v1=xh_u*yh_u*zh1
		v2=xh_w*yh_v*zh1
		v3=xh_w*yh_u*zh2
                wose(jjj,iii)=wose(jjj,iii)+2.*real(uh(kk)*v1)
                vose(jjj,iii)=vose(jjj,iii)+2.*real(vh(kk)*v2)
                uose(jjj,iii)=uose(jjj,iii)+2.*real(wh(kk)*v3)
             enddo
        enddo
      enddo  
c
      if (myrank.eq.0) then 
          do k=1,mysize-1
            if (k*nxloc+2.lt.nx) then  
	     call MPI_RECV(uose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,0,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(vose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,1,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(wose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,2,MPI_COMM_HYDRO,status,ierr)
            endif
  	  enddo
      else
	  if (myrank*nxloc+2.lt.nx) then  
 	     call MPI_SEND(uose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,0,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(vose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,1,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(wose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,2,MPI_COMM_HYDRO,ierr)
          endif
      endif         
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine waves_write(jcounter,kx,ky,kz,uh,vh,wh)
      implicit none
      include 'headers/ose.h'
      include 'headers/common.h'
      
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      
      integer jcounter 
      integer k
      
      if (myrank.eq.0) then        
       open(101,file=NOISEFILE,form='formatted',status='unknown')
       write(101,'(I5)') jcounter,kose
         do k=1,jcounter 
           write(101,'(2 I5,9G25.15)') k,0,kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
         enddo
        close(101)
      endif	
      return
      end
cccccccccccccccccccccccccccccccccccccccccccc
      subroutine waves_write_old(jcounter,kx,ky,kz,count,is,uh,vh,wh)
      implicit none
      include 'headers/ose.h'
      include 'headers/common.h'
      
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      
      integer jcounter, count(nmodes), is(nmodes,kose)
      integer k, ipos
      
      if (myrank.eq.0) then        
       open(101,file=NOISEFILE,form='formatted',status='unknown')
       write(101,'(I5)') jcounter,kose
         do k=1,jcounter 
           write(101,'(2 I5,9G25.15)') k,count(k),kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
           do ipos=1,count(k)
             write(101,'(I5)') is(k,ipos)
           enddo	  	       
         enddo
        close(101)
      endif ! myrank
      return
      end
	