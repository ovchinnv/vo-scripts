       integer nb2, loadnum
       parameter (nb2=224)
       character*80 box2name
       data box2name/'b2/fstle2h.new'/
       data loadnum/0/
       
       
!
      if (box2) call box2_setup(a_u,a_w,b_u,b_w,ps_u,ps_w,tlevel,xu,xcf,index,type)
!
       
       
       
       
       
       
       
       
      subroutine box2_setup(a_u,a_w,b_u,b_w,ps_u,ps_w,tlevel,xu,xcf,idx,type)
      implicit none
     
      include 'headers/common.h'
      include 'headers/dimension.h' 
      include 'headers/box2.h'   
      include 'mpif.h'
cccccc
c
      character*3 idx(0:999)
      integer nx_b2, ny_b2
      real dt_b2, x_b2(nb2), xc_b2(nb2), y_b2(ny), t_b2, u_b2(ny,nb2),v_b2(ny,nb2),w_b2(ny,nb2)
      integer i,j,k,type,ierr
      real tlevel, xu(nx), xcf(nx)
c     INTERPOLATION DATA
      real a_w(nx),b_w(nx),a_u(nx),b_u(nx)
      integer ps_w(nx),ps_u(nx)
c             
      integer nlen
      nlen=index(box2name,' ')-1
      if(nlen.le.0) nlen=80
c     
      if (myrank.eq.0) then 
	  open(63, file=box2name(1:nlen)//idx(loadnum),form='binary',status='old')
	  read(63) t_b2,nx_b2,ny_b2,(x_b2(i),i=1,nx_b2),(y_b2(j),j=1,ny_b2)
          tlevel=t_b2
	  if (lcfl) write(0,*) 'WARNING: box2 with cfl ON'
	  if (ny_b2.ne.ny) write(0,*) 'WARNING: box2 and box1 spanwise domains do not match'
	  if (nx_b2.ne.nb2) write(0,*) 'WARNING: box2 and box1 transverse domains do not match'
	  do i=2,nx_b2
	   xc_b2(i)=0.5*(x_b2(i)+x_b2(i-1))
	  enddo
	  xc_b2(1)=-xc_b2(2)

c         offset grids -- throw away near-wall region in the disturbance
          do i=1,nx_b2
	   x_b2(i)=x_b2(i)-1.3
	   xc_b2(i)=xc_b2(i)-1.3
	  enddo  

	  call setinterp_ose( x_b2, xu(2),nx_b2,nx-1,a_w(2),b_w(2),ps_w(2),1.) !w-interpolation (streamwise) 
	  call setinterp_ose(xc_b2,xcf(2),nx_b2,nx-1,a_u(2),b_u(2),ps_u(2),1.) !u,v-interpolation (span & transverse)
	 endif !myrank
	 call MPI_BCAST(tlevel,1,type,0,MPI_COMM_HYDRO,ierr)
ccccccccccccc skip until start index
         if (myrank.eq.0) then 
          do k=1,mod(itini,500)  
          read(63) dt_b2
	   do i=2,nx_b2-1
	    read(63) (u_b2(j,i), v_b2(j,i), w_b2(j,i),j=2,ny-1)
           enddo !i
          enddo !k
      endif !myrank 	   
      end ! subroutine
cccccccccccc

      subroutine box2dt(dtm1,icycle,idx,type)
      implicit none

      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/box2.h'
      include 'mpif.h'
      
      character*3 idx(0:999)
      real dtm1
      integer nx_b2, ny_b2
      real dt_b2, x_b2(nb2), xc_b2(nb2), y_b2(ny), t_b2
      integer i,j,k,type,ierr,icycle
c      
      save loadnum      
c        
      integer nlen
      nlen=index(box2name,' ')-1
      if(nlen.le.0) nlen=80
c       
      if (myrank.eq.0) then 
          if ((icycle-itini).gt.1.and.mod(icycle,500).eq.1) then 
 	   loadnum=loadnum+1
	   open(63, file=box2name(1:nlen)//idx(loadnum),form='binary',status='old')
	   read(63) t_b2,nx_b2,ny_b2,(x_b2(i),i=1,nx_b2),(y_b2(j),j=1,ny_b2)
	   if (lcfl) write(0,*) 'WARNING: box2 with cfl ON'
	   if (ny_b2.ne.ny) write(0,*) 'WARNING: box2 and box1 spanwise domains do not match'
	   if (nx_b2.ne.nb2) write(0,*) 'WARNING: box2 and box1 transverse domains do not match'
          endif !icycle
	  read(63) dt_b2
	  dtm1=dt_b2
      endif !myrank
      call MPI_BCAST(dtm1,1,type,0,MPI_COMM_HYDRO,ierr)
      end !subroutine
cccccccccccccccccccccccc
      subroutine box2read(ustar,vstar,wstar,a_u,a_w,b_u,b_w,ps_u,ps_w)
      implicit none
c
      include 'headers/dimension.h'
      include 'headers/common.h'
      include 'headers/box2.h'
       
c     INTERPOLATION DATA
      real a_w(nx),b_w(nx),a_u(nx),b_u(nx)
      integer ps_w(nx),ps_u(nx)
c             
      real ustar(nx,ny,nz),vstar(nx,ny,nz),wstar(nx,ny,nz)
      integer nx_b2, ny_b2, i_cen, i_cor,i,j,k
      real u_b2(ny,nb2),v_b2(ny,nb2),w_b2(ny,nb2)
cccccccccccc
      if (myrank.eq.0) then 
	  do i=2,nb2-1
	   read(63) (u_b2(j,i), v_b2(j,i), w_b2(j,i),j=2,ny-1)
          enddo !i
c   bcs   
          do j=2,ny-1 
           u_b2(j,1)=-u_b2(j,2)
	   v_b2(j,1)=-v_b2(j,2)	       
	   w_b2(j,1)=0.	       
           u_b2(j,nb2)=u_b2(j,nb2-1)
	   v_b2(j,nb2)=v_b2(j,nb2-1)	       
	   w_b2(j,nb2)=w_b2(j,nb2-1)
c   interpolate
	  enddo !j
	  do i=2,nx-1
	   do j=2,ny-1
	    i_cen=ps_u(i)
	    i_cor=ps_w(i)
	    ustar(i,j,1)=a_w(i)*u_b2(j,i_cor)+b_w(i)*u_b2(j,i_cor+1)
	    vstar(i,j,1)=a_u(i)*v_b2(j,i_cen)+b_u(i)*v_b2(j,i_cen+1)
	    wstar(i,j,1)=a_u(i)*w_b2(j,i_cen)+b_u(i)*w_b2(j,i_cen+1)
	   enddo !j
	  enddo !i 	       
	 endif !myrank 
      return
      end
c
      subroutine box2bc(ustar,vstar,wstar,uo,vo,wo)
      implicit none
      include 'headers/dimension.h'
      include 'headers/common.h'
c
      real ustar(nx,ny,nz),vstar(nx,ny,nz),wstar(nx,ny,nz)
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      integer i,j,k
c      
      if (myrank.eq.0) then 
       do i=1,nx
        do j=1,ny
	 uo(i,j,1)=ustar(i,j,1)
	 vo(i,j,1)=vstar(i,j,1)
	 wo(i,j,1)=wstar(i,j,1)
        enddo
       enddo
      endif 	
      return
      end
c
      subroutine box2close()
      implicit none
      include 'headers/common.h'
      
          if (myrank.eq.0) close(63)
      end
