










!ccccccccccccc READ INITIAL DATA ccccccccccccccccccc      
      call readini(uo,vo,wo,p,uab1,vab1,wab1,tv,
#ifdef LES
     &             ilm,imm,
#endif
     &             xu,zw,xcf,ycf,tlevel)
!cccccccccccccccccccccccccccccccccccccccccccccccccccc       

















      subroutine readini(uo,vo,wo,p,uab1,vab1,wab1,tv,
#ifdef LES
     &             ilm,imm,
#endif LES
     &             xu,zw,xcf,ycf,tlevel)
c     
      use IO
      use size, only: nx, ny, nz, 
c
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
c       
      integer i,ii,j,k,kk
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      real uab1(nx,ny,nz),vab1(nx,ny,nz),wab1(nx,ny,nz),tv(nx,ny,nz)
      real xu(nx),zw(nz),xcf(nx),ycf(ny),tlevel
      integer nlen1, nlen3
#ifdef LES
      real ilm(nx,ny,nz),imm(nx,ny,nz)
#endif            
c
      real xloc,dstar
      real u1(1:nx),w1(1:nx)
      real u2(1:nx),w2(1:nx)
      integer idir
      real pi,dummy1,dummy2

      if(ifield.EQ.0) THEN
c     
      pi=4.0*atan(1.)
         do i=1,ix2+1
	    dummy1=xcf(i)/xu(ix2)
	    dummy2=4.*1./ru1/2.*dummy1*(1-dummy1)
            do j=1,jy2+1
               do k=1,kz2+1
                  uo(i,j,k)=0.0!+0.1/ru1*sin(2.*pi/2*xu(i))
                  vo(i,j,k)=0.0!+0.1/ru1*sin(2.*pi/2*ycf(j))
                  wo(i,j,k)=uinfty!dummy2!+sin(2*2.*pi/2*ycf(j))
c                  wo(i,j,k)=1./ru1/8*.66!uinfty+0.2*uinfty*(ran1(idum)-0.5)
               enddo
            enddo
         enddo
         tlevel=0.
c
      else if (ifIELD.EQ.4) THEN !initialize to Blasius profile
c     
         dstar=0.0              !displacement thickness at inflow 
                                !x_ref--location of inflow (in the input file)
c--------begin bl generation
c--------first point
         k=1
         xloc=zw(k)          !u-location
         call blasius(xu, u1, w1, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar) 
         xloc=0.5*(3*zw(k)-zw(k+1)) !w-location
         u1(1)=-u1(2)
c     
         call blasius(xu, u2, w2, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar)  
         do i=ix1, ix2+1
           do j=jy1+1, jy2
               wo(i,j,k)=u1(i)
               uo(i,j,k)=w2(i)
               vo(i,j,k)=0.
           enddo
         enddo
ccccccccccccccccc
         do k=kz1+1, kz2+1
            xloc=zw(k)          !u-location
            call blasius(xu, u1, w1, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar) 
            xloc=0.5*(zw(k)+zw(k-1)) !w-location
c     
            u1(1)=-u1(2)
c     
            call blasius(xu, u2, w2, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar)  
            do i=ix1, ix2+1
               do j=jy1+1, jy2
                  wo(i,j,k)=u1(i)
                  uo(i,j,k)=w2(i)
                  vo(i,j,k)=0.
               enddo
            enddo
         enddo
c     --------Boundary conditions
            do j=jy1+1,jy2
               do k=kz1+1,kz2
                  uo(ix1,j,k)=0. !wall
                  wo(ix1,j,k)=-wo(ix1+1,j,k) !wall
                  vo(ix1,j,k)=-vo(ix1+1,j,k) !wall
               enddo
            enddo
           do i=ix1+1, ix2
              do k=kz1+1, kz2
                 uo(i,jy1,k)=uo(i,jy2,k)
                 uo(i,jy2+1,k)=uo(i,jy1+1,k)
                 vo(i,jy1,k)=vo(i,jy2,k)
                 wo(i,jy1,k)=wo(i,jy2,k)
                 wo(i,jy2+1,k)=wo(i,jy1+1,k)
              enddo
           enddo
           tlevel=0.
c         
      else if (ifield.eq.2) then
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOSCALAR(STR4,P,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOVECTOR(STR3,UAB1,VAB1,WAB1,NX,NY,NZ,NZG,idir,TLEVEL)
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then
	  if (read_vt_tv) then  
           call ioscalar(vt_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(tv_file,tv,NX,NY,NZ,NZG,idir,TLEVEL)
	  else 
	   do k=1,nz
	    do j=1,ny
	     do i=1,nx
	      Ilm(i,j,k)=ru1
	      tv(i,j,k)=ru1
	     enddo
	    enddo
	   enddo
	  endif  
         endif
#endif

      else if (ifield.eq.22) then	 
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOSCALAR(STR4,P,NX,NY,NZ,NZG,idir,TLEVEL)
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then
	  if (read_vt_tv) then  
           call ioscalar(vt_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(tv_file,tv,NX,NY,NZ,NZG,idir,TLEVEL)
	  else 
	   do k=1,nz
	    do j=1,ny
	     do i=1,nx
	      Ilm(i,j,k)=ru1
	      tv(i,j,k)=ru1
	     enddo
	    enddo
	   enddo
	  endif  
         endif
#endif
c
      else if (ifield.eq.5) then	 
         nlen1=index(STR1,' ')-1
         nlen3=index(STR3,' ')-1
         idir  = -1
         call IOSCALAR(STR1(1:nlen1)//'.u_f',uo,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOSCALAR(STR1(1:nlen1)//'.v_f',vo,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOSCALAR(STR1(1:nlen1)//'.w_f',wo,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOSCALAR(STR3(1:nlen3)//'.uab_f',uab1,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOSCALAR(STR3(1:nlen3)//'.vab_f',vab1,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOSCALAR(STR3(1:nlen3)//'.wab_f',wab1,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOSCALAR(STR4,P,NX,NY,NZ,NZG,idir,TLEVEL)
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then
	  if (read_vt_tv) then  
           call ioscalar(vt_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(tv_file,tv,NX,NY,NZ,NZG,idir,TLEVEL)
	  else 
	   do k=1,nz
	    do j=1,ny
	     do i=1,nx
	      Ilm(i,j,k)=ru1
	      tv(i,j,k)=ru1
	     enddo
	    enddo
	   enddo
	  endif  
         endif
#endif
c
      else if (ifield.eq.33) then	 
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
c        set pressure to 0 since it's not available
         do k=1,nz
	   do j=1,ny
	     do i=1,nx
	       p(i,j,k)=0.
             enddo  
	   enddo
	 enddo   
c hardwiring below!!!!!!!!	 
c         i=1
c         do while (i.le.nx.and.xu(i).le.40.)
c	  i=i+1
c	 enddo
c	 do k=1,nz
c	  do j=1,ny
c	   do ii=1,i
c	    wo(ii,j,k)=1.
c	    vo(ii,j,k)=0.
c	    uo(ii,j,k)=0.
c	   enddo
c	  enddo
c	 enddo             
cccccccccccccc hardwiring above	        
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then
	  if (read_vt_tv) then  
           call ioscalar(vt_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(tv_file,tv,NX,NY,NZ,NZG,idir,TLEVEL)
	  else 
	   do k=1,nz
	    do j=1,ny
	     do i=1,nx
	      Ilm(i,j,k)=ru1
	      tv(i,j,k)=ru1
	     enddo
	    enddo
	   enddo
	  endif  
         endif
#endif
c
      else if (ifield.eq.44) then
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOVECTOR(STR3,UAB1,VAB1,WAB1,NX,NY,NZ,NZG,idir,TLEVEL)
c        set pressure to 0 since it's not available
         do k=1,nz
	   do j=1,ny
	     do i=1,nx
	       p(i,j,k)=0.
             enddo  
	   enddo
	 enddo   
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then
	  if (read_vt_tv) then  
           call ioscalar(vt_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(tv_file,tv,NX,NY,NZ,NZG,idir,TLEVEL)
	  else 
	   do k=1,nz
	    do j=1,ny
	     do i=1,nx
	      Ilm(i,j,k)=ru1
	      tv(i,j,k)=ru1
	     enddo
	    enddo
	   enddo
	  endif  
         endif
#endif
ccccccccc hardwiring to read 2D field
      else if (ifield.eq.55) then	 
         idir  = -1
         call IOVECTOR2D(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
c        set pressure to 0 since it's not available
         do k=1,nz
	   do j=1,ny
	     do i=1,nx
	       p(i,j,k)=0.
             enddo  
	   enddo
	 enddo   
cccccccc hardwiring	 
      endif
      return
      end
      