      module LES
! subroutines for large-eddy simulation
       use size
       
       
       
       
!.....sxx,sxy,sxz,syy,syz,szz:	strain rates
       
       
       !     LES DATA 
      real dxdydz(nx,ny,nz),g(nx,ny,nz)
      real sxx(nx,ny,nz),syy(nx,ny,nz),szz(nx,ny,nz)
      real sxy(nx,ny,nz),sxz(nx,ny,nz),syz(nx,ny,nz)
      real uoc(nx,ny,nz),voc(nx,ny,nz),woc(nx,ny,nz)
      real gbar(nx,ny,nz),lm(nx,ny,nz),mm(nx,ny,nz)
      real Ilm(nx,ny,nz),Imm(nx,ny,nz)
      real wrk1(nx,ny,nz),wrk2(nx,ny,nz),wrk3(nx,ny,nz)
      real wrk4(nx,ny,nz)
!      stats
      real tvavg(nx,nz)!,tvavgg(nx,nzg),tvi




!          part of LES initialization
!
#ifdef LES
      CALL filtersize(xu,yv,zw,dxdydz,nx,ny,nz)
      CALL REFRESHBC(dxdydz,nx*ny,(nz-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      if (myrank.eq.0) then 
        do i=2,ix2
	 do j=2,jy2
	  dxdydz(i,j,1)=dxdydz(i,j,2)
	 enddo
	enddo
      endif	  
      CALL weights(xc,yc,zc,nx,ny,nz)
#endif
