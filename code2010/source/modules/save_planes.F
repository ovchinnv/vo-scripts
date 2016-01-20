


      integer nyplanes,nzplanes,nxplanes
      parameter (nyplanes=1, nxplanes=1, nzplanes=1)
      real nyploc(nyplanes), nxploc(nxplanes), nzploc(nzplanes)

      character*80 planesname
      data planesname/'planes/t3b'/
   
      data nxploc/1.5/
      data nzploc/55/
      data nyploc/25/
      


!ccccccc variables for 2D planes output


      include 'headers/planes.h'

      integer cpuxy(nzplanes),kxy(nzplanes)
!      real u2dxz(nx,nz,2),v2dxz(nx,nz,2),w2dxz(nx,nz,2),p2dxz(nx,nz,2)
!      real u2dxy(nx,ny),v2dxy(nx,ny),w2dxy(nx,ny),p2dxy(nx,ny)
!      real u2dyz(ny,nz,2),v2dyz(ny,nz,2),w2dyz(ny,nz,2),p2dyz(ny,nz,2)
!ccccccc



!cccccccccccc setup for output planes cccccccccccccc      
      if (planes) call planes_setup(cpuxy,kxy,zwg)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccc
      subroutine planes_setup(cpuxy,kxy,zwg)
      implicit none
      include 'headers/dimension.h'
      include 'headers/planes.h'
      include 'headers/common.h'
 
      integer cpuxy(nzplanes),kxy(nzplanes)
      integer k,n
      real zwg(nzg)
 
c      find save locations for xy planes       
       do n=1,nzplanes
         kxy(n)=1
	 cpuxy(n)=0 
         k=1
         do while (zwg(k).lt.(nzploc(n)))
            k=k+1
         enddo
	 if (k.gt.1) then 
	   cpuxy(n)=(k-2)/(nz-2)
	   kxy(n)=mod(k-2,nz-2)+1
         endif
	 if (cpuxy(n).gt.mysize-1) then 
	   cpuxy(n)=cpuxy(n)-1
	   kxy(n)=nz
	 endif  
	enddo !n
      return
      end 
ccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_planes(uo,vo,wo,p,
#ifdef LES
     &                        tv,
#endif
     &  xu,xcf,yv,ycf,zcf,zcg,cpuxy,kxy,type, real2d, cross2d,idx,idx2,icycle)
      
      implicit none
      include 'headers/planes.h'
      include 'headers/dimension.h'
      include 'headers/common.h'
      include 'mpif.h'
      
c
      integer cpuxy(nzplanes),kxy(nzplanes)
      real u2dxz(nx,nz,2),v2dxz(nx,nz,2),w2dxz(nx,nz,2),p2dxz(nx,nz,2)
      real u2dxy(nx,ny),v2dxy(nx,ny),w2dxy(nx,ny),p2dxy(nx,ny)
      real u2dyz(ny,nz,2),v2dyz(ny,nz,2),w2dyz(ny,nz,2),p2dyz(ny,nz,2)
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),p(nx,ny,nz)
      real xu(nx),xcf(nx),yv(ny),ycf(ny),zcf(nz),zcg(nzg)
      
      integer ierr, status(MPI_STATUS_SIZE)
      integer type, real2d, cross2d
      integer i,j,k,m,n,nn,icycle
      character*3 idx(0:999) 
      character*5 idx2(0:99999) 

      integer nlen
      nlen=index(planesname,' ')-1
      if (nlen.le.0) nlen=80

cccccccc xz
      do m=1,nyplanes
        j=1
	do while (yv(j).le.nyploc(m).and.j.le.jy2)
         j=j+1
	enddo
         do k=2,kz2
	  do i=2,ix2
           u2dxz(i,k,2)=-0.25*((wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1)+
     &                    (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)-
     &                    (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1)-
     &                    (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k))  !omx
           v2dxz(i,k,2)=-0.25*((uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1)+
     &                    (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)-
     &                    (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1)-
     &		          (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i))  !omy
           w2dxz(i,k,2)=-0.25*((vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1)+
     &		          (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)-
     &                    (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1)-
     &                    (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j))  ! omz
          enddo
	 enddo 
c
         if (myrank.eq.0) then !1st processor writes datafile
           open(62,file=planesname(1:nlen)//'.'//idx(m)//'.xz-'//idx2(icycle),status='new',form='unformatted')
           write(62) nx-2,nzg-2  
           write(62) (xcf(i),i=2,ix2)		 
           write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (0.5*(uo(i,j,k)+uo(i-1,j,k)),
     &		             vo(i,j,k),
     &			     0.5*(wo(i,j,k)+wo(i,j,k-1)),
     &			     p(i,j,k),                   
     &		             u2dxz(i,k,2),v2dxz(i,k,2),w2dxz(i,k,2),                
#ifdef LES
     &                       tv(i,j,k),                 
#endif
     &                       i=2,ix2)
              enddo
c              receive velocity from other processors
           do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(u2dxz(1,1,1),nx*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dxz(1,1,1),nx*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dxz(1,1,1),nx*nz,type,n,nn+2,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dxz(1,1,1),nx*nz,type,n,nn+3,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(u2dxz(1,1,2),nx*nz,type,n,nn+4,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dxz(1,1,2),nx*nz,type,n,nn+5,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dxz(1,1,2),nx*nz,type,n,nn+6,mpi_comm_hydro,status,ierr)
#ifdef LES
                  call MPI_Recv(p2dxz(1,1,2),nx*nz,type,n,nn+7,mpi_comm_hydro,status,ierr)
#endif
c     
                  do k=2,nz-1
                    write(62) (0.5*(u2dxz(i,k,1)+u2dxz(i-1,k,1)),
     &		               v2dxz(i,k,1),
     &			       0.5*(w2dxz(i,k,1)+w2dxz(i,k-1,1)),
     &			       p2dxz(i,k,1),
     &                         u2dxz(i,k,2),v2dxz(i,k,2),w2dxz(i,k,2),
#ifdef LES
     &		               p2dxz(i,k,2),
#endif     
     &		    i=2,nx-1)
                  enddo
            enddo		           	       
	    write(62) ycf(j) ! write position last		           	       
            close(62)
cc
         else
           nn=myrank*10
           call MPI_Send(uo(1,j,1),1,real2d,0,nn,  MPI_COMM_HYDRO,ierr)
           call MPI_Send(vo(1,j,1),1,real2d,0,nn+1,MPI_COMM_HYDRO,ierr)
           call MPI_Send(wo(1,j,1),1,real2d,0,nn+2,MPI_COMM_HYDRO,ierr)
           call MPI_Send(p (1,j,1),1,real2d,0,nn+3,MPI_COMM_HYDRO,ierr)	      
           call MPI_Send(u2dxz(1,1,2),nx*nz,type,0,nn+4,MPI_COMM_HYDRO,ierr)
           call MPI_Send(v2dxz(1,1,2),nx*nz,type,0,nn+5,MPI_COMM_HYDRO,ierr)
           call MPI_Send(w2dxz(1,1,2),nx*nz,type,0,nn+6,MPI_COMM_HYDRO,ierr)
#ifdef LES
           call MPI_Send(tv (1,j,1),1,real2d,0,nn+7,MPI_COMM_HYDRO,ierr)	      
#endif
         endif !(myrank.eq.0) 
      enddo !m
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccc yz
      do m=1,nxplanes
        i=1
	do while (xu(i).le.nxploc(m).and.i.le.ix2)
         i=i+1
	enddo
      	 
         do k=2,kz2
	  do j=2,jy2
           u2dyz(j,k,2)=-0.25*((wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1)+
     &                    (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)-
     &                    (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1)-
     &                    (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k))  !omx
           v2dyz(j,k,2)=-0.25*((uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1)+
     &                    (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)-
     &                    (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1)-
     &		          (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i))  !omy
           w2dyz(j,k,2)=-0.25*((vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1)+
     &		          (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)-
     &                    (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1)-
     &                    (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j))  ! omz
          enddo
	 enddo 
ccccccccccccccccccccccc
ccccccccccccccccccccccc
         if (myrank.eq.0) then !1st processor writes datafile
           open(62,file=planesname(1:nlen)//'.'//idx(m)//'.yz-'//idx2(icycle),status='new',form='unformatted')
           write(62) ny-2,nzg-2  
           write(62) (ycf(j),j=2,jy2)		 
           write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (uo(i,j,k),0.5*(vo(i,j,k)+vo(i,j-1,k)),0.5*(wo(i,j,k)+wo(i,j,k-1)),p(i,j,k),                   
     &		             u2dyz(j,k,2),v2dyz(j,k,2),w2dyz(j,k,2),                
#ifdef LES
     &                       tv(i,j,k),                 
#endif
     &                                j=2,jy2)
              enddo
c              receive velocity from other processors
           do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(u2dyz(1,1,1),ny*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dyz(1,1,1),ny*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dyz(1,1,1),ny*nz,type,n,nn+2,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dyz(1,1,1),ny*nz,type,n,nn+3,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(u2dyz(1,1,2),ny*nz,type,n,nn+4,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dyz(1,1,2),ny*nz,type,n,nn+5,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dyz(1,1,2),ny*nz,type,n,nn+6,mpi_comm_hydro,status,ierr)
#ifdef LES
                  call MPI_Recv(p2dyz(1,1,2),ny*nz,type,n,nn+7,mpi_comm_hydro,status,ierr)
#endif
c     
                  do k=2,nz-1
                    write(62) (u2dyz(j,k,1),0.5*(v2dyz(j,k,1)+v2dyz(j-1,k,1)),0.5*(w2dyz(j,k,1)+w2dyz(j,k-1,1)),p2dyz(j,k,1),
     &                         u2dyz(j,k,2),v2dyz(j,k,2),w2dyz(j,k,2),
#ifdef LES
     &		    p2dxz(j,k,2),
#endif     
     &		    j=2,jy2)
                  enddo
            enddo
	    write(62) xcf(i) ! write position last		           	       
            close(62)
         else
           nn=myrank*10
           call MPI_Send(uo(i,1,1),1,cross2d,0,nn,  MPI_COMM_HYDRO,ierr)
           call MPI_Send(vo(i,1,1),1,cross2d,0,nn+1,MPI_COMM_HYDRO,ierr)
           call MPI_Send(wo(i,1,1),1,cross2d,0,nn+2,MPI_COMM_HYDRO,ierr)
           call MPI_Send(p (i,1,1),1,cross2d,0,nn+3,MPI_COMM_HYDRO,ierr)	      
           call MPI_Send(u2dyz(1,1,2),ny*nz,type,0,nn+4,MPI_COMM_HYDRO,ierr)
           call MPI_Send(v2dyz(1,1,2),ny*nz,type,0,nn+5,MPI_COMM_HYDRO,ierr)
           call MPI_Send(w2dyz(1,1,2),ny*nz,type,0,nn+6,MPI_COMM_HYDRO,ierr)
#ifdef LES
           call MPI_Send(tv (i,1,1),1,cross2d,0,nn+7,MPI_COMM_HYDRO,ierr)	      
#endif
         endif !(myrank.eq.0) 
      enddo !m
      
cccccccc surface pressure/skin friction
c
         i=2
ccccccccccccccccccccccc
         if (myrank.eq.0) then !1st processor writes datafile
           open(62,file=planesname(1:nlen)//'.pcf0-'//idx2(icycle),status='new',form='unformatted')
           write(62) ny-2,nzg-2  
           write(62) (ycf(j),j=2,jy2)		 
           write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (p(i,j,k),ru1*0.5*(wo(i,j,k)+wo(i,j,k-1))*au(1)*2.*2./uinfty**2,j=2,jy2) !p & cf
               enddo
c              receive velocity from other processors
               do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(w2dyz(1,1,1),ny*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dyz(1,1,1),ny*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
c     
                  do k=2,nz-1
                    write(62) (p2dyz(j,k,1),ru1*0.5*(w2dyz(j,k,1)+w2dyz(j,k-1,1))*au(1)*2.*2./uinfty**2,j=2,jy2)
                  enddo
            enddo		           	       
            write(62) xcf(i)
            close(62)
         else
           nn=myrank*10
           call MPI_Send(wo(i,1,1),1,cross2d,0,nn,  MPI_COMM_HYDRO,ierr)
           call MPI_Send(p (i,1,1),1,cross2d,0,nn+1,MPI_COMM_HYDRO,ierr)	      
         endif !(myrank.eq.0) 
ccccccccccccccccccccccccccccc
cccccccccccccc xy ccccccccccc
      do m=1,nzplanes
       if (myrank.eq.cpuxy(m)) then 
        k=kxy(m)
c
         do j=2,jy2
	  do i=2,ix2
           u2dxy(i,j)=-0.25*((wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1)+
     &                    (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)-
     &                    (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1)-
     &                    (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k))  !omx
           v2dxy(i,j)=-0.25*((uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1)+
     &                    (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)-
     &                    (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1)-
     &		          (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i))  !omy
           w2dxy(i,j)=-0.25*((vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1)+
     &		          (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)-
     &                    (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1)-
     &                    (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j))  ! omz
          enddo
	 enddo 
ccccccccccccccccccccccc
ccccccccccccccccccccccc
           open(62,file=planesname(1:nlen)//'.'//idx(m)//'.xy-'//idx2(icycle),status='new',form='unformatted')
           write(62) nx-2,ny-2  
           write(62) (xcf(i),i=2,ix2)		 
           write(62) (ycf(j),j=2,ny-1)
               do j=2,jy2
                  write(62) (0.5*(uo(i,j,k)+uo(i-1,j,k)),0.5*(vo(i,j,k)+vo(i,j-1,k)),wo(i,j,k),p(i,j,k),                   
     &		             u2dxy(i,j),v2dxy(i,j),w2dxy(i,j),                
#ifdef LES
     &                       tv(i,j,k),                 
#endif
     &                                i=2,ix2)
               enddo
               write(62) zcf(kxy(m))
               close(62)
       endif !myrank
      enddo !m
      return
      end 
cccccccc write 1 xz-plane ccccccccccccccccc
      subroutine plane_xz(uo,vo,wo,p,xcf,ycf,zcg,fname,j,real2d)
      implicit none
      include 'headers/dimension.h'
      include 'headers/common.h'
      include 'mpif.h'
      
      real uo(nx,ny,nz), vo(nx,ny,nz), wo(nx,ny,nz), p(nx,ny,nz)
      real xcf(nx), ycf(ny), zcg(nzg)
      character*80 fname
      integer i,j,k,n,nn,type,ierr,status(MPI_STATUS_SIZE)
      real u2dxz(nx,nz),v2dxz(nx,nz),w2dxz(nx,nz),p2dxz(nx,nz)
      integer real2d

      integer nlen
      nlen=index(fname,' ')-1
      if (nlen.le.0) nlen=80

#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
      

      if (myrank.eq.0) then 
c     write(60,*) 'variables="x","z","u","v","w"'
c     write(60,*) 'zone i=',nx-2,',j=',nz-2,'f=point'
       open(62,file=fname(1:nlen),status='unknown',form='unformatted')
       write(62) nx-2,nzg-2 
       write(62) (xcf(i),i=2,ix2)		 
       write(62) (zcg(k),k=2,nzg-1)
       do k=2,kz2
          write(62) (   !xcf(i),zcf(k),
     &                 0.5*(uo(i,j,k)+uo(i-1,j,k)),
     &                 vo(i,j,k),
     &                 0.5*(wo(i,j,k)+wo(i,j,k-1)),
     &                 p(i,j,k),
     &                 i=2,ix2)
       enddo
       do n=1,mysize-1
          nn=n*10
          call MPI_Recv(u2dxz(1,1),nx*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
          call MPI_Recv(v2dxz(1,1),nx*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
          call MPI_Recv(w2dxz(1,1),nx*nz,type,n,nn+2,mpi_comm_hydro,status,ierr)
          call MPI_Recv(p2dxz(1,1),nx*nz,type,n,nn+3,mpi_comm_hydro,status,ierr)
c     
          do k=2,nz-1
             write(62) (0.5*(u2dxz(i,k)+u2dxz(i-1,k)),
     &                          v2dxz(i,k),
     &                          0.5*(w2dxz(i,k)+w2dxz(i,k-1)),
     &                          p2dxz(i,k),
     &                          i=2,nx-1)
          enddo    
       enddo		           	       
       write(62) ycf(j)
       close(62)
      else
       nn=myrank*10
       call MPI_Send(uo(1,j,1),1,real2d,0,nn,  MPI_COMM_HYDRO,ierr)
       call MPI_Send(vo(1,j,1),1,real2d,0,nn+1,MPI_COMM_HYDRO,ierr)
       call MPI_Send(wo(1,j,1),1,real2d,0,nn+2,MPI_COMM_HYDRO,ierr)
       call MPI_Send(p (1,j,1),1,real2d,0,nn+3,MPI_COMM_HYDRO,ierr)	      
      endif !myrank
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
