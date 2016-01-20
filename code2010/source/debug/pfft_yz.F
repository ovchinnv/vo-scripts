      subroutine pfft_yz(p,dp,dummyt,dummy,dummy2,au,ap,dy,dz,me,comm,datatype)
c     V. Ovchinnikov     8/06/03
c     DO NOT MODIFY!!!
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
      return
      end
      










