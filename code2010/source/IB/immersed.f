      subroutine IBInput(str11,
     $     nibmax,nibu,nibv,nibw,nibp,
     $     ijku,ijkv,ijkw,ijkp,coef,ibheight)
C
C-----------------------------------------------------------------------
C
c      IMPLICIT NONE
c
      CHARACTER*80 str11
      INTEGER i,j,k,n,nd
c
      INTEGER nibu,nibv,nibw,nibp,nibmax
      INTEGER ijku(nibmax,3,2),ijkv(nibmax,3,2)
      INTEGER ijkw(nibmax,3,2),ijkp(nibmax,3,2)
      REAL    coef(nibmax,4,2)
      REAL    ibheight
c
c     Open file
c
      OPEN(UNIT=22,FILE=str11,STATUS='old' )
c
      read(22,*) ibheight
c
      read(22,*) nibu
      if(nibu .gt. nibmax) then
         write(6,*) ' Error on the immersed boundary stuff',nibu,nibmax
      endif
      write(6,*) ' u '
      do n=1,nibu
         read(22,1210) nd,((ijku(n,i,j),i=1,3),j=1,2),
     $        (coef(n,1,j),j=1,2)
      enddo
c
      write(6,*) ' v '
      read(22,1200) nibv
      if(nibv .gt. nibmax) then
         write(6,*) ' Error on the immersed boundary stuff',nibv,nibmax
         stop
      endif
      do n=1,nibv
         read(22,1210) nd,((ijkv(n,i,j),i=1,3),j=1,2),
     $        (coef(n,2,j),j=1,2)
      enddo
c
      write(6,*) ' w '
      read(22,1200) nibw
      if(nibw .gt. nibmax) then
         write(6,*) ' Error on the immersed boundary stuff',nibw,nibmax
         stop
      endif
c....my lines......
      do n=1,nibw
         read(22,1210) nd,((ijkw(n,i,j),i=1,3),j=1,2),
     $        (coef(n,3,j),j=1,2)
      enddo
c..............................
c
      write(6,*) ' p '
      read(22,1200) nibp
      if(nibp .gt. nibmax) then
         write(6,*) ' Error on the immersed boundary stuff',nibp,nibmax
         stop
      endif
      do n=1,nibp
         read(22,1210) nd,((ijkp(n,i,j),i=1,3),j=1,2),
     $        (coef(n,4,j),j=1,2)
      enddo
c
      CLOSE(22)
c$$$
c$$$      write(6,*) nibu,nibv,nibw
c$$$      write(6,*) ((ijku(1,i,j),i=1,3),j=1,2)
c$$$      write(6,*) ((ijkv(1,i,j),i=1,3),j=1,2)
c$$$      write(6,*) ((ijkw(1,i,j),i=1,3),j=1,2)
c$$$      write(6,*) ((coef(1,i,j),i=1,3),j=1,2)
c$$$      stop
c$$$
c     
      return
c
 1200 format (i12)
 1210 format (7i5,1p2e25.15)
c
      end 
c ---------------------------------------------------------------
      subroutine ibboundary(nibmax,nibu,nibv,nibw,nibp,
     $        ijku,ijkv,ijkw,ijkp,
     $        coef,
     $        USTAR,VSTAR,WSTAR,nx,ny,nz)
c
c     subroutine to set the velocity at the points above the 
c     immersed boundary to zero.
c
      IMPLICIT NONE
      include 'headers/common.h'
c
      INTEGER ii,jj,kk,iip,jp,kp,n,nd
      INTEGER nx,ny,nz
      INTEGER nibu,nibv,nibw,nibp,nibmax
      INTEGER ijku(nibmax,3,2),ijkv(nibmax,3,2)
      INTEGER ijkw(nibmax,3,2),ijkp(nibmax,3,2)
      REAL    coef(nibmax,4,2)
      REAL    USTAR(nx,ny,nz),VSTAR(nx,ny,nz),WSTAR(nx,ny,nz)
c
c     assign points for the u velocity
c      
      do n=1,nibu
         ii=ijku(n,1,1)
         jj=ijku(n,2,1)
         kk=ijku(n,3,1)
         iip=ijku(n,1,2)
         jp=ijku(n,2,2)
         kp=ijku(n,3,2)
         USTAR(ii,jj,kk) = coef(n,1,1)*USTAR(iip,jp,kp) + coef(n,1,2)
      enddo
c
c     assign points for the v velocity
c      
      do n=1,nibv
         ii=ijkv(n,1,1)
         jj=ijkv(n,2,1)
         kk=ijkv(n,3,1)
         iip=ijkv(n,1,2)
         jp=ijkv(n,2,2)
         kp=ijkv(n,3,2)
         VSTAR(ii,jj,kk) = coef(n,2,1)*VSTAR(iip,jp,kp) + coef(n,2,2)
      enddo
c
c     assign points for the w velocity
c      
      do n=1,nibw
         ii=ijkw(n,1,1)
         jj=ijkw(n,2,1)
         kk=ijkw(n,3,1)
         iip=ijkw(n,1,2)
         jp=ijkw(n,2,2)
         kp=ijkw(n,3,2)
         WSTAR(ii,jj,kk) = coef(n,3,1)*WSTAR(iip,jp,kp) + coef(n,3,2)
      enddo
c
      return
      end
c
c ---------------------------------------------------------------
      subroutine iblag(nibmax,nibp,ijkp,coef,
     $                 Ilm,Imm,nx,ny,nz)
c
c     subroutine which modifies numerator and denominator of
c     coefficient C of the lagrandian dynamic model
c      - Ilm is set to zero
c      - Imm is extrapolated
c
      IMPLICIT NONE
      include 'headers/common.h'
c
      INTEGER ii,jj,kk,iip,jp,kp,n,nd
      INTEGER nx,ny,nz
      INTEGER nibp,nibmax
      INTEGER ijkp(nibmax,3,2)
      REAL    coef(nibmax,4,2)
      REAL    Ilm(nx,ny,nz),Imm(nx,ny,nz)
c
c     assign points for the Ilm and Imm
c      
      do n=1,nibp
         ii=ijkp(n,1,1)
         jj=ijkp(n,2,1)
         kk=ijkp(n,3,1)
         iip=ijkp(n,1,2)
         jp=ijkp(n,2,2)
         kp=ijkp(n,3,2)
         Ilm(ii,jj,kk) = coef(n,4,1)*Ilm(iip,jp,kp)
         Imm(ii,jj,kk) = Imm(iip,jp,kp)
      enddo
c
      return
      end

c
c ---------------------------------------------------------------
      subroutine ibvt(nibmax,nibp,ijkp,coef,vt,nx,ny,nz)
c
c     subroutine which modifies turbulent viscosity
c     in a way that its value on the boundary is zero
c
      IMPLICIT NONE
      include 'headers/common.h'
c
      INTEGER ii,jj,kk,iip,jp,kp,n,nd
      INTEGER nx,ny,nz
      INTEGER nibp,nibmax
      INTEGER ijkp(nibmax,3,2)
      REAL    coef(nibmax,4,2)
      REAL    vt(nx,ny,nz)
c
c     assign points for the Ilm and Imm
c      
      do n=1,nibp
         ii=ijkp(n,1,1)
         jj=ijkp(n,2,1)
         kk=ijkp(n,3,1)
         iip=ijkp(n,1,2)
         jp=ijkp(n,2,2)
         kp=ijkp(n,3,2)
         vt(ii,jj,kk) = coef(n,4,1)*vt(iip,jp,kp)
      enddo
c
      return
      end
