



















C-----SUBROUTINE-Rhsfrg-------------------------P. Flohr--08/07/1994----
C
      SUBROUTINE Rhsfrg(UNI,VNI,WNI,rhsu,rhsv,rhsw,wfrg,
     &                  nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - modifies RHS (add fringe terms)
C                                                                       
C     PARAMETERS: - velocity field, RHS, fringe velocities
C                                                                       
C-----------------------------------------------------------------------
C
c      IMPLICIT NONE
      INCLUDE"mpif.h"
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      INTEGER nx,ny,nz
      REAL    UNI(nx,ny,nz),VNI(nx,ny,nz),
     &        WNI(nx,ny,nz),
     &        rhsu(nx,ny,nz),rhsv(nx,ny,nz),
     &        rhsw(nx,ny,nz),
     &        wfrg(nx,ny,nz)
C
      INTEGER i,j,k
      REAL    wfrgxpzm, uzp,uzm,wxppzm,wzm,wfrgypzm,vzp,vzm,wyppzm,
     &        rhs_mod,rhsmodav
C 
C-----------------------------------------------------------------------
C                                                          u - component
C-----------------------------------------------------------------------
C

c      open (66,status='unknown',form='unformatted')

c      write(66) IBU,IEU,ix1,ix2,kz1,kz2
      DO I=IBU,IEU
        DO K=KZ1+1,KZ2-1

          rhsmodav=0.

          DO J=JY1+1,JY2
            WFRGXPZM = 0.5 * (0.5 * (wfrg(I+1,J,K) + wfrg(I,J,K))
     &               +   0.5 * (wfrg(I+1,J,K-1) + wfrg(I,J,K-1))) 

            UZP    = (UNI(I,J,K+1) + UNI(I,J,K)) * 0.5
            UZM    = (UNI(I,J,K) + UNI(I,J,K-1)) * 0.5
            WXPPZM = (WNI(I+1,J,K) + WNI(I+1,J,K-1)) * 0.5
            WZM    = (WNI(I,J,K) + WNI(I,J,K-1)) * 0.5

            rhs_mod = - WFRGXPZM * ( CU(K) * (UZP - UZM) / DZ
     &                             - AU(I) * (WXPPZM - WZM) / DX )
            rhsu(I,J,K) = rhsu(I,J,K) + rhs_mod
            rhsmodav=rhsmodav+rhs_mod




          END DO

c          write(66) i,k,rhsmodav/(JY2-JY1)

        END DO
      END DO
C 
C-----------------------------------------------------------------------
C                                                          v - component
C-----------------------------------------------------------------------
C
      DO I=IX1+1,IX2
        DO K=KZ1+1,KZ2-1

          rhsmodav=0.

          DO J=JBV,JEV
            WFRGYPZM = 0.5 * (0.5 * (wfrg(I,J+1,K) + wfrg(I,J,K))
     &               +   0.5 * (wfrg(I,J+1,K-1) + wfrg(I,J,K-1))) 

            VZP    = (VNI(I,J,K+1) + VNI(I,J,K)) * 0.5
            VZM    = (VNI(I,J,K) + VNI(I,J,K-1)) * 0.5
            WYPPZM = (WNI(I,J+1,K) + WNI(I,J+1,K-1)) * 0.5
            WZM    = (WNI(I,J,K) + WNI(I,J,K-1)) * 0.5


            rhs_mod = WFRGYPZM * ( BV(J) * (WYPPZM - WZM) / DY
     &                           - CV(K) * (VZP - VZM) / DZ )
c             if(i.eq.3 .and. j.eq.10) then
c                 write(126,*) k,rhsv(i,j,k),rhs_mod
c             endif


            rhsv(I,J,K) = rhsv(I,J,K) + rhs_mod

            rhsmodav=rhsmodav+rhs_mod

          END DO

c          write(66) i,k,rhsmodav/(JEV-JBV+1)

        END DO
      END DO
C 
C-----------------------------------------------------------------------
C                                                          w - component
C-----------------------------------------------------------------------
C
C.....no modification in w
C


c      close (66)

      RETURN 
      END 






C
C
C-----SUBROUTINE-Frginit------------------------P. Flohr--08/07/1994----
C
      SUBROUTINE Frginit_flohr(wfrg,xu,zw,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - impose fringe velocities
C                                                                       
C     PARAMETERS: - fringe velocity, domain size
C                                                                       
C-----------------------------------------------------------------------
C
c      IMPLICIT NONE
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      INTEGER nx,ny,nz
      REAL    wfrg(nx,ny,nz),xu(nx),zw(nz)
C
      INTEGER i,j,k
      REAL    x,xmax,z,zmax,wmax,func
      real zfrg, xfrg, sfrg

C 
C-----------------------------------------------------------------------
C
      xmax = 0.5 * (xu(ix2+1) + xu(ix2))
      zmax = zw(kz2+1)
C
      DO 10 I=IX1+1, IX2+1
        x = 0.5 * (xu(i) + xu(i-1))
c        write(6,*) 'x, xu, i : ', x, xu(i), i
        DO 10 J=JY1, JY2+1
          DO 10 K=KZ1, KZ2+1
            IF ((x.LT.(2.*xfrg)) .OR. ((xmax-x).LT.(2.*xfrg))) THEN
C.....Spalart
c              wfrg(I,J,K) = ( zw(k) * zw(k) / (zw(k) + zfrg) )
c     &                    * ( EXP(- (x / xfrg)**2 )
c     &                    +   EXP(- ((xmax - x) / xfrg )**2 ) )
C.....Spalart mod.
c              wfrg(I,J,K) = ( zw(k) * zw(k) / (zw(k) + zfrg) )
c     &                    * ( EXP(- 0.15 * zw(k) ) )
c     &                    * ( EXP(- (x / xfrg)**2 )
c     &                    +   EXP(- ((xmax - x) / xfrg )**2 ) )
C.....new 1
c              wfrg(I,J,K) = ( 1 - EXP(- (2. * zw(k)/zfrg)**2) )
c     &                    * ( 1 - EXP(- (2. * (zw(k)-zmax)/zfrg)**2) )
c     &                    * ( EXP(- ((zw(k)-zfrg)/zfrg)**2) )
c     &                    * ( EXP(- (x / xfrg)**2 )
c     &                    +   EXP(- ((xmax - x) / xfrg )**2 ) )
C.....new 2
              IF (zw(k) .LE. zfrg) THEN
                z = zw(k)
                func = 0.5*(1.-TANH(2.*(1.-2.*z/zfrg)))
              ELSE IF (zw(k) .LE. (zfrg + (zmax-zfrg)) ) THEN
                z = zw(k) - zfrg
                func = 0.5*(1.+TANH(2.*(1.-2.*z/(zmax-zfrg))))
              ELSE
                func = 0.
              END IF
              wfrg(I,J,K) = func
     &                    * ( EXP(- (x / xfrg)**2. )
     &                    +   EXP(- ((xmax - x) / xfrg )**2. ) )
C...........
            ELSE
              wfrg(I,J,K) = 0.
            END IF
 10   CONTINUE
C
      DO 15 J=JY1, JY2+1
        DO 15 K=KZ1, KZ2+1
          wfrg(ix1,J,K) = wfrg(ix2+1,j,k)
 15   CONTINUE
C
      wmax = 0.
      DO 20 I=IX1, IX2+1
        DO 20 J=JY1, JY2+1
          DO 20 K=KZ1, KZ2+1
            IF (wmax .LT. wfrg(i,j,k)) wmax = wfrg(i,j,k)
 20   CONTINUE
      IF (wmax .EQ. 0.) WRITE(6,*) '*..error in fringe setup'
      DO 30 I=IX1, IX2+1
        DO 30 J=JY1, JY2+1
          DO 30 K=KZ1, KZ2+1
            wfrg(i,j,k) =  - SFRG * wfrg(i,j,k) / wmax
c            wfrg(i,j,k) =  - SFRG * wfrg(i,j,k)
 30   CONTINUE
C
c      do i =ix1, ix2+1
c        write(6,*) xu(i),wfrg(i,15,8),wfrg(i,15,10),wfrg(i,15,12)
c      end do
C
        
         do i=1,ix2+1
           write(112,*) i,xu(i),wfrg(i,5,5)
         enddo

      RETURN
      END

C
C
C-----SUBROUTINE-Frginit------------------------E. Balaras--28/05/1998----
C
      SUBROUTINE Frginit_ori(wfrg,xu,zw,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - impose fringe velocities
C                                                                       
C     PARAMETERS: - fringe velocity, domain size
C                                                                       
C-----------------------------------------------------------------------
C
c      IMPLICIT NONE
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      INTEGER nx,ny,nz
      REAL    wfrg(nx,ny,nz),xu(nx),zw(nz)
C
      INTEGER i,j,k
      REAL    x,xmax,z,zmax,wmax,func
      real zfrg, sfrg, xfrg
C 
C-----------------------------------------------------------------------
C
c      xmax = 0.5 * (xu(ix2+1) + xu(ix2))
      xmax=xu(ix2)
      zmax = zw(kz2+1)
      DO 10 I=IX1+1, IX2
        x = 0.5 * (xu(i) + xu(i-1))
c        write(6,*) 'x, xu, i : ', x, xu(i), i
        DO 10 J=JY1, JY2+1
          DO 10 K=KZ1+1, KZ2+1
            IF ((x.LT.(2.*xfrg)) .OR. ((xmax-x).LT.(2.*xfrg))) THEN
C.....Spalart
c              wfrg(I,J,K) = ( zw(k) * zw(k) / (zw(k) + zfrg) )
c     &                    * ( EXP(- (x / xfrg)**2 )
c     &                    +   EXP(- ((xmax - x) / xfrg )**2 ) )
C.....Spalart mod.
c              wfrg(I,J,K) = ( zw(k) * zw(k) / (zw(k) + zfrg) )
c     &                    * ( EXP(- 0.15 * zw(k) ) )
c     &                    * ( EXP(- (x / xfrg)**2 )
c     &                    +   EXP(- ((xmax - x) / xfrg )**2 ) )
C.....new 1
c              wfrg(I,J,K) = ( 1 - EXP(- (2. * zw(k)/zfrg)**2) )
c     &                    * ( 1 - EXP(- (2. * (zw(k)-zmax)/zfrg)**2) )
c     &                    * ( EXP(- ((zw(k)-zfrg)/zfrg)**2) )
c     &                    * ( EXP(- (x / xfrg)**2 )
c     &                    +   EXP(- ((xmax - x) / xfrg )**2 ) )
C.....new 2
              IF (zw(k) .LE. zfrg) THEN
                z = zw(k)
                func = 0.5*(1.-TANH(2.5*(1.-2.*z/zfrg)))
              ELSE IF (zw(k) .LE. (zfrg + (zmax-zfrg)) ) THEN
                z = zw(k) - zfrg
                func = 0.5*(1.+TANH(2.5*(1.-2.*z/(zmax-zfrg))))
              ELSE
                func = 0.
              END IF
c
c......new3
c              IF (zw(k) .LE. zfrg) THEN
c                z = zw(k)
c                func = 0.5*(1.-TANH(4.*(1.-2.*z/zfrg)))
c              ELSE IF (zw(k) .LE. (zfrg + (zmax-zfrg)) ) THEN
c                z = zw(k) - zfrg
c                func = 0.5*(1.+TANH(4.*(1.-2.*z/(zmax-zfrg))))
c              ELSE
c                func = 0.
c              END IF
              wfrg(I,J,K) = func
     &                    * ( EXP(- (x / xfrg)**2. )
     &                    +   EXP(- ((xmax - x) / xfrg )**2. ) )

C...........
            ELSE
              wfrg(I,J,K) = 0.
            END IF
 10   CONTINUE
C
c  bc for fringe velocity
      DO 15 J=JY1, JY2+1
        DO 15 K=KZ1, KZ2+1
          wfrg(ix1,J,K) = wfrg(ix2,j,k)
          wfrg(ix2+1,J,K) = wfrg(ix1+1,j,k)
 15   CONTINUE
C
      DO 16 J=JY1, JY2+1
        DO 16 I=IX1, IX2+1
          wfrg(I,J,KZ1) = 0.0
          wfrg(I,J,KZ1) = 0.0
 16    CONTINUE


      wmax = 0.
      DO 20 I=IX1, IX2+1
        DO 20 J=JY1, JY2+1
          DO 20 K=KZ1, KZ2+1
            IF (wmax .LT. wfrg(i,j,k)) wmax = wfrg(i,j,k)
 20   CONTINUE

         write(6,*) 'wmax=',wmax
      IF (wmax .EQ. 0.) WRITE(6,*) '*..error in fringe setup'
      DO 30 I=IX1, IX2+1
        DO 30 J=JY1, JY2+1
          DO 30 K=KZ1, KZ2+1
            wfrg(i,j,k) =  - SFRG * wfrg(i,j,k) / wmax
cc            wfrg(i,j,k) =  - SFRG * wfrg(i,j,k)
 30   CONTINUE
C


        do k=1,kz2+1
          write(123,*) k,zw(k),wfrg(2,20,k),wfrg(30,20,k),wfrg(50,20,k)
        enddo

        do i=1,ix2+1
          write(124,*) xu(i),wfrg(i,20,20),wfrg(i,20,30),wfrg(i,20,40)
        enddo

C
      RETURN
      END

C
C
C-----SUBROUTINE-Frginit------------------------E. Balaras--04/06/1999----
C
      SUBROUTINE Frginit(wfrg,xu,zw,nx,ny,nz)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - impose fringe velocities
C                                                                       
C     PARAMETERS: - fringe velocity, domain size
C                                                                       
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      INTEGER nx,ny,nz
      REAL    wfrg(nx,ny,nz),xu(nx),zw(nz)
      REAL    frgmax(mmz)
C
      INTEGER i,j,k
      REAL    x,xmax,zmax,wmax,func,alpha
      real zfrg, sfrg, xfrg, radius
C 
C-----------------------------------------------------------------------
C
c...alpha controls how fast TANH decays. Increase alpha to increase decay
      alpha =1.7      ! No2 1.8
c      alpha =2.5
      radius=19.6    ! No2 19.8
      xmax =xu(ix2)
      zmax =zw(kz2+1)
c...create function in z (its the maximum at i=2)
      do k=kz1,kz2+1
c        frgmax(k)= 1. -tanh(alpha*(zw(k)-zfrg)**2)
        if(abs((zw(k)-zfrg)).le.radius) then
         frgmax(k)= 1. 
     &        -(tanh(alpha*((zw(k)-zfrg)/radius)**2))/tanh(alpha)
        else
         frgmax(k)= 0.0
        endif
      enddo
c...calculate x distribution
      do i=ix1+1,ix2
       x = 0.5 * (xu(i) + xu(i-1))
       do j=jy1,jy2+1
        do k=kz1,kz2+1
         IF ((x.LT.(2.*xfrg)) .OR. ((xmax-x).LT.(2.*xfrg))) THEN
           wfrg(I,J,K) = - SFRG*frgmax(k) 
     &                   * ( EXP(- (x / xfrg)**2. )
     &                   +   EXP(- ((xmax - x) / xfrg )**2. ) )
         ELSE
           wfrg(I,J,K) = 0.
         ENDIF
        enddo
       enddo
      enddo

c...make periodic in x
      do J=JY1, JY2+1
        do K=KZ1, KZ2+1
          wfrg(ix1,J,K) = wfrg(ix2,j,k)
          wfrg(ix2+1,J,K) = wfrg(ix1+1,j,k)
        enddo
      enddo
c...zero on walls
      do J=JY1, JY2+1
        do I=IX1, IX2+1
          wfrg(I,J,KZ1) = 0.0
          wfrg(I,J,KZ1) = 0.0
        enddo
      enddo

c...write information
        do k=1,kz2+1
          write(123,*) k,zw(k),wfrg(2,20,k),wfrg(30,20,k),wfrg(50,20,k)
        enddo

        do i=1,ix2+1
          write(124,*) xu(i),wfrg(i,20,20),wfrg(i,20,30),wfrg(i,20,40)
        enddo

C
      RETURN
      END
