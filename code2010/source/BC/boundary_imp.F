C
C-----------------------------------------------------------------------
C                 ***************************                         
C                 *        boun_imp.F       *                        
C                 ***************************                       
C----------------------------------------------------------------------- 
c        THESE ROUTINES INTENDED FOR FLIPPED CHANNEL ONLY
C
C        Boundary conditions for implicit code
C
C-----------------------------------------------------------------------
C
C-----SUBROUTINE-BOUNDARY-----------------------E. Balaras--23/05/1998----
c      Recoding for flipped channel by V. Ovchinnikov 8/11/03 
C
      SUBROUTINE BOUN_IMP(botbc,topbc,imax,jmax,kmax)
      implicit none
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - calls all corresponding routines
C                   variables are defined at lower level routines
C                              
C-----------------------------------------------------------------------
C 
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
      integer imax,jmax,kmax,ibnd
      real  botbc(3,jmax,kmax), topbc(3,jmax,kmax)
C-----------------------------------------------------------------------
C                                  select appropriate boundary condition
C-----------------------------------------------------------------------
C
      DO 10 ibnd=1,2
C
C-----------------------------------------------------------------------
C                                                          non-slip wall
C-----------------------------------------------------------------------
C
        IF( itype(ibnd).eq.50 ) THEN
         CALL BND_IMP050(botbc,topbc,imax,jmax,kmax,ibnd)
c
        ENDIF
C
 10   CONTINUE
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C

*//
*//
*//////// NON-SLIP WALL BOUNDARY CONDITIONS  /////////////////////
*//
*//

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BND_IMP050(botbc,topbc,im,jm,km,ibnd)
*                                                             *
*     Non-slip boundary conditions                            *
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      REAL  botbc(3,jm,km),topbc(3,jm,km)

      integer ibx, iex, jby, jey, kbz, kez, ibnd, im, jm, km, i,j,k
*
****  Starting and ending indexes for ibnd 
*
      ibx=ibegin(ibnd)
      iex=iend(ibnd)
      jby=jbegin(ibnd)
      jey=jend(ibnd)
      kbz=kbegin(ibnd)
      kez=kend(ibnd)
*
**** Main loop
*
c...set all velocity components on the boundary

      if(ibnd.eq.1) then
       do k=kz1,kz2+1
        do j=jy1,jy2+1
          botbc(1,j,k)=0.0 !u
          botbc(2,j,k)=0.0 !v
          botbc(3,j,k)=0.0 !w
        enddo
       enddo
      elseif (ibnd.eq.2) then 
       do k=kz1,kz2+1
        do j=jy1,jy2+1
          topbc(1,j,k)=0.0 !u
          topbc(2,j,k)=0.0 !v
          topbc(3,j,k)=0.0 !w
        enddo
       enddo
      endif

c...set flags for Dirichlet (1) or Neumann (0) 

      ibc(1,ibnd)=1.!u
      ibc(2,ibnd)=1.!v      
      ibc(3,ibnd)=1.!w
c
      RETURN
      END
