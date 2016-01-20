
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE FORTERM(UUAV,YF,ZF,QQO,QQOLD,DDPDX,DTM,JMAX,KMAX)
*                                                                 *
*       Computes the forcing term DPDX which represents           *
*       the overall presure gradient in the streamwise            *
*        direction                                                *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      real UUAV(JMAX,KMAX), YF(JMAX), ZF(KMAX)
      real qstar, qqo, qqold, ddpdx, tmp, dtm
      integer i,j,k, jmax, kmax
*
**** Flow rate
*
      QSTAR = 0.0
      DO 200 J=JY1+1,JY2
      DO 200 K=KZ1+1,KZ2
       QSTAR = QSTAR + UUAV(J,K)*(YF(J)-YF(J-1))*(ZF(K)-ZF(K-1))
200   CONTINUE  
*
**** Cross sectional area
*
       TMP = ( YF(JY2)-YF(JY1) ) * ( ZF(KZ2)-ZF(KZ1) )
c       AREA = 1./TMP
*
**** Average streamwise presure gradient
*
c       write(6,*) 'qref=',qqo
c       write(6,*) 'qold=',qqold 
c       write(6,*) 'qnew=',qstar
       DDPDX=DDPDX + ( 2.*(QSTAR-QQO)/dtm-1.*(QQOLD-QQO)/dtm )/tmp
       QQOLD = QSTAR
       print*,'- DDPDX',DDPDX

       RETURN
       END




C 
C-----SUBROUTINE-FORTERMNEW---------------------P. FLOHR--30/10/1993----
C
      SUBROUTINE FORTERMNEW(UAV,YF,ZF,QQO,QQOLD,DDPDX,DTM,
     &                      imax,JMAX,KMAX)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - Computes the forcing term DPDX which represents           *
C                   the overall presure gradient in the streamwise
C                   direction
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C

      real UAV(imax,JMAX,KMAX), YF(JMAX), ZF(KMAX)
      real qstar, qqo, qqold, ddpdx, tmp, dtm
      integer i,j,k, imax, jmax, kmax
*
**** Flow rate
*
      QSTAR = 0.0
      DO 200 i = ix1+1, ix2
        DO 200 J=JY1+1,JY2
          DO 200 K=KZ1+1,KZ2
            QSTAR = QSTAR + UAV(i,J,K)*(YF(J)-YF(J-1))*(ZF(K)-ZF(K-1))
200   CONTINUE  
      qstar = qstar / float(ix2 - ix1)
*
**** Cross sectional area
*
       TMP = ( YF(JY2)-YF(JY1) ) * ( ZF(KZ2)-ZF(KZ1) )
c       AREA = 1./TMP
*
**** Average streamwise presure gradient
*
c       write(6,*) 'qref=',qqo
c       write(6,*) 'qold=',qqold 
c       write(6,*) 'qnew=',qstar
       DDPDX=DDPDX + ( 2.*(QSTAR-QQO)/dtm-1.*(QQOLD-QQO)/dtm )/tmp
       QQOLD = QSTAR
       print*,'NEW- DDPDX',DDPDX

       RETURN
       END



