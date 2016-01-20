* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
*            A D A M S  -  B A S H F O R T H   S T E P            *
*                                                                 *
*    All nesecery routines to perform an A-B step are included    *
*                                                                 *
*           Predictor step  :   USTARAB                           *
*           Divergence      :   DIVAB                             *
*           Velocity corect.:   CORAB                             *
*           Store prev. val.:   STOREAB                           *
*                                                                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      SUBROUTINE predictor(uo,vo,wo,p,uab1,vab1,wab1,
     &                     uab2,vab2,wab2,ucn,vcn,wcn,tv,
     &                     us,vs,ws,zf,dpdx,dpdy,dpdz,dtm,a1,a2,botbc,topbc,
     &                     im,jm,km                               )
c
      include 'headers/common.h'
c
      REAL   uo(im,jm,km),vo(im,jm,km),wo(im,jm,km),
     &       p(im,jm,km),
     &       uab1(im,jm,km),vab1(im,jm,km),wab1(im,jm,km),
     &       uab2(im,jm,km),vab2(im,jm,km),wab2(im,jm,km),
     &       ucn(im,jm,km),vcn(im,jm,km),wcn(im,jm,km),
     &       tv(im,jm,km),
     &       us(im,jm,km),vs(im,jm,km),ws(im,jm,km)
      REAL   zf(km)
c
      REAL   a(mmz),b(mmz),c(mmz),r(mmz),u1d(mmz)
c
      double precision bd(mmz),cd(mmz),ad(mmz),rd(mmz)
c
      REAL   dtm,dpdx,dpdy,dpdz,rflag,a1,a2
      REAL   coef
c
      REAL   botbc(im,jm,3),topbc(im,jm,3)
c
      rflag=real(impfla)
*
**** Compute provisional values for each momentum eq.
*
      coef=   rflag*0.5*dtm/dz**2

c------------------------------------------------------------------
c                                                    compute ustar
c------------------------------------------------------------------
      
      do k=kz1+1,kz2
         do j=jy1+1,jy2
            do i=ibu,ieu
               us(i,j,k)=uo(i,j,k)
     &              +(a2*uab2(i,j,k)-a1*uab1(i,j,k)
     &              -au(i)*(p(i+1,j,k)-p(i,j,k))/dx
     &              -dpdx)*dtm
            enddo
         enddo
      enddo



c------------------------------------------------------------------
c                                                    compute vstar
c------------------------------------------------------------------

      do k=kz1+1,kz2
         do j=jbv,jev
            do i=ix1+1,ix2
               vs(i,j,k)=vo(i,j,k)
     &              +(a2*vab2(i,j,k)-a1*vab1(i,j,k)
     &              -bv(j)*(p(i,j+1,k)-p(i,j,k))/dy-dpdy)*dtm
            enddo
         enddo
      enddo
c------------------------------------------------------------------
c                                                    compute wstar
c------------------------------------------------------------------
      do k=kbw,kew

         do j=jy1+1,jy2
            do i=ix1+1,ix2

               ws(i,j,k)=wo(i,j,k)+
     &              (a2*wab2(i,j,k)-a1*wab1(i,j,k)
     &              -cw(k)*(p(i,j,k+1)-p(i,j,k))/dz-dpdz)*dtm
            enddo
         enddo
      enddo

      RETURN
      END    

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE DIVERGENCE(US,VS,WS,IM,JM,KM,DTM,DIV)
*                                                                 *
*          Computes the divergence for an AB step                 *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      DIMENSION US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),
     &          DIV(IM,JM,KM)
c
      REAL      rflag,coef
*
**** Compute divergence
*
      rflag=real(impfla)
      coef = 1.+ rflag
c
      do 10 i=ix1+1,ix2
       do 20 j=jy1+1,jy2
        do 30 k=kz1+1,kz2
         div(i,j,k)=coef*(ap(i)*(us(i,j,k)-us(i-1,j,k))/dx
     &                   +bp(j)*(vs(i,j,k)-vs(i,j-1,k))/dy
     &                   +cp(k)*(ws(i,j,k)-ws(i,j,k-1))/dz)
     &                   /(dtm)
30      continue
20     continue
10    continue


      RETURN
      END

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE CORRECT(US,VS,WS,DPP,DTM,IM,JM,KM,UUO,VVO,WWO,
     &                 IT,RSU,RSV,RSW)
*                                                                 *
*                Performs an AB corector step                     *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      DIMENSION UUO(IM,JM,KM),VVO(IM,JM,KM),WWO(IM,JM,KM),
     &          US(IM,JM,KM),VS(IM,JM,KM),WS(IM,JM,KM),
     &          DPP(IM,JM,KM)
c
      REAL      rflag,coef
c
      rflag=real(impfla)
      coef = 1./(rflag+1.)
*
**** Update u values
*
        do  k=kz1+1,kz2
           do  j=jy1+1,jy2
              do  i=ibu,ieu


                 uuo(i,j,k)=us(i,j,k)-au(i)*( dpp(i+1,j,k)-dpp(i,j,k) )
     &                             *(coef*dtm/dx )
              enddo
           enddo
        enddo
*
**** Update v values
*
                do  k=kz1+1,kz2
                   do  j=jbv,jev
                      do  i=ix1+1,ix2


                         vvo(i,j,k)=vs(i,j,k)-bv(j)*
     &                        ( dpp(i,j+1,k)-dpp(i,j,k) )
     &                        *( coef*dtm/dy )
                      enddo
                   enddo
                enddo

*
**** Update w values
*

        do k=kbw,kew
           do  j=jy1+1,jy2
              do i=ix1+1,ix2


                 wwo(i,j,k)=ws(i,j,k)-cw(k)*( dpp(i,j,k+1)-dpp(i,j,k) )
     &                             *( coef*dtm/dz )
         
      enddo
      enddo
      enddo

      RETURN
      END


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                 *
      SUBROUTINE STOREN(UUAB1,VVAB1,WWAB1,UUAB2,VVAB2,WWAB2,
     &                   tmeth,IM,JM,KM                      )
*                                                                 *
*              Stores RHS + Presure gradient                      *                
*            to be used in the predictor step                     *
*                                                                 *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'

      DIMENSION UUAB1(IM,JM,KM),VVAB1(IM,JM,KM),WWAB1(IM,JM,KM),
     &          UUAB2(IM,JM,KM),VVAB2(IM,JM,KM),WWAB2(IM,JM,KM)
c
       REAL     tmeth
*
**** Store RHS + presure gradient for 3 momentum equations
*
        do  k=kz1+1,kz2
           do  j=jy1+1,jy2
              do  i=ibu,ieu


                uuab1(i,j,k)= (1.-tmeth)*uuab1(i,j,k)+tmeth*uuab2(i,j,k)
              enddo
           enddo
        enddo

*


        do  k=kz1+1,kz2
           do  j=jbv,jev
              do  i=ix1+1,ix2
                 vvab1(i,j,k)=(1.-tmeth)*vvab1(i,j,k)+tmeth*vvab2(i,j,k)
              enddo
           enddo
        enddo
*
        do  k=kbw,kew
           do  j=jy1+1,jy2
               do  i=ix1+1,ix2
                 wwab1(i,j,k)=(1.-tmeth)*wwab1(i,j,k)+tmeth*wwab2(i,j,k)
               enddo
            enddo
         enddo
      


      RETURN
      END




