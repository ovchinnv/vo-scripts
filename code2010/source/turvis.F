c
c-----------------------------------------------------------------------
c                 ***************************                         
c                 *         turvis.f        *                        
c                 ***************************                       
c----------------------------------------------------------------------- 
c
c	- Turvis:       calls the different eddy viscosity routines
c
c----------------------------------------------------------------------- 
c
c
c-----SUBROUTINE-Turvis------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE Turvis(uo,vo,wo,dxdydz,xf,yf,zf,urun,
     &                  sxx,syy,szz,sxy,sxz,syz,g,
     &                  uoc,voc,woc,upr,vpr,wpr,gbar,
     &                  lm,mm,Ilm,Imm,xc,yc,zc,
     &                  wrk1,wrk2,wrk3,wrk4,
     &                  tv,dt,icycle,im,jm,km         )
c       Driver Routine
c       Calls the different eddy viscosity routines
c
c-----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c-----------------------------------------------------------------------
c
      integer im,jm,km
      real tv(im,jm,km)
      real lm(im,jm,km),mm(im,jm,km),Ilm(im,jm,km),Imm(im,jm,km)
      real uo(im,jm,km),vo(im,jm,km),wo(im,jm,km),dxdydz(im,jm,km),
     &     sxx(im,jm,km),syy(im,jm,km),szz(im,jm,km),sxy(im,jm,km),
     &     sxz(im,jm,km),syz(im,jm,km),g(im,jm,km),uoc(im,jm,km),
     &     voc(im,jm,km),woc(im,jm,km),upr(im,jm,km),vpr(im,jm,km),
     &     wpr(im,jm,km),gbar(im,jm,km),xc(im),yc(jm),zc(km),
     &     xf(im),yf(jm),zf(km),
     &     wrk1(im,jm,km),wrk2(im,jm,km),wrk3(im,jm,km),wrk4(im,jm,km),
     &     dt, urun(im,km)
      INTEGER   ibnd
c     for Spalart-Allmaras model  
      real dtsa
      integer nitsa, nit, ier, icycle
c
c-----------------------------------------------------------------------
c

c----------------------------------------------------------------------
c                        compute turbulent viscosity at cell centers
c----------------------------------------------------------------------

c.....Laminar flow
      IF (    Isgsmodel .EQ. 0 ) THEN
        CALL vtzero(tv,im,jm,km)
c.....Smagorinsky
      ELSEIF (Isgsmodel .EQ. 1 ) THEN
        CALL vtsmag(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     &              dxdydz,im,jm,km                       )
      CALL REFRESHBC(tv,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
c.....Smagorinsky with van Driest damping (one side)
      ELSEIF (Isgsmodel .EQ. 2 ) THEN
         IF(MYRANK.EQ.0) THEN
            print*,'Van Driest damping not implemented 
     &           (use dynamic instead!!!)'
            CALL MPI_FINALIZE(IER)
            STOP
         ENDIF
         
c.....Smagorinsky with van Driest damping (two sides)
      ELSEIF (Isgsmodel .EQ. 3 ) THEN
         IF(MYRANK.EQ.0) THEN
            print*,'Van Driest damping not implmented 
     &           (use dynamic instead!!!)'
            CALL MPI_FINALIZE(IER)
            STOP
         ENDIF
c.....Dynamic Plane average
      ELSEIF (Isgsmodel .EQ. 10) THEN     
        CALL vtdnplane(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     &                 uoc,voc,woc,upr,vpr,wpr,gbar,lm,mm,
     &                 wrk1,wrk2,wrk2,wrk4,
     &                 dxdydz,im,jm,km                       )
        CALL REFRESHBC(tv,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
c.....Dynamic model with lagrangian average
      ELSEIF (Isgsmodel .EQ. 20) THEN  
         
       CALL vtdnlagr(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     &                uoc,voc,woc,upr,vpr,wpr,gbar,
     &                lm,mm,Ilm,Imm,xc,yc,zc,
     &                wrk1,wrk2,wrk3,wrk4,
     &                dxdydz,dt,icycle,im,jm,km              )
        CALL REFRESHBC(tv,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
        CALL REFRESHBC(Ilm,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
        CALL REFRESHBC(Imm,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
c
      ELSEIF (Isgsmodel.eq.30) then ! Spalart-Allmaras model 
        nitsa=1
	dtsa=dt/nitsa
	do nit=1,nitsa
	  call vtsaimp(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,Ilm,Imm,tv,
     &                 mm,xc,yc,zc,xf,yf,zf,dxdydz,dtsa,icycle,
     &                 im,jm,km)
        enddo
c     
c      return
      END IF

c----------------------------------------------------------------------
c                        set boundary values for turbulent viscosity
c----------------------------------------------------------------------
      
      
      DO 10 ibnd=1,6
       IF(       itype(ibnd).eq.50 ) THEN
          CALL BOUNDWALL( tv,im,jm,km,ibnd)
       
       ELSEIF( itype(ibnd).eq.500 ) THEN
        CALL BOUNDPERIOD(tv,im,jm,km,ibnd)

       ELSEIF(itype(ibnd).GT.0) then !do not touch the internal ghost layers!
          CALL BOUNDGRAD(tv,im,jm,km,ibnd)
       ENDIF
10    CONTINUE
c
      
      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE vtzero(tv,im,jm,km)

c          Sets turbulent viscosity to zero for 
c          LAMINAR FLOW calculations
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      REAL      tv(im,jm,km)
      INTEGER   i, j, k, im, jm, km
c
c----------------------------------------------------------------------
c
        do k=kz1,kz2+1
           do j=jy1,jy2+1
              do i=ix1,ix2+1


         tv(i,j,k)=0.0
        enddo
       enddo
      enddo
c
      RETURN 
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE vtsmag(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     %                  dxdydz,imax,jmax,kmax                 )


c          Computes turbulent viscosity at the cell center 
c          using SMAGORISKY Model         
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
       integer imax, jmax, kmax
       REAL  uo(imax,jmax,kmax),vo(imax,jmax,kmax),
     %       wo(imax,jmax,kmax)
       REAL  sxx(imax,jmax,kmax),syy(imax,jmax,kmax),
     %       szz(imax,jmax,kmax),sxy(imax,jmax,kmax),
     %       sxz(imax,jmax,kmax),syz(imax,jmax,kmax)
       REAL  g(imax,jmax,kmax),tv(imax,jmax,kmax)
       REAL  dxdydz(imax,jmax,kmax)
c
      INTEGER   i, j, k
c
c----------------------------------------------------------------------
c

c----------------------------------------------------------------------
c                                     calculate Strain Rates and |S|
c----------------------------------------------------------------------
      CALL Strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,
     %            imax,jmax,kmax                     ) 

c----------------------------------------------------------------------
c                                calculate turbulent viscosity
c                                 nut=(cs^2)*(Del^2)*|S| --> vt 
c----------------------------------------------------------------------
        do k=2,kz2
           do j=2,jy2
              do i=2,ix2


                 tv(i,j,k)=css*css*dxdydz(i,j,k)*g(i,j,k)
              enddo
           enddo
        enddo
c
      RETURN 
      END
c
c----------------------------------------------------------------------
c
      SUBROUTINE vtsmagd(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     *                   dxdydz,zf,urun,imax,jmax,kmax              )


c          Computes turbulent viscosity at the cell center 
c          using SMAGORISKY Model with Van Driest damping at the
c          bottom wall         
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
       integer imax, jmax, kmax
       REAL  uo(imax,jmax,kmax),vo(imax,jmax,kmax),
     %       wo(imax,jmax,kmax)
       REAL  sxx(imax,jmax,kmax),syy(imax,jmax,kmax),
     %       szz(imax,jmax,kmax),sxy(imax,jmax,kmax),
     %       sxz(imax,jmax,kmax),syz(imax,jmax,kmax)
       REAL  g(imax,jmax,kmax),tv(imax,jmax,kmax)
       REAL  dxdydz(imax,jmax,kmax)
       REAL  zf(kmax), urun(imax,kmax)
c
      INTEGER   i, j, k
      REAL      alpha, aplus, coef1, coef2, h1, h2, utau
      REAL      zeta, zplus, dudz
c
      INTRINSIC exp
c
c----------------------------------------------------------------------
c

c----------------------------------------------------------------------
c                                     calculate Strain Rates and |S|
c----------------------------------------------------------------------
      CALL Strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,
     %            imax,jmax,kmax                     ) 

c----------------------------------------------------------------------
c                                calculate turbulent viscosity
c                    nut=(cs^2)*(Del^2)*((1-exp(-z+/A+))^2)*|S| --> vt 
c----------------------------------------------------------------------
      aplus = 25.0
c
c     Compute compute friction velocity 
c
      h1 = 0.5*(zf(2) + zf(1))
      h2 = 0.5*(zf(3) + zf(2))
c
      coef1 = h2/(h1*(h2-h1))
      coef2 = h1/(h2*(h2-h1))
c
c
c
      do k=2, kz2
        do j=2,jy2
          do i=2,ix2
            dudz  = abs(coef1*urun(i,2) - coef2*urun(i,3)) 
            utau  = sqrt(ru1*dudz)
            zeta  = 0.5*(zf(k)+zf(k-1))
            zplus = utau*zeta/ru1
            alpha = (1.0 - exp(-zplus/aplus))**2
            tv(i,j,k) = css*css*dxdydz(i,j,k)*alpha*g(i,j,k)
          end do
        end do
      end do
c
      RETURN 
      END
c
c
c----------------------------------------------------------------------
c
      SUBROUTINE vtsmagd2(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     *                   dxdydz,zf,imax,jmax,kmax              )


c          Computes turbulent viscosity at the cell center 
c          using SMAGORISKY Model with Van Driest damping at the
c          bottom & top walls         
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
       integer imax, kmax, jmax
       REAL  uo(imax,jmax,kmax),vo(imax,jmax,kmax),
     %       wo(imax,jmax,kmax)
       REAL  sxx(imax,jmax,kmax),syy(imax,jmax,kmax),
     %       szz(imax,jmax,kmax),sxy(imax,jmax,kmax),
     %       sxz(imax,jmax,kmax),syz(imax,jmax,kmax)
       REAL  g(imax,jmax,kmax),tv(imax,jmax,kmax)
       REAL  dxdydz(imax,jmax,kmax)
       REAL  zf(kmax)
       REAL  udum(mmz)
c
      INTEGER   i, j, k
      REAL      alpha, aplus, chahw, coef1, coef2, coef3, coef4 
      REAL      dudzb, dudzt, h1, h2, h3, h4, utaub, utaut, utau
      REAL      zbottom, zeta, zplus, ztop
      REAL      fjp
c
      INTRINSIC exp
c
c----------------------------------------------------------------------
c

c----------------------------------------------------------------------
c                                     calculate Strain Rates and |S|
c----------------------------------------------------------------------
      CALL Strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,
     %            imax,jmax,kmax                     ) 

c----------------------------------------------------------------------
c                                calculate turbulent viscosity
c                    nut=(cs^2)*(Del^2)*((1-exp(-z+/A+))^2)*|S| --> vt 
c----------------------------------------------------------------------
      zbottom  = zf(1)
      ztop     = zf(kz2)
      chahw    = 0.5*(ztop-zbottom)
      aplus    = 25.0
c
c     plane averaged streamwise velocity
c
        fjp  = 1.0/real((ix2-1)*(jy2-1))
        do k = 2, kz2
         udum(k) = 0.0
          do i = ix1+1, ix2
            do j = jy1+1, jy2
              udum(k) = udum(k) + uo(i,j,k)
            end do
          end do
            udum(k) = udum(k)*fjp
        end do
c
c     Compute compute friction velocity 
c
      h1 = 0.5*(zf(2) + zf(1))
      h2 = 0.5*(zf(3) + zf(2))
c
      coef1 = h2/(h1*(h2-h1))
      coef2 = h1/(h2*(h2-h1))
c
      h3 =  0.5*(zf(kz2)   + zf(kz2-1))  - ztop
      h4 =  0.5*(zf(kz2-1) + zf(kz2-2))  - ztop
c
      coef3 = h4/(h3*(h4-h3))
      coef4 = h3/(h4*(h4-h3))
c
c
      do k=2, kz2
        do j=2,jy2
          do i=2,ix2
            dudzb =   coef1*udum(2)   - coef2*udum(3) 
            dudzt = - coef3*udum(kz2) + coef4*udum(kz2-1)
            utaub = sqrt(ru1*dudzb)
            utaut = sqrt(ru1*dudzt)
            utau  = 0.5*(utaub + utaut)
            zeta  = 0.5*(zf(k)+zf(k-1))
            zplus = utau*(1.0-abs(zeta - chahw)/chahw)/ru1
            alpha = (1.0 - exp(-zplus/aplus))**2
            tv(i,j,k) = css*css*dxdydz(i,j,k)*alpha*g(i,j,k)
          end do
        end do
      end do
c
      RETURN 
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE vtdnplane(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     &                     uoc,voc,woc,upr,vpr,wpr,gbar,lm,mm,
     &                     wrk1,wrk2,wrk3,wrk4,
     &                     dxdydz,im,jm,km                       )

c          Computes turbulent viscosity at the cell center 
c          using SMAGORISKY Model         
c          FLIPPED CHANNEL VERSION ( V.OVCHINNIKOV)
c----------------------------------------------------------------------
c
      include 'mpif.h'
      include 'headers/common.h'
      include 'headers/dimension.h'
c
c----------------------------------------------------------------------
c
      REAL      sxx(im,jm,km),syy(im,jm,km),    
     &          szz(im,jm,km),sxy(im,jm,km),
     &          sxz(im,jm,km),syz(im,jm,km),
     &          tv(im,jm,km),g(im,jm,km),
     &          dxdydz(im,jm,km),
     &          uo(im,jm,km),vo(im,jm,km),
     &          wo(im,jm,km)
c
      REAL      uoc(im,jm,km),voc(im,jm,km),
     &          woc(im,jm,km),upr(im,jm,km),
     &          vpr(im,jm,km),wpr(im,jm,km),
     &          gbar(im,jm,km),lm(im,jm,km),
     &          mm(im,jm,km),
     &          wrk1(im,jm,km),wrk2(im,jm,km),
     &          wrk3(im,jm,km),wrk4(im,jm,km)
c
      INTEGER   i, j, k, ivt, im, km, jm
      REAL      r,cut, sc
      REAL      lij,mij
      REAL      lmav(mmx),mmav(mmx),csav(mmx)
      real      lmavg(mmx,ncpu),mmavg(mmx,ncpu),lmvg(ncpu*mmx),
     &          mmvg(ncpu*mmx),dummy1,dummy2
      
      integer type,ierr
c
      equivalence (lmavg,lmvg), (mmavg,mmvg)
c      
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
#ifdef SINGLE
      type=MPI_REAL
#endif
            
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c                                                    initialize arrays
c----------------------------------------------------------------------
      do k=1,kz2+1
       do j=1,jy2+1
        do i=1,ix2+1
         lm(i,j,k)=0.0
         mm(i,j,k)=0.0
        enddo
       enddo
      enddo 
c
      ivt= 0 
      cut= 1.0e-12    

c----------------------------------------------------------------------
c                      square of the ratio betwenn test and grid filter
c----------------------------------------------------------------------
      r=6. 

c----------------------------------------------------------------------
c                             calculate velocity field at cell centers
c----------------------------------------------------------------------
c      cpub=mpi_wtime()
      CALL center(uo,vo,wo,uoc,voc,woc,im,jm,km)

c----------------------------------------------------------------------
c                       calculate Strain Rates and |S| at cell centers
c----------------------------------------------------------------------
      CALL Strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,
     %            im,jm,km                            ) 
C NEED TO UPDATE THE BOUNDARY POINTS!
       CALL REFRESHBC(uoc,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(voc,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(woc,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)

       CALL REFRESHBC(sxx,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(syy,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(szz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(sxy,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(sxz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(syz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(g,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)

c       cpue=mpi_wtime()-cpub
c       if(myrank.eq.0) then
c          print*,'*******************',cpue,'************'
c       endif
c----------------------------------------------------------------------
c                                   calculate numerator and denominator
c----------------------------------------------------------------------
c.....Component XX
c...uu --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=uoc(i,j,k)*uoc(i,j,k)

        enddo
       enddo
      enddo

c...\hat{uu} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(u) --> vpr
      CALL filter3d(uoc,vpr,wrk1,wrk2,im,jm,km,2)
c...\hat{|S|} --> gbar
      CALL filter3d(g,gbar,wrk1,wrk2,im,jm,km,2) 
c-----
c... |S|Sxx --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*sxx(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Sxx} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxx}  --> wrk4
      CALL filter3d(sxx,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*vpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k)
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lij*mij
         mm(i,j,k)= mij*mij
        enddo
       enddo
      enddo 
c.....Component XY
c...uv --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=uoc(i,j,k)*voc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{uv} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(v) --> vpr
      CALL filter3d(voc,wpr,wrk1,wrk2,im,jm,km,2)
c... |S|Sxy --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*sxy(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Sxy} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxy}  --> wrk4
      CALL filter3d(sxy,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+2.*lij*mij
         mm(i,j,k)= mm(i,j,k)+2.*mij*mij
        enddo
       enddo
      enddo 
c.....Component XZ
c...uw --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=uoc(i,j,k)*woc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{uw} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(w) --> vpr
      CALL filter3d(woc,wpr,wrk1,wrk2,im,jm,km,2)
c... |S|Sxz --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*sxz(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Sxz} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxz}  --> wrk4
      CALL filter3d(sxz,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+2.*lij*mij
         mm(i,j,k)= mm(i,j,k)+2.*mij*mij
        enddo
       enddo
      enddo 
c.....Component ZZ
c...ww --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=woc(i,j,k)*woc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{ww} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c... |S|Szz --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*szz(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Szz} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Szz}  --> wrk4
      CALL filter3d(szz,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-wpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+lij*mij
         mm(i,j,k)= mm(i,j,k)+mij*mij
        enddo
       enddo
      enddo 
c.....Component YZ
c...vw --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=voc(i,j,k)*woc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{vw} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(v) --> vpr
      CALL filter3d(voc,vpr,wrk1,wrk2,im,jm,km,2)
c... |S|Syz --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*syz(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Syz} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Syz}  --> wrk4
      CALL filter3d(syz,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+2.*lij*mij
         mm(i,j,k)= mm(i,j,k)+2.*mij*mij
        enddo
       enddo
      enddo 
c.....Component YY
c...vv --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=voc(i,j,k)*voc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{vv} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c... |S|Syy --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*syy(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Syy} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxz}  --> wrk4
      CALL filter3d(syy,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*vpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+lij*mij
         mm(i,j,k)= mm(i,j,k)+mij*mij
        enddo
       enddo
      enddo
 
c----------------------------------------------------------------------
c                                                       calculate cs^2
c----------------------------------------------------------------------
      sc=1./float( (kz2-1)*(jy2-1) )
      do i=2,ix2
       lmav(i)=0.0
       mmav(i)=0.0     
       do k=2,kz2
        do j=2,jy2
         lmav(i)= lmav(i)+lm(i,j,k)
         mmav(i)= mmav(i)+mm(i,j,k)
        enddo
       enddo
        lmav(i)=lmav(i)*sc
        mmav(i)=mmav(i)*sc
      enddo
c     for flipped channel
      call MPI_ALLGATHER(lmav,nx,type,lmvg,nx,type, MPI_COMM_HYDRO, ierr)
      call MPI_ALLGATHER(mmav,nx,type,mmvg,nx,type, MPI_COMM_HYDRO, ierr)
      do i=1,nx
       dummy1=0.
       dummy2=0.
c
       do k=1,ncpu
        dummy1=dummy1+lmavg(i,k)
        dummy2=dummy2+mmavg(i,k)
       enddo
       lmav(i)=dummy1/ncpu	
       mmav(i)=dummy2/ncpu	
      enddo       
      
c
      do i=2,ix2
        csav(i)=-0.5*(lmav(i)/(mmav(i)+cut))
c        write(6,*) k,lmav(k),mmav(k)
      enddo

c----------------------------------------------------------------------
c                                         calculate turbulent viscosity
c----------------------------------------------------------------------
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         tv(i,j,k)=csav(i)*dxdydz(i,j,k)*g(i,j,k)
         wrk4(i,j,k)=csav(i)
        enddo
       enddo
      enddo

c...check for negative total viscosities
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         if( (tv(i,j,k)+ru1) .lt. 0.) then
           tv(i,j,k) = 0.0
           ivt=ivt+1
         endif
        enddo
       enddo
      enddo
c
      if(ivt.ne.0) then 
       write(6,*) 'negative total viscosity values = ',ivt
      endif
c
      RETURN 
      END
c
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE vtdnlagr(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,tv,
     &                    uoc,voc,woc,upr,vpr,wpr,gbar,
     &                    lm,mm,Ilm,Imm,xc,yc,zc,
     &                    wrk1,wrk2,wrk3,wrk4,
     &                    dxdydz,dt,icycle,im,jm,km              )

c          Computes turbulent viscosity at the cell center 
c          using Lagrangian Dynamic model         
c
c----------------------------------------------------------------------
c
      include 'mpif.h'
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      REAL      sxx(im,jm,km),syy(im,jm,km),    
     &          szz(im,jm,km),sxy(im,jm,km),
     &          sxz(im,jm,km),syz(im,jm,km),
     &          tv(im,jm,km),g(im,jm,km),
     &          dxdydz(im,jm,km),
     &          uo(im,jm,km),vo(im,jm,km),
     &          wo(im,jm,km)
c
      REAL      uoc(im,jm,km),voc(im,jm,km),
     &          woc(im,jm,km),upr(im,jm,km),
     &          vpr(im,jm,km),wpr(im,jm,km),
     &          gbar(im,jm,km)

      REAL      lm(im,jm,km),mm(im,jm,km),
     &          Ilm(im,jm,km),Imm(im,jm,km)

      REAL      wrk1(im,jm,km),wrk2(im,jm,km),
     &          wrk3(im,jm,km),wrk4(im,jm,km)
      
      REAL      xc(im),yc(jm),zc(km)
c
      INTEGER   i, j, k, ivt, im, jm, km, ii, ibase, ibnd
      INTEGER   icycle
      REAL      dt
      REAL      r,base,eps,coeff,dtt,cut,cslag
      REAL      lij,mij
      REAL      lmav(mmx,mmz),mmav(mmx,mmz),csav(mmz),tvav(mmz)
c...immersed boundary arrays
      real dum, tconst, sc
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c                                                    initialize arrays
c----------------------------------------------------------------------
c
      do k=1,kz2+1
       do j=1,jy2+1
        do i=1,ix2+1
         lm(i,j,k)=0.0
         mm(i,j,k)=0.0
        enddo
       enddo
      enddo 
c
      ivt=0     

c----------------------------------------------------------------------
c                      square of the ratio betwenn test and grid filter
c----------------------------------------------------------------------
      r  =6. 
      cut=1.e-12

c----------------------------------------------------------------------
c                             calculate velocity field at cell centers
c----------------------------------------------------------------------
      CALL center(uo,vo,wo,uoc,voc,woc,im,jm,km)
c----------------------------------------------------------------------
c                       calculate Strain Rates and |S| at cell centers
c----------------------------------------------------------------------
      CALL Strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,
     %            im,jm,km                            ) 

      CALL REFRESHBC(uoc,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(voc,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(woc,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      
      CALL REFRESHBC(sxx,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(syy,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(szz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(sxy,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(sxz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(syz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(g,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)


      if(icycle.eq.0.and.read_ilm_imm) goto 123
c----------------------------------------------------------------------
c                                   calculate numerator and denominator
c----------------------------------------------------------------------
c.....Component XX
c...uu --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=uoc(i,j,k)*uoc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{uu} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(u) --> vpr
      CALL filter3d(uoc,vpr,wrk1,wrk2,im,jm,km,2)
c...\hat{|S|} --> gbar
      CALL filter3d(g,gbar,wrk1,wrk2,im,jm,km,2) 
c... |S|Sxx --> wrk3

      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*sxx(i,j,k)
        enddo

       enddo
      enddo
      ii=ix2/2
c...\hat{|S|Sxx} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxx}  --> wrk4
      CALL filter3d(sxx,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*vpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k)
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lij*mij
         mm(i,j,k)= mij*mij
        enddo
c            write(2000+myrank,'(2I,6G15.5)') k,j,upr(ii,j,k),vpr(ii,j,k),dxdydz(ii,j,k),gbar(ii,j,k),wrk4(ii,j,k),wrk3(ii,j,k)
       enddo
      enddo 
c.....Component XY
c...uv --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=uoc(i,j,k)*voc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{uv} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(v) --> vpr
      CALL filter3d(voc,wpr,wrk1,wrk2,im,jm,km,2)
c... |S|Sxy --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*sxy(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Sxy} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxy}  --> wrk4
      CALL filter3d(sxy,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+2.*lij*mij
         mm(i,j,k)= mm(i,j,k)+2.*mij*mij
        enddo
       enddo
      enddo 
c.....Component XZ
c...uw --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=uoc(i,j,k)*woc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{uw} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(w) --> vpr
      CALL filter3d(woc,wpr,wrk1,wrk2,im,jm,km,2)
c... |S|Sxz --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*sxz(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Sxz} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxz}  --> wrk4
      CALL filter3d(sxz,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+2.*lij*mij
         mm(i,j,k)= mm(i,j,k)+2.*mij*mij
        enddo
       enddo
      enddo 
c.....Component ZZ
c...ww --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=woc(i,j,k)*woc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{ww} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c... |S|Szz --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*szz(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Szz} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Szz}  --> wrk4
      CALL filter3d(szz,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-wpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+lij*mij
         mm(i,j,k)= mm(i,j,k)+mij*mij
        enddo
       enddo
      enddo 
c.....Component YZ
c...vw --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=voc(i,j,k)*woc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{vw} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c...\hat(v) --> vpr
      CALL filter3d(voc,vpr,wrk1,wrk2,im,jm,km,2)
c... |S|Syz --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*syz(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Syz} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Syz}  --> wrk4
      CALL filter3d(syz,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*wpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+2.*lij*mij
         mm(i,j,k)= mm(i,j,k)+2.*mij*mij
        enddo
       enddo
      enddo 
c.....Component YY
c...vv --> upr
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         upr(i,j,k)=voc(i,j,k)*voc(i,j,k)
        enddo
       enddo
      enddo
c...\hat{vv} --> upr
      CALL filter3d(upr,upr,wrk1,wrk2,im,jm,km,2)
c... |S|Syy --> wrk3
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         wrk3(i,j,k)=g(i,j,k)*syy(i,j,k)
        enddo
       enddo
      enddo
c...\hat{|S|Syy} --> wrk3
      CALL filter3d(wrk3,wrk3,wrk1,wrk2,im,jm,km,2) 
c... \hat{Sxz}  --> wrk4
      CALL filter3d(syy,wrk4,wrk1,wrk2,im,jm,km,2)
c...accumulate into numerator and denominator
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         lij=upr(i,j,k)-vpr(i,j,k)*vpr(i,j,k)
         mij=dxdydz(i,j,k)*(r*gbar(i,j,k)*wrk4(i,j,k) 
     %                      -wrk3(i,j,k)             )
         lm(i,j,k)= lm(i,j,k)+lij*mij
         mm(i,j,k)= mm(i,j,k)+mij*mij
        enddo
       enddo
      enddo

c----------------------------------------------------------------------
c                                                    Lagrangian Average
c----------------------------------------------------------------------
      if(icycle.eq.0) then
c...get average numerator and denominator
          sc=1./(jy2-1)
          do k=1,kz2+1
           do i=2,ix2
            lmav(i,k)=0.0
            mmav(i,k)=0.0     
            do j=2,jy2
             lmav(i,k)= lmav(i,k)+lm(i,j,k)
             mmav(i,k)= mmav(i,k)+mm(i,j,k)
            enddo
            lmav(i,k)= lmav(i,k)*sc
            mmav(i,k)= mmav(i,k)*sc
           enddo
          enddo
c...get values of integrals
           do k=1,kz2+1
            do j=2,jy2
             do i=2,ix2
              wrk1(i,j,k)=lmav(i,k)
              wrk2(i,j,k)=mmav(i,k)
             enddo
            enddo
           enddo
c     initialize ghost points
           do i=1,ix2+1
	    do j=1,jy2+1
	     Ilm(i,j,1)=0.
	     Imm(i,j,1)=0.
	     Ilm(i,j,kz2+1)=0.
	     Imm(i,j,kz2+1)=0.
	    enddo
	   enddo     
           do k=1,kz2+1
	    do j=1,jy2+1
	     Ilm(1,j,k)=0.
	     Imm(1,j,k)=0.
	     Ilm(ix2+1,j,k)=0.
	     Imm(ix2+1,j,k)=0.
	    enddo
	   enddo     
           do i=1,ix2+1
	    do k=1,kz2+1
	     Ilm(i,1,k)=0.
	     Imm(i,1,k)=0.
	     Ilm(i,jy2+1,k)=0.
	     Imm(i,jy2+1,k)=0.
	    enddo
	   enddo     
c     (initialize ghost points)	   
      else
c...Interplolate numerator and denominator
       CALL lagrange(Ilm,Imm,wrk1,wrk2,
     &               uoc,voc,woc,xc,yc,zc,dt,im,jm,km)

CCCCCCCCCCC NEED TO UPDATE BOUNDARY VALUES OF wrk1 and wrk2
	 

       CALL REFRESHBC(WRK1,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
       CALL REFRESHBC(WRK2,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
c       boundary conditions on Wrk1, Wrk2
       if (myrank.eq.0) then 
        k=1
        do i=2,ix2+1
	 do j=2,jy2+1
	  wrk1(i,j,k)=wrk1(i,j,k+1)
	  wrk2(i,j,k)=wrk2(i,j,k+1)
	 enddo
	enddo  
       elseif (myrank.eq.mysize-1) then 
        k=kz2+1	
        do i=2,ix2+1
	 do j=2,jy2+1
	  wrk1(i,j,k)=wrk1(i,j,k-1)
	  wrk2(i,j,k)=wrk2(i,j,k-1)
	 enddo
	enddo 
      endif	 

c...Calculate T^n, eps
       ibase=0
       do k=1,kz2+1
        do j=2,jy2
         do i=2,ix2
          base=-8.*Ilm(i,j,k)*Imm(i,j,k)
          if(base.le.cut) then
           ibase=ibase+1
           eps=0.0
          else
           coeff=1.5*sqrt(dxdydz(i,j,k))
           tconst=coeff/base**0.125
           dtt   =dt/tconst
           eps   =dtt/(1.+dtt)
          endif
           wrk1(i,j,k)=eps*lm(i,j,k)+(1.-eps)*wrk1(i,j,k)
           wrk2(i,j,k)=eps*mm(i,j,k)+(1.-eps)*wrk2(i,j,k)
         enddo
        enddo
       enddo
c
c          write(6,*) 'eps=0',ibase
c
      endif 

c...update the value of the integral (truncate Ilm<0)
       do k=1,kz2+1
        do j=2,jy2
         do i=2,ix2
           Ilm(i,j,k)=0.5*( wrk1(i,j,k)-abs(wrk1(i,j,k)) )
           Imm(i,j,k)=wrk2(i,j,k)
         enddo
        enddo
       enddo

************************************************************************
****  NEED TO CHECK THIS PART  **************************************
************************************************************************
       DO 10 ibnd=1,6
          IF(       itype(ibnd).eq.50.or.itype(ibnd).eq.60 ) THEN
             CALL BOUNDZERO( Ilm,im,jm,km,ibnd)
             CALL BOUNDGRAD( Imm,im,jm,km,ibnd)
             
          ELSEIF( itype(ibnd).eq.500 ) THEN

             CALL BOUNDPERIOD(Ilm,im,jm,km,ibnd)
             CALL BOUNDPERIOD(Imm,im,jm,km,ibnd)

          ELSEIF(itype(ibnd).GT.0) then !do not touch the internal ghost layers!

             CALL BOUNDGRAD(Ilm,im,jm,km,ibnd)
             CALL BOUNDGRAD(Imm,im,jm,km,ibnd)
          ENDIF
 10    ENDDO

c...  boundary values for Inm, Imm
c   (only valid for time-developing mixing layer)
c..periodic
c       do k=1,kz2+1
c        do j=1,jy2+1
c           Ilm(1,j,k)=Ilm(ix2,j,k)
c           Imm(1,j,k)=Imm(ix2,j,k)
c           Ilm(ix2+1,j,k)=Ilm(2,j,k)
c           Imm(ix2+1,j,k)=Imm(2,j,k)
c        enddo
c       enddo
c...modification for bl
c       do k=1,kz2+1
c        do j=1,jy2+1
c           Ilm(1,j,k)=Ilm(2,j,k)
c           Imm(1,j,k)=Imm(2,j,k)
c           Ilm(ix2+1,j,k)=Ilm(ix2,j,k)
c           Imm(ix2+1,j,k)=Imm(ix2,j,k)
c        enddo
c       enddo

c       do j=1,jy2+1
c        do i=1,ix2+1
c           Ilm(i,j,1)=0.0
c           Imm(i,j,1)=Imm(i,j,2)
c           Ilm(i,j,kz2+1)=0.0
c           Imm(i,j,kz2+1)=Imm(i,j,kz2)
c        enddo
c       enddo

c       do k=1,kz2+1
c        do i=1,ix2+1
c           Ilm(i,1,k)=Ilm(i,jy2,k)
c           Imm(i,1,k)=Imm(i,jy2,k)
c           Ilm(i,jy2+1,k)=Ilm(i,2,k)
c           Imm(i,jy2+1,k)=Imm(i,2,k)
c        enddo
c       enddo

*********************************************************************************
*********************************************************************************
c      return

c...calculate turbulent viscosity
 123  continue
      do k=1,kz2+1
       do j=2,jy2
        do i=2,ix2
         cslag=-0.5*(Ilm(i,j,k)/(Imm(i,j,k)+cut))
         wrk4(i,j,k)=cslag
         tv(i,j,k)  =cslag*dxdydz(i,j,k)*g(i,j,k)
        enddo
       enddo
      enddo

c----------------------------------------------------------------------
c                                 truncate negative turbulent viscosity
c----------------------------------------------------------------------
c      do k=1,kz2+1
c       do j=2,jy2
c        do i=2,ix2
c         if( (tv(i,j,k)+ru1) .lt. 0.) then
c         write(6,*) i,j,k,tv(i,j,k),-1.*Ilm(i,j,k),Imm(i,j,k)
c           tv(i,j,k) = 0.0
c           ivt=ivt+1
c         endif
c        enddo
c       enddo
c      enddo
c
c      if(ivt.ne.0) then 
c       write(6,*) 'negative total viscosity values = ',ivt
c      endif
c
      RETURN 
      END
c
c     ------------------------------------------ U Piomelli 5/10/03 -------
c
c     Spalart-Allmaras eddy viscosity 
c
c     ---------------------------------------------------------------------
c
      subroutine vtsaimp(uo,vo,wo,
     $     sxx,syy,szz,sxy,sxz,syz,sabs,
     $     vt,vtnew,tv,rhsnuo,
     &     xc,yc,zc,xf,yf,zf,
     &     dxdydz,dt,icycle,im,jm,km)
c
c----------------------------------------------------------------------
c     
      implicit none
      include 'headers/common.h'
      include 'mpif.h'
      include 'headers/dimension.h'
c     
c----------------------------------------------------------------------
c     
      integer i,j,k,im,jm,km,icycle,type,ierr,status(MPI_STATUS_SIZE)
c
      REAL uo (im,jm,km),vo (im,jm,km),wo (im,jm,km)
      REAL vt (im,jm,km),vtnew(im,jm,km),tv(im,jm,km),rhsnuo(im,jm,km)
      REAL sxx(im,jm,km),syy(im,jm,km)
      REAL szz(im,jm,km),sxy(im,jm,km)
      REAL sxz(im,jm,km),syz(im,jm,km)
      REAL sabs(im,jm,km)
      REAL dxdydz(im,jm,km)
c      
      REAL xc(im),yc(jm),zc(km),xf(im),yf(jm),zf(km)
c
c     
      real flag,coflhs,abash1,abash2,sigma,oosig,ck,cke2,cw1,
     &     cb1,cb2,cw2,cw3,ct3,ct4,cv1,cw3e6,cv1e3,dt,cdes
      real oodxcp(mmx),oodxce(mmx),oodxp(mmx)
      real oodycp(mmy),oodycf(mmy),oodyp(mmy)
      real oodzcp(mmz),oodzcn(mmz),oodzp(mmz)
      
      real dist1,dist2,dist,vte,vtw,vtf,vtb,vtn,vts,
     &     chi,fv1,fv2,ft2,stilde,r,g,fw,
     &     bprod,bdest,flyf,flyb,flzn,flzs,vdiffexp,
     &     avdiff,cvdiff,bvdiff,coef,diffexp,adiff,cdiff,
     &     advexp,hn,hnm1
c     
       
      logical impsa,init
      real xnu(mmx),anu(mmx),bnu(mmx),cnu(mmx),fnu(mmx)

      save oodxcp,oodxce,oodxp,oodycp,oodycf,oodyp,oodzcp,oodzcn,oodzp,init
      data init /.true./
c
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
#ifdef SINGLE
      type=MPI_REAL
#endif
c
      if(icycle .lt. 1) return ! tv and vt are assumed available for 1st timestep
c
c     impsa = .false.
      impsa = .true.
      if (impsa) then
         flag   = 0.
         coflhs = 0.5*dt
      else
         flag   = 1.
         coflhs = 0.
      endif
c
      call strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,sabs,im,jm,km)
c
      CALL REFRESHBC(sxx,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(syy,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(szz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(sxy,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(sxz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(syz,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(sabs,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)

c     begin the time-advancement for nutilde
c
c      write(0,*) myrank,' inside SA...'
c     debug      
      if(icycle .eq. 1) then 
         abash1   = 0.
         abash2   = 1.
         do k=1,kz2+1
            do j=1,jy2+1
               do i=1,ix2+1
                  rhsnuo(i,j,k) = 0.
               enddo
            enddo
         enddo
c      boundary conditions
c      write(0,*) 'icycle is 1 -- Euler',icycle
c      goto 998 ! skip boundary conditions
c      
        call REFRESHBC(tv,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO) 
        call REFRESHBC(vt,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO) ! vt
        if (myrank.eq.0) then 
	  call MPI_RECV(vt(1,1,1),im*jm,type,mysize-1,0,MPI_COMM_HYDRO,status,ierr)
	  call MPI_SEND(vt(1,1,2),im*jm,type,mysize-1,1,MPI_COMM_HYDRO,ierr)
c
	  call MPI_RECV(tv(1,1,1),im*jm,type,mysize-1,2,MPI_COMM_HYDRO,status,ierr)
	  call MPI_SEND(tv(1,1,2),im*jm,type,mysize-1,3,MPI_COMM_HYDRO,ierr)

	elseif (myrank.eq.mysize-1) then 

	  call MPI_SEND(vt(1,1,kz2),im*jm,type,0,  0,MPI_COMM_HYDRO,ierr)
	  call MPI_RECV(vt(1,1,kz2+1),im*jm,type,0,1,MPI_COMM_HYDRO,status,ierr)
c
	  call MPI_SEND(tv(1,1,kz2),im*jm,type,0,  2,MPI_COMM_HYDRO,ierr)
	  call MPI_RECV(tv(1,1,kz2+1),im*jm,type,0,3,MPI_COMM_HYDRO,status,ierr)
	endif  
c       write(0,*) myrank,' after 1st call'
c
         do k=1,kz2+1
            do j=1,jy2+1
               vt(1,j,k)     = -vt(2,j,k)
               vt(ix2+1,j,k) = -vt(ix2,j,k)
            enddo
         enddo
         do k=1,kz2+1
            do i=1,ix2+1
               tv(i,1    ,k) = tv(i,jy2,k)
               vt(i,1    ,k) = vt(i,jy2,k)
               tv(i,jy2+1,k) = tv(i,2,k)
               vt(i,jy2+1,k) = vt(i,2,k)
            enddo
         enddo
  998  continue 
      else
         abash1 =-0.5
         abash2 = 1.5
      endif      
c
c     coefficients of the model
c     
      cw2=0.3
      cw3=2.
      ct3=1.2
      ct4=0.5
      cv1=7.1
      cb1=0.135
      cb2=0.622
c
      sigma = 2./3.
      oosig = 1./sigma
      ck    = 0.41
      cke2  = ck**2
      cw1   = cb1/cke2 + (1.+ cb2)/sigma
      cw3e6 = cw3**6
      cv1e3 = cv1**3
      cdes=0.65
c
c      calculate velocity field at cell centers
c
c      CALL center(uo,vo,wo,uoc,voc,woc,im,jm,km)
c
c     start i - j loop
c
c   debug
c       j=ny/2
c       i=nx/2
c        do k=2,kz2
c	 write(myrank+400,'(I, 2G20.10)') k,0.5*(vt(i,j,k)+vt(i,j,k-1)),0.5*(vt(i,j,k)+vt(i,j,k+1))
c	enddo
c	close(myrank+400) 
c   debug
      if (init) then
       init=.false. 
       do i=2,ix2
         oodxcp(i) = 1./(xc(i+1)-xc(i-1))
         oodxce(i) = 1./(xc(i+1)-xc(i))
         oodxp(i)  = 1./(xf(i)-xf(i-1))
       enddo
       oodxce(1) = 1./(xc(2)-xc(1))
       do j=2,jy2
         oodycp(j) = 1./(yc(j+1)-yc(j-1))
         oodycf(j) = 1./(yc(j+1)-yc(j))
         oodyp(j)  = 1./(yf(j)-yf(j-1))
       enddo
       oodycf(1) = 1./(yc(2)-yc(1))
       do k=2,kz2
         oodzcp(k)   = 1./(zc(k+1)-zc(k-1))
         oodzcn(k)   = 1./(zc(k+1)-zc(k))
         oodzp(k)    = 1./(zf(k)-zf(k-1))
       enddo
       oodzcn(1) = 1./(zc(2)-zc(1))
c       write(300+myrank,'(I,5G20.10)') (k,oodzcp(k),oodzcn(k),oodzp(k),zc(k),zf(k),k=1,kz2+1)
c       close(300+myrank)
      endif ! init
c 
      do k=2,kz2
         do j=2,jy2
      	    do i=2,ix2
c      	    
               dist1    = MIN(ABS(xc(i)-xf(1)),ABS(xf(ix2)-xc(i)))
               dist2    = cdes*ABS( dxdydz(i,j,k) )
               dist     = MIN(dist1,dist2)
c     shortcuts
               vte      = 0.5*(vt(i,j,k)+vt(i+1,j,k))
               vtw      = 0.5*(vt(i,j,k)+vt(i-1,j,k))
               vtf      = 0.5*(vt(i,j,k)+vt(i,j+1,k))
               vtb      = 0.5*(vt(i,j,k)+vt(i,j-1,k))
               vtn      = 0.5*(vt(i,j,k)+vt(i,j,k+1))
               vts      = 0.5*(vt(i,j,k)+vt(i,j,k-1))
c     functions and parameters ..........................................
               chi      = vt(i,j,k)/ru1
               fv1      = chi**3/(chi**3+cv1e3)
               fv2      = 1. - chi/(1.+fv1*chi)
               ft2      = ct3/EXP(ct4*chi*chi)
               stilde   = sabs(i,j,k) + vt(i,j,k)*( fv2/cke2/dist**2 )
               r        = MIN(vt(i,j,k)/stilde/cke2/dist**2,10.)
               g        = r +cw2*(r**6-r)
               fw       = g*( (1+cw3e6)/(g**6+cw3e6) )**0.167
c     production coefficient ............................................
               bprod    = cb1*stilde
c     destruction coefficient ...........................................
               bdest    = -cw1*fw*vt(i,j,k)/dist**2
c     viscous diffusion .................................................
c     explicit part .....................................................
               flyf     = (vtf+ru1)*(vt(i,j+1,k)-vt(i,j,k))*oodycf(j)
               flyb     = (vtb+ru1)*(vt(i,j,k)-vt(i,j-1,k))*oodycf(j-1)
               flzn     = (vtn+ru1)*(vt(i,j,k+1)-vt(i,j,k))*oodzcn(k)
               flzs     = (vts+ru1)*(vt(i,j,k)-vt(i,j,k-1))*oodzcn(k-1)
               vdiffexp = ( (flzn-flzs)*oodzp(k)
     $              + (flyf-flyb)*oodyp(j) )*oosig
c     implicit cooefficient ..............................................
               avdiff   = (vtw+ru1)*oodxce(i-1)*oodxp(i)*oosig
               cvdiff   = (vte+ru1)*oodxce(i)*oodxp(i)*oosig
               bvdiff   = -avdiff-cvdiff
c     diffusion ........................................................
               coef     = cb2*oosig
c     explicit part .....................................................
               diffexp  = coef*( ((vt(i,j,k+1)-vt(i,j,k-1))*oodzcp(k))**2 
     $                          +((vt(i,j+1,k)-vt(i,j-1,k))*oodycp(j))**2 )
c     implicit coefficients .............................................
               adiff    = -coef*(vt(i+1,j,k)-vt(i-1,j,k))*oodxcp(i)**2 
               cdiff    = -adiff
c     advection .........................................................
               advexp   =  oodxp(i)*( vte*uo(i,j,k) - vtw*uo(i-1,j,k) )
     $                    +oodyp(j)*( vtf*vo(i,j,k) - vtb*vo(i,j-1,k) )
     $                    +oodzp(k)*( vtn*wo(i,j,k) - vts*wo(i,j,k-1) ) 
c     put stuff together ................................................
               hn       = - advexp + vdiffexp  + diffexp
     $              +flag*( (avdiff+adiff)*vt(i-1,j,k)
     $                     +(bprod+bdest+bvdiff)*vt(i,j,k)
     $                     +(cvdiff+cdiff)*vt(i+1,j,k) )
               hnm1     = rhsnuo(i,j,k)
               anu(i)   =   -coflhs*(                adiff + avdiff)
               bnu(i)   = 1.-coflhs*(bprod + bdest +         bvdiff)
               cnu(i)   =   -coflhs*(                cdiff + cvdiff)
               fnu(i)   = - anu(i)*vt(i-1,j,k)
     $              + (2.-bnu(i))*vt(i,j,k)
     $              - cnu(i)*vt(i+1,j,k)
     $              + abash2*dt*hn 
     $              + abash1*dt*hnm1
c     update h^{n-1} ....................................................
               rhsnuo(i,j,k) = hn
            enddo
c     
c     wall b.c. ..........................................................
c     
c     0 at wall
      bnu(2)=bnu(2)-anu(2)
      bnu(ix2)=bnu(ix2)-cnu(ix2)

c            anu(1) = 0.
c            bnu(1) = 1.
c            cnu(1) = 1.
c            fnu(1) = 0.
c     
c            anu(ix2+1) = 1.
c            bnu(ix2+1) = 1.
c            cnu(ix2+1) = 0.
c            fnu(ix2+1) = 0.
c     
c     invert system ......................................................
c     
            call tridag(anu(2),bnu(2),cnu(2),fnu(2),xnu(2),ix2-1,mmx-1)
c
            do i=2,ix2
               vtnew(i,j,k) = MAX(0.,xnu(i))
               chi          = vtnew(i,j,k)/ru1
               fv1          = chi**3/(chi**3+cv1e3)
               tv(i,j,k)    = vtnew(i,j,k)*fv1
            enddo
c     
c     wall b.c. ..........................................................
c     
            vtnew(1,j,k)    =-vtnew(2,j,k)
            vtnew(ix2+1,j,k)=-vtnew(ix2,j,k)
            tv   (1,j,k)    =-tv   (2,j,k)
            tv   (ix2+1,j,k)=-tv   (ix2,j,k)
c
         enddo
      enddo
      
c     dump vtnew into vt

      do k=2,kz2
         do j=2,jy2
            do i=1,ix2+1
               vt(i,j,k) = vtnew(i,j,k)
            enddo
         enddo
      enddo

c     periodic bc 
996   continue
c        goto 997
        call REFRESHBC(tv,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO) 
        call REFRESHBC(vt,im*jm,(km-2),MYSIZE,MYRANK,MPI_COMM_HYDRO) ! vt
        if (myrank.eq.0) then 
	  call MPI_RECV(vt(1,1,1),im*jm,type,mysize-1,0,MPI_COMM_HYDRO,status,ierr)
	  call MPI_SEND(vt(1,1,2),im*jm,type,mysize-1,1,MPI_COMM_HYDRO,ierr)
	  call MPI_RECV(tv(1,1,1),im*jm,type,mysize-1,2,MPI_COMM_HYDRO,status,ierr)
	  call MPI_SEND(tv(1,1,2),im*jm,type,mysize-1,3,MPI_COMM_HYDRO,ierr)
	elseif (myrank.eq.mysize-1) then 
	  call MPI_SEND(vt(1,1,kz2),im*jm,type,  0,0,MPI_COMM_HYDRO,ierr)
	  call MPI_RECV(vt(1,1,kz2+1),im*jm,type,0,1,MPI_COMM_HYDRO,status,ierr)
	  call MPI_SEND(tv(1,1,kz2),im*jm,type,  0,2,MPI_COMM_HYDRO,ierr)
	  call MPI_RECV(tv(1,1,kz2+1),im*jm,type,0,3,MPI_COMM_HYDRO,status,ierr)
	endif  
  997 continue

      do k=1,kz2+1
         do i=1,ix2+1
            tv(i,1    ,k) = tv(i,jy2,k)
            vt(i,1    ,k) = vt(i,jy2,k)
            tv(i,jy2+1,k) = tv(i,2,k)
            vt(i,jy2+1,k) = vt(i,2,k)
         enddo
      enddo
c
c      write(0,*) myrank,' vtsaimp completed, iteration ',icycle
      return
 1000 format(' Imp ',i5,1p10e25.16)
 1956 format(f15.5,1p10e15.5)
      end
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDWALL(vt,im,jm,km,ibnd)

c       Sets the boundary values of nu_t for the case 
c       of a solid wall 
c
c           i.e.  nu_t(i,j,1)=-nu_t(i,j,2) 
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibx,iex,jby,jey,kbz,kez
      INTEGER isx,jsy,ksz,ippx,jpy,kpz
      INTEGER icx,jcy,kcz
      INTEGER i,j,k,ii,jj,kk
      INTEGER ibnd
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
****  Find orientation of boundary : 
*      i=const (isx=1,jsy=0,ksz=0)
*      j=const (isx=0,jsy=1,ksz=0)
*      k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez
*
****  Find possition of boundary
*     ippx=1 : begining, ippx=-1 : end, etc.
*     position indexes for the other surfaces are 0
*
      ippx=isign( 1,(ix2-ibx) ) * isx
      jpy=isign( 1,(jy2-jby) ) * jsy
      kpz=isign( 1,(kz2-kbz) ) * ksz
*
****  Indexes pointing the first inner values
*
      icx=( ((ix2-2)*ibx+ix2+2)/ix2 - ibx ) * isx
      jcy=( ((jy2-2)*jby+jy2+2)/jy2 - jby ) * jsy
      kcz=( ((kz2-2)*kbz+kz2+2)/kz2 - kbz ) * ksz
*
**** Main loop
*
      DO 10 i=ibx,iex
         ii=i+icx
       DO 20 j=jby,jey
          jj=j+jcy
        DO 30 k=kbz,kez
           kk=k+kcz
           vt(i,j,k) = -1.*vt(ii,jj,kk) 
30      CONTINUE
20     CONTINUE
10    CONTINUE

      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDGRAD(vt,im,jm,km,ibnd)

c       Sets the boundary values of nu_t for the case 
c       of a freestream boundaries 
c
c           i.e.  nu_t(i,j,1)= nu_t(i,j,2) 
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibx,iex,jby,jey,kbz,kez
      INTEGER isx,jsy,ksz,ippx,jpy,kpz
      INTEGER icx,jcy,kcz
      INTEGER i,j,k,ii,jj,kk
      INTEGER ibnd
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
****  Find orientation of boundary : 
*      i=const (isx=1,jsy=0,ksz=0)
*      j=const (isx=0,jsy=1,ksz=0)
*      k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez
*
****  Find possition of boundary
*     ippx=1 : begining, ippx=-1 : end, etc.
*     position indexes for the other surfaces are 0
*
      ippx=isign( 1,(ix2-ibx) ) * isx
      jpy=isign( 1,(jy2-jby) ) * jsy
      kpz=isign( 1,(kz2-kbz) ) * ksz
*
****  Indexes pointing the first inner values
*
      icx=( ((ix2-2)*ibx+ix2+2)/ix2 - ibx ) * isx
      jcy=( ((jy2-2)*jby+jy2+2)/jy2 - jby ) * jsy
      kcz=( ((kz2-2)*kbz+kz2+2)/kz2 - kbz ) * ksz
*
**** Main loop
*
      DO 10 i=ibx,iex
         ii=i+icx
       DO 20 j=jby,jey
          jj=j+jcy
        DO 30 k=kbz,kez
           kk=k+kcz
           vt(i,j,k) = vt(ii,jj,kk) 
30      CONTINUE
20     CONTINUE
10    CONTINUE

      RETURN
      END
c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDPERIOD(vt,im,jm,km,ibnd)

c      Sets boundary values for nu_t for the case of
c      periodic boundaries
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibx,iex,jby,jey,kbz,kez
      INTEGER isx,jsy,ksz,ippx,jpy,kpz
      INTEGER icx,jcy,kcz
      INTEGER i,j,k,ii,jj,kk
      INTEGER ibnd
*
****  Starting and ending indexes for ibnd 
*
      if (mysize.gt.1.and.(ibnd.eq.5.or.ibnd.eq.6)) then 
       call boundperiod_par(vt,im,jm,km,ibnd)
       return
      endif  

      ibx=ibegin(ibnd)
      iex=iend(ibnd)
      jby=jbegin(ibnd)
      jey=jend(ibnd)
      kbz=kbegin(ibnd)
      kez=kend(ibnd)
*
****  Find orientation of boundary : 
*      i=const (isx=1,jsy=0,ksz=0)
*      j=const (isx=0,jsy=1,ksz=0)
*      k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez     
*
****  Indexes pointing the coresponding inner values
*
      icx=isx*( (ix1+1-ix2)*(ibx-ix1)/ix2 + ix2 - ibx )
      jcy=jsy*( (jy1+1-jy2)*(jby-jy1)/jy2 + jy2 - jby )
      kcz=ksz*( (kz1+1-kz2)*(kbz-kz1)/kz2 + kz2 - kbz )
*
****  Shift the indexes
*
       ibx=ibx-1*(jsy+ksz)
       iex=iex+1*(jsy+ksz)
       jby=jby-1*(isx+ksz)
       jey=jey+1*(isx+ksz)
       kbz=kbz-1*(isx+jsy)
       kez=kez+1*(isx+jsy)

*
****  Main loop
*
      DO 10 i=ibx,iex
         ii=i+icx
       DO 20 j=jby,jey
          jj=j+jcy

        DO 30 k=kbz,kez
           kk=k+kcz

           vt(i,j,k) = vt(ii,jj,kk)

30      CONTINUE
20     CONTINUE
10    CONTINUE
c
      RETURN
      END
c
c...........................................
c
      subroutine boundperiod_par(vt,imax,jmax,kmax,ibnd)
c     this is the parallel version
c      V. Ovchinnikov 8/10/03. 
      implicit none  
      include 'headers/common.h'
      include 'mpif.h'
      integer imax,jmax,kmax,ibnd
      integer type,comm,ierr,status(MPI_STATUS_SIZE)
      real vt(imax,jmax,kmax)

      comm=MPI_COMM_HYDRO
c      
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
#ifndef DOUBLE 
      type=MPI_REAL
#endif
      if (ibnd.eq.5) then
        if (myrank.eq.0) then            
          call MPI_RECV(vt(1,1,1),imax*jmax,type,mysize-1,1,comm,status,ierr)
        else if (myrank.eq.mysize-1) then
          call MPI_SEND(vt(1,1,kz2),imax*jmax,type,0,1,comm,ierr)
	endif
      elseif (ibnd.eq.6) then 
        if (myrank.eq.mysize-1) then            
          call MPI_RECV(vt(1,1,kz2+1),imax*jmax,type,0,4,comm,status,ierr)
	else if (myrank.eq.0) then 
          call MPI_SEND(vt(1,1,2),imax*jmax,type,mysize-1,4,comm,ierr)
	endif
      endif
      return
      end

c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE BOUNDZERO(vt,im,jm,km,ibnd)

c       Sets the boundary values of vt to zero
c
c           i.e.  nu_t(i,j,1)=0.0
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    vt(im,jm,km)
c
      INTEGER ibx,iex,jby,jey,kbz,kez
      INTEGER isx,jsy,ksz,ippx,jpy,kpz
      INTEGER icx,jcy,kcz
      INTEGER i,j,k,ii,jj,kk
      INTEGER ibnd
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
****  Find orientation of boundary : 
*      i=const (isx=1,jsy=0,ksz=0)
*      j=const (isx=0,jsy=1,ksz=0)
*      k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez
*
****  Find possition of boundary
*     ippx=1 : begining, ippx=-1 : end, etc.
*     position indexes for the other surfaces are 0
*
      ippx=isign( 1,(ix2-ibx) ) * isx
      jpy=isign( 1,(jy2-jby) ) * jsy
      kpz=isign( 1,(kz2-kbz) ) * ksz
*
****  Indexes pointing the first inner values
*
      icx=( ((ix2-2)*ibx+ix2+2)/ix2 - ibx ) * isx
      jcy=( ((jy2-2)*jby+jy2+2)/jy2 - jby ) * jsy
      kcz=( ((kz2-2)*kbz+kz2+2)/kz2 - kbz ) * ksz
*
**** Main loop
*
      DO 10 i=ibx,iex
         ii=i+icx
       DO 20 j=jby,jey
          jj=j+jcy
        DO 30 k=kbz,kez
           kk=k+kcz
           vt(i,j,k) = 0.0
30      CONTINUE
20     CONTINUE
10    CONTINUE

      RETURN
      END

c-----SUBROUTINE-Strain---------------------- E. Balaras 7/12/98 ----
c
      SUBROUTINE Strain(uo,vo,wo,sxx,syy,szz,sxy,sxz,syz,g,im,jm,km)
c
c       Calculate strain rates at cell centers
c
c
c-----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c-----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    uo(im,jm,km),vo(im,jm,km),wo(im,jm,km),
     &        sxx(im,jm,km),syy(im,jm,km),szz(im,jm,km),
     &        sxy(im,jm,km),sxz(im,jm,km),syz(im,jm,km),
     &        g( im,jm,km)
c
      INTEGER i,j,k
      REAL    dx1,dy1,dz1
c
c-----------------------------------------------------------------------
c
      dx1=1./dx
      dy1=1./dy
      dz1=1./dz

c...calc du/dy -->sxx, dvdx -->syy
      do k=2,kz2
       do j=1,jy2
        do i=1,ix2
         sxx(i,j,k)=bv(j)*(uo(i,j+1,k)-uo(i,j,k))*dy1
         syy(i,j,k)=au(i)*(vo(i+1,j,k)-vo(i,j,k))*dx1
        enddo
       enddo
      enddo
c...Sxy at cell center
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         sxy(i,j,k)= 0.125*( sxx(i,j  ,k)+sxx(i-1,j  ,k)
     %                      +syy(i,j  ,k)+syy(i-1,j  ,k)
     %                      +sxx(i,j-1,k)+sxx(i-1,j-1,k)
     %                      +syy(i,j-1,k)+syy(i-1,j-1,k) )
        enddo
       enddo
      enddo
c...calc du/dz -->sxx, dwdx -->syy
      do k=1,kz2
       do j=2,jy2
        do i=1,ix2
         sxx(i,j,k)=cw(k)*(uo(i,j,k+1)-uo(i,j,k))*dz1
         syy(i,j,k)=au(i)*(wo(i+1,j,k)-wo(i,j,k))*dx1
        enddo
       enddo
      enddo
c...Sxz at cell center
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         sxz(i,j,k)= 0.125*( sxx(i,j  ,k)+sxx(i-1,j  ,k)
     %                      +syy(i,j  ,k)+syy(i-1,j  ,k)
     %                      +sxx(i,j,k-1)+sxx(i-1,j,k-1)
     %                      +syy(i,j,k-1)+syy(i-1,j,k-1) )         
        enddo
       enddo
      enddo    
c...calc dv/dz -->sxx, dwdy -->syy
      do k=1,kz2
       do j=1,jy2
        do i=2,ix2
         sxx(i,j,k)=cw(k)*(vo(i,j,k+1)-vo(i,j,k))*dz1
         syy(i,j,k)=bv(j)*(wo(i,j+1,k)-wo(i,j,k))*dy1
        enddo
       enddo
      enddo
c...Syz at cell center
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         syz(i,j,k)= 0.125*( sxx(i,j  ,k)+sxx(i,j-1,k  )
     %                      +syy(i,j  ,k)+syy(i,j-1,k  )
     %                      +sxx(i,j,k-1)+sxx(i,j-1,k-1)
     %                      +syy(i,j,k-1)+syy(i,j-1,k-1) )         
        enddo
       enddo
      enddo 
c...Sxx,Syy,Szz at cell center
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         sxx(i,j,k)=app(i)*(uo(i,j,k)-uo(i-1,j,k))*dx1
         syy(i,j,k)=bpp(j)*(vo(i,j,k)-vo(i,j-1,k))*dy1
         szz(i,j,k)=cpp(k)*(wo(i,j,k)-wo(i,j,k-1))*dz1
        enddo
       enddo
      enddo
c...|S| at cell center
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         g(i,j,k)=(  2.*sxx(i,j,k)*sxx(i,j,k)
     %              +2.*syy(i,j,k)*syy(i,j,k)
     %              +2.*szz(i,j,k)*szz(i,j,k)
     %              +4.*sxy(i,j,k)*sxy(i,j,k)
     %              +4.*sxz(i,j,k)*sxz(i,j,k)
     %              +4.*syz(i,j,k)*syz(i,j,k) )**0.5
        enddo
       enddo
      enddo
C
      RETURN
      END

c
c
c---------------------------------------------E. Balaras 7/12/98-------
c
      SUBROUTINE FILTERSIZE(xf,yf,zf,dxdydz,im,jm,km)

c       Calculates filter size on the non-uniform grid 
c----------------------------------------------------------------------
c
      implicit none
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      INTEGER im,jm,km
      REAL    xf(im),yf(jm),zf(km)
      REAL    dxdydz(im,jm,km)
      REAL    ppp,filt,filty,filtz,sq6
      INTEGER i,j,k,filttype
c
c     filter types:
c
c     0--geometric mean SQUARED -- use for dynamic models
c     1--geometric mean * sqrt(6)
c     2--cell diagonal
c     3--maximum cell side length      

      filttype=0

      if (Isgsmodel.eq.30) filttype=3
      if (Isgsmodel.eq.10) filttype=0
      if (Isgsmodel.eq.20) filttype=0
c
      if(filttype.eq.0) then 
       ppp=2./3.
       do i=2,ix2
        do j=2,jy2
         do k=2,kz2+1
          dxdydz(i,j,k)=( (xf(i)-xf(i-1))*
     %                    (yf(j)-yf(j-1))*
     %                    (zf(k)-zf(k-1)) )**ppp     
         enddo
        enddo
       enddo

      elseif(filttype .eq. 1) then 
         sq6=sqrt(6.)
         ppp=1./3.
         do i=2,ix2
            do j=2,jy2
               do k=2,kz2
                  dxdydz(i,j,k)=sq6* ( (xf(i)-xf(i-1))*
     $                 (yf(j)-yf(j-1))*
     $                 (zf(k)-zf(k-1)) )**ppp     
               enddo
            enddo
         enddo
      elseif(filttype .eq. 2) then 
         do i=2,ix2
            do j=2,jy2
               do k=2,kz2
                  dxdydz(i,j,k)=SQRT( (xf(i)-xf(i-1))**2
     $                 +(yf(j)-yf(j-1))**2
     $                 +(zf(k)-zf(k-1))**2 )
               enddo
            enddo
         enddo
      elseif(filttype .eq. 3) then 
         do i=2,ix2
            do j=2,jy2
               do k=2,kz2
                  filt = xf(i)-xf(i-1)
                  filty = yf(j)-yf(j-1)
                  filtz = zf(k)-zf(k-1)
                  if(filt .lt. filty) filt = filty
                  if(filt .lt. filtz) filt = filtz
                  dxdydz(i,j,k)=filt
               enddo
            enddo
         enddo
      endif    

      RETURN
      END


