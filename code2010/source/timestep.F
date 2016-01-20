C
C-----------------------------------------------------------------------
C                 ***************************                         
C                 *         timestep.f      *                        
C                 ***************************                       
C----------------------------------------------------------------------- 
C
C       - TSTEPCEN: computes allowable timestep
C
C-----------------------------------------------------------------------
C
C 
C-----SUBROUTINE-TSTEPCEN-----------------------E. BALARAS--07/12/1998----
C
      SUBROUTINE TSTEPCEN( XF,YF,ZF,UUO,VVO,WWO,TV,IMAX,
     &                     JMAX,KMAX,DTM                )
C
C-----------------------------------------------------------------------
C
      include 'mpif.h'
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
C     ...Parameters...
      REAL EPS
      PARAMETER (EPS = 1.0E-6)
C     ...Intrinsic Functions...
      INTRINSIC MIN, REAL
C
      real XF(IMAX),YF(JMAX),ZF(KMAX),
     &          UUO(IMAX,JMAX,KMAX),VVO(IMAX,JMAX,KMAX),
     &          WWO(IMAX,JMAX,KMAX),TV(IMAX,JMAX,KMAX)
C
      INTEGER  imax, jmax, kmax,type
      REAL     cfl, cflv, dtm
      REAL     dtconv, dtvis, dtvisx, dtvisy, dtvisz, dtvijk, totvX, dtvis1, dtvis2
      REAL     dtconvout,dtvisout
      REAL     dti, dtj, dtk
      REAL     dxl, dx2, dyl, dy2, dzl, dz2
      REAL     rflag, uabs, vabs, wabs
      INTEGER  i, j, k, ier
      real tvmax,tvmin

#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

C 
C-----------------------------------------------------------------------
C
      cfl   = 0.35
      cflv  = 0.25
c
      dtconv  = 1000.
      dtvis   = 1000.
      tvmax=0.
      tvmin=1000.
c
      rflag = real(impfla)
C 
C
      do 10 k = kz1+1, kz2
        dzl = zf(k)-zf(k-1)
        dz2 = dzl*dzl
        do 11 j = jy1+1, jy2
          dyl = yf(j)-yf(j-1) 
          dy2 = dyl*dyl
          do 12 i = ix1+1, ix2
            dxl = xf(i)-xf(i-1) 
            dx2 = dxl*dxl
c
c           Viscous constraint 
c
            totvX = ru1 + tv(i,j,k)
c            write(6,*) ru1, tv(i,j,k)
c
            dtvisx = dx2/totvX
            dtvisy = dy2/totvX
            dtvisz = dz2/totvX
c
            dtvis1 = min(dtvisx,dtvisy,dtvisz)
            dtvis2 = min(dtvisy,dtvisz) ! NOTE: for flipped channel z replaced by x
            dtvijk = dtvis1*(1.0-rflag) + dtvis2*rflag
c
            dtvis = min(dtvis,dtvijk)
c
c           Convective constraint 
c
            uabs = eps + abs(UUO(i,j,k))
            vabs = eps + abs(VVO(i,j,k))
            wabs = eps + abs(WWO(i,j,k))
c
            dti  = dxl/uabs
            dtj  = dyl/vabs
            dtk  = dzl/wabs
c
            dtconv = min(dtconv,dti,dtj,dtk)
c
   12     continue
   11   continue
   10 continue
       do k=1,kmax
        do j=1,jmax
	 do i=1,imax
            tvmax=MAX(tv(i,j,k),tvmax)
            tvmin=MIN(tv(i,j,k),tvmin)
         enddo
	enddo
       enddo	 
c 
       tvmax=tvmax+ru1
       tvmin=tvmin+ru1    
c      
      dtvis   = cflv*dtvis
      dtconv  = cfl *dtconv
C
      dtm = min(dtvis,dtconv)

      call mpi_allreduce(dtvis,dtvisout,1,type,MPI_MIN,
     &     MPI_COMM_HYDRO,ier)
      call mpi_allreduce(dtconv,dtconvout,1,type,
     &     MPI_MIN,MPI_COMM_HYDRO,ier)
      
c      write(6,*) '1/Re=',ru1
c      if(myrank.eq.0) then
c         write(6,*) 'dtvis=',dtvisout,'dtconv=',dtconvout
c      endif
      
      dtm = min(dtvisout,dtconvout)

C
      RETURN 
      END 


C 
C-----SUBROUTINE-TSTEPCEN-----------------------E. BALARAS--07/12/1998----
C
      SUBROUTINE TSTEPCEN_ori( XF,YF,ZF,UUO,VVO,WWO,TV,IMAX,
     &                     JMAX,KMAX,DTM                )
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C 
C-----------------------------------------------------------------------
C
      real XF(IMAX),YF(JMAX),ZF(KMAX),
     &          UUO(IMAX,JMAX,KMAX),VVO(IMAX,JMAX,KMAX),
     &          WWO(IMAX,JMAX,KMAX),TV(IMAX,JMAX,KMAX)
C
      INTEGER  imax, jmax, kmax
      REAL     dtm, totvX, dtvisz, dtvisy, dtvisx, dtvis
C
      REAL     dxl, dyl, dzl, dt1, dt2, dt3, uabs, vabs, wabs
      INTEGER  i, j, k
C 
C-----------------------------------------------------------------------
C
      dt1=1000.
C 
C-----------------------------------------------------------------------
C                                                        DIFFUSIVE TERMS
C-----------------------------------------------------------------------
C
      do 10 i = IX1+1, IX2
        dxl = (XF(i+1)-XF(i-1)) * 0.5
        do 10 k = kz1+1, kz2
          dzl = (ZF(k+1)-ZF(k-1)) * 0.5
          do 10 j=jy1+1,jy2
            dyl = (YF(j+1)-YF(j-1)) * 0.5		
            totvX = ru1 + abs(TV(i,j,k))
C
            dtvisx = dxl*dxl/totvX
            dtvisy = dyl*dyl/totvX
            dtvisz = dzl*dzl/totvX
C
            dtvis = 0.1 * amin1(dtvisx,dtvisy,dtvisz)
            dt1   = amin1(dt1,dtvis)
 10   continue
      dtm = dt1
C 
C-----------------------------------------------------------------------
C                                                        ADVECTIVE TERMS
C-----------------------------------------------------------------------
C
      do 20 i = ix1+1, ix2
        do 20 j = jy1+1, jy2
          do 20 k = kz1+1, kz2
            dxl = (XF(i+1)-XF(i-1)) * 0.5
            dyl = (YF(j+1)-YF(j-1)) * 0.5
            dzl = (ZF(k+1)-ZF(k-1)) * 0.5
            uabs = 0.0001 + abs(UUO(i,j,k))
            vabs = 0.0001 + abs(VVO(i,j,k))
            wabs = 0.0001 + abs(WWO(i,j,k))
            dt2=0.35/((uabs/dxl + vabs/dyl + wabs/dzl))
            dt3 = amin1(dt1,dt2)
            DTM = amin1(DTM, 0.9*dt3)
 20   continue
C
      RETURN 
      END 
