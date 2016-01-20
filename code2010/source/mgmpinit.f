


      SUBROUTINE mgmpi_start(xu,yv,zw,nnx,nny,nnz)
      
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/mgmpi.h'
      
c **************************************
c The following subroutine initializes the 
c poisson solver based on mgmpi
c AS 12/9/01

c **************************************


      
c************************************************************
c**** The folowing parameters may need to be changed depending
c**** on the problem to be solved. In genral, use the smallest
c**** values possible
c*************************************************************

c     ---------------------------------------------------
c     array sizes
c     --------------------------------------------------

c     =====================================================
c     DECLARATIONS
c     =====================================================

      integer nnx,nny,nnz
      REAL xu(nnx),yv(nny),zw(nnz)


c============================================================
c     bigin code
c============================================================

c-------------------------------------------------------------
c     initialize MPI
c-------------------------------------------------------------

      ip=0
      print*,'in mginit'
c-----------------------------------------------------------
c     set initial parameters
c----------------------------------------------------------

      npx=1 ! # of processor
      npy=1
      npz=1
      ifine_solver=1
      maxit_fine=2000
      rtol_fine=1.0D-8
c      rtol_fine=1.0D-6
      icoarse_solver=1
      maxit_coarse=200
      rtol_coarse=1.0D-6
      m=1
      nu_pre=8
      nu_post=16
      iSmoothIt=20
      iSmoother=1
      wgt=0.88
      ic1=1
      ic2=1
      ic3=1
      is1=1
      is2=1
      is3=1
      iRestrict=1
      iProlong=1
      iCoarse=1
      iOutput=MM_NO_OUTPUT
      gtsDC=1.0
      iscollapsed=0
      iscaling=1



c-----------------------------------------------------
c     Determine processor indices (ipx,ipy,ipz) and ip
c-----------------------------------------------------

      ipx=0
      ipy=0
      ipz=0
      
c------------------------------------------------------
c     determine local <----> global index conversions
c------------------------------------------------------

c------------------------------------------------------
c     check that enough memory is allocated
c------------------------------------------------------
      if (nnx.gt.nxlim .or. nny.gt.nylim .or. nnz.gt.nzlim) then
       if (nnx.gt.nxlim) print*, ip, ':NXLIM must at least ',nnx
       if (nny.gt.nylim) print*, ip, ':NYLIM must at least ',nny
       if (nnz.gt.nzlim) print*, ip, ':NZLIM must at least ',nnz

        stop
      end if


c------------------------------------------------------
c     computr mesh coordinates
c------------------------------------------------------
 
      do ix = 2,nnx-1
         gtsxp(ix)=(DBLE(XU(ix))+DBLE(XU(ix-1)))/2.
         gtshx(ix)=gtsxp(ix)-gtsxp(ix-1)
      enddo
      do iy = 2,nny-1
         gtsyp(iy)=(DBLE(YV(iy))+DBLE(YV(iy-1)))/2.
         gtshy(iy)=gtsyp(iy)-gtsyp(iy-1)
      enddo
      do iz = 2,nnz-1
         gtszp(iz)=(DBLE(ZW(iz))+DBLE(ZW(iz-1)))/2.
         gtshz(iz)=gtszp(iz)-gtszp(iz-1)
      enddo

      gtsxp(1)=DBLE(XU(1))-(gtsxp(2)-DBLE(XU(1)))
      gtshx(1)=gtshx(2)
      gtsxp(nnx)=DBLE(XU(nnx-1))+(DBLE(XU(nnx-1))-gtsxp(nnx-1))
      gtshx(nnx)=gtshx(nnx-1)

      gtsyp(1)=DBLE(YV(1))-(gtsyp(2)-DBLE(YV(1)))
      gtshy(1)=gtshy(2)
      gtsyp(nny)=DBLE(YV(nny-1))+(DBLE(YV(nny-1))-gtsyp(nny-1))
      gtshy(nny)=gtshy(nny-1)

      gtszp(1)=DBLE(ZW(1))-(gtszp(2)-DBLE(ZW(1)))
      gtshz(1)=gtshz(2)
      gtszp(nnz)=DBLE(ZW(nnz-1))+(DBLE(ZW(nnz-1))-gtszp(nnz-1))
      gtshz(nnz)=gtshz(nnz-1)      

      xx1=gtsxp(1)
      yy1=gtsyp(1)
      zz1=gtszp(1)

c-------------------------------------------------
c     Neumann Z-faces
c-------------------------------------------------

         do iy= 1, nny

            do ix = 1, nnx

               gtsZMB(ix,iy) = 0.0D0
               gtsZPB(ix,iy) = 0.0D0

            enddo

         enddo


      RETURN
 
      END









