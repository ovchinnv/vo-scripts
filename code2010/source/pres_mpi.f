


      SUBROUTINE MPI_SOLVER(DIV,DPP,DPO,NNX,NNY,NNZ)
      IMPLICIT NONE      
      include 'headers/common.h'

      include 'headers/dimension.h'
      include 'headers/mgmpi.h'
      
c ********************************************************
c Solves the poisson equation using a variety of 
c direct and multigrid methods
c
c     INPUT: DIV (Divergence)
c                  DPP(Guess from previous step)
c     OUTPUT: DPP (Pressure correction)
c     OUTPUT: DPO (PRessure corrected)
c
c     Requires: compile with libDmgNAMEOFCOMPILER.a
c     By AS and YSCHANG 12/10/01

      
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

C**************************************************************
      integer ix,iy,iz,nnx,nny,nnz
      real div(nnx,nny,nnz),dpp(nnx,nny,nnz),dpo(nnx,nny,nnz)
      real*8 resfine,resredfine,time,xmean
      integer numitfine,istat
      
c============================================================
c     begin code
c============================================================

      



c-------------------------------------------------------
c     initialize initial guess X and right-hand side F
c-------------------------------------------------------
      
c   ***   User Dependent   ***

      do iz = kz1+1,kz2

         do iy = jy1+1,jy2

            do ix = ix1+1,ix2
               gtsX(ix,iy,iz)=DBLE(dpp(ix,iy,iz))
               gtsF(ix,iy,iz)=DBLE(DIV(ix,iy,iz))
            enddo


         enddo

      enddo   


c **************************
c call the poisson solver
c *************************





         CALL mg3p(
     &        gtsX,gtsF,MM_POISSON,MM_CARTESIAN,
     &        MM_PERIODIC,MM_PERIODIC,MM_PERIODIC,MM_PERIODIC,
     &        MM_NEUMANN,MM_NEUMANN,
     &        gtsXMB,gtsXPB,gtsYMB,gtsYPB,gtsZMB,gtsZPB,
     &        NXLIM,NYLIM,NZLIM,
     &        MM_IS_UNIFORM,MM_IS_UNIFORM,MM_IS_NONUNIFORM,
     &        gtsxp,gtsyp,gtszp,
     &        npx,npy,npz,ipx,ipy,ipz,ip,
     &        nnx,nny,nnz,
     &        comm,nbr,
     &        maxit_fine,maxit_coarse,
     &        rtol_fine,rtol_coarse,
     &        m,wgt,nu_pre,nu_post,iSmoothit,
     &        ifine_solver,icoarse_solver,iSmoother,
     &        iOutput,
     &        numItFine,resFine,resRedFine,istat,time,
     &        gtsZ,NWLIM)



C-----------------------------------------------------
C     Update Pressure
C-----------------------------------------------------

         xmean=0.
      do iz = kz1+1,kz2
        do iy = jy1+1,jy2
           do ix = ix1+1,ix2
             DPP(ix,iy,iz)=SNGL(gtsX(ix,iy,iz))
             DPO(ix,iy,iz)=DPO(ix,iy,iz)+DPP(ix,iy,iz)
             xmean=xmean+dpp(ix,iy,iz)
           enddo
         enddo
      enddo

      do iz = kz1+1,kz2
        do iy = jy1+1,jy2
           do ix = ix1+1,ix2
             DPP(ix,iy,iz)=dpp(ix,iy,iz)-
     *             xmean/((kz2-kz1)*(jy2-jy1)*(ix2-ix1))
           enddo
         enddo
      enddo


C-------------------------------------------------------
C      Boundary Conditions   
C------------------------------------------------------

C
C     CYCLIC B.C. IN X DIRECTION
C
      DO iz=1,KZ2+1
       DO iy=1,JY2+1
        DPO(1,iy,iz)=DPO(IX2,iy,iz)
        DPO(IX2+1,iy,iz)=DPO(2,iy,iz)
        DPP(1,iy,iz)=DPP(IX2,iy,iz)
        DPP(IX2+1,iy,iz)=DPP(2,iy,iz)
       ENDDO
      ENDDO

C
C     CYCLIC B.C. IN Y DIRECTION
C
      DO ix=1,IX2+1
       DO iz=1,KZ2+1
        DPO(ix,1,iz)=DPO(ix,JY2,iz)
        DPO(ix,JY2+1,iz)=DPO(ix,2,iz)
        DPP(ix,1,iz)=DPP(ix,JY2,iz)
        DPP(ix,JY2+1,iz)=DPP(ix,2,iz)
       ENDDO
      ENDDO

C
C  NEUMAN B.C. IN Z DIRECTION
C
      DO ix=1,IX2+1
       DO iy=1,JY2+1
        DPO(ix,iy,1)=DPO(ix,iy,2)
        DPO(ix,iy,KZ2+1)=DPO(ix,iy,KZ2)
        DPP(ix,iy,1)=DPP(ix,iy,2)
        DPP(ix,iy,KZ2+1)=DPP(ix,iy,KZ2)
       ENDDO
      ENDDO



      RETURN
 
      END









