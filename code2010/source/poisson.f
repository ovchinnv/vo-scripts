C
C-----------------------------------------------------------------------
C                 ***************************                         
C                 *        poisson.f        *                        
C                 ***************************                       
C----------------------------------------------------------------------- 
C
C	- Poisson:	calls the different Poisson solvers
C	- PresDP: 	direct Poisson solver (pressure correction)
C                       X dir periodic, Y dir periodic
C       - Multi3D       3D multigrid (any dir periodic or dp/dn=0)
C
C----------------------------------------------------------------------- 
C
C 
C-----SUBROUTINE-Poisson-----------------------E. Balaras--17/01/99----
C
      SUBROUTINE Poisson(divv,dp,p,xu,yv,zw,dummy,dtm1,
     &                   jxc,jyc,jzc,nlevel,
     &                   in_dx1,in_dx2,in_dx3,in_dx4,
     &                   in_sn1,in_sn2,in_sn3,in_sn4,
     &                   in_sp1,in_sp2,in_sp3,in_sp4,
     &                   in_st1,in_st2,in_st3,in_st4,
     &                   in_av1,in_av2,in_av3,in_av4,
     &                   in_in1,in_in2,in_in3,in_in4,
     &                   f,fi,rhsn,rhs2v,rhs2n,rhs3v,rhs3n,rhs4v,
     &                   pr2,pr3,pr4,
     &                   g11,g22,g33,g11_2,g22_2,g33_2,
     &                   g11_3,g22_3,g33_3,g11_4,g22_4,g33_4,
     &                   giac,
     &                   n1,n2,n3,n12,n22,n32,n13,n23,n33,
     &                   n14,n24,n34,icycle                        
     &                                                            )
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - different Poisson solvers call
C                                                                       
C-----------------------------------------------------------------------
C
c      IMPLICIT NONE
      include 'headers/dimension.h'
      include 'headers/common.h'
C
      REAL     dtar1(2),DTIME,deltatime
      EXTERNAL DTIME
      real divv(mmx,mmy,mmz),dp(mmx,mmy,mmz),p(mmx,mmy,mmz)
      real xu(mmx),yv(mmy),zw(mmz)


C 
C-----------------------------------------------------------------------
C
c      deltatime = DTIME(DTAR1)
c      WRITE(6,*)'*..dtime before Poisson: ',deltatime
C
c...Direct solver, 2 periodic directions (X,Y), streching only Z
      IF (Ipresmethod .EQ. 1) THEN
        CALL PresDP2P(divv,dp,p,xu,yv,zw,dummy,nx,ny,nz)
c...Direct solver, 1 periodic direction (Y), streching Z
      ELSE IF(Ipresmethod .EQ. 2) THEN
        CALL PresDP1P(divv,dp,p,xu,yv,zw,dummy,nx,ny,nz)
c...Multigrid solver, streching in all directions
      ELSE IF(Ipresmethod .EQ. 10) THEN
        CALL multig3d(jxc,jyc,jzc,
     &                in_dx1,in_dx2,in_dx3,in_dx4,
     &                in_sn1,in_sn2,in_sn3,in_sn4,
     &                in_sp1,in_sp2,in_sp3,in_sp4,
     &                in_st1,in_st2,in_st3,in_st4,
     &                in_av1,in_av2,in_av3,in_av4,
     &                in_in1,in_in2,in_in3,in_in4,
     &                f,fi,rhsn,rhs2v,rhs2n,rhs3v,rhs3n,rhs4v,
     &                pr2,pr3,pr4,
     &                g11,g22,g33,g11_2,g22_2,g33_2,
     &                g11_3,g22_3,g33_3,g11_4,g22_4,g33_4,
     &                giac,
     &                divv,dp,p,xu,yv,zw,
     &                nx,ny,nz,dtm1,nlevel,                      
     &                n1,n2,n3,n12,n22,n32,n13,n23,n33,
     &                n14,n24,n34,icycle                        )
c...Direct Solver, 1 periodic direction (Y), streching in X and Z 
      ELSE IF(Ipresmethod .EQ. 22) THEN
        CALL PresDirect(divv,dp,p,xu,yv,zw,au,cw,icycle,nx,ny,nz)

c...MPI Solver, 2 periodic direction (X,Y), stretching in z.
c...modified by Y.S. Chang.  Nov.14, 2001
      ELSE IF(Ipresmethod .eq. 33) THEN
        CALL MPI_SOLVER(divv,dp,p,nx,ny,nz)

      ELSE
c
        STOP 'Invalid solver assignement' 
c
      END IF
c
      RETURN
      END
C
C 
C-----SUBROUTINE-PresDP-------------------------P. FLOHR--01/02/1994----
C
      SUBROUTINE PresDP2P(DIV,DPP,PPO,XF,YF,ZF,DUM,IM,JM,KM  )

      include 'headers/common.h'

      DIMENSION DIV(IM,JM,KM),DPP(IM,JM,KM),DUM(IM,JM,KM),
     &          PPO(IM,JM,KM),
     &          XF(IM),YF(JM),ZF(KM)
*
**** Local arrays
*
      DIMENSION WORK(30000),XLMB(mmx),ZH(mmz)
*
**** Transform array
*
      do 7 i=1,ix2+1
      do 7 j=1,jy2
      do 7 k=1,kz2
       dum(i,j,k)=0.0
7     continue

      DO 5 I=IX1+1,IX2
      DO 5 J=JY1,JY2-1
      DO 5 K=KZ1,KZ2-1
       DUM(I,J,K)=DIV(I,J+1,K+1)
5     CONTINUE
*
*
**** Compute indexes
*
      IY=IX2+1
      JY=JY2+1
      KY=KZ2+1

      IIM=IX2-1
      JJM=JY2-1
      KKM=KZ2-1
*
**** Mesh spacing
*
      DXP=XF(IX2)/FLOAT(IX2-1)
      DYP=YF(JY2)/FLOAT(JY2-1)
      DO 10 K=KZ1,KZ2
       ZH(K)=ZF(K)
10    CONTINUE
*
**** Dimension of working array
*
      LW = (IIM+2)*JMAX0(KKM,JJM) +
     &     JMAX0(2*JJM,3*IIM/2+1) +
     &     2*(IIM+2)              +
     &     JJM                    +
     &     2*KKM+(KKM+1)+(KKM+2)  +
     &     14
c      write(6,*)'dim. of work array', lw
*
**** Initialize past arrays
*
      DO 20 J=1,KKM
       XLMB(J)=0.0
20    CONTINUE

      DO 30 I=1,LW
       WORK(I)=0.
30    CONTINUE
*
**** BC: NEUMANN=1,CYCLIC=0
*
c      NBI=1
c      NBJ=1
c      NBK=1

      NBI=0
      NBJ=0
      NBK=1

*
**** Call direct solver
*
C

      CALL H3DCY2(DUM,IY,JY,IIM,JJM,KKM,NBI,NBJ,NBK,XLMB,
     &            DXP,DYP,ZH,WORK,LW,IERROR)

C      WRITE(6,*) 'IERROR= ',IERROR,'LW= ',LW
*
**** Transfer back
*
      DO 90 I=IX1+1,IX2
       DO 80 J=JY1+1,JY2
        DO 70 K=KZ1+1,KZ2
         DPP(I,J,K)=DUM(I,J-1,K-1)
70      CONTINUE
80     CONTINUE      
90    CONTINUE
C

*
**** update pressure field
*
      do i=ix1+1,ix2
       do j=jy1+1,jy2
        do k=kz1+1,kz2
         ppo(i,j,k)=ppo(i,j,k)+dpp(i,j,k)
        enddo
       enddo
      enddo
*
**** B.C. 
*
C
C  CYCLIC B.C. IN X DIRECTION
C
      DO 175 K=1,KZ2+1
       DO 185 J=1,JY2+1
        PPO(1,J,K)=PPO(IX2,J,K)
        PPO(IX2+1,J,K)=PPO(2,J,K)
        DPP(1,J,K)=DPP(IX2,J,K)
        DPP(IX2+1,J,K)=DPP(2,J,K)
185    CONTINUE
175   CONTINUE

C
C  CYCLIC B.C. IN Y DIRECTION
C
      DO 130 I=1,IX2+1
       DO 140 K=1,KZ2+1
        PPO(I,1,K)=PPO(I,JY2,K)
        PPO(I,JY2+1,K)=PPO(I,2,K)
        DPP(I,1,K)=DPP(I,JY2,K)
        DPP(I,JY2+1,K)=DPP(I,2,K)
140    CONTINUE
130    CONTINUE
C
C  NEUMAN B.C. IN Z DIRECTION
C
      DO 150 I=1,IX2+1
       DO 160 J=1,JY2+1
        PPO(I,J,1)=PPO(I,J,2)
        PPO(I,J,KZ2+1)=PPO(I,J,KZ2)
        DPP(I,J,1)=DPP(I,J,2)
        DPP(I,J,KZ2+1)=DPP(I,J,KZ2)
160    CONTINUE
150   CONTINUE

      RETURN
      END

c
c
c---------------------------------------------------------------------
c
      SUBROUTINE PresDP1P(div,dp,po,xf,yf,zf,dum,imax,jmax,kmax)
c
c      solves poisson equation on a staggered grid
c      j is the periodic direction
c      i is NS/NS and uniform grid
c      k is the stretched grid direction
c

c---------------------------------------------------------------------
      include 'headers/common.h'
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      REAL  div(imax,jmax,kmax),dp(imax,jmax,kmax),
     &      po(imax,jmax,kmax),dum(imax,jmax,kmax)
      REAL  xf(imax),yf(jmax),zf(kmax)
c
      REAL  work(100000),xlmb(mmz),zh(mmz)
c      REAL  work(mmx*mmy*mmz),xlmb(mmz),zh(mmz)
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c                                                    procces variables
c---------------------------------------------------------------------
c..dimension of the arrays
      IY=IX2+1
      JY=JY2+1
      KY=KZ2+1

c..number of unknows
      IM=IX2-1
      JM=JY2-1
      KM=KZ2-1

c..bc set up (NEUMANN=1,CYCLIC=0)
      NBI=1        ! poissx
      NBJ=0
      NBK=1

c      NBI=5         ! poistg  
c      NBJ=0
c      NBK=1

c..grid spacing
      DXP=XF(IX2)/FLOAT(IX2-1)
      DYP=YF(JY2)/FLOAT(JY2-1)
      do K=KZ1,KZ2
       ZH(K)=ZF(K)
      enddo

c..dimension of working array
c      LW=50000
c      LW = (IM+2)*JMAX0(KM,JM) +
c     &     JMAX0(2*JM,3*IM/2+1) +
c     &     2*(IM+2)              +
c     &     JM                    +
c     &     2*KM+(KM+1)+(KM+2)  +
c     &     14


      LWFFT = 3*JM/2+1 + IM*(JM+1)
      LWPOI = 4*IM + (IM+INT(LOG(float(IM)))+10)*KM  !FOR NBI=1,..,4 
c      LWPOI = INT(4*IM + (IM+LOG(float(IM))+12)*KM)  !FOR NBI=5,..,8
      LW = 2*KM+4 + MAX0(LWFFT,LWPOI)          
C                                                                       
C                                   
C      write(6,*) 'lw=',lw

c---------------------------------------------------------------------
c                                              Initialize past arrays
c---------------------------------------------------------------------
      DO K=1,MMZ
       XLMB(K)=0.0
      ENDDO
c
      DO I=1,100000 
c      DO I=1,mmx*mmy*mmz
       WORK(I)=0.
      ENDDO

c---------------------------------------------------------------------
c                                                     Transform array
c---------------------------------------------------------------------
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         DUM(I,J,K)=DIV(I,J,K)
        enddo
       enddo
      enddo

c---------------------------------------------------------------------
c                                                   Call direct solver
c---------------------------------------------------------------------
      CALL H3DNDC(DUM,IY,JY,IM,JM,KM,NBI,NBJ,NBK,XLMB,
     &            DXP,DYP,ZH,WORK,LW,IERROR)  

c         write(6,*) 'LW=',LW
c         write(6,*) 'IERROR=',IERROR

c---------------------------------------------------------------------
c                                                       Transfer back
c---------------------------------------------------------------------
      do k=2,kz2
       do j=2,jy2
        do i=2,ix2
         dp(I,J,K)=dum(I,J,K)
        enddo
       enddo
      enddo

c---------------------------------------------------------------------
c                                                update pressure field      
c---------------------------------------------------------------------
      do i=ix1+1,ix2
       do j=jy1+1,jy2
        do k=kz1+1,kz2
         po(i,j,k)=po(i,j,k)+dp(i,j,k)
        enddo
       enddo
      enddo
c---------------------------------------------------------------------
c                                                          pressure bc     
c---------------------------------------------------------------------
c...NEUMANN B.C. IN X DIRECTION
      DO K=1,KZ2+1
       DO J=1,JY2+1
        PO(1,J,K)    =PO(2,J,K)
        PO(IX2+1,J,K)=PO(IX2,J,K)
        DP(1,J,K)    =DP(2,J,K)
        DP(IX2+1,J,K)=DP(IX2,J,K)
       ENDDO
      ENDDO

c...CYCLIC B.C. IN Y DIRECTION
      DO I=1,IX2+1
       DO K=1,KZ2+1
        PO(I,1,K)    =PO(I,JY2,K)
        PO(I,JY2+1,K)=PO(I,2,K)
        DP(I,1,K)    =DP(I,JY2,K)
        DP(I,JY2+1,K)=DP(I,2,K)
       ENDDO
      ENDDO       

c...NEUMAN B.C. IN Z DIRECTION
      DO I=1,IX2+1
       DO J=1,JY2+1
        PO(I,J,1)    =PO(I,J,2)
        PO(I,J,KZ2+1)=PO(I,J,KZ2)
        DP(I,J,1)    =DP(I,J,2)
        DP(I,J,KZ2+1)=DP(I,J,KZ2)
       ENDDO
      ENDDO  
c
      RETURN
      END
