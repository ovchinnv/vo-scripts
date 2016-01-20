
      include"mpif.h"
      real xu(nx),yv(ny),zw(nz)
      real p(nx,ny,(nz-2)/ncpu+2),dp(nx,ny,(nz-2)/ncpu+2)
     &  ,div(nx,ny,(nz-2)/ncpu+2)
     &  ,dum(nx,ny,(nz-2)/ncpu)
c      real xu(66),yv(66),zw(34)
c      real p(66,66,32/2+2),dp(66,66,32/2+2),div(66,66,32/2+2)
c      real dum(66,64,32/2)
      real time(3)

      open(8,file='input_size')

      read(8,*) nx,ny,nz

      close(8)

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,mysize,ierr)

      np=mysize

      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
      do i=1,nx
         xu(i)=(i-1)/float((nx-2))
      enddo
      do i=1,ny
         yv(i)=(i-1)/float((ny-2))
      enddo
      do i=1,nz
         zw(i)=(i-1)/float((nz-2))
      enddo

      idum=-1
      r=ran1(idum)

      do k=1,(nz-2)/np+2
         do j=1,ny
            do i=1,nx
               div(i,j,k)=ran1(idum)
               p(i,j,k)=ran1(idum)
            enddo
         enddo
      enddo


      

c      call h3dumm(P,PT,DZ,DZH,B,A,C,Z,
c     &     DX,DY,66,64,32,XLAM,ifirst,aln,alm,
c     &     ifax,trigs,66,64,1,ifirst,ierr,mysize,myrank)

      lw=4693+66*64*32/mysize
      nbi=0
      nbj=0
      nbk=1



c      call h3dcy2(P,66,64,64,64,32,nbi,nbj,nbk,xlam,dx,dy,zh,work,
c     &    lw,ierr,mysize,myrank,time)
      call presdp2p(div,dp,p,xu,yv,zw,dum,nx,ny,nz,mysize,myrank,TIME)
      cpub=mpi_wtime()
      t1=0
      t2=0
      t3=0
      do jjj=1,10
        call presdp2p(div,dp,p,xu,yv,zw,dum,nx,ny,nz,mysize,myrank,TIME)
        t1=t1+time(1)
        t2=t2+time(2)
        t3=t3+time(3)
      enddo
      
      cpue=(mpi_wtime()-cpub)/10
      t1=t1/10
      t2=t2/10
      t3=t3/10
      IF (MYRANK.EQ.0) THEN
         write(6,1000) nx,ny,nz,np
         PRINT*,' FFT TIME = ',t1
         XFLOPS=5*NX*NY*aLOG(float(NX*NY))*
     &        (NZ-2)/NP/T1/1.E6/alog(2.)
         PRINT*,' MFLOPS OF FFT (LOCAL)',XFLOPS
         PRINT*,'TRIDIAG TIME = ',T2
         XFLOPS=8*NX*NY*(NZ-2)/NP/T2/1.E6
         PRINT*,'MFLOPS OF TRIDIAG (LOCAL)',XFLOPS
         PRINT*,'MPI TIME = ',T3
         print*,'time TOT= ',cpue
         xflops=(5*alog(float(nx*ny))/alog(2.)+9)*nx*ny*nz
         print*,'CUMULATIVE MFLOPS = ',xflops/cpue/1.e6

         print*,'sample of solution ',dp(10,10,10)
      ENDIF
 1000 format('dimension ',i3,' x ',i3,' x ',i3, 
     &     ' and number of processors ',i2)
         call mpi_finalize(ierr)
      stop
      end

C 
C-----SUBROUTINE-PresDP-------------------------P. FLOHR--01/02/1994----
C
      SUBROUTINE PresDP2P(DIV,DPP,PPO,XF,YF,ZF,DUM,
     &     ,myrank,TIME)


      include "mpif.h"
      DIMENSION DIV(NX,NY,*),DPP(NX,NY,*),DUM(NX,NY-2,*),
     &          PPO(NX,NY,*),
     &          XF(NX),YF(NY),ZF(NZG)
      real TIME(3)

**** Dimension of working array
*
     parameter ( LW = (NX)*(MAX0(NZG,NY)-2) +
     &     MAX0(2*NY-4,3*NX/2) +
     &     2*(NX)              +
     &     NY-2                    +
     &     2*(NZG-2)+(NZG-1)+(NZG)  +
     &     14+(NY-2)*NX*(NZG-2)/NCPU)

*
**** Local arrays
*
      real WORK(lw),XLMB(nzg),ZH(nzg)
*
**** Transform array
*
C    print*,'remember to check beta in pressure definition!'
      IX1=1
      JY1=1
      KZ1=1
      IX2=NX-1
      JY2=NY-1
      KZ2=NZG-1
      do 7 k=1,(km-2)/mysize
      do 7 j=1,jm-2
      do 7 i=1,ix2+1


       dum(i,j,k)=0.0
7     continue

      DO 5 K=KZ1,(KM-2)/MYSIZE
      DO 5 J=JY1,JY2-1
      
#ifdef CRAY
         DO 5 I=IX1+1,IX2
       DUM(I,J,K)=DIV(I,J+1,K+1)
#endif
#ifdef FFTW
       do 5 i=ix1,ix2-1
          DUM(I,J,K)=DIV(I+1,J+1,K+1)
#endif
#ifdef SGI
       do 5 i=ix1,ix2-1
          DUM(I,J,K)=DIV(I+1,J+1,K+1)
#endif

5     CONTINUE
*
*
**** Compute indexes
*
      IY=IX2+1
      JY=JY2-1


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

c      if (LW.GT.100000) then
c         if(MYRANK.EQ.0) then
c            print*,'LW = ',LW,' IS LARGER THAN 100000, adjust value in 
c     &           presdp2p'
c         endif
c         call mpi_finalize(ierr)
c         stop
c      endif


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
     &            DXP,DYP,ZH,WORK,LW,IERROR,MYSIZE,MYRANK,TIME)
 

*
**** Transfer back
*
      DO  K=KZ1+1,(KM-2)/MYSIZE+1

         DO  J=JY1+1,JY2
            DO  I=IX1+1,IX2
#ifdef CRAY

               DPP(I,J,K)=DUM(I,J-1,K-1)
#endif
#ifdef FFTW
               DPP(I,J,K)=DUM(I-1,J-1,K-1)
#endif               
#ifdef SGI
               DPP(I,J,K)=DUM(I-1,J-1,K-1)
#endif               
            enddo
         enddo
      enddo
C     
**** B.C. 
*
C
C  CYCLIC B.C. IN X DIRECTION
C
      DO 175 K=1,(KM-2)/MYSIZE+2
       DO 185 J=1,JY2+1
        DPP(1,J,K)=DPP(IX2,J,K)
        DPP(IX2+1,J,K)=DPP(2,J,K)
185    CONTINUE
175   CONTINUE

C
C  CYCLIC B.C. IN Y DIRECTION
C
      DO 140 K=1,(KM-2)/MYSIZE+2      
         DO 130 I=1,IX2+1
      
            DPP(I,1,K)=DPP(I,JY2,K)
            DPP(I,JY2+1,K)=DPP(I,2,K)
130    CONTINUE
140    CONTINUE



C
C  NEUMAN B.C. IN Z DIRECTION
C
       CALL REFRESHBC(DPP,IY*(JY+2),(KM-2)/MYSIZE,MYSIZE,MYRANK)

       IF(MYRANK.EQ.0) THEN
           DO J=1,JY+2  
              DO I=1,IY
  
                DPP(I,J,1)=DPP(I,J,2)
             ENDDO
          ENDDO
       ELSEIF(MYRANK.EQ.MYSIZE-1) THEN
          DO J=1,JY+2
             DO I=1,IY
                
                DPP(I,J,(KM-2)/MYSIZE+2)=DPP(I,J,(KM-2)/MYSIZE+1)
             ENDDO
          ENDDO
       ENDIF
          


*
**** update pressure field
*     
       do k=kz1,(km-2)/mysize+2
          do j=jy1,jy2+1
             do i=ix1,ix2+1
                
                
                ppo(i,j,k)=ppo(i,j,k)+dpp(i,j,k)
             enddo
          enddo
       enddo
*     





      RETURN
      END
c

