C 
C-----SUBROUTINE-PresDP-------------------------A. SCOTTI 1/27/02
C
      SUBROUTINE PresDP2P(DIV,DPP,PPO,XF,YF,ZF,DUM,
     &     myrank,mpi_comm_hydro,TIME)
c
      include 'headers/dimension.h'
      include 'mpif.h' 
c
      REAL DIV(NX,NY,*),DPP(NX,NY,*),DUM(NX,NY-2,*),
     &          PPO(NX,NY,*),
     &          XF(NX),YF(NY),ZF(NZG)
      real TIME(3)
      integer mpi_comm_hydro
      integer IM,JM,KM, NBI, NBJ, NBK, IIM, JJM, KKM, IY, JY, KY
      real dxp, dyp

**** Dimension of working array
*
c     parameter ( LW = (NX)*(MAX0(NZG,NY)-2) +
c     &     MAX0(2*NY-4,3*NX/2) +
c     &     2*(NX)              +
c     &     NY-2                    +
c     &     2*(NZG-2)+(NZG-1)+(NZG)  +
c     &     14+(NY-2)*NX*(NZG-2)/NCPU)

*************I AM NOT SURE ABOUT THE CALCULATION BELOW LW MAY BE TOO SMALL. VO 03 12 03
      integer LW, mysize, myrank, ierror
      integer ix1, jy1, kz1, ix2, jy2, kz2, i,j,k
      parameter ( LW = (NX)*(NZG-2) +
     &     3*NX/2 + 2*NY-4+
     &     2*(NX)              +
     &     NY-2                    +
     &     2*(NZG-2)+(NZG-1)+(NZG)  +
     &     14+(NY-2)*NX*(NZG-2)/NCPU)
      parameter (mysize=ncpu)
* 
**** Local arrays
*
      real WORK(lw),XLMB(nzg),ZH(nzg)
*
**** Transform array
*
c     
c      write(0,*) 'LW=',lw
      IM=NX
      JM=NY
      KM=NZG
c
      IX1=1
      JY1=1
      KZ1=1
      IX2=IM-1
      JY2=JM-1
      KZ2=KM-1
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
**** Dimension of working array
*
c      LW = (IIM+2)*MAX0(KKM,JJM) +
c     &     MAX0(2*JJM,3*IIM/2+1) +
c     &     2*(IIM+2)              +
c     &     JJM                    +
c     &     2*KKM+(KKM+1)+(KKM+2)  +
c     &     14+JY*IY*KKM/MYSIZE
c
c
c      allocate(work(lw))

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
     &            DXP,DYP,ZH,WORK,LW,IERROR,MYSIZE,MYRANK,
     &     MPI_COMM_HYDRO,
     &  TIME)


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
       CALL REFRESHBC(DPP,IY*(JY+2),(KM-2)/MYSIZE,MYSIZE,MYRANK,
     &      MPI_COMM_HYDRO)
c
       IF(MYRANK.EQ.0) THEN
           DO J=1,JY+2  
              DO I=1,IY
c 
                DPP(I,J,1)=DPP(I,J,2)
             ENDDO
          ENDDO
       ELSEIF(MYRANK.EQ.MYSIZE-1) THEN
          DO J=1,JY+2
             DO I=1,IY
c                
                DPP(I,J,(KM-2)/MYSIZE+2)=DPP(I,J,(KM-2)/MYSIZE+1)
             ENDDO
          ENDDO
       ENDIF
c          


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

c       deallocate(zh)
c       deallocate(xlmb)
c       deallocate(work)



      RETURN
      END
c

