
      SUBROUTINE H3DCY2(Y,IY,JY
     *                 ,IM,JM,KM
     *                 ,NBI,NBJ,NBK,XLMB
     *                 ,DX,DY,ZH
     *                 ,WORK,LW,IERROR,MYSIZE,MYRANK,MPI_COMM_HYDRO,
     *     TIME)
c                                                                       
C                                                                       
C     * * * * * * * * *  PURPOSE    * * * * * * * * * * * * * * * * * * 
C                                                                       
C     THIS ROUTINE CALCULATES THE SOLUTION OF A HELMHOLTZ EQUATION      
C     IN A CARTESIAN COORDINATE-SYSTEM WITH PRESCRIBED BOUNDARY
C     CONDITIONS (CYCLIC IN I - AND J - DIRECTION, NEUMANN OR
C     DIRICHLET IN K - DIRECTION) FOR A THREE DIMENSIONAL DOMAIN,
C     DISCRETIZED ON A STAGGERED GRID. THE PROBLEM IS SOLVED IN PARALLEL
C     USING MYSIZE NUMBER OF PROCESSORS.
C                                                                      
C     THE 3D-PROBLEM IS REDUCED TO 'IMP2*JM' 1D-PROBLEMS VIA A FAST-
C     FOURIER-TRANSFORMATION (FFT). THE 1D-PROBLEMS ARE SOLVED BY USE
C     OF A ROUTINE 'H3DUMM', CONTAINING A TRIDIAGONAL ALGORITHM, WHICH
C     IS OPTIMIZED FOR CRAY COMPUTERS BUT ALSO APPLICABLE TO OTHER
C     MACHINES. THE 2D-FFT ARE PERFORMED BY THE ROUTINES FFTFAX,FFT99,
C     CFTFAX AND CFFT99 (BY C. TEMPERTON), WHICH ARE OPTIMIZED
C     FOR CRAY-COMPUTERS VIA SPECIAL ROUTINES IN 'CRAY ASSEMBLER
C     LANGUAGE (CAL)'. FOR THE ORDINARY FORTRAN VERSION SEE FILE
C     'PA54.JOB.CNTL(FFT99)'.
C                                                                       
C                                                                       
C                                                                       
C     SPECIFICALLY SUBROUTINE 'H3DCY2' SOLVES THE LINEAR SYSTEM OF      
C     EQUATIONS                                                         
C                                                                      
C       (X(I-1,J,K) - 2.*X(I,J,K)   + X(I+1,J,K)) / (DX*DX) +         
C       (X(I,J-1,K) - 2.*X(I,J,K)   + X(I,J+1,K)) / (DY*DY) +           
C       A(K)*X(I,J,K-1) + B(K)*X(I,J,K) + C(K)*X(I,J,K+1)   +           
C                       - XLMB(K)*X(I,J,K)                = Y(I,J,K)    
C                                                                      
C       FOR I=2,..,IM+1;  J=1,..,JM  AND K=1,..,KM                      
C                                                                       
C     WITH SPECIFIC BOUNDARY CONDITIONS FOR EACH CO-ORDINATE. THE       
C     K-DEPENDENT FACTORS 'A', 'B' AND 'C' ALLOW FOR NON-EQUIDISTANT    
C     GRID SPACING IN THE K-DIRECTION; THEY NEED NOT BE PREDEFINED,     
C     BUT ARE CALCULATED FROM THE "MESH-BOUNDARY POSITIONS" STORED      
C     IN 'ZH' (SEE DESCRIPTION OF PARAMETERS FOR DETAILS). 
C
C     ************ DOMAIN DECOMPOSITION   ******************
C     THE DATA IS DISTRIBUTED AMONG THE PROCESSOR USING A SLAB
C     DECOMPOSITION. SPECIFICALLY, THE PROCESSOR MYRANK OWNS
C     THE FOLLOWING SLICE OF DATA
C     
C     I=1,...,IM+2; J=1,...,JM AND K=MYRANK*KM/MYSIZE+1,...(MYRANK+1)*KM/MYSIZE
C
C     NOTICE THAT THERE IS NO OVERLAP. ALSO, THERE ARE NO GHOSTLAYERS
C     IN THE Z DIRECTION. 
C     INFORMATION ON THE GRID SPACING IS NOT DISTRIBUTED, THAT IS EACH 
C     PROCESSOR OWNS ZH(1,...,KM+1). 
C
C     ************* ALGORITHM DETAILS *********************
C     THE DATA IN Y ARE FFTED LOCALLY ALONG THE X AND Y DIRECTION. 
C     THEN, Y IS SWAPPED AMONG PROCESSORS USING SWAP_XZ, SO THAT NOW
C     EACH PROCESSOR OWNS A VERTICAL SLICE
C     
C     I=1,...,IM+2; J=MYRANK*JM/MYSIZE+1,...,(MYRANK+1)*JM/MYSIZE AND K=1,...,KM
C     
C     THE TRIDIAGONAL SYSTEM IS SOLVED AND FINALLY THE DATA IS SWAPPED BACK
C     INTO THE ORIGINAL LAYOUT. IT IS IMPORTANT THAT BOTH JM AND KM BE EXACTLY
C     MULTIPLES OF MYSIZE. ALL COMMUNICATION IS HANDLED IN SWAP_XX ROUTINES.
C                                                                       
C                                                                       
C                                                                       
C     * * * * * * * *    PARAMETER DESCRIPTION     * * * * * * * * * *  
C                                                                       
C             * * * * * *   ON INPUT    * * * * * *                     
C                                                                       
C   Y                                                                   
C     THREE-DIMENSIONAL ARRAY WHICH CONTAINS THE RIGHT HAND SIDE OF     
C     THE LINEAR SYSTEM ABOVE AT LOCATIONS I=2,..,IM+1; J=1,..,JM       
C     AND K=1,..,KM. POSITIONS Y(IM+2,*,*) AND Y(1,*,*) MUST NOT BE     
C     SPECIFIED IN THE CALLING PROGRAM, BUT ARE COMPUTED BY SUBROU-     
C     TINE CYCLIC SO, THAT Y(1   ,*,*) = Y(IM+1,*,*) AND
C                          Y(IM+2,*,*) = Y(1   ,*,*)
C     'Y' MUST BE DIMENSIONED AS Y(IY,JY,KY).
C                                                                       
C   IY                                                                 
C     FIRST DIMENSION OF 'Y' EXACTLY AS IN THE CALLING PROGRAMME;     
C     'IY' MUST AT LEAST EQUAL IM+2.                                 
C                                                                   
C   JY                                                                 
C     SECOND DIMENSION OF 'Y' EXACTLY AS IN THE CALLING PROGRAMME;    
C     'JY' MUST BE EQUAL JM. JY/MYSIZE HAS TO BE AN INTEGER.                                   
C                                                                   
C   IM                                                                 
C     THE NUMBER OF UNKNOWNS IN THE I-DIRECTION; IM MUST BE EVEN AND  
C     MAY ONLY CONTAIN PRIME FACTORS 2, 3 OR 5; IM >= 4
C                                                                      
C   JM                                                                
C     THE NUMBER OF UNKNOWNS IN THE J-DIRECTION; JM MUST BE EVEN AND 
C     MAY ONLY CONTAIN PRIME FACTORS 2, 3 OR 5; IM >= 4
C                                                                      
C   KM                                                                 
C     THE NUMBER OF UNKNOWNS IN THE K-DIRECTION; KM >= 3. KM/MYSIZE HAS TO
C     BE AN INTEGER.
C
C
C   NBI
C     BECAUSE OF THE CYCLIC BOUNDARY CONDITION IN I-DIRECTION NBI
C     IS INTRODUCED FOR REASONS OF CONSISTANCY WITH H3DNDC; FURTHER-
C     MORE NBI IS USED AS AN INFORMATIVE PARAMETER:
C     NBI=0   INVOKES A SPECIAL OUTPUT OPTION FOR CPU-TIMES FOR THE  
C             FFT-PART AND THE GAUSS-ALGORITHM, RESPECTIVELY.        
C                                                                    
C   NBJ                                                              
C     SPECIFICATION OF BOUNDARY CONDITIONS IN J-DIRECTION.          
C     BECAUSE OF THE CYCLIC BOUNDARY CONDITION IN THE J -          
C     DIRECTION NBJ IS ONLY INTRODUCED FOR REASONS OF
C     CONSISTANCY WITH H3DNDC
C
C                                                                    
C   NBK                                                             
C     SPECIFICATION OF BOUNDARY CONDITIONS IN K-DIRECTION          
C                                                                 
C     NBK = 1 ==> NEUMANN   - NEUMANN   - BOUNDARY-CONDITION
C     NBK = 2 ==> DIRICHLET - DIRICHLET - BOUNDARY-CONDITION
C     NBK = 3 ==> NEUMANN   - DIRICHLET - BOUNDARY-CONDITION
C     NBK = 4 ==> DIRICHLET - NEUMANN   - BOUNDARY-CONDITION
C                                                                      
C                                                                      
C   XLMB                                                               
C     K-DEPENDENT HELMHOLTZ PARAMETER; MUST BE DIMENSIONED AS XLMB(KM) 
C                                                                      
C   DX                                                                
C     CONSTANT MESH-SIZE IN I-DIRECTION                               
C                                                                     
C   DY                                                                
C     CONSTANT MESH-SIZE IN J-DIRECTION                               
C                                                                     
C   ZH                                                                
C     NOT NECESSARILY EQUIDISTANT POSITIONS OF MESH-BOUNDARIES IN     
C     THE K-DIRECTION; IF THIS IS THE VERTICAL, 'ZH(1)' STANDS FOR    
C     DOMAIN'S GROUND LEVEL AND 'ZH(KM+1)' FOR ITS TOP LEVEL.         
C     'ZH' MUST BE DIMENSIONED AS ZH(KM+1).                            
C                                                                     
C     NOTE: OTHER GEOMETRIC QUANTITIES DERIVED FROM 'ZH' ARE COM-     
C           COMPUTED IN SUBROUTINE 'HOEFEL'; THEY ARE USED FOR        
C           DETERMINATION OF THE FACTORS 'A', 'B' AND 'C' IN THE      
C           LINEAR SYSTEM ABOVE.                                      
C                                                                     
C                                                                     
C   WORK                                                              
C     WORK ARRAY, AT LEAST OF LENGTH 'LW'                             
C                                                                     
C                                                                     
C   LW                                                                
C     LENGTH OF 'WORK'; LW = S1 + S2 + S3 + S4 + S5 + S6 +S7              
C                                                                     
C        S1 = (IM+2)*MAX(KM,JM)                                       
C        S2 = (MAX(2*JM,3*IM/2+1)
C        S3 = 2*(IM+2)
C        S4 = JM
C        S5 = 2*KM + (KM+1) + (KM+2)
C        S6 = IY*JY*KM/MYSIZE
C        S7 = 13

C     MYSIZE
C        THE TOTAL NUMBER OF PROCESSES IN THE COMMUNICATOR, AS REPORTED
C        BY A CALL TO MPI_COMM_SIZE(COMM,MYSIZE,IERR). 
C
C     MYRANK
C        THE RANK OF THE PROCESS IN THE COMMUNICATOR, AS REPORTED BY A
C        CALL TO MPI_COMM_RANK(COMM,MYRANK,IERR). 
C                                                                     
C                                                                     
C                                                                    
C             * * * * * *   ON OUTPUT     * * * * * *                
C                                                                    
C   Y                                                                
C     CONTAINS THE SOLUTION X ON THE INPUT SLAB
C     WITH X(1,*,*)    = X(IM+1,*,*) AND  X(IM+2,*,*) = X(2,*,*)                      
C                                                                    
C   IERROR                                                           
C     AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.  EXCEPT 
C     FOR NUMBER ZERO, A SOLUTION IS NOT ATTEMPTED.                  
C                                                                    
C     IERROR = 0:  NO ERROR                                          
C     IERROR = 1:  IY<IM+2  OR  JY<JM                                
C     IERROR = 2:  IM<4     OR  JM<4 OR KM<3                         
C     IERROR = 3:  IM OR JM CONTAIN FACTORS OTHER THAN 2,3,5         
C     IERROR = 4:  LW TOO SMALL                                      
C     IERROR = 5:  NBK OUT OF RANGE (1,..,4)
C     IERROR = 6:  JY OR KM NOT DIVISIBLE BY MYSIZE
C                                                                    
C     TIME
C       CONTAINS THE TIME SPENT IN DIFFERENT PARTS OF 
C       THE SUBROUTINE BY THE CALLING PROCESS.
C       THE TOTAL TIME IS OBTAINED SUMMING THE TIMES REPORTED BY ALL 
C       THE PROCESSES. TIME HAS DIMENSION 3. TIME(1) IS THE TIME SPENT IN FFT
C       TIME(2) FOR THE TRIDIAGONAL PART AND TIME(3) FOR MESSAGE PASSING
C     
C                                                                    
C                                                                    
C                                                                    
C     * * * * * * *   PROGRAM SPECIFICATIONS    * * * * * * * * * * * *
C                                                                      
C     DIMENSION OF   Y(IY,JY,KM), XLMB(KM), ZH(KM+1)                   
C     ARGUMENTS      WORK(LW)  (SEE ARGUMENT LIST FOR DETAILS)         
C                                                                      
C     LATEST         JUNE 1984  (SERIAL CODE)  BY UNKNOWN                                      
C     REVISION     JANUARY 2002 (PARALLEL VERSION) BY A. SCOTTI. SEND
C     QUESTIONS, BUG REPORTS AND COMMENTS TO ASCOTTI@UNC.EDU    
C
C     SUBPROGRAMS  - 'FFT99','CFFT99' AND LOWER LEVEL
C     REQUIRED       SUBROUTINES FOR FAST FOURIER TRANSFORMATION:
C                    FFT99
C                    CFFT99
C                    FFTFAX
C                    CFTFAX
C                    ALL THESE SUBROUTINES ARE CONTAINED IN 'FFT99.F' .
C     
C     COMMUNICATION BETWEEN PROCESSORS IS HANDLED BY CALLS TO MPI
C     ROUTINES (MPI_ALLTOALLV), CONTAINED IN SWAP.F 
C     TIMING IS CALCULATED USING MPI_WTIME(), A PORTABLE TIMER 
C     SPECIFIED BY THE MPI STANDARD.
C
C     SPECIAL           REQUIRES 'MPIF.H' AND LIBMPI                                   
C     CONDITIONS                                                       
C                                                                      
C     COMMON         NONE                                              
C     BLOCKS                                                           
C                                                                      
C     I/O            FOR NBI=0 (OUTPUT OF CPU-TIMES) AND IF            
C                    IERROR NOT EQUALS ZERO.
C                                                                      
C     PRECISION      BECAUSE OF THE WAY MESSAGES ARE ASSEMBLED IN MPI
C                              IT IS NECESSARY TO SPECIFY -DSINGLE WHEN COMPILING
C                              IN SINGLE OR -DDOUBLE (AND -r8 OR EQUIVALENT FLAG) FOR
C                              DOUBLE PRECISION.
C                                                                      
C     LANGUAGE       FORTRAN                                           
C                                                                      
C                                                                     
C     REQUIRED       ATAN
C     RESIDENT                                                        
C     ROUTINES                                                        
C                                                                     
C     REFERENCE      PARALLEL THREE DIMENSIONAL, DIRECT AND VECTORIZED ELLIPTIC
C                    SOLVERS FOR VARIOUS BOUNDARY CONDITIONS 
C                    DFVLR-MITTEILUNG
C                             
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C                             
C                             
      include"mpif.h"
      include"fftw_f77.i"
      DIMENSION Y(IY,JY,*),WORK(*),XLMB(*),ZH(*)
      REAL TIME(3)
C      DIMENSION Y(IY,JY,1),WORK(1),XLMB(1),ZH(1)
C
      integer MPI_COMM_HYDRO
      IMP2=IM+2
      KMP2=KM+2
      IFIRST=1
      IERROR=0
C
C   COMPUTATION OF START-ADRESSES FOR WORK-FIELD IN THE CALL OF
C   SUBROUTINES HOEFEL AND H3DUMM
C

      NZ=1
      NDZ=NZ+KM+2
      NDZH=NDZ+KM+2
      NB=NDZH
      NA=NB+IMP2*MAX0(KM,JM)
      NC=NA+KM
      ND=NC+KM
      NALN=ND+IMP2
      NALM=NALN+JM
      NFAX=NALM+IMP2
      NYT=NFAX+13
      NTRIGS=NYT+IY*JY*KM/mysize

C
      NTOTAL=NTRIGS+MAX0(2*JM,3*IM/2+1)-1
C
      IF (NBI.EQ.0)  THEN                                              
      NBI  = 1                                                         
      NTIME= 1                                                         
      ELSE                                                             
      NTIME= 0                                                         
      ENDIF                                                            
C
      IF((MOD(KM,MYSIZE).GT.0).OR.(MOD(JY,MYSIZE).GT.0))GOTO 997
      IF(NTOTAL.GT.LW)GOTO 999

C
      kml=km/mysize

#ifdef CRAY
      CALL CYCLIC(Y,IY,JY,IM,JM,KML)
#endif

      CALL HOEFEL1(ZH,WORK(NZ),WORK(NDZ),WORK(NDZH),KMP2)


      CALL H3DUMM
     *         (Y,WORK(NYT),WORK(NDZ),WORK(NDZH),WORK(NB),WORK(NA)
     *         ,WORK(NC),WORK(ND),DX,DY,IMP2,JM,KM
     *         ,XLMB,IFIRST,WORK(NALN),WORK(NALM),WORK(NFAX)
     *         ,WORK(NTRIGS),IY,JY,NBK,NTIME,IERROR,mysize,myrank,
     *     MPI_COMM_HYDRO,TIME)
C


      RETURN
C
 997  IERROR=6
      IF(MYRANK.EQ.0) THEN
         WRITE(6,998) IERROR,KM,JY,MYSIZE
      ENDIF
 998  FORMAT('IERROR =',I2,'; ',I3,' OR', I3,' NOT DIVISIBLE BY',I3)
      CALL MPI_FINALIZE()
      STOP
         RETURN

999   IERROR=4
      IF(MYRANK.EQ.0) THEN
         WRITE(6,9990)IERROR
         
         WRITE(6,9991)LW

         WRITE(6,9992)NTOTAL
      ENDIF
      CALL MPI_FINALIZE()
      STOP
9990  FORMAT('IERROR =',I2,' LW TOO SMALL')    
9991  FORMAT(' ALLOCATED LW : ',I5)
9992  FORMAT(' MINIMUM LW REQUIRED: ',I5)
      RETURN
C
      END
C
C
C *************************************************************
C
C

      SUBROUTINE H3DUMM
     +(P,PT,DZ,DZH,B,A,C,Z,DX,DY,IMP2,JM,KM
     +,XLAM,IFIRST,ALN,ALM
     +,IFAX,TRIGS,IY,JY,NBZ,NTIME,IERROR,mysize,myrank,
     +     MPI_COMM_HYDRO,TIME)


C
CB    PArallel solver for Helmholtz problem. Assume slab decomposition.
C
CP    P(IMP2,JM,KM/mysize) , INPUT= DIVERGENCE ON SLAB, OUTPUT=POTENTIAL
CP    XLAM(KM) = HELMHOLTZ-PARAMETER
C     NBZ = 1 ==> NEUMANN-NEUMANN - BOUNDARY-CONTINTION
C     NBZ = 2 ==> DIRICHLET-DIRICHLET - BOUNDARY-CONTINTION
C     NBZ = 3 ==> NEUMANN-DIRICHLET - BOUNDARY-CONTINTION
C     NBZ = 4 ==> DIRICHLET-NEUMANN - BOUNDARY-CONTINTION
CP    ARBEITSFELDER
CP    A,C  (KM) , B(IMP2,MAX(KM,JM)),PT(IMP2,JM/MYSZE,KM)
CP    ALN (JM)
CP    ALM,Z (IMP2)
CP    IFAX(13)
C
      include"mpif.h"
      include"fftw_f77.i"
      DIMENSION C(*),A(*),ALN(*),B(IMP2,*),Z(*)
     *         ,ALM(*),P(IY,JY,*),PT(IY,JY/mysize,*),DZ(*),DZH(*),
     *     XLAM(*)
      DIMENSION IFAX(*),TRIGS(*)
      REAL TIME(3),coeff(1000)
      logical flag
      integer mpi_comm_hydro
#ifdef SGIPOINTER
      integer*8 planf,planb
#else
      integer planf,planb
#endif
      data flag /.true./
      save flag
      save planf,planb
      save coeff
      
C
C      DIMENSION C(1),A(1),ALN(1),B(IMP2,1),Z(1)
C     *         ,ALM(1),P(IY,JY,1),DZ(1),DZH(1),XLAM(1)
C      DIMENSION IFAX(1),TRIGS(1)
C
C      TGAUS1= secnds(0.)
c	TGAUS1=TIME(0.)
      IM=IMP2-2
      IMP1=IM+1
      IYH=IY/2
      IMH=IMP2/2
      KMM1=KM-1
      PI2=ATAN(1.)*8.
      PI=PI2/2.
#ifdef CRAY      
      Q =JM
#endif
#ifdef FFTW
      Q=JM*IM
#endif
#ifdef SGI
      Q=JM*IM
#endif

C

      IF(IMP2.GT.IY.OR.JM.GT.JY)        GOTO 991
      IF(IM.LT.4.OR.JM.LT.4.OR.KM.LT.3) GOTO 992
      IF(NBZ.GT.4.OR.NBZ.LT.1)          GOTO 995
C
c      write(*,*) IFIRST
      IF(IFIRST.EQ.0) GOTO 10
      IFIRST=0
C
      DO 1 K = 1,KM
    1 A(K)=Q/(DZ(K+1)*DZH(K))
C
C       SUBDIAGONAL COEFFICIENT
C
      DO 2 K=1,KM
    2 C(K)=Q/(DZ(K+1)*DZH(K+1))
C
C       SUPERDIAGONAL COEFFICIENT
C
      IF((NBZ.EQ.1).OR.(NBZ.EQ.3))A(1)=0.
      IF((NBZ.EQ.1).OR.(NBZ.EQ.4))C(KM)=0.
C
      HM=2*Q/(DX*DX)
      HN=2*Q/(DY*DY)
C
      DO 3 N=1,JM
    3 ALN(N)=HN*(1.-COS(FLOAT(N-1)*PI2/FLOAT(JM)))
C
      DO 4 M=1,IMP1,2
    4 ALM(M)=HM*(1.-COS(FLOAT(M-1)*PI /FLOAT(IM)))
C
   10 CONTINUE
C     CALL DRUCK(P,JM,IM,IMP2,KM,IY,JY)
C      CPU= secnds(0.)
c	CPU=TIME(0.)
C
C  FFT - SYNTHESIS
C

#ifdef CRAY

      CALL FFTFAX(IM,IFAX,TRIGS)
      IF(IFAX(1).EQ.-99) GOTO 993

      XL=0.
      do k=1,km
         XL=XL+XLAM(K)
      enddo
      cpub=mpi_wtime()
      DO 100 K=1,KM/mysize


      CALL FFT99(P(1,1,K),B,TRIGS,IFAX,1,IY  ,IM,JM,-1)

  100 CONTINUE
      CALL CFTFAX(JM,IFAX,TRIGS)
      DO 101 K=1,KM/mysize
      CALL CFFT99(P(1,1,K),B,TRIGS,IFAX,IYH,1,JM,IMH,-1)
  101 CONTINUE
      TIME(1)=-CPUB+MPI_WTIME()

#endif
#ifdef FFTW

      if(flag) then
         flag=.false.
         call rfftw2d_f77_create_plan(planf,im,jm,FFTW_REAL_TO_COMPLEX,
     &        FFTW_MEASURE+FFTW_IN_PLACE)
         call rfftw2d_f77_create_plan(planb,im,jm,FFTW_COMPLEX_TO_REAL,
     &        FFTW_MEASURE+FFTW_IN_PLACE)
      endif
      cpub=mpi_wtime()

      call rfftwnd_f77_real_to_complex(planf,km/mysize,p,1,(im+2)*jm,b,
     &     1,(im+2)*jm/2)

      TIME(1)=-CPUB+MPI_WTIME()

#endif
#ifdef SGI

      if(flag) then
         flag=.false.
         call scfft2dui(im,jm,coeff)
      endif
      cpub=mpi_wtime()

      do k=1,km/mysize
         call scfft2du(-1,im,jm,p(1,1,k),im+2,coeff)
      enddo

      TIME(1)=-CPUB+MPI_WTIME()
      
#endif
      

      


C     now we swap the arrays 
      CPUB=MPI_WTIME()
      call swap_xz(P,PT,IY*JY,KM,mysize,myrank,mpi_comm_hydro)
      TIME(3)=-CPUB+MPI_WTIME()
C

C BEGIN OF TRIDIAGONAL GAUSS-ALGORITHM
C
      cpub=mpi_wtime()
      DO 300 NC=1,JM/mysize
         N=JM/mysize*myrank+NC
      X=C(1)+ALN(N)+XLAM(1)*Q+A(1)
CVD$  NODEPCHK
      DO 201 M=1,IMP1,2
      Z(M)=-1./(X+ALM(M))
  201 Z(M+1)=-1./(X+ALM(M))
C
      DO 203 M=1,IMP2
      PT(M,NC,1)=PT(M,NC,1)*Z(M)

  203 B(M,1)=C(1)*Z(M)
C
      DO 250 K=2,KMM1
      X=-ALN(N)-A(K)-C(K)-XLAM(K)*Q
CVD$  NODEPCHK
      DO 212 M=1,IMP1,2
      Z(M)=1./(X-ALM(M)-A(K)*B(M,K-1))
  212 Z(M+1)=1./(X-ALM(M)-A(K)*B(M,K-1))
C
      DO 214 M=1,IMP2
      PT(M,NC,K)=(PT(M,NC,K)-A(K)*PT(M,NC,K-1))*Z(M)
  214 B(M,K)=C(K)*Z(M)
C
  250 CONTINUE
C
      X=-ALN(N)-A(KM)-XLAM(KM)*Q-C(KM)
      DO 221 M=1,IMP1,2
  221 B(M,KM)=X-ALM(M)
C
C     VERAENDERUNG DES DIAGONALELEMENTS FUER DEN (0,0)-FOURIER-
C     MODE BEI SINGULAEREM GLEICHUNGSSYSTEM.
C
      IF((N.EQ.1).AND.(XL.EQ.0.).AND.NBZ.EQ.1)
     *B(1,KM)=B(1,KM)-C(KMM1)
CVD$  NODEPCHK
      DO 222 M=1,IMP1,2
      Z(M)=1./(B(M,KM)-A(KM)*B(M,KMM1))
  222 Z(M+1)=1./(B(M,KM)-A(KM)*B(M,KMM1))
C
      DO 224 M=1,IMP2

      PT(M,NC,KM)=(PT(M,NC,KM)-A(KM)*PT(M,NC,KMM1))*Z(M)

  224 CONTINUE
C
      DO 280 K=KMM1,1,-1
      DO 270 M=1,IMP2
  270 PT(M,NC,K)=PT(M,NC,K)-B(M,K)*PT(M,NC,K+1)

  280 CONTINUE
C
  300 CONTINUE
C
      TIME(2)=MPI_WTIME()-CPUB

C     now we swap back

      CPUB=MPI_WTIME()

      call swap_zx(PT,P,IY*JY,KM,mysize,myrank,mpi_comm_hydro)
      TIME(3)=TIME(3)-CPUB+MPI_WTIME()
C      CPU1= secnds(0.)
c	CPU1=TIME(0.)
C
C FFT - ANALYSIS
C
#ifdef CRAY
      cpub=mpi_wtime()
      DO 502 K=1,KM/mysize
      CALL CFFT99(P(1,1,K),B,TRIGS,IFAX,IYH,1,JM,IMH,+1)
  502 CONTINUE
      CALL FFTFAX(IM,IFAX,TRIGS)
      DO 503K=1,KM/mysize
      CALL FFT99(P(1,1,K),B,TRIGS,IFAX,1,IY  ,IM,JM,+1)
  503 CONTINUE
      TIME(1)=TIME(1)-CPUB+MPI_WTIME()

#endif
#ifdef FFTW
      cpub=mpi_wtime()
      call rfftwnd_f77_complex_to_real(planb,km/mysize,p,1,(im+2)*jm/2,
     &     b,1,(im+2)*jm)
      TIME(1)=TIME(1)-CPUB+MPI_WTIME()

#endif

#ifdef SGI
      cpub=mpi_wtime()
      do k=1,km/mysize
         call csfft2du(1,im,jm,p(1,1,k),im+2,coeff)
      enddo

      TIME(1)=TIME(1)-CPUB+MPI_WTIME()

#endif

C
C     CALL DRUCK(P,JM,IM,IMP2,KM,IY,JY)
C
C      CPU2= secnds(0.)
c	CPU2=TIME(0.)
      TFFT=CPU+CPU2-CPU1
C      TGAUSS= secnds(0.)-TFFT-TGAUS1
c	TGAUSS=TIME(0.)-TFFT-TGAUS1
      IF (NTIME.EQ.1) THEN                                           
      TGES = TFFT + TGAUSS                                           
      TRFF = TFFT/TGES * 100.                                        
      TRGA = TGAUSS/TGES * 100.                                      
c      WRITE (6,1000)   TGES,TFFT,TRFF,TGAUSS,TRGA                   
1000  FORMAT (1X,'CPU-TIMES IN S OR %: TOT, FFT, FFT/TOT, GAU, '     
     2,          'PAU/TOT ',2E12.3,F8.2,E12.3,F8.2,/)                
      ENDIF                                                          
C                                                                    
      RETURN
C                                                                   
991   IERROR = 1                                                   
      IF (MYRANK.EQ.0) THEN
         WRITE(6,9910)IERROR
      ENDIF
      CALL MPI_FINALIZE()
      STOP

9910  FORMAT('IERROR = ',I2,' IY<IM+2  OR  JY<JM')                
      RETURN                                                     
992   IERROR = 2                                                
      IF(MYRANK.EQ.0) THEN
         WRITE(6,9920)IERROR
      ENDIF
      CALL MPI_FINALIZE()
      STOP
9920  FORMAT('IERROR = ',I2,' IM<4 OR JM<4 OR KM<3')                  
      RETURN                                                          
993   IERROR = 3                                                      
      IF(MYRANK.EQ.0) THEN
         WRITE(6,9930)IERROR
      ENDIF
      CALL MPI_FINALIZE()
      STOP
9930  FORMAT('IERROR= ',I2,'IM OR JM CONTAIN FACTORS OTHER THAN 2,3,5')
      RETURN                                                          
995   IERROR = 5
      IF(MYRANK.EQ.0) THEN
         WRITE(6,9950)IERROR
      ENDIF
      CALL MPI_FINALIZE()
      STOP
9950  FORMAT('IERROR = ',I2,' NBK OUT OF RANGE (1,...,4)')          
      RETURN                                                       
C                                                                 
      END
C                                                                
C                                                               
C ********************************************************************* 00064900
C                                                                       00065000
C                                                                       00065100
          SUBROUTINE HOEFEL1 (ZH,Z,DZ,DZH,KMP2)                
C                                                                       00010000
C **********************************************************************00020000
C *                                                                    *00030000
C *                DEFINIERT EIN STAGGERED GRID                        *00040000
C *                                                                    *00050000
C *    EINGABE    KMP2 : ANZAHL DER MASCHEN + 2                        *00060000
C *               ZH   : POSITION DER MASCHENRAENDER; DIMENSION KMP1   *00070000
C *                      (KMP1 = KMP2-1)                               *00080000
C *    AUSGABE    Z    : POSITION DER MASCHENMITTEN; DIMENSION KMP2    *00080000
C *                      (JEWEILS EINE POSITION AUSSERHALB DES IN IN-  *00080000
C *                      TEGRATIONSGEBIETS IST ZUSAETZLICH DEFINIERT)  *00080000
C *               DZ   : ABSTAND ZWISCHEN BENACHBARTEN POSITIONEN Z    *00080000
C *                      DIMENSION KMP2                                *00080000
C *               DZH  : ABSTAND ZWISCHEN BENACHBARTEN POSITIONEN ZH   *00080000
C *                      DIMENSION KMP1                                *00080000
C *                                                                    *00080000
C **********************************************************************00090000
C                                                                      
          DIMENSION ZH(1),DZ(1),DZH(1),Z(1)                           
          
C                                                                    
          KMP1 = KMP2-1                                             
C                                                                  
      DO  10010  K=2,KMP1                                         
          Z(K)  = (ZH(K)+ZH(K-1))/2.                             
          DZ(K) = ZH(K)-ZH(K-1)                                 
10010 CONTINUE                                                 
C                                                             
          Z(1)  = ZH(1) - (Z(2)-ZH(1))                       
          DZ(1) = DZ(2)                                     
          Z(KMP2)  = ZH(KMP1) + (ZH(KMP1)-Z(KMP1))         
          DZ(KMP2) = DZ(KMP1)                             
C                                                        
      DO 10020  K=1,KMP1                                
          DZH(K) = Z(K+1)-Z(K)                         
10020 CONTINUE
C                                                     
      RETURN                                         
      END                                           
C                                                  
C                                                 
C ********************************************************************* 00059300
C                                                                       00059400
C                                                                       00059500
      SUBROUTINE CYCLIC(Y,IY,JY,IM,JM,KM)
C                                                                       00010000
C **********************************************************************00020000
C *                                                                    *00030000
C * CYCLICAL COMPLETION OF INPUT ARRAY "Y" IN THE I - DIRECTION. AS    *00040000
C * RESULT "Y" IS CONTAINING SEQUENCES                                 *00050000
C *         Y(IM,*,*),Y(1,*,*),Y(2,*,*),....,Y(IM,*,*),Y(1,*,*)        *00060000
C *                                                                    *00080000
C **********************************************************************00090000
C                                                                       00100000
      DIMENSION Y(IY,JY,*)
C
C      DIMENSION Y(IY,JY,1)
C
      IMP1 = IM + 1
      IMP2 = IM + 2
C
      DO 10 K=1,KM
      DO 20 J=1,JM
          Y(1   ,J,K) = Y(IMP1,J,K)
          Y(IMP2,J,K) = Y(2   ,J,K)
20    CONTINUE
10    CONTINUE
C
      RETURN
      END
