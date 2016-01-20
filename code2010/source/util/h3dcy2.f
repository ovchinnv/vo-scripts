
      SUBROUTINE H3DCY2(Y,IY,JY
     *                 ,IM,JM,KM
     *                 ,NBI,NBJ,NBK,XLMB
     *                 ,DX,DY,ZH
     *                 ,WORK,LW,IERROR)
c                                                                       
C                                                                       
C     * * * * * * * * *  PURPOSE    * * * * * * * * * * * * * * * * * * 
C                                                                       
C     THIS ROUTINE CALCULATES THE SOLUTION OF A HELMHOLTZ EQUATION      
C     IN A CARTESIAN COORDINATE-SYSTEM WITH PRESCRIBED BOUNDARY
C     CONDITIONS (CYCLIC IN I - AND J - DIRECTION, NEUMANN OR
C     DIRICHLET IN K - DIRECTION) FOR A THREE DIMENSIONAL DOMAIN,
C     DISCRETIZED ON A STAGGERED GRID.
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
C     'JY' MUST AT LEAST EQUAL JM.                                   
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
C     THE NUMBER OF UNKNOWNS IN THE K-DIRECTION; KM >= 3              
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
C     LENGTH OF 'WORK'; LW = S1 + S2 + S3 + S4 + S5 + S6              
C                                                                     
C        S1 = (IM+2)*MAX(KM,JM)                                       
C        S2 = (MAX(2*JM,3*IM/2+1)
C        S3 = 2*(IM+2)
C        S4 = JM
C        S5 = 2*KM + (KM+1) + (KM+2)
C        S6 = 13
C                                                                     
C                                                                     
C                                                                     
C                                                                    
C             * * * * * *   ON OUTPUT     * * * * * *                   00014600
C                                                                       00014700
C   Y                                                                   00014800
C     CONTAINS THE SOLUTION X AT LOCATIONS I=2,..,IM+1; J=1,..,JM       00014900
C     AND K=1,..,KM WITH X(1,*,*)    = X(IM+1,*,*) AND                  00015000
C                        X(IM+2,*,*) = X(2,*,*)                         00015100
C                                                                       00015100
C   IERROR                                                              00015200
C     AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.  EXCEPT    00015300
C     FOR NUMBER ZERO, A SOLUTION IS NOT ATTEMPTED.                     00015400
C                                                                       00015500
C     IERROR = 0:  NO ERROR                                             00015600
C     IERROR = 1:  IY<IM+2  OR  JY<JM                                   00015700
C     IERROR = 2:  IM<4     OR  JM<4 OR KM<3                            00015800
C     IERROR = 3:  IM OR JM CONTAIN FACTORS OTHER THAN 2,3,5            00015900
C     IERROR = 4:  LW TOO SMALL                                         00016000
C     IERROR = 5:  NBK OUT OF RANGE (1,..,4)                            00016100
C                                                                       00017000
C   LW                                                                  00017100
C     LW CONTAINS THE REQUIRED LENGTH OF 'WORK'                         00017200
C                                                                       00017300
C                                                                       00017400
C                                                                       00017500
C     * * * * * * *   PROGRAM SPECIFICATIONS    * * * * * * * * * * * * 00017600
C                                                                       00017700
C     DIMENSION OF   Y(IY,JY,KM), XLMB(KM), ZH(KM+1)                    00017800
C     ARGUMENTS      WORK(LW)  (SEE ARGUMENT LIST FOR DETAILS)          00017900
C                                                                       00018000
C     LATEST         JUNE 1984                                          00018100
C     REVISION                                                          00018200
C
C     SUBPROGRAMS  - 'FFT99','CFFT99' AND LOWER LEVEL
C     REQUIRED       SUBROUTINES FOR FAST FOURIER TRANSFORMATION:
C                    FFT99
C                    CFFT99
C                    FFTFAX
C                    CFTFAX
C                    ALL THESE SUBROUTINES ARE CONTAINED IN 'FFT99' ON
C                    TAPE.
C
C     SPECIAL        NONE                                               00019100
C     CONDITIONS                                                        00019200
C                                                                       00019300
C     COMMON         NONE                                               00019400
C     BLOCKS                                                            00019500
C                                                                       00019600
C     I/O            FOR NBI=0 (OUTPUT OF CPU-TIMES) AND IF             00019700
C                    IERROR NOT EQUALS ZERO.
C                                                                       00019800
C     PRECISION      SINGLE                                             00019900
C                                                                       00020000
C     LANGUAGE       FORTRAN                                            00020400
C                                                                       00020500
C                                                                       00020600
C     TIMING AND     SOME TYPICAL VALUES FOR THE EXECUTION TIME ON THE  00020700
C     ACCURACY       DFVLR CRAY-1 COMPUTER ARE LISTED BELOW.            00020800
C                                                                       00021100
C                    TO MEASURE THE ROUTINE'S ACCURACY A RANDOM SOLUTION00021200
C                    ARRAY 'X' IS  CREATED BY USE OF THE CRAY INTRINSIC 00021300
C                    FUNCTION 'RANF'. ACCORDING TO THE LINEAR SYSTEM    00021400
C                    DESCRIBED ABOVE THE RIGHT HAND SIDE 'Y' IS COM-    00021500
C                    PUTED. WITH THIS AS INPUT 'H3DCY2' IS CALLED TO    00021600
C                    PRODUCE AN APPROPRIATE SOLUTION 'Z', WHICH EXHIBITS00021700
C                    A RELAITVE ERROR, DEFINED AS                       00021800
C                                                                       00021900
C                     E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))00022000
C                                                                       00022100
C                    WHERE THE TWO MAXIMA ARE TAKEN OVER ALL I=2,..,IM+100022200
C                    J=1,..,JM AND K=1,..,KM.                           00022300
C                                                                       00022400
C                    TYPICAL VALUES OF 'E' AND 'T' ARE GIVEN IN THE     00022500
C                    TABLE BELOW FOR SOME TYPICAL VALUES OF 'IM', 'JM', 00022600
C                    'KM'.                                              00022700
C                                                                       00022800
C                                                                       00022900
C                    IM  JM  KM     XLMB   T(MSECS)      E              00023000
C                    ----------    ------   -----   --------            00023100
C                                                                       00023200
C                     8   8   8      0         3.0     2.E-13           00023300
C                    16  16  16      0        12.3     4.E-13           00023700
C                    32  32  32      0        77.2     8.E-13           00024100
C                    64  64  64      0       522.0     5.E-10           00024600
C
C                     8   8   8      1        11.4     8.E-14
C                    16  16  16      1        12.3     9.E-14
C                    32  32  32      1        77.2     9.E-14
C                    64  64  64      1       521.0     1.E-13
C                                                                       00024900
C     REQUIRED       ATAN
C     RESIDENT                                                          00025100
C     ROUTINES                                                          00025200
C                                                                       00025800
C     REFERENCE      THREE DIMENSIONAL, DIRECT AND VECTORIZED ELLIPTIC
C                    SOLVERS FOR VARIOUS BOUNDARY CONDITIONS
C                    DFVLR-MITTEILUNG
C                                                                       00025800
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 00025900
C                                                                       00026000
C                                                                       00026100
      DIMENSION Y(IY,JY,*),WORK(*),XLMB(*),ZH(*)
C      DIMENSION Y(IY,JY,1),WORK(1),XLMB(1),ZH(1)
C
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
      NTRIGS=NFAX+13
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
      IF(NTOTAL.GT.LW)GOTO 999
C
      CALL CYCLIC(Y,IY,JY,IM,JM,KM)
      CALL HOEFEL1(ZH,WORK(NZ),WORK(NDZ),WORK(NDZH),KMP2)
c      do i=1,30
c      write(6,*) ZH(i), work(i), work(ndz-1+i)
c      enddo
c      write(*,*) '???'
c      write(6,*) ndz
      CALL H3DUMM
     *         (Y,WORK(NDZ),WORK(NDZH),WORK(NB),WORK(NA)
     *         ,WORK(NC),WORK(ND),DX,DY,IMP2,JM,KM
     *         ,XLMB,IFIRST,WORK(NALN),WORK(NALM),WORK(NFAX)
     *         ,WORK(NTRIGS),IY,JY,NBK,NTIME,IERROR)
C
      RETURN
C
999   IERROR=4
      WRITE(6,9990)IERROR
9990  FORMAT('IERROR =',I2,' LW TOO SMALL')                           
      WRITE(6,9991)LW
9991  FORMAT(' ANGEGEBENES LW : ',I5)
      WRITE(6,9992)NTOTAL
9992  FORMAT(' MINDESTENS ERFORDERLICHES LW : ',I5)
      RETURN
C
      END
C
C
C *************************************************************
C
C
      SUBROUTINE H3DUMM
     +(P,DZ,DZH,B,A,C,Z,DX,DY,IMP2,JM,KM
     +,XLAM,IFIRST,ALN,ALM
     +,IFAX,TRIGS,IY,JY,NBZ,NTIME,IERROR)
C
CB    LOESUNG DER HELMHOLTZ-GLEICHUNG MIT FFT
C     VOLL VEKTORISIERT
C
CP    P(IMP2,JM,KM) , INPUT= DIVERGENZ, OUTPUT=POTENTIAL
CP    XLAM(KM) = HELMHOLTZ-PARAMETER
C     NBZ = 1 ==> NEUMANN-NEUMANN - BOUNDARY-CONTINTION
C     NBZ = 2 ==> DIRICHLET-DIRICHLET - BOUNDARY-CONTINTION
C     NBZ = 3 ==> NEUMANN-DIRICHLET - BOUNDARY-CONTINTION
C     NBZ = 4 ==> DIRICHLET-NEUMANN - BOUNDARY-CONTINTION
CP    ARBEITSFELDER
CP    A,C  (KM) , B(IMP2,MAX(KM,JM))
CP    ALN (JM)
CP    ALM,Z (IMP2)
CP    IFAX(13)
C
      DIMENSION C(*),A(*),ALN(*),B(IMP2,*),Z(*)
     *         ,ALM(*),P(IY,JY,*),DZ(*),DZH(*),XLAM(*)
      DIMENSION IFAX(*),TRIGS(*)
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
      Q =JM
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
      CALL FFTFAX(IM,IFAX,TRIGS)
      IF(IFAX(1).EQ.-99) GOTO 993
      XL=0.
      DO 100 K=1,KM
      XL=XL+XLAM(K)
c      write(6,*) p(1,1,k)
      CALL FFT99(P(1,1,K),B,TRIGS,IFAX,1,IY  ,IM,JM,-1)
c      write(6,*) b
  100 CONTINUE
      CALL CFTFAX(JM,IFAX,TRIGS)
      DO 101 K=1,KM
      CALL CFFT99(P(1,1,K),B,TRIGS,IFAX,IYH,1,JM,IMH,-1)
  101 CONTINUE
C
C      CPU=-CPU+ secnds(0.)
c	CPU=-CPU+TIME(0.)
C     CALL DRUCK(P,JM,IM,IMP2,KM,IY,JY)
C
C BEGIN OF TRIDIAGONAL GAUSS-ALGORITHM
C
      DO 300 N=1,JM
      X=C(1)+ALN(N)+XLAM(1)*Q+A(1)
CVD$  NODEPCHK
      DO 201 M=1,IMP1,2
      Z(M)=-1./(X+ALM(M))
  201 Z(M+1)=-1./(X+ALM(M))
C
      DO 203 M=1,IMP2
      P(M,N,1)=P(M,N,1)*Z(M)
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
      P(M,N,K)=(P(M,N,K)-A(K)*P(M,N,K-1))*Z(M)
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
      P(M,N,KM)=(P(M,N,KM)-A(KM)*P(M,N,KMM1))*Z(M)
  224 CONTINUE
C
      DO 280 K=KMM1,1,-1
      DO 270 M=1,IMP2
  270 P(M,N,K)=P(M,N,K)-B(M,K)*P(M,N,K+1)
  280 CONTINUE
C
  300 CONTINUE
C
C      CPU1= secnds(0.)
c	CPU1=TIME(0.)
C
C FFT - ANALYSIS
C
      DO 502 K=1,KM
      CALL CFFT99(P(1,1,K),B,TRIGS,IFAX,IYH,1,JM,IMH,+1)
  502 CONTINUE
      CALL FFTFAX(IM,IFAX,TRIGS)
      DO 503K=1,KM
      CALL FFT99(P(1,1,K),B,TRIGS,IFAX,1,IY  ,IM,JM,+1)
  503 CONTINUE
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
      WRITE(6,9910)IERROR
9910  FORMAT('IERROR = ',I2,' IY<IM+2  OR  JY<JM')                
      RETURN                                                     
992   IERROR = 2                                                
      WRITE(6,9920)IERROR
9920  FORMAT('IERROR = ',I2,' IM<4 OR JM<4 OR KM<3')                  
      RETURN                                                          
993   IERROR = 3                                                      
      WRITE(6,9930)IERROR
9930  FORMAT('IERROR= ',I2,'IM OR JM CONTAIN FACTORS OTHER THAN 2,3,5')
      RETURN                                                          
995   IERROR = 5                                                     
      WRITE(6,9950)IERROR
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
