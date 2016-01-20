C                                                                       
      SUBROUTINE H3DNDC (Y,IY,JY,IM,JM,KM,NBI,NBJ,NBK,XLMB,
     &                   DX,DY,ZH,WORK,LW,IERROR)                       
C                                                                       
C                                                                       
C     * * * * * * * * *  PURPOSE    * * * * * * * * * * * * * * * * * * 
C                                                                       
C     THIS ROUTINE CALCULATES THE SOLUTION OF A HELMHOLTZ EQUATION      
C     WITH PRESCRIBED BOUNDARY CONDITONS (NEUMANN, DIRICHLET OR         
C     CYCLIC) FOR A THREE DIMENSIONAL DOMAIN, DISCRETIZED ON A          
C     STAGGERED GRID.                                                   
C                                                                       
C     THE 3D-PROBLEM IS REDUCED TO 'JM' 2D-PROBLEMS VIA A FAST-FOURIER- 
C     TRANSFORMATION (FFT). THE 2D-PROBLEMS ARE SOLVED BY USE OF        
C     ROUTINES 'POISSX' (BY U. SCHUMANN) OR 'POISTG' (BY R. SWEET)      
C     WHICH USE DIFFERENT DEFINITIONS FOR THE DOMAIN'S BOUNDARIES       
C     THE DIRICHLET CASE. THE FFT IS PERFORMED USING ROUTINE 'FFT99'    
C     (BY C. TEMPERTON) WHICH IS OPTIMIZED FOR CRAY COMPUTERS VIA       
C     SPECIAL ROUTINES IN 'CRAY ASSEMBLER LANGUAGE (CAL)'. FOR THE      
C     ORDINARY FORTRAN VERSION SEE FILE 'PA54.JOB.CNTL(FFT99)'.         
C                                                                       
C                                                                       
C                                                                       
C     SPECIFICALLY SUBROUTINE 'H3DNDC' SOLVES THE LINEAR SYSTEM OF      
C     EQUATIONS                                                         
C                                                                       
C       (X(I-1,J,K) - 2.*X(I,J,K)   + X(I+1,J,K)) / (DX*DX) +           
C       (X(I,J-1,K) - 2.*X(I,J,K)   + X(I,J+1,K)) / (DY*DY) +           
C       A(K)*X(I,J,K-1) + B(K)*X(I,J,K) + C(K)*X(I,J,K+1)   +           
C                       - XLMB(K-1)*X(I,J,K)                = Y(I,J,K)  
C                                                                       
C       FOR I=2,..,IM+1;  J=2,..,JM+1  AND K=2,..,KM+1                  
C                                                                       
C     WITH SPECIFIC BOUNDARY CONDITIONS FOR EACH CO-ORDINATE. THE       
C     K-DEPENDENT FACTORS 'A', 'B' AND 'C' ALLOW FOR NON-EQUIDISTANT    
C     GRID SPACING IN THE K-DIRECTION; THEY NEED NOT BE PREDEFINED,     
C     BUT ARE CALCULATED FROM THE "MESH-BOUNDARY POSITIONS" STORED      
C     IN 'ZH' (SEE DESCRIPTION OF PARAMETERS FOR DETAILS)               
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
C     THE LINEAR SYSTEM ABOVE AT LOCATIONS I=2,..,IM+1; J=2,..,JM+1     
C     AND K=2,..,KM+1. THE CONTENTS OF LOCATIONS WITH I=1, I=IM+2,      
C     J=1, JM=JM+2, K=1 OR K=KM+2 IS IRRELEVANT. 'Y' MUST BE            
C     DIMENSIONED AS Y(IY,JY,(KM+2))                                    
C                                                                       
C   IY                                                                  
C     FIRST DIMENSION OF 'Y' EXACTLY AS IN THE CALLING PROGRAMME;       
C     'IY' MUST AT LEAST EQUAL IM+2.                                    
C                                                                       
C   JY                                                                  
C     SECOND DIMENSION OF 'Y' EXACTLY AS IN THE CALLING PROGRAMME;      
C     'JY' MUST AT LEAST EQUAL JM+2.                                    
C                                                                       
C   IM                                                                  
C     THE NUMBER OF UNKNOWNS IN THE I-DIRECTION; 2 < IM                 
C                                                                       
C   JM                                                                  
C     THE NUMBER OF UNKNOWNS IN THE J-DIRECTION; JM/2 MAY ONLY CONTAIN  
C     PRIME FACTORS 2, 3 OR 5                                           
C                                                                       
C   KM                                                                  
C     THE NUMBER OF UNKNOWNS IN THE K-DIRECTION; 2 < KM                 
C                                                                       
C   NBI                                                                 
C     SPECIFICATION OF BOUNDARY CONDITIONS IN I-DIRECTION               
C            (NS: NEUMANN-BOUNDARY ON A STAGGERD GRID;                  
C            (DS: DIRICHLET-BOUNDARY ON A STAGGERD GRID;                
C            (D : DIRICHLET-BOUNDARY ON A REGULAR GRID)                 
C             NBI=1,..,4 SOLUTION OF 2D-PROBLEM WITH 'POISSX'           
C             NBI=5,..,8 SOLUTION OF 2D-PROBLEM WITH 'POISTG')          
C                                                                       
C     NBI=1:  X(1,J,K)= X(2,J,K)  AND  X(IM+2,J,K)= X(IM+1,J,K)  (NS/NS)
C     NBI=2:  X(1,J,K)= 0.        AND  X(IM+2,J,K)= 0.           (D /D )
C     NBI=3:  X(1,J,K)= X(2,J,K)  AND  X(IM+2,J,K)= 0.           (NS/D )
C     NBI=4:  X(1,J,K)= 0.        AND  X(IM+2,J,K)= X(IM+1,J,K)  (D /NS)
C     NBI=5:  X(1,J,K)= X(2,J,K)  AND  X(IM+2,J,K)= X(IM+1,J,K)  (NS/NS)
C     NBI=6:  X(1,J,K)=-X(2,J,K)  AND  X(IM+2,J,K)=-X(IM+1,J,K)  (DS/DS)
C     NBI=7:  X(1,J,K)= X(2,J,K)  AND  X(IM+2,J,K)=-X(IM+1,J,K)  (NS/DS)
C     NBI=8:  X(1,J,K)=-X(2,J,K)  AND  X(IM+2,J,K)= X(IM+1,J,K)  (DS/NS)
C                                                                       
C                                                                       
C   NBJ                                                                 
C     SPECIFICATION OF BOUNDARY CONDITIONS IN J-DIRECTION               
C            (NS: NEUMANN-BOUNDARY ON A STAGGERD GRID;                  
C             C: CYCLIC-BOUNDARY)                                       
C                                                                       
C     NBJ=0:  X(I,1,K)=X(I,JM+1,K)  AND  X(I,JM+2,K)=X(I,2,K)    (C)    
C     NBJ=1:  X(I,1,K)=X(I,2,K)     AND  X(I,JM+2,K)=X(I,JM+1,K) (NS/NS)
C                                                                       
C                                                                       
C   NBK                                                                 
C     SPECIFICATION OF BOUNDARY CONDITIONS IN K-DIRECTION               
C     NBK=1:  X(I,J,1)=X(I,J,2)     AND  X(I,J,KM+2)=X(I,J,KM+1) (NS/NS)
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
C     'ZH' MUST BE DIMENSIONED AS ZH(KM+1). SPECIFY                     
C                                                                       
C        ZH(K) = (K-1) * H / KM      K=1,..,KM+1                        
C                                                                       
C     FOR MODEL HEIGHT 'H' AND EQUAL MESH SIZE.                         
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
C     LENGTH OF 'WORK'; LW = 2*KM+4 + MAX(LWFFT,LWPOI)    WITH          
C                                                                       
C        LWFFT = 3*JM/2+1 + IM*(JM+1)                                   
C        LWPOI = 4*IM     + (IM+LOG2(IM)+10)*KM    FOR NBI=1,..,4  AND  
C        LWPOI = 4*IM     + (IM+LOG2(IM)+12)*KM    FOR NBI=5,..,8       
C                                                                       
C     'LWPOI' IS A SUFFICIENT ESTIMATE FOR THE WORK SPACE NEEDED BY     
C     THE POISSON PART OF 'H3DNDC'. AS THE PRECISE VALUE FOR 'LW' IS    
C     RETURNED A DUMMY-CALL IS RECOMMENDED IN CASES WHERE STORAGE       
C     SPACE IS CRUCIAL.                                                 
C                                                                       
C                                                                       
C                                                                       
C             * * * * * *   ON OUTPUT     * * * * * *                   
C                                                                       
C   Y                                                                   
C     CONTAINS THE SOLUTION X AT LOCATIONS I=2,..,IM+1; J=2,..,JM+2     
C     AND K=2,..,KM+1 AND IN THE BOUNDARY-PLANES ACCORDING TO THE       
C     BOUNDARY CONDITION USED. THE GRID-CELLS AT THE DOMAIN'S EDGES     
C     CONTAIN STILL THE VALUES THEY HAD WHEN 'H3DNDC' WAS ENTERED.      
C                                                                       
C     NOTE! FOR THE POISSON CASE (MIN(XLMB)=0) WITH ONLY NEUMANN OR     
C           CYCLIC BOUNDARY CONDITIONS 'Y' CONTAINS AN ARBITRARY MEAN   
C           VALUE. NO ADJUSTMENTS ARE MADE IN 'H3DNDC'.                 
C                                                                       
C   IERROR                                                              
C     AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.  EXCEPT    
C     FOR NUMBER ZERO, A SOLUTION IS NOT ATTEMPTED.                     
C                                                                       
C     IERROR = 0:  NO ERROR                                             
C     IERROR = 1:  IY<IM+2  OR  JY<JM+2                                 
C     IERROR = 2:  IM<3     OR  KM<3                                    
C     IERROR = 3:  JM/2 CONTAINS FACTORS OTHER THAN 2,3,5               
C     IERROR = 4:  LW TOO SMALL                                         
C     IERROR = 5:  NBI OUT OF RANGE (1,..,8)                            
C     IERROR = 6:  NBJ OUT OF RANGE (0,1)                               
C                                                                       
C     IF ONE OF THE FOLLOWING ERROR FLAGS APPEARS A THOROUGH CHECK OF   
C     ROUTINE 'H3DNDC' IS ADVISABLE                                     
C                                                                       
C     IERROR = 10,..,60:   ERROR 1,..,6 IN ROUTINE 'POISSX'             
C     IERROR = 100,..,700: ERROR 1,..,7 IN ROUTINE 'POISTG'             
C                                                                       
C                                                                       
C   LW                                                                  
C     LW CONTAINS THE REQUIRED LENGTH OF 'WORK', EVEN IF IERROR=4       
C                                                                       
C                                                                       
C                                                                       
C     * * * * * * *   PROGRAM SPECIFICATIONS    * * * * * * * * * * * * 
C                                                                       
C     DIMENSION OF   Y(IY,JY,KM+2), XLMB(KM), ZH(KM+1)                  
C     ARGUMENTS      WORK(LW)  (SEE ARGUMENT LIST FOR DETAILS)          
C                                                                       
C     LATEST         MARCH 1984                                         
C     REVISION                                                          
C                                                                       
C     SUBPROGRAMS  - 'FFT99', 'FFTFAX' AND LOWER LEVEL SUBROUTINES FOR  
C     REQUIRED        FAST FOURIER TRANSFORMATION (FROM 'M06LIB')       
C                  - 'POISSX' AND LOWER LEVEL ROUTINES AS "POISSON-SOL- 
C                     VER"; FOR NBI=1,..,4; FROM 'PA08.SRC.FORT(POISSX)'
C                  - 'POISTG' AND LOWER LEVEL ROUTINES AS "POISSON-SOL- 
C                     VER"; FOR NBI=5,..,8; FROM 'PA08.SRC.FORT(POISTG)'
C                                                                       
C     SPECIAL        NONE                                               
C     CONDITIONS                                                        
C                                                                       
C     COMMON         /TEST/ LTE(20) - LOGICAL ARRAY FOR TEST-OUTPUT     
C     BLOCKS         HERE ONLY: LTE(1)=.T. --> OUTPUT OF CPU-TIMES      
C                                                                       
C     I/O            NONE, EXCEPT LTE(1)=.T.  (OUTPUT OF CPU-TIMES)     
C                                                                       
C     PRECISION      SINGLE                                             
C                                                                       
C     ORIGINATOR     HANS VOLKERT, INSTITUT FUER PHYSIK DER ATMOSPHAERE,
C                    DFVLR, D-8031 WESSLING                             
C                                                                       
C     LANGUAGE       FORTRAN 77                                         
C                                                                       
C     TIMING AND     SOME TYPICAL VALUES FOR THE EXECUTION TIME ON THE  
C     ACCURACY       DFVLR CRAY-1 COMPUTER ARE LISTED BELOW. THE RELA-  
C                    TIVE CONTRIBUTIONS OF THE FAST FOURIER TRANSFORMA- 
C                    TION DECREASE WEAKLY WITH INCREASING DIMENSIONS;   
C                    THEY LIE IN THE ORDER OF 3% FOR NBJ=0 AND IN THE   
C                    ORDER OF 10% FOR NBJ=1.                            
C                                                                       
C                    TO MEASURE THE ROUTINE'S ACCURACY A RANDOM SOLUTION
C                    ARRAY 'X' IS  CREATED BY USE OF THE CRAY INTRINSIC 
C                    FUNCTION 'RANF'. ACCORDING TO THE LINEAR SYSTEM    
C                    DESCRIBED ABOVE THE RIGHT HAND SIDE 'Y' IS COM-    
C                    PUTED. WITH THIS AS INPUT 'H3DNDC' IS CALLED TO    
C                    PRODUCE AN APPROPRIATE SOLUTION 'Z', WHICH EXHIBITS
C                    A RELAITVE ERROR, DEFINED AS                       
C                                                                       
C                     E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))
C                                                                       
C                    WHERE THE TWO MAXIMA ARE TAKEN OVER ALL I=2,..,IM+1
C                    J=2,..,JM+1 AND K=2,..,KM+1.                       
C                                                                       
C                    TYPICAL VALUES OF 'E' AND 'T' ARE GIVEN IN THE     
C                    TABLE BELOW FOR SPECIFIC VALUES OF 'IM', 'JM', 'KM'
C                    AND BOUNDARY PARAMETERS 'NBI' AND 'NBJ'. FOR       
C                    THESE TESTS A UNIFORM MESHSIZE AND VANISHING       
C                    HELMHOLTZ-PARAMETERS WERE CHOSEN; I.E.             
C                    DX=DY=1; ZH(K)=0,..,KM; XLMB=0.                    
C                    NOTE: MACHINE ACCURACY OF CRAY-1S = 7.015 E-15     
C                                                                       
C                                                                       
C                    IM  JM  KM     NBI      NBJ     T (S)        E     
C                    ----------    -----    -----   -------     ------  
C                                                                       
C                     9   8   8      1        0       .008      9.E-14  
C                     9   8   8      5        0       .010      9.E-14  
C                     9   8   8      1        1       .010      9.E-14  
C                     9   8   8      5        1       .012      9.E-14  
C                                                                       
C                    17  16  16      1        0       .055      2.E-13  
C                    17  16  16      5        0       .065      2.E-13  
C                    17  16  16      1        1       .061      2.E-13  
C                    17  16  16      5        1       .072      2.E-13  
C                                                                       
C                    33  32  32      1        0       .42       3.E-13  
C                    33  32  32      5        0       .51       3.E-13  
C                    33  32  32      1        1       .45       3.E-13  
C                    33  32  32      5        1       .54       3.E-13  
C                                                                       
C                    65  64  64      1        0      3.5        4.E-13  
C                    65  64  64      5        0      4.4        4.E-13  
C                    65  64  64      1        1      3.7        5.E-13  
C                    65  64  64      5        1      4.5        5.E-13  
C                                                                       
C     REQUIRED       SIN, COS, ATAN, ALOG, MAX0, SECOND                 
C     RESIDENT                                                          
C     ROUTINES                                                          
C                                                                       
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C                                                                       
C                                                                       
      DIMENSION Y(IY,JY,KM+2),ZH(*),WORK(*),XLMB(*)                        
C                                                                       
C
      NDZ  = 1                                                          
      NDZH = NDZ  + KM+2                                                
      NZ   = NDZH + KM+2                                                
      NTRI = NZ                                                         
      NWFF = NTRI + 3*JM/2+1                                            
      NHI  = NZ                                                         
      NSI  = NHI  + JM                                                  
      NCO  = NSI  + JM                                                  
      NARG = NZ                                                         
      NA   = NARG + JM                                                  
      NB   = NA   + KM                                                  
      NC   = NB   + KM                                                  
      NQE  = NC   + KM                                                  
      NWPO = NQE  + IM*KM                                          
C
      SOMMA = NDZ+NDZH+NA+NB+NC + NTRI+NWFF+NHI+
     &       NSI+NCO+NARG+NQE+NWPO+NZ
C
      NTOTAL = NTRIGS + MAX0(2*JM,3*IM/2+1) - 1
C                                                                       
      CALL HDXX   (Y,IY,JY,IM,JM,KM,NBI,NBJ,NBK,XLMB,DX,DY,ZH,LW,IERROR 
     2,            WORK(NDZ),WORK(NDZH),WORK(NZ),WORK(NTRI),WORK(NWFF)  
C** ERKLAERUNG:        DZ        DZH        Z       TRIGS      WFFT     
     3,            WORK(NHI),WORK(NSI),WORK(NCO),WORK(NARG),WORK(NA)    
C** ERKLAERUNG:        HI        SI        CO        ARGJ       A       
     4,            WORK(NB),WORK(NC),WORK(NQE),WORK(NWPO))              
C** ERKLAERUNG:        B        C        QE        WPOI                 
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE HDXX   (Y,IY,JY,IM,JM,KM,NBI,NBJ,NBK,XLMB,DX,DY,ZH,LW  
     2,                  IERROR,DZ,DZH,Z,TRIGS,WFFT,HI,SI,CO,ARGJ,A,B,C 
     3,                  QE,WPOI)                                       
C                                                                       
      LOGICAL   LTE(20)                                                 
      DIMENSION Y(IY,JY,*),ZH(*),XLMB(*)                                
      DIMENSION DZ(*),DZH(*),Z(*),TRIGS(*),WFFT(*),A(*),B(*),C(*)       
      DIMENSION QE(*),WPOI(*),HI(*),SI(*),CO(*),ARGJ(*)                 
      DIMENSION NBITRA(8),IFAX(13)                                      
      COMMON /TEST/ LTE                                                 
      PI  = 4.*ATAN(1.)                                                 
      IERROR = 0                                                        
C                                                                       
C
      IMP1= IM+1                                                        
      IMP2= IM+2                                                        
C                                                                       
      JMM1= JM-1                                                        
      JMP1= JM+1                                                        
      JMP2= JM+2                                                        
      FAKNS = PI / (2.*FLOAT(JM))                                       
C                                                                       
      KMP1= KM+1                                                        
      KMP2= KM+2                                                        
C                                                                       
      NBITRA(1) = 1                                                     
      NBITRA(2) = 2                                                     
      NBITRA(3) = 3                                                     
      NBITRA(4) = 4                                                     
      NBITRA(5) = 3                                                     
      NBITRA(6) = 1                                                     
      NBITRA(7) = 4                                                     
      NBITRA(8) = 2                                                     
C                                                                       
      TVN  = -1.E10                                                     
      TFFT = -1.E10                                                     
      TPOI = -1.E10                                                     
C                                                                       
      IF (IY.LT.IMP2 .OR. JY.LT.JMP2)                      GOTO 991     
      IF (IM.LT.3    .OR. KM.LT.3)                         GOTO 992     
C                                                                       
      LWFFT= 3*JM/2+1 + IM*(JM+1)                                       
      LWA  = LW                                                         
C                                                                       
      IF (NBI.LE.0   .OR. NBI.GT.8)                        GOTO 995     
      IF (NBJ.LT.0   .OR. NBJ.GT.1)                        GOTO 996     
      IF (NBK.LT.1   .OR. NBK.GT.1)                        GOTO 997     
C                                                                       
C                                                                       
      DXQ = DX*DX                                                       
      DYQ = DY*DY                                                       
      XLM = 1.E10                                                       
      DO 222 K=1,KM                                                     
      XLM     = AMIN1(XLM,XLMB(K))                                      
222   CONTINUE                                                          
      CALL HOEFEL (ZH,Z,DZ,DZH,KMP2)                                    
C                                                                       
C                                                                       
C      CPU1 = secnds(0.)                                                 
C                                                                       
C== ES FOLGT VORBEREITUNG FUER "HIN-TRANSFORMATION" =================== 
C                                                                       
C==       A)  ZYKLISCHE RAENDER IN J-RICHTUNG       =================== 
      IF (NBJ.EQ.0) THEN                                                
      DO 111   I=2,IMP1                                                 
      DO 111   K=2,KMP1                                                 
      Y(I,JMP2,K) = Y(I,2,K)                                            
      Y(I,1,K)    = Y(I,JMP1,K)                                         
111   CONTINUE                                                          
C                                                                       
C==       B)  NEUMANN-RAENDER IN J-RICHTUNG         =================== 
      ELSE                                                              
      DO 122   I=2,IMP1                                                 
      DO 122   K=2,KMP1                                                 
      Y(I,1,K) = Y(I,2,K)                                               
      Y(I,2,K) = 0.                                                     
      Y(I,JMP2,K) = 0.                                                  
122   CONTINUE                                                          
      DO 133   I=2,IMP1                                                 
      DO 133   K=2,KMP1                                                 
      DO 133   J=3,JMM1,2                                               
      Y(I,J,K) = (Y(I,J,K) + Y(I,J+1,K)) * 0.5                          
133   CONTINUE                                                          
      DO 155   I=2,IMP1                                                 
      DO 155   K=2,KMP1                                                 
CDIR$ IVDEP                                                             
      DO 155   J=4,JM,2                                                 
      Y(I,J,K) =  Y(I,J,K) - Y(I,J-1,K)                                 
155   CONTINUE                                                          
      ENDIF                                                             
C      CPU2 = secnds(0.)                                                 
      TVN  = CPU2 - CPU1                                                
C                                                                       
C                                                                       
C== ES FOLGT "HIN-TRANSFORMATION" MIT 'FFT99' (SCHLEIFE '11') ========= 
C                                                                       
C      CPU1 = secnds(0.)                                                 
      IF (NBJ.EQ.0) THEN                                                
      ISIGN = -1                                                        
      ELSE                                                              
      ISIGN = 1                                                         
      ENDIF                                                             
C                                                                       
      JUMP  =  1                                                        
      CALL FFTFAX (JM,IFAX,TRIGS)                                       
      IF (IFAX(1).EQ.-99)                                  GOTO 993     
      DO  11   K=2,KMP1                                                 
      CALL FFT99  (Y(2,1,K),WFFT,TRIGS,IFAX,IY,JUMP,JM,IM,ISIGN)        
11    CONTINUE                                                          
C
      TFFT = CPU2 - CPU1                                                
C                                                                       
C                                                                       
C== ES FOLGT NACHBEHANDLUNG DER "HIN-FFT" (NEUM.-RAENDER IN J-RICHT.) = 
C      CPU1 =  secnds(0.)                                                
      IF (NBJ.EQ.1) THEN                                                
      DO 166   I=2,IMP1                                                 
      DO 166   K=2,KMP1                                                 
      Y(I,1,K) = Y(I,2,K)                                               
166   CONTINUE                                                          
      DO 200   NS=2,JM                                                  
      ARG = (NS-1)*FAKNS                                                
      SI(NS) = SIN(ARG)                                                 
      CO(NS) = COS(ARG)                                                 
200   CONTINUE                                                          
      DO 199   I=2,IMP1                                                 
      DO 199   K=2,KMP1                                                 
      DO 177   NS=2,JM                                                  
      HI(NS) = Y(I,NS+1,K)                                              
177   CONTINUE                                                          
      DO 188   NS=2,JM                                                  
      Y(I,NS+1,K) = ((SI(NS)+CO(NS)) * HI(NS) -                         
     2               (SI(NS)-CO(NS)) * HI(JM-NS+2)) * 0.5               
188   CONTINUE                                                          
199   CONTINUE                                                          
      ENDIF                                                             
C      CPU2 =  secnds(0.)                                                
      TVN  = CPU2 - CPU1 + TVN                                          
C                                                                       
C                                                                       
C      CPU1 =  secnds(0.)                                                
C== ES FOLGT FESTLEGUNG DER "WELLENZAHLEN" ============================ 
      IF (NBJ.EQ.0)  THEN                                               
      DO  77   J=1,JM                                                   
      ARGJ(J) = FLOAT(J/2) * 2.*PI/FLOAT(JM)                            
77    CONTINUE                                                          
      ELSE                                                              
      DO  88   J=1,JM                                                   
      ARGJ(J) = FLOAT(J-1) * PI/FLOAT(JM)                               
88    CONTINUE                                                          
      ENDIF                                                             
C                                                                       
C== ES FOLGT FESTLEGUNG DER K-ABHAENGIGEN "AUSSENFAKTOREN" ============ 
      DO 244   K=1,KM                                                   
      A(K) = DXQ/(DZ(K+1)*DZH(K))                                       
      C(K) = DXQ/(DZ(K+1)*DZH(K+1))                                     
244   CONTINUE                                                          
C                                                                       
C                                                                       
C                                                                       
C== ES FOLGT SCHLEIFE '22' UEBER DIE MODEN ============================ 
C                                                                       
C                                                                       
      JP = 0                                                            
      DO  22   J=1,JM                                                   
C                                                                       
      IK = 0                                                            
      DO  33   I=1,IM                                                   
      DO  33   K=1,KM                                                   
      IK = IK+1                                                         
      QE(IK)= Y(I+1,J+JP,K+1) * DXQ                                     
33    CONTINUE                                                          
C                                                                       
C                                                                       
C== ES FOLGT FESTLEGUNG DER K- UND MODE-ABHAENGIGEN "ZENTRALFAKTOREN" = 
      A(1)  = DXQ/(DZ(2)*DZH(1))                                        
      C(KM) = DXQ/(DZ(KM+1)*DZH(KM+1))                                  
      DO  44   K=1,KM                                                   
      B(K)  = -A(K)-C(K) - XLMB(K)*DXQ + DXQ/DYQ*(2.*COS(ARGJ(J))-2.)   
44    CONTINUE                                                          
      B(1)  = B(1) + A(1)                                               
      B(KM) = B(KM)+ C(KM)                                              
      A(1)  = 0.                                                        
      C(KM) = 0.                                                        
C                                                                       
C                                                                       
C                                                                       
C== ES FOLGT LOESUNG DER POISSON-GLEICHUNG FUER MODE 'J/2' BZW. 'J-1' = 
      IF (NBI.LE.4)  THEN                                               
C                                                                       
C== VERSION A: MIT 'POISSX' (NBI=1,..,4) ============================== 
C                                                                       
      MPEROD = 0                                                        
      SING   = 1.                                                   
      IF (NBI.EQ.1.AND.JP.EQ.0.AND.NBK.EQ.1.AND.XLM.EQ.0.) SING = 0.    
      CALL POISSX (IM,NBI,MPEROD,SING,KM,A,B,C,KM,QE,WPOI,LW,IER)       
      LWPOI = LW + (IM+3)*KM + JM                                       
      LW    = 2*KM+4 + MAX0(LWFFT,LWPOI)
      IF (IER.EQ.5)                                        GOTO 994     
      IERROR = IER*10                                                   
      ELSE                                                              
C                                                                       
C== VERSION B: MIT 'POISTG' (NBI=5,..,8) ============================== 
C                                                                       
52    MPEROD = 1                                                        
      NPEROD = NBITRA(NBI)                                              
C ++++ 
      CALL POISTG (NPEROD,IM,MPEROD,KM,A,B,C,KM,QE,IERR,WPOI)           
      LWPOI = JM + (IM+3)*KM  + INT(WPOI(1))                            
      LW    = 2*KM+4 + MAX0(LWFFT,LWPOI)                                
      ENDIF                                                             
      IF (LWA.LT.LW)                                       GOTO 994     
      IERROR = IERR*100                                                 
C                                                                       
C                                                                       
C                                                                       
      IK = 0                                                            
      DO  55   I=1,IM                                                   
      DO  55   K=1,KM                                                   
      IK = IK+1                                                         
      Y(I+1,J+JP,K+1) = QE(IK)                                          
55    CONTINUE                                                          
C                                                                       
      JP = 1                                                            
22    CONTINUE                                                          
C      CPU2 =  secnds(0.)                                                
      TPOI = CPU2 - CPU1                                                
C                                                                       
C== ENDE SCHLEIFE '22' UEBER DIE MODEN ================================ 
C                                                                       
C                                                                       
C                                                                       
C== ES FOLGT VORBEHANDLUNG FUER "RUECK-FFT" (NEUM.-RAENDER IN J-RICH.) =
C      CPU1 =  secnds(0.)                                                
      IF (NBJ.EQ.1) THEN                                                
      DO 311   I=2,IMP1                                                 
      DO 311   K=2,KMP1                                                 
      Y(I,2,K) = Y(I,1,K)                                               
311   CONTINUE                                                          
      DO 299   NS=2,JM                                                  
      ARG = (NS-1)*FAKNS                                                
      SI(NS) = SIN(ARG)                                                 
      CO(NS) = COS(ARG)                                                 
299   CONTINUE                                                          
      DO 300   I=2,IMP1                                                 
      DO 300   K=2,KMP1                                                 
      DO 322   NS=2,JM                                                  
      HI(NS) = Y(I,NS+1,K)                                              
322   CONTINUE                                                          
      DO 333   NS=2,JM                                                  
      ARG = (NS-1)*FAKNS                                                
      Y(I,NS+1,K) = (SI(NS)+CO(NS)) * HI(NS) +                          
     2              (SI(NS)-CO(NS)) * HI(JM-NS+2)                       
333   CONTINUE                                                          
300   CONTINUE                                                          
      DO 399   I=2,IMP1                                                 
      DO 399   K=2,KMP1                                                 
      Y(I,JMP2,K) = Y(I,2,K)                                            
      Y(I,1,K)    = Y(I,JMP1,K)                                         
399   CONTINUE                                                          
      ENDIF                                                             
C      CPU2 =  secnds(0.)                                                
      TVN  = CPU2 - CPU1 + TVN                                          
C                                                                       
C                                                                       
C== ES FOLGT "RUECK-TRANSFORMATION" MIT 'FFT99' (SCHLEIFE '66') ======= 
C      CPU1 =  secnds(0.)                                                
      IF (NBJ.EQ.0) THEN                                                
      ISIGN = 1                                                         
      ELSE                                                              
      ISIGN =-1                                                         
      ENDIF                                                             
C                                                                       
      CALL FFTFAX (JM,IFAX,TRIGS)                                       
      DO  66   K=2,KMP1                                                 
      CALL FFT99  (Y(2,1,K),WFFT,TRIGS,IFAX,IY,JUMP,JM,IM,ISIGN)        
66    CONTINUE                                                          
      TFFT = CPU2 - CPU1 + TFFT                                         
C                                                                       
C                                                                       
C      CPU1 =  secnds(0.)                                                
C== ES FOLGT NACHBEHANDLUNG DER "RUECK-FFT" (NEU.-RAEND. IN J-RICHT.) = 
      IF (NBJ.EQ.1) THEN                                                
      DO 344   I=2,IMP1                                                 
      DO 344   K=2,KMP1                                                 
      Y(I,2,K) = Y(I,1,K)                                               
344   CONTINUE                                                          
      DO 355   I=2,IMP1                                                 
      DO 355   K=2,KMP1                                                 
      DO 355   J=3,JMM1,2                                               
      Y(I,J,K) = Y(I,J,K) - Y(I,J+1,K)                                  
355   CONTINUE                                                          
      DO 366   I=2,IMP1                                                 
      DO 366   K=2,KMP1                                                 
CDIR$ IVDEP                                                             
      DO 366   J=4,JM,2                                                 
      Y(I,J,K) = 2.*Y(I,J,K) + Y(I,J-1,K)                               
366   CONTINUE                                                          
      ENDIF                                                             
C                                                                       
C                                                                       
C== ES FOLGT VERLAENGERUNG VON 'Y' IN RANDMASCHEN GEMAESS RANDBED. ===  
C====== A)  I- RICHTUNG            ==================================== 
      IF (NBI.EQ.1 .OR. NBI.EQ.5)  THEN                                 
      DO 411   K=2,KMP1                                                 
      DO 411   J=2,JMP1                                                 
      Y(1,J,K)    = Y(2,J,K)                                            
      Y(IMP2,J,K) = Y(IMP1,J,K)                                         
411   CONTINUE                                                          
      ELSE IF (NBI.EQ.2)  THEN                                          
      DO 422   K=2,KMP1                                                 
      DO 422   J=2,JMP1                                                 
      Y(1,J,K)    = 0.                                                  
      Y(IMP2,J,K) = 0.                                                  
422   CONTINUE                                                          
      ELSE IF (NBI.EQ.3)  THEN                                          
      DO 433   K=2,KMP1                                                 
      DO 433   J=2,JMP1                                                 
      Y(1,J,K)    = Y(2,J,K)                                            
      Y(IMP2,J,K) = 0.                                                  
433   CONTINUE                                                          
      ELSE IF (NBI.EQ.4)  THEN                                          
      DO 444   K=2,KMP1                                                 
      DO 444   J=2,JMP1                                                 
      Y(1,J,K)    = 0.                                                  
      Y(IMP2,J,K) = Y(IMP1,J,K)                                         
444   CONTINUE                                                          
      ELSE IF (NBI.EQ.6)  THEN                                          
      DO 466   K=2,KMP1                                                 
      DO 466   J=2,JMP1                                                 
      Y(1,J,K)    = -Y(2,J,K)                                           
      Y(IMP2,J,K) = -Y(IMP1,J,K)                                        
466   CONTINUE                                                          
      ELSE IF (NBI.EQ.7)  THEN                                          
      DO 477   K=2,KMP1                                                 
      DO 477   J=2,JMP1                                                 
      Y(1,J,K)    =  Y(2,J,K)                                           
      Y(IMP2,J,K) = -Y(IMP1,J,K)                                        
477   CONTINUE                                                          
      ELSE IF (NBI.EQ.8)  THEN                                          
      DO 488   K=2,KMP1                                                 
      DO 488   J=2,JMP1                                                 
      Y(1,J,K)    = -Y(2,J,K)                                           
      Y(IMP2,J,K) =  Y(IMP1,J,K)                                        
488   CONTINUE                                                          
      ELSE IF (NBI.EQ.9)  THEN                                          
      DO 498   K=2,KMP1                                                 
      DO 498   J=2,JMP1                                                 
      Y(IMP2,J,K) =  Y(1,J,K)                                        
      Y(1,J,K)    =  Y(IMP22,J,K)                                           
498   CONTINUE                                                          
      ENDIF                                                             
C                                                                       
C====== B)  J- RICHTUNG            ==================================== 
      IF (NBJ.EQ.0)  THEN                                               
      DO 511   I=2,IMP1                                                 
      DO 511   K=2,KMP1                                                 
      Y(I,1,K)    = Y(I,JMP1,K)                                         
      Y(I,JMP2,K) = Y(I,2,K)                                            
511   CONTINUE                                                          
      ELSE IF (NBJ.EQ.1)  THEN                                          
      DO 522   I=2,IMP1                                                 
      DO 522   K=2,KMP1                                                 
      Y(I,1,K)    = Y(I,2,K)                                            
      Y(I,JMP2,K) = Y(I,JMP1,K)                                         
522   CONTINUE                                                          
      ENDIF                                                             
C                                                                       
C====== C)  K- RICHTUNG            ==================================== 
      IF (NBK.EQ.1)  THEN                                               
      DO 611   I=2,IMP1                                                 
      DO 611   J=2,JMP1                                                 
      Y(I,J,1)    = Y(I,J,2)                                            
      Y(I,J,KMP2) = Y(I,J,KMP1)                                         
611   CONTINUE                                                          
      ENDIF                                                             
C      CPU2 =  secnds(0.)                                                
      TVN  = CPU2 - CPU1 + TVN                                          
C                                                                       
C                                                                       
C====== CPU-ZEITEN FALLS LTE(1)=.TRUE.  =============================== 
      IF (LTE(1))  THEN                                                 
      TGES = TVN + TFFT + TPOI                                          
      TRVN = TVN /TGES * 100.                                           
      TRFF = TFFT/TGES * 100.                                           
      TRPO = TPOI/TGES * 100.                                           
      WRITE (6,1000)   TGES,TVN,TFFT,TPOI,TRVN,TRFF,TRPO                
1000  FORMAT (1X,'CPU-TIMES (IN S OR %) - TOT,VN,FFT,POI,VN/TOT,FFT/TOT'
     2,         ' POI/TOT:',4E11.3,3F7.2)                               
      ENDIF                                                             
C                                                                       
C                                                                       
C                                                                       
      RETURN                                                            
C                                                                       
C                                                                       
991   IERROR = 1                                                        
      RETURN                                                            
992   IERROR = 2                                                        
      RETURN                                                            
993   IERROR = 3                                                        
      RETURN                                                            
994   IERROR = 4                                                        
      RETURN                                                            
995   IERROR = 5                                                        
      RETURN                                                            
996   IERROR = 6                                                        
      RETURN                                                            
997   IERROR = 7                                                        
      RETURN                                                            
C                                                                       
      END                                                               
C                                                                       
C                                                                       
C ********************************************************************* 
C                                                                       
C                                                                       
      SUBROUTINE HOEFEL (ZH,Z,DZ,DZH,KMP2)                              
      DIMENSION ZH(1),DZ(1),DZH(1),Z(1)                                 
      KMP1 = KMP2-1                                                     
C                                                                       
      DO  11   K=2,KMP1                                                 
      DZ(K) = ZH(K)-ZH(K-1)                                             
      Z(K)  = (ZH(K)+ZH(K-1))/2.                                        
11    CONTINUE                                                          
      DZ(1) = DZ(2)                                                     
      Z(1)  = ZH(1) - 0.5*DZ(1)                                         
      DZ(KMP2) = DZ(KMP1)                                               
      Z(KMP2)  = ZH(KMP1) + 0.5*DZ(KMP2)                                
C                                                                       
      DO  22   K=1,KMP1                                                 
22    DZH(K) = Z(K+1)-Z(K)                                              
C                                                                       
      RETURN                                                            
      END                                                               
C                                                                       
C                                                                       
C ********************************************************************* 
C                                                                       
C                                                                       
c      SUBROUTINE FFTFAX(N,IFAX,TRIGS)                                  
c      DIMENSION IFAX(13),TRIGS(1)                                      
c      DATA MODE /3/                                                    
c      CALL FAX (IFAX, N, MODE)                                         
c      I = IFAX(1)                                                      
c      IF (IFAX(I+1) .GT. 5)   IFAX(1) = -99                            
c      IF (IFAX(1) .GT. 0)     CALL FFTRIG (TRIGS,N,MODE)               
c      RETURN                                                           
c      END                                                              
