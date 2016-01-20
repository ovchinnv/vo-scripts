      SUBROUTINE POISSX(N,ND,MPEROD,SING,M,A,B,C,IDIMY,Y,W,LW,IERROR)   
C                                                                       
C POISSX                 SOLUTION OF POISSON(S) EQUATION WITH NEUMANN   
C                        OR DIRICHLET                                   
C                        BOUNDARY CONDITION ON A STAGGERED GRID.        
C                                                                       
C                        NOVEMBER 1975                                  
C                        CHANGED NOV. 9, 83 TO TO ENSURE CORRECT EXECU- 
C                                        TION ON CRAY AND IBM COMPUTERS 
C                                        (EXACT POSITIONS ARE MARKED)   
C                                                                       
C PURPOSE                 POISSX SOLVES THE LINEAR SYSTEM OF EQUATIONS  
C                                                                       
C                            A(I)*X(I-1,J)+B(I)*X(I,J)+C(I)*X(I+1,J)+   
C                            X(I,J-1)-2*X(I,J)+X(I,J+1) = Y(I,J)        
C                                I = 1,2,...,M, J = 1,2,...,N.          
C                                                                       
C                        HERE, THE INDICES I+1 AND I-1 ARE EVALUATED    
C                        MODULO M; I.E., X(0,J) = X(M,J) AND            
C                        X(M-1,J) = X(1,J).  MOREOVER, WE ASSUME        
C                        EITHER DIRICHLET OR NEUMANN BOUNDARY CONDITIONS
C                        (THE LATTER ON A STAGGERED GRID), E.G. FOR J=0:
C                        NEUMANN   X(I,0)=X(I,1), DIRICHLET  X(I,0)=0   
C                                                                       
C                        PARAMETER DESCRIPTION AND PROGRAM SPECIFICATION
C                        -----------------------------------------------
C                                                                       
C                          THE ROUTINE IS USED BY THE STATEMENT         
C                                                                       
C                        CALL POISSX(N,ND,MPEROD,SING,M,A,B,C,IDIMY,Y   
C                                    ,W,LW,IERROR)                      
C                                                                       
C                          WHERE THE ARGUMENTS ARE DEFINED AS           
C                                                                       
C ON INPUT               N                                              
C                        THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.  N  
C                        MAY BE ANY INTEGER GREATER THAN 1  .           
C                                                                       
C                                                                       
C                        ND             J=0    -     J=N+1              
C                          = 1       NEUMANN   -   NEUMANN              
C                          = 2       DIRICHLET -   DIRICHLET            
C                          = 3       NEUMANN   -   DIRICHLET            
C                          = 4       DIRICHLET -   NEUMANN              
C                        MPEROD                                         
C                          = 0  IF A(1) AND C(M) ARE BOTH ZERO.         
C                               E.G. NON-PERIODIC BOUNDARY CONDITION    
C                          = 1  IF A(1) AND C(M) ARE NOT ZERO.          
C                               E.G. PERIODIC BOUNDARY CONDITION IN X   
C                                                                       
C                        SING                                           
C                          SING CONTROLS THE INTERNAL TEST FOR          
C                          SINGULARITY OF THE GIVEN EQUATION-SYSTEM     
C                          SING = 0. IF SINGULAR SYSTEM                 
C                            IN THIS CASE ONE EQUATION IS REPLACED      
C                            TO MAKE THE SYSTEM DEFINITE                
C                          SING .EQ. 1., IT IS ASSUMED THAT THE         
C                            SYSTEM IS EITHER  DEFINITE OR THE          
C                            SINGULARITY CAN BE DETECTED BY THE         
C                            NUMERICAL SCHEME IN SPITE OF ROUND-OFF     
C                            ERRORS. IN THE LATTER CASE ONE EQUATION    
C                            IS REPLACED AS FOR SING = 0.               
C                                                                       
C                        M                                              
C                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.  M
C                          MAY BE ANY INTEGER GREATER THAN 1  .         
C                                                                       
C                        A,B,C                                          
C                          ONE-DIMENSIONAL ARRAYS OF LENGTH M WHICH     
C                          SPECIFY THE COEFFICIENTS IN THE LINEAR       
C                          EQUATIONS GIVEN ABOVE.                       
C                                                                       
C                        IDIMY                                          
C                          THE ROW (OR FIRST) DIMENSION OF THE          
C                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS IN THE 
C                          PROGRAM CALLING POISSN.  THIS PARAMETER IS   
C                          USED TO SPECIFY THE VARIABLE DIMENSION OF Y. 
C                          IDIMY MUST BE AT LEAST M.                    
C                                                                       
C                        Y                                              
C                          A TWO-DIMENSIONAL ARRAY WHICH SPECIFIES THE  
C                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
C                          OF EQUATIONS GIVEN ABOVE.  Y MUST BE         
C                          DIMENSIONED IDIMY*N.                         
C                                                                       
C                        W                                              
C                          A ONE-DIMENSIONAL ARRAY WHICH MUST BE        
C                          PROVIDED BY THE USER FOR WORK SPACE.  IF     
C                          STORAGE SPACE IS NOT CRUCIAL, TAKE           
C                          4*N+(ALOG(N)/ALOG(2.)+7)*M                   
C                          AS THE LENGTH OF W.                          
C                          OTHERWISE, THE LENGTH OF W HAS TO BE         
C                          L1*N+(L2+MPEROD)*M-L3.  IF N2 IS THE LARGEST 
C                          POWER OF 2 LESS THAN N (E.G., FOR            
C                          N = 9,10,11,...,16     N2 = 8) THEN FOR      
C                          N = N2+1     L1 = 1, L2 = 5, L3 = 2; FOR     
C                          N2+1 .LT. N .LE. 3*N2/2     L1 = 3.75,       
C                          L2=ALOG(N)/ALOG(2.) + 7, L3 = 4; FOR         
C                          3*(N2/2).LT.N.LE.2*N2     L1 = 2,            
C                            L2 = ALOG(N)/ALOG(2.)+5, L3 = 3  .         
C                        IF((ND.GT.2).AND.(N.EQ.N2+1))                  
C                           L1=3.75,L2=ALOG(N)/ALOG(2.), L3=4           
C                          IF(ND.NE.1)  L3=L3+1                         
C                                                                       
C                          (IT IS SUFFICIENT TO TAKE THE INTEGER PART   
C                           OF L1*N.)                                   
C                                                                       
C                                                                       
C                        LW                                             
C                          ACTUAL LENGTH OF LW IN CALLING ROUTINE       
C                        IERROR                                         
C                          NOT USED ON INPUT.                           
C                                                                       
C ON OUTPUT              ALL PARAMETERS ARE UNCHANGED EXCEPT            
C                                                                       
C                        Y                                              
C                          CONTAINS THE SOLUTION X(I,J), I = 1,2,...,M, 
C                                                        J = 1,2,...,N. 
C                                                                       
C                        LW                                             
C                          CONTAINS THE MINIMUM LENGTH OF W             
C                          NECESSARY FOR THESE VALUES OF M AND N        
C                        IERROR                                         
C                          AN ERROR FLAG WHICH INDICATES INVALID INPUT  
C                          PARAMETERS.  IF IERROR IS NOT ZERO, A        
C                          SOLUTION IS NOT ATTEMPTED.                   
C                               = 0  NO ERROR.                          
C                               = 1  M LESS THAN 2  .                   
C                               = 2  N LESS THAN 2  .                   
C                               = 3  IDIMY LESS THAN M.                 
C                               = 4  MPEROD IS NEITHER 0 NOR 1  .       
C                               = 5  LW IS TOO SMALL                    
C                               = 6  UNVALID  ND                        
C ENTRY POINTS           POISXX, PSXSUB, TRIX, TRIXP, TRI3, TRI3P,      
C                        TRIXS, TRISPS, TRI3S, TRI3PS                   
C                                                                       
C SPECIAL CONDITIONS     NONE                                           
C                                                                       
C COMMON BLOCKS          NONE                                           
C                                                                       
C I/O                    NONE                                           
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C ORIGINATOR             ULRICH SCHUMANN AND ROLAND SWEET               
C                                                                       
C SPACE REQUIRED         ABOUT 5000 DECIMAL LOCATIONS ON THE NCAR       
C                        CONTROL DATA 7600 AND ABOUT  31400 BYTES       
C                        ON THE  IBM 370/168 USING FORTRAN-H-EXTENDED.  
C                                                                       
C TIMING AND             THE EXECUTION TIME IS PROPORTIONAL TO MN LOG N.
C ACCURACY               HOWEVER, THE TIME IS UP TO 20 PERCENT SMALLER  
C                        FOR SOME VALUES OF N WHICH ARE OF THE FORM     
C                        2**L+1 THAN FOR VALUES OF N WHICH DIFFER FROM  
C                        THIS BY +1 OR -1  .  SOME TYPICAL VALUES ARE   
C                        TABULATED BELOW.  THE ACCURACY WITH RESPECT TO 
C                        ROUNDOFF ERRORS FOR M = N = 65 IS OF THE ORDER 
C                        10**(-12) ON THE CDC 7600 AND  10**(-4) ON THE 
C                        IBM 370 SYSTEM. IT IS ROUGHLY PROPORTIONAL M*M 
C                        THE ROUNDOFF ERROR SHOULD NOT BE CONFUSED WITH 
C                        THE DISCRETIZATION ERROR WHICH IS PROPORTIONAL 
C                        TO M**(-2)+N**(-2).                            
C                                                                       
C                                       EXECUTION TIME                  
C                                       (MSEC ON THE    (MSEC ON THE    
C                                       NCAR-CDC 7600)  IBM 370/165)    
C                       M = N                (FOR MPEROD=0)             
C                       -----           -------------   ------------    
C                                                                       
C                         32                 44             90          
C                         33                 35             87          
C                         34                 45            107          
C                         50                103            257          
C                         64                182            413          
C                         65                144            333          
C                         66                178            427          
C                         80                274            673          
C                        100                              1000          
C                                                                       
C PORTABILITY            AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN  
C                        WITH ONE MACHINE-DEPENDENT CONSTANT            
C                        (PI IN PSXSUB).                                
C                                                                       
C REQUIRED ROUTINES      COS                                            
C FURTHER DESCRIPTION    THE THEORY USED FOR THE PRESENT ROUTINE IS     
C                        DESCRIBED IN U. SCHUMANN AND R.A. SWEET, A     
C                        DIRECT METHOD FOR THE SOLUTION OF POISSON(S)   
C                        EQUATION WITH NEUMANN BOUNDARY CONDITIONS ON A 
C                        STAGGERED GRID OF ARBITRARY SIZE.              
C                        J. COMPUTATIONAL PHYSICS VOL.20 (1976) 171-182 
C                        THE EXTENSION  TO OTHER THAN  NEUMANN BOUNDARY 
C                        CONDITIONS IS GIVEN BY U.SCHUMANN IN           
C                        REPORT KFK 2645 (1979), ANHANG 5.              
C                                                                       
C RECOMMENDATIONS        IF ONLY SPECIFIC VALUES OF N OR MPEROD ARE     
C                        USED, THE ROUTINE CAN BE SHORTENED.  REFER TO  
C                        THE LISTING FOR DETAILS.                       
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C     THE FOLLOWING CARD IS NEEDED IF MPEROD .EQ. 0                     
      EXTERNAL TRIX,TRIXS,TRI3S,TRI3                                    
C     THE FOLLOWING CARD IS NEEDED IF MPEROD.NE.0                       
      EXTERNAL TRIXP,TRIXPS,TRI3PS,TRI3P                                
      DIMENSION W(1),B(1),A(1),C(1),Y(IDIMY,1)                          
      IF(M.LE.1) GO TO 401                                              
      IF(N.LE.1) GO TO 402                                              
      IF(M.GT.IDIMY) GO TO 403                                          
      IF((ND.LT.1).OR.(ND.GT.4)) GOTO 406                               
      IERROR=0                                                          
      DO 50 I=1,M                                                       
      A(I)=-A(I)                                                        
      B(I)=-B(I)+2.                                                     
   50 C(I)=-C(I)                                                        
      LN2=0                                                             
      N2=4                                                              
    1 IF(N2.GT.N) GO TO 2                                               
      N2=N2*2                                                           
      LN2=LN2+1                                                         
      GO TO 1                                                           
    2 N2=N2/2                                                           
      IF((N.EQ.(N2+1)).AND.(ND.LE.2)) GOTO 11                           
      IF(N.LE.3*(N2/2)) GO TO 12                                        
      L2=2*N-2                                                          
      GO TO 13                                                          
   11 L2=N-1                                                            
      GO TO 13                                                          
   12 L2=3.75*N-4                                                       
   13 IF (ND.NE.1 ) L2=L2+1                                             
      L3=L2+M                                                           
      L4=L3+M                                                           
      L5=L4+M                                                           
      L6=L5+M                                                           
      IF(((N.EQ.(N2+1)).AND.(ND.LE.2)).OR.(N.GT.3*(N2/2))) GOTO 21      
      L7=L6+M                                                           
      L8=L7+M                                                           
      GO TO 22                                                          
   21 L7=L6                                                             
      L8=L7                                                             
   22 L9=L8+M                                                           
      IF(MPEROD.EQ.1) GO TO 200                                         
      IF(MPEROD.NE.0) GO TO 404                                         
      L10=L9                                                            
      I=L10+LN2*M-1                                                     
      IF(LW.LT.I) GOTO 405                                              
      LW=I                                                              
      CALL PSXSUB(N,ND,M,IDIMY,A,B,C,Y,W,W(L2),W(L3),W(L4),W(L5),W(L6)  
     1,W(L7),W(L8),W(L9),W(L10),TRIX,TRIXS,TRI3,TRI3S,SING)             
      GO TO 300                                                         
  200 L10=L9+M                                                          
      I=L10+LN2*M-1                                                     
      IF(LW.LT.I) GOTO 405                                              
      LW=I                                                              
      CALL PSXSUB(N,ND,M,IDIMY,A,B,C,Y,W,W(L2),W(L3),W(L4),W(L5),W(L6)  
     1,W(L7),W(L8),W(L9),W(L10),TRIXP,TRIXPS,TRI3P,TRI3PS,SING)         
  300 DO 301 J=1,N                                                      
      DO 301 I=1,M                                                      
  301 Y(I,J)=-Y(I,J)                                                    
  400 CONTINUE                                                          
      DO 51 I=1,M                                                       
      A(I)=-A(I)                                                        
      B(I)=2.-B(I)                                                      
   51 C(I)=-C(I)                                                        
  500 CONTINUE                                                          
      RETURN                                                            
  401 IERROR=1                                                          
      GO TO 500                                                         
  402 IERROR=2                                                          
      GO TO 500                                                         
  403 IERROR=3                                                          
      GO TO 500                                                         
  404 IERROR=4                                                          
      GO TO 400                                                         
  405 IERROR=5                                                          
      LW=I                                                              
      GOTO  400                                                         
  406 IERROR=6                                                          
      GOTO 500                                                          
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE PSXSUB(N,ND,M,IDIMQ,BA,BB,BC,Q,TCOS,D,B,B2,B3,W,W2,W3,S
     1,P,TRIX,TRIXS,TRI3,TRI3S,SING)                                    
      DIMENSION Q(IDIMQ,1),TCOS(1),B(1),D(1),W(1),S(1)                  
     1,BA(1),BB(1),BC(1)                                                
     2,P(M,1)                                                           
      DIMENSION KKK(4) ,B2(1),B3(1),W2(1),W3(1)                         
      INTEGER H                                                         
      DATA PI/3.14159265358979/                                         
C                                                                       
C     PSXSUB SOLVES THE POISSON-EQUATION WITH NEUMANN- OR DIRICHLET-    
C     BOUNDARY-CONDITIONS IN FINITE-DIFFERENCE FORM                     
C     ON A STAGGERED GRID                                               
C                                                                       
C     NO RESTRICTION ON N AND M    EXCEPT N.GE.2 , M .GE.2              
C                                                                       
C                                                                       
C      N   NUMBER OF POINTS IN Y-DIRECTION                              
C     ND  =    1    2    3    4                                         
C             N-N  D-D  N-D  D-N     BOUNDARY-COND. AT J =  0-(N+I)     
C      M   NUMBER OF POINTS IN X-DIRECTION                              
C      BA,BB,BC ARE ARRAYS OF LENGTH M CONTAINING THE LOWER, MIDDLE     
C                                      AND UPPER DIAGONAL OF A          
C       IDIMQ FIRST DIMENSION OF Q IN CALLING ROUTINE                   
C      Q  INPUT - RIGHT HAND SIDE                                       
C        OUTPUT - SOLUTION                                              
C     TCOS   WORK FIELD OF LENGTH 2*N-3                                 
C       B,D,W ARE WORK FIELDS OF LENGTH M                               
C     P      WORK FIELD OF LENGTH FLOOR(LOG(N)/LOG(2)-EPS)  * M         
C       FLOOR = ROUND-OFF FUNCTION, EPS IS THE SMALLEST POSITVE NUMBER  
C                                   FOR WHICH (1+EPS).GT. 1             
C     S WORK-SPACE OF LENGTH M, NEEDED IF PERIODIC BOUNDARY COND. ONLY  
C      TRIX,TRIXS ARE EXTERNAL ROUTINES CONTAINING THE PROPER           
C            METHOD TO SOLVE TRIDIAGONAL SYSTEMS                        
C     SING   SINGULARITY INDICATOR   =0. FOR SINGULAR SYSTEMS           
C                                                                       
C                                                                       
C     SPECIAL CASE N=3                                                  
      IF(N.NE.3) GOTO 19                                                
      IF(ND.GT.2) GOTO 19                                               
C     ROOT OF A(0) IS 0.                                                
      DO 14 I=1,M                                                       
   14 B(I)=Q(I,2)                                                       
      CALL TRIX(1,0,M,BA,BB,BC,B,0.,D,W,S)                              
      DO 10 I=1,M                                                       
      B(I)=Q(I,1)+Q(I,3)+2.*B(I)                                        
   10 Q(I,2)=B(I)                                                       
      IF(ND.EQ.2) GOTO 15                                               
C     ROOTS OF M(0) ARE (-1.,2.), SINGULAR                              
      CALL TRIXS(1,0,M,BA,BB,BC,B,-1.,D,W,SING,S)                       
      TCOS(1)=1.                                                        
      GOTO 16                                                           
C    D-D                                                                
   15 TCOS(1)=2.*COS(PI/4.)                                             
      TCOS(2)=-TCOS(1)                                                  
      CALL TRIX(2,0,M,BA,BB,BC,B,TCOS,D,W,S)                            
      TCOS(1)=0.                                                        
   16 DO 11 I=1,M                                                       
      Q(I,2)=B(I)+.5*(Q(I,2)-Q(I,1)-Q(I,3))                             
   11 B(I)=Q(I,1)+Q(I,2)                                                
C     ROOT OF K(0) IS 1.                                                
      CALL TRIX(1,0,M,BA,BB,BC,B,TCOS,D,W,S)                            
      DO 12 I=1,M                                                       
      Q(I,1)=B(I)                                                       
   12 B(I)=Q(I,2)+Q(I,3)                                                
      CALL TRIX(1,0,M,BA,BB,BC,B,TCOS,D,W,S)                            
      DO 13 I=1,M                                                       
   13 Q(I,3)=B(I)                                                       
      GOTO 900                                                          
   19 CONTINUE                                                          
C     END SPECIAL N=3- CASE                                             
C                                                                       
C     SET INTEGERS FOR R=0                                              
C     NL=2**R                                                           
      NL=1                                                              
C      H=NL/2 =2**(R-1)                                                 
      H=0                                                               
C     NL2=2**(R+1)                                                      
      NL2=2                                                             
C     NR=N(R)                                                           
      NR=N                                                              
C     NA=LAST ROW AT STEP R                                             
      NA=N                                                              
C     DEGREE OF B(R) IS KR  , THAT OF C(R) IS LR                        
      KR=1                                                              
      LR=0                                                              
      NODD=MOD(N,2)                                                     
      NRJ=NA-2+NODD                                                     
      IRP=1                                                             
      KKK(2)=0                                                          
      KKK(3)=0                                                          
      KKK(4)=0                                                          
C     BOUNDARY CONDITION INDICATORS  NNL = 1  IF NEUMANN AT J=0         
C                                    NNR = 1  IF NEUMANN AT J=N+1       
C                                        = 0  OTHERWISE                 
      NNR=1                                                             
      IF((ND.EQ.2).OR.(ND.EQ.3)) NNR=0                                  
      NNL=1                                                             
      TCOS(1)=1.                                                        
      IF((ND.EQ.1).OR.(ND.EQ.3)) GOTO 22                                
      NNL=0                                                             
      TCOS(1)=0.                                                        
C      REDUCTION FOR R=0                                                
C     ROOTS OF K(0) ARE TCOS(1)=1.                                      
   22 DO 20 I=1,M                                                       
   20 B(I)=Q(I,1)                                                       
      CALL TRIX(1,0,M,BA,BB,BC,B,TCOS,D,W,S)                            
      DO 21 I=1,M                                                       
   21 Q(I,1)=Q(I,2)+B(I)                                                
      IF(NR.EQ.2) GOTO 901                                              
C     ROOTS OF A(0) ARE TCOS(1)=0.                                      
      IF(3.GT.NRJ) GOTO 40                                              
      KKK(1)=1                                                          
      J0=3                                                              
   70 J1=J0+4                                                           
      IF(J1.GT.NRJ) GOTO 71                                             
      J=J0+2                                                            
      DO 72  I=1,M                                                      
      B (I)=Q(I,J0)                                                     
      B2(I)=Q(I,J )                                                     
   72 B3(I)=Q(I,J1)                                                     
      CALL TRI3(M,BA,BB,BC,KKK,B,B2,B3, 0. ,D,W,W2,W3,S)                
      DO 73 I=1,M                                                       
      Q(I,J0)=Q(I,J0-1)+Q(I,J0+1)+2.*B (I)                              
      Q(I,J )=Q(I,J -1)+Q(I,J +1)+2.*B2(I)                              
   73 Q(I,J1)=Q(I,J1-1)+Q(I,J1+1)+2.*B3(I)                              
      J0=J1+2                                                           
      IF(J0-NRJ) 70,71,40                                               
C      70 TEST FOR ANOTHER TRIPLE, 71 ONE LEFT, 40 FINISHED             
   71 DO 30 J=J0,NRJ,2                                                  
      DO 31 I=1,M                                                       
   31 B(I)=Q(I,J)                                                       
      CALL TRIX(1,0,M,BA,BB,BC,B,0.,D,W,S)                              
      DO 32 I=1,M                                                       
   32 Q(I,J)=Q(I,J-1)+Q(I,J+1)+2.*B(I)                                  
   30 CONTINUE                                                          
   40 CONTINUE                                                          
      IF(NODD) 51,51,41                                                 
C     41 = ODD                                                          
C     ROOTS OF B(0) ARE TCOS(1)=1., OF C(0) ARE EMPTY                   
   41 TCOS(1)=1.                                                        
      IF(NNR.EQ.0) TCOS(1)=0.                                           
      DO 46 I=1,M                                                       
   46 B(I)=Q(I,NA)                                                      
      CALL TRIX(1,0,M,BA,BB,BC,B,TCOS,D,W,S)                            
      DO 47 I=1,M                                                       
   47 Q(I,NA)=Q(I,NA-1)+B(I)                                            
      NR=(NR+1)/2                                                       
      IF((NR.EQ.3).AND.(ND.LE.2)) GOTO 800                              
      KR=2                                                              
      NODOLD=1                                                          
      IRP=0                                                             
      GOTO 140                                                          
C      51 EVEN                                                          
C     ROOTS OF A(0) ARE TCOS(1)=0                                       
   51 DO 52 I=1,M                                                       
   52 B(I)=Q(I,NA-1)                                                    
      CALL TRIX(1,0,M,BA,BB,BC,B,0.,D,W,S)                              
      DO 53 I=1,M                                                       
      P(I,1)=B(I)                                                       
   53 B(I)=Q(I,NA)+B(I)                                                 
C     ROOTS OF B(0) = TCOS(1)=1., ROOTS OF A(0) ARE TCOS(2)=0,          
C                   THOSE OF C(0) ARE EMPTY                             
      IF(NNR.EQ.0) GOTO 59                                              
      TCOS(1)=1.                                                        
      TCOS(2)=0.                                                        
      CALL TRIX(1,1,M,BA,BB,BC,B,TCOS,D,W,S)                            
   59 DO 58 I=1,M                                                       
   58 Q(I,NA-1)=Q(I,NA-2)+B(I)+P(I,1)                                   
      LR=1                                                              
      KR=3                                                              
      NR=NR/2                                                           
      NA=NA-1                                                           
      NODOLD=0                                                          
      GOTO 140                                                          
C      END OF RDUCTION FOR R=0                                          
C     REDUCTION LOOP STARTS AT 200                                      
  200 CONTINUE                                                          
      NODD=MOD(NR,2)                                                    
C     NODD=1  IF NR IS ODD, OTHERWISE NODD=0                            
C     NRJ IS THE LIMIT FOR THE INNER REDUCTION LOOP                     
      NRJ=NA-NL2+NODD                                                   
C     J=1                                                               
      IF((1+NL).EQ.NA) GOTO 104                                         
      DO 101 I=1,M                                                      
  101 B(I)=Q(I,1)+.5*(Q(I,1+NL)-Q(I,1+H)-Q(I,1+NL+H))                   
      GOTO 105                                                          
  104 IF(NODOLD.NE.0) GOTO 107                                          
      DO 106 I=1,M                                                      
  106 B(I)=Q(I,1)+P(I,IR)                                               
      GOTO 105                                                          
  107 DO 108 I=1,M                                                      
  108 B(I)=Q(I,1)-Q(I,NA-H)+Q(I,NA)                                     
C     ROOTS OF K(R)                                                     
  105 PIX= PI/(NL2+2-NNL)                                               
      DO 102 I=1,NL                                                     
  102 TCOS(I)= 2.* COS((2*I-NNL)*PIX)                                   
      CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      DO 103 I=1,M                                                      
  103 Q(I,1)=Q(I,1+NL)+Q(I,1)-Q(I,1+H)+B(I)                             
      IF(NR.EQ.2) GO TO 901                                             
C     ROOTS OF  A(R)                                                    
      DO 109 I=1,NL                                                     
  109 TCOS(I)=2.*COS((2*I-1)*PI/NL2)                                    
      J0=1+NL2                                                          
      IF(J0.GT.NRJ) GO TO 120                                           
C     INNER REDUCTION LOOP                                              
      KKK(1)=NL                                                         
  150 J1=J0+2*NL2                                                       
      IF(J1.GT.NRJ) GOTO 159                                            
      J=J0+NL2                                                          
      DO 151 I=1,M                                                      
      B (I)=Q(I,J0-NL)-Q(I,J0-NL-H)-Q(I,J0-H)+2.*Q(I,J0)-Q(I,J0+H)      
     1+Q(I,J0+NL)-Q(I,J0+NL+H)                                          
      B2(I)=Q(I,J -NL)-Q(I,J -NL-H)-Q(I,J -H)+2.*Q(I,J )-Q(I,J +H)      
     1+Q(I,J +NL)-Q(I,J +NL+H)                                          
  151 B3(I)=Q(I,J1-NL)-Q(I,J1-NL-H)-Q(I,J1-H)+2.*Q(I,J1)-Q(I,J1+H)      
     1+Q(I,J1+NL)-Q(I,J1+NL+H)                                          
      CALL TRI3(M,BA,BB,BC,KKK,B,B2,B3,TCOS,D,W,W2,W3,S)                
      DO 152 I=1,M                                                      
      Q(I,J0)=Q(I,J0-NL)-Q(I,J0-H)+Q(I,J0)-Q(I,J0+H)+Q(I,J0+NL)+B (I)   
      Q(I,J )=Q(I,J -NL)-Q(I,J -H)+Q(I,J )-Q(I,J +H)+Q(I,J +NL)+B2(I)   
  152 Q(I,J1)=Q(I,J1-NL)-Q(I,J1-H)+Q(I,J1)-Q(I,J1+H)+Q(I,J1+NL)+B3(I)   
      J0=J1+NL2                                                         
      IF(J0-NRJ) 150,159,120                                            
C     150 TEST FOR ANOTHER TRIPLE, 159 ONE LEFT, 120 FINISHED           
  159 DO 110 J=J0,NRJ,NL2                                               
      DO 111 I=1,M                                                      
  111 B(I)=Q(I,J-NL)-Q(I,J-NL-H)-Q(I,J-H)+2.*Q(I,J)-Q(I,J+H)+Q(I,J+NL)  
     1-Q(I,J+NL+H)                                                      
      CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      DO 112 I=1,M                                                      
  112 Q(I,J)=Q(I,J-NL)-Q(I,J-H)+Q(I,J)-Q(I,J+H)+Q(I,J+NL)+B(I)          
  110 CONTINUE                                                          
  120 CONTINUE                                                          
C     REDUCTION FOR LAST ROW                                            
      IF(NODD) 131,131,121                                              
C     121=ODD                                                           
C     ROOTS OF B(R)                                                     
  121 PIX= PI/(2*KR+2-NNR)                                              
      DO  123 I=1,KR                                                    
  123 TCOS(I)=2.*COS((2*I-NNR)*PIX)                                     
      IF(LR.EQ.0) GO TO 125                                             
      PIX = PI/(2*LR+2-NNR)                                             
      DO 124 I=1,LR                                                     
  124 TCOS(KR+I) = 2.*COS((2*I-NNR)*PIX)                                
  125 DO 126 I=1,M                                                      
  126 B(I)=Q(I,NA)+.5*(Q(I,NA-NL)-Q(I,NA-NL-H)-Q(I,NA-H))               
      CALL TRIX(KR,LR,M,BA,BB,BC,B,TCOS,D,W,S)                          
      IF(NODOLD.EQ.1) GOTO 122                                          
      DO 127 I=1,M                                                      
  127 Q(I,NA)=Q(I,NA-NL)+P(I,IR)+B(I)                                   
      IRP=IR-1                                                          
      GOTO 128                                                          
  122 DO 129 I=1,M                                                      
  129 Q(I,NA)=Q(I,NA-NL)-Q(I,NA-H  )+Q(I,NA)+B(I)                       
      IRP=IR                                                            
  128 NR=(NR+1)/2                                                       
      IF((NR.EQ.3).AND.(LR.EQ.0).AND.(ND.LE.2)) GOTO 800                
C     LR=LR , NA=NA                                                     
      KR=KR+NL                                                          
      NODOLD=1                                                          
      GO TO 140                                                         
C     131=EVEN                                                          
C     TCOS CONTAINS STILL THE ROOTS OF A(R)                             
  131 IF(NODOLD.EQ.0) GOTO 139                                          
      DO 141 I=1,M                                                      
  141 B(I)=Q(I,NA-NL)+Q(I,NA)+.5*(Q(I,NA-NL2)-Q(I,NA-NL2-H)-Q(I,NA-NL2+H
     1)) -Q(I,NA-H)                                                     
      GOTO 142                                                          
  139 DO 132 I=1,M                                                      
  132 B(I)=Q(I,NA-NL)+P(I,IR)+.5*(Q(I,NA-NL2)-Q(I,NA-NL2-H)             
     1-Q(I,NA-NL2+H))                                                   
  142 CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      DO 133 I=1,M                                                      
      P(I,IRP)=B(I)+.5*(Q(I,NA-NL)-Q(I,NA-NL-H)-Q(I,NA-H))              
  133 B(I)=Q(I,NA)+P(I,IRP)                                             
C     ROOTS OF B(R)                                                     
      PIX = PI/(2*KR+2-NNR)                                             
      DO 134 I=1,KR                                                     
  134 TCOS(I)=2.*COS((2*I-NNR)*PIX)                                     
      KK=KR                                                             
C     ROOTS OF A(R) AND  C(R)  IN  MERGED SEQUENCE                      
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      J=(NL+1.)/(LR+1.)                                                 
      X =  PI/NL2                                                       
      J0=1                                                              
      PIX= PI/(2*LR+2-NNR)                                              
      DO 135 I=1,NL                                                     
      KK=KK+1                                                           
      TCOS(KK) = 2.*COS((2*I-1)*X)                                      
      IF(I.NE.J) GOTO 135                                               
      KK=KK+1                                                           
      TCOS(KK) = 2.*COS((2*J0-NNR)*PIX)                                 
      J0=J0+1                                                           
      J= (J0*(NL+1.))/(LR+1.)                                           
  135 CONTINUE                                                          
      CALL TRIX(KR,NL+LR,M,BA,BB,BC,B,TCOS,D,W,S)                       
      DO 138 I=1,M                                                      
  138 Q(I,NA-NL)=Q(I,NA-NL2)+B(I)+P(I,IRP)                              
      LR=KR                                                             
      KR=KR+NL2                                                         
      NR=NR/2                                                           
      NA=NA-NL                                                          
      NODOLD=0                                                          
C     END OF STANDARD REDUCTION                                         
  140 H=NL                                                              
      NL=NL2                                                            
      NL2=NL2*2                                                         
      IF(NR.EQ.3) GO TO 700                                             
      IR=IRP                                                            
      IRP=IRP+1                                                         
C     LOOP                                                              
      GO TO 200                                                         
C     SPECIAL FINAL SITUATION                                           
C     (NR.EQ.3).AND.(LR.EQ.0).AND.(N.GE.5).AND.(ND.LE.2)                
C     HAPPENS IF N.EQ. 2**(R-1)+1  ONLY                                 
  800 H=NL                                                              
      NL=NL2                                                            
      NL2=2*NL                                                          
      J=1+NL                                                            
C     ROOTS OF A(R+1)                                                   
      DO 819 I=1,NL                                                     
  819 TCOS(I)=2.*COS(((2*I-1)*PI)/NL2)                                  
      DO 820 I=1,M                                                      
  820 B(I)=Q(I,1)-Q(I,1+H)+Q(I,J)-Q(I,NA-H)+Q(I,NA)                     
      CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      DO 821 I=1,M                                                      
      B(I)=2.*B(I)+Q(I,1)-Q(I,1+H)+Q(I,J)-Q(I,NA-H)+Q(I,NA)             
  821 Q(I,J)=B(I)                                                       
C     ROOTS OF M(R)                                                     
      IF(ND.EQ.2) GOTO 810                                              
      DO 802 I=1,NL                                                     
  802 TCOS(I)=2.*COS(((2*I)*PI)/(NL2+1))                                
      J0=NL-1                                                           
      DO 803 I=1,J0                                                     
  803 TCOS(I+NL)=2.*COS((I*PI)/NL)                                      
      CALL TRIXS(NL+J0,0,M,BA,BB,BC,B,TCOS,D,W,SING,S)                  
      GOTO 811                                                          
  810 J0=NL+1                                                           
      PIX =PI/(2*J0)                                                    
      DO 812 I=1,J0                                                     
  812 TCOS(I)=2.*COS((2*I-1)*PIX)                                       
      KK= NL-1                                                          
      PIX = PI/NL                                                       
      DO 813 I=1,KK                                                     
  813 TCOS(J0+I)=2.*COS(I*PIX)                                          
      CALL TRIX(NL2,0,M,BA,BB,BC,B,TCOS,D,W,S)                          
  811 DO 806 I=1,M                                                      
      Q(I,J   )=B(I)+0.5*(Q(I,J   )-Q(I,1)-Q(I,NA))                     
  806 B(I)=Q(I,1)+Q(I,J)                                                
C     ROOTS OF K(R)                                                     
      PIX=PI/(NL2+2-NNL)                                                
      DO 807 I=1,NL                                                     
  807 TCOS(I)=2.*COS((2*I-NNL)*PIX)                                     
      CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      DO 808 I=1,M                                                      
      Q(I,1)=B(I)+Q(I,1)-Q(I,1+H)                                       
  808 B(I)=Q(I,NA)+Q(I,J)                                               
      CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      DO 809 I=1,M                                                      
  809 Q(I,NA)=B(I)+Q(I,NA)-Q(I,NA-H)                                    
C     END OF CASE N.EQ.2**(R-1)+1                                       
      GOTO 300                                                          
  700 CONTINUE                                                          
C     SPECIAL FINAL SITUATION (NR.EQ.3) .AND. (LR.GT.0).OR.(ND.GT.2)    
C                                                                       
C     ROOTS OF  D   IN  MERGED SEQUENCE                                 
      J1=KR+NL+1-NNL*NNR                                                
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      J=  (J1+1.)/FLOAT(NL)                                             
      J0=1                                                              
      X=PI/FLOAT(NL)                                                    
      PIX = PI/(2*KR+NL2+4-NNL-NNR)                                     
      KK = 0                                                            
      DO  701 I=1,J1                                                    
      KK=KK+1                                                           
      TCOS(KK) = 2.*COS((2*I-IABS(NNL-NNR))*PIX)                        
      IF (I.NE.J) GOTO 701                                              
      KK=KK+1                                                           
      TCOS(KK) = 2.*COS( J0*X)                                          
      J0=J0+1                                                           
      J=  (J0*(J1+1.))/FLOAT(NL)                                        
  701 CONTINUE                                                          
      KKK(1)=KK                                                         
      KKK(2)=NL+LR                                                      
      KKK(3)=KR                                                         
      KKK(4)=KR+NL                                                      
C     ROOTS OF K(R) AND C(R) IN MERGED SEQUENCE                         
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      J=(NL+1.)/(LR+1.)                                                 
      X=PI/(2*LR+2-NNR)                                                 
      PIX=PI/(NL2+2-NNL)                                                
      J0=1                                                              
      DO 702 I=1,NL                                                     
      KK=KK+1                                                           
      TCOS(KK)=2.*COS((2*I-NNL)*PIX)                                    
      IF(I.NE.J) GOTO 702                                               
      KK=KK+1                                                           
      TCOS(KK)=2.*COS((2*J0-NNR)*X)                                     
      J0=J0+1                                                           
      J=(J0*(NL+1.))/(LR+1.)                                            
  702 CONTINUE                                                          
C     ROOTS OF B(R)                                                     
      PIX=PI/(2*KR+2-NNR)                                               
      DO 703 I=1,KR                                                     
      KK=KK+1                                                           
  703 TCOS(KK)=2.*COS((2*I-NNR)*PIX)                                    
C     ROOTS OF B(R) AND K(R) IN MERGED SEQUENCE                         
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      J=(KR+1.)/(NL+1.)                                                 
      X=PI/(NL2+2-NNL)                                                  
      J0=1                                                              
      DO 704 I=1,KR                                                     
      KK=KK+1                                                           
      TCOS(KK)=2.*COS((2*I-NNR)*PIX)                                    
      IF(I.NE.J) GOTO 704                                               
      KK=KK+1                                                           
      TCOS(KK)=2.*COS((2*J0-NNL)*X)                                     
      J0=J0+1                                                           
      J=(J0*(KR+1.))/(NL+1.)                                            
  704 CONTINUE                                                          
C                                                                       
C      K1=KKK(1) , K2=KKK(1)+KKK(2), K3=KKK(1)+KKK(2)+KKK(3),           
C     K4=KKK(1)+KKK(2)+KKK(3)+KKK(4)                                    
C                                                                       
C      TCOS CONTAINS                                                    
C     PART 1                                                            
C                                      KR+2*NL-1                        
C      FROM 1 TO K1                    ROOTS OF D                       
C     PART 2                                                            
C                                      NL+LR                            
C     FROM K1+1 TO K2   ROOTS OF K(R) AND C(R) MERGED                   
C                                      KR                               
C     PART 3                                                            
C      FROM K2+1 TO K3                 ROOTS OF B(R)                    
C     PART 4                                                            
C                                      NL+KR                            
C     FROM K3=1  TO   K4               ROOTS OF  B(R)  AND  K(R)        
C                                      IN  A  MERGED SEQUENCE           
C                                                                       
      J=1+NL                                                            
      IF(NODOLD.EQ.0) GOTO 707                                          
      DO 706 I=1,M                                                      
      X=.5*(Q(I,J)-Q(I,J-H)-Q(I,J+H))                                   
      B(I)=X+Q(I,NA)                                                    
      B2(I)=X+Q(I,1)                                                    
  706 B3(I)=Q(I,1)-Q(I,1+H)+Q(I,J)-Q(I,NA-H)+Q(I,NA)                    
      GO TO 709                                                         
  707 DO 708 I=1,M                                                      
      X=.5*(Q(I,J)-Q(I,J-H)-Q(I,J+H))                                   
      B(I)=X + Q(I,NA)                                                  
      B2(I)=X+Q(I,1)                                                    
  708 B3(I)=Q(I,1)-Q(I,1+H)+Q(I,J)+P(I,IRP)                             
  709 IF(ND.NE.1) GOTO 722                                              
      CALL TRI3S(M,BA,BB,BC,KKK,B,B2,B3,TCOS,D,W,W2,W3,SING,S)          
      GOTO 724                                                          
  722 CONTINUE                                                          
      CALL TRI3  (M,BA,BB,BC,KKK,B,B2,B3,TCOS,D,W,W2,W3,S)              
      DO 723 I=1,M                                                      
  723 B(I)=B(I)+B2(I)+B3(I)                                             
  724 DO 710 I=1,M                                                      
      X=B(I)+.5*(Q(I,J)-Q(I,J-H)-Q(I,J+H))                              
      Q(I,J)=X                                                          
      B(I)=Q(I,1)+X                                                     
  710 B2(I)=Q(I,NA)+X                                                   
C     ROOTS OF K(R)                                                     
      PIX=PI/(NL2+2-NNL)                                                
      DO 705 I=1,NL                                                     
  705 TCOS(I)=2.*COS((2*I-NNL)*PIX)                                     
      CALL TRIX(NL,0,M,BA,BB,BC,B ,TCOS,D,W,S)                          
C     COPY ROOTS OF B(R) TO THE BEGINNING OF TCOS                       
      J=KKK(1)+KKK(2)                                                   
      DO 713 I=1,KR                                                     
      J=J+1                                                             
  713 TCOS(I)=TCOS(J)                                                   
C     ROOTS OF C(R)                                                     
      IF(LR.EQ.0) GOTO 712                                              
      PIX=PI/(2*LR+2-NNR)                                               
      KK=KR                                                             
      DO 711 I=1,LR                                                     
      KK=KK+1                                                           
  711 TCOS(KK)=2.*COS((2*I-NNR)*PIX)                                    
  712 CALL TRIX (KR,LR,M,BA,BB,BC,B2,TCOS,D,W,S)                        
      KKK(2)=0                                                          
      KKK(3)=0                                                          
      KKK(4)=0                                                          
      IF(NODOLD.EQ.0) GO TO 715                                         
      DO 714 I=1,M                                                      
      Q(I,1)=B(I) + Q(I,1)-Q(I,1+H)                                     
  714 Q(I,NA)=B2(I)-Q(I,NA-H)+Q(I,NA)                                   
      GOTO 303                                                          
  715 DO 716 I=1,M                                                      
      Q(I,1)=B(I)+Q(I,1)-Q(I,1+H)                                       
  716 Q(I,NA)=B2(I)+P(I,IRP)                                            
      IRP=IRP-1                                                         
      GOTO 303                                                          
C                                                                       
C     END OF SPECIAL CASE ((NR.EQ.3).AND.(LR.GT.0))                     
C                                                                       
  901 CONTINUE                                                          
C     STANDARD-END-CASE, NR.EQ.2                                        
C     SOLVE FOR X(1)                                                    
      DO 902 I=1,M                                                      
  902 B(I)=Q(I,1)                                                       
C     ROOTS OF  E , C(R+1)                                              
      KK=KR+1-NNL*NNR                                                   
      PIX= PI/(2*KR+4-NNL-NNR)                                          
      DO 903 I=1,KK                                                     
  903 TCOS(I)= 2.*COS( (2*I-IABS(NNL-NNR))*PIX)                         
      J0=NL-1                                                           
      IF(J0.EQ.0) GO TO 905                                             
      PIX =PI/NL                                                        
      DO 904 I=1,J0                                                     
  904 TCOS (KK+I)=2.*COS(I*PIX)                                         
  905 KK=KK+J0                                                          
      IF(LR.EQ.0) GOTO 907                                              
      PIX = PI/(2*LR+2-NNR)                                             
      DO 906 I=1,LR                                                     
      KK=KK+1                                                           
  906 TCOS (KK) = 2.*COS((2*I-NNR)*PIX)                                 
      KK=KK-LR                                                          
  907 IF(ND.EQ.1) CALL TRIXS(KK,LR,M,BA,BB,BC,B,TCOS,D,W,SING,S)        
      IF(ND.NE.1) CALL TRIX (KK,LR,M,BA,BB,BC,B,TCOS,D,W,S)             
      DO 908 I=1,M                                                      
  908 Q(I,1)=B(I)+Q(I,1)-Q(I,1+NL)                                      
      IRP=IRP-1                                                         
C     NA=1+NL                                                           
      GOTO 400                                                          
C     END OF STANDARD-END-CASE                                          
  300 IF(NL.EQ.1) GO TO 900                                             
  303 NL2=NL                                                            
      NL=H                                                              
      H=H/2                                                             
      J=NA+NL                                                           
      IF(J.GT.N)GO TO 301                                               
C     EVEN                                                              
      NA=J                                                              
      KR=KR-NL2                                                         
      LR=KR-NL                                                          
C     SOLVE FOR LAST ROW                                                
C     ROOTS OF B                                                        
  400 PIX = PI/(2*KR+2-NNR)                                             
      DO 401 I = 1,KR                                                   
  401 TCOS(I)=2.*COS((2*I-NNR)*PIX)                                     
C     ROOTS OF C                                                        
      IF(LR.EQ.0) GO TO 403                                             
      PIX =PI /(2*LR+2-NNR)                                             
      DO 402 I=1,LR                                                     
  402 TCOS(KR+I)=2.*COS((2*I-NNR)*PIX)                                  
  403 DO 404 I=1,M                                                      
  404 B(I)=Q(I,NA)+Q(I,NA-NL)                                           
      CALL TRIX(KR,LR,M,BA,BB,BC,B,TCOS,D,W,S)                          
      IF(NL.EQ.1) GOTO 406                                              
      IF((NA+H).GT.N) GOTO 409                                          
      DO 405 I=1,M                                                      
  405 Q(I,NA)=B(I)+P(I,IRP)                                             
      IRP=IRP-1                                                         
      GOTO 407                                                          
  406 DO 408 I=1,M                                                      
  408 Q(I,NA)=B(I)                                                      
      GOTO 407                                                          
C     R .GT. 0, ODD NR                                                  
  409 DO 410 I=1,M                                                      
  410 Q(I,NA)=B(I)+Q(I,NA)-Q(I,NA-H)                                    
C      END OF ELIMINATION FOR LAST ROW                                  
  407 NRJ=NA-NL2                                                        
      GO TO 302                                                         
  301 KR=KR-NL                                                          
      LR=KR-NL                                                          
      NRJ=NA-NL                                                         
  302 J0=1+NL                                                           
      IF(J0.GT.NRJ) GO TO 300                                           
C     ROOTS OF A(R)                                                     
      DO 309 I=1,NL                                                     
  309 TCOS(I)=2.*COS((2*I-1)*PI/NL2)                                    
C     ELIMINATION OF INNER VARIABLES                                    
      J0=1+NL                                                           
      NRJ=NA-NL                                                         
      KKK(1)=NL                                                         
  340 J1=J0+2*NL2                                                       
      IF(J1.GT.NRJ) GOTO 350                                            
      J=J0+NL2                                                          
      DO 320 I=1,M                                                      
      B(I)=Q(I,J0)+Q(I,J0-NL)+Q(I,J0+NL)                                
      B2(I)= Q(I,J )+Q(I,J -NL)+Q(I,J +NL)                              
  320 B3(I)= Q(I,J1)+Q(I,J1-NL)+Q(I,J1+NL)                              
      CALL TRI3(M,BA,BB,BC,KKK,B,B2,B3,TCOS,D,W,W2,W3,S)                
      IF(NL.EQ.1) GOTO 323                                              
      DO 321 I=1,M                                                      
      Q(I,J0)=B (I)+0.5*(Q(I,J0)-Q(I,J0-H)-Q(I,J0+H))                   
      Q(I,J )=B2(I)+0.5*(Q(I,J )-Q(I,J -H)-Q(I,J +H))                   
  321 Q(I,J1)=B3(I)+0.5*(Q(I,J1)-Q(I,J1-H)-Q(I,J1+H))                   
      GOTO 324                                                          
  323 DO 326 I=1,M                                                      
      Q(I,J0)=B (I)                                                     
      Q(I,J )=B2(I)                                                     
  326 Q(I,J1)=B3(I)                                                     
  324 J0=J1+NL2                                                         
      IF(J0-NRJ) 340,350,300                                            
C     340 TEST FOR ANOTHER TRIPLE, 350 ONE LEFT, 300 FINISHED           
  350 DO 312 J=J0,NRJ,NL2                                               
      DO 311 I=1,M                                                      
  311 B(I)=Q(I,J)+Q(I,J-NL)+Q(I,J+NL)                                   
      CALL TRIX(NL,0,M,BA,BB,BC,B,TCOS,D,W,S)                           
      IF(NL.EQ.1) GOTO 313                                              
      DO 310 I=1,M                                                      
  310 Q(I,J)=B(I)+0.5*(Q(I,J)-Q(I,J-H)-Q(I,J+H))                        
      GOTO 314                                                          
  313 DO 316 I=1,M                                                      
  316 Q(I,J)=B(I)                                                       
  314 CONTINUE                                                          
  312 CONTINUE                                                          
C     END OF INNER ELIMINATION                                          
C     END OF ELIMINATION AT STEP R                                      
  315 GO TO 300                                                         
C     FINISHED                                                          
  900 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRIX(IDEGBR,IDEGCR,M,A,B,C,Y,TCOS,D,W,S)               
      DIMENSION A(1),B(1),C(1),Y(1),TCOS(1),D(1),W(1),S(1)              
C                                                                       
C     SUBROUTINE TO SOLVE TRIDIAGONAL SYSTEMS                           
C                                                                       
C      SOLVES MM.X=NN.Y                                                 
C      WHERE Y IS THE INPUT VECTOR OF LENGTH M                          
C            X  IS THE SOLUTION WHICH IS RETURNED IN  Y                 
C      THE MATRICES MM AND NN ARE GIVEN AS PRODUCTS                     
C      MM= PROD(K=1,IDEGBR) OF (A-LAMBDA(K) I)                          
C      NN = PROD (L=1,IDEGCR) OF (A-MU(L)   I)                          
C      THE MATRIX A IS  TRI-DIAGONAL THE LOWER,MAIN, AND UPPER          
C      DIAGONALS ARE GIVEN IN A,B,C                                     
C      THE ROOTS LAMBDA AND MU ARE STORED IN THE FIELD TCOS             
C      LAMBDA(K)=TCOS(K), K=1 TO IDEGBR                                 
C      MU(L) =TCOS(IDEGBR+L) FOR L=1 TO IDEGCR                          
C      0 .LE. IDEGCR .LT. IDEGBR                                        
C       M .GT. 1                                                        
C      D,W ARE WORKING FIELDS OF LENGTH M                               
      MM1 = M - 1                                                       
C ********** THE NEXT 3 CARDS WERE CHANGED ON NOV. 9, 1983 ************ 
      FB=    IDEGBR+1                                                   
      FC=    IDEGCR+1                                                   
      L=  FB/FC                                                         
      LINT=1                                                            
      DO 300 K=1,IDEGBR                                                 
      X =  TCOS(K)                                                      
      IF(K.NE.L) GOTO 90                                                
      XX=X-TCOS(IDEGBR+LINT)                                            
      DO 20 I=1,M                                                       
      W(I) = Y(I)                                                       
   20 Y(I) = XX*Y(I)                                                    
   90 CONTINUE                                                          
      Z=1./(B(1)-X)                                                     
      D(1)=C(1)*Z                                                       
      Y(1)=Y(1)*Z                                                       
      IM1=1                                                             
      DO 100 I=2,M                                                      
      Z=1./(B(I)-X-A(I)*D(IM1))                                         
      D(I)=C(I)*Z                                                       
      Y(I)=(Y(I)-A(I)*Y(IM1))*Z                                         
  100 IM1=I                                                             
      DO 200 IP=1,MM1                                                   
      I=M-IP                                                            
  200 Y(I)=Y(I)-D(I)*Y(I+1)                                             
      IF(K.NE.L) GOTO 300                                               
      LINT=LINT+1                                                       
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L= (FLOAT(LINT)*FB)/FC                                            
      DO 250 I=1,M                                                      
  250 Y(I) = Y(I) + W(I)                                                
  300 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRIXP(KR,LR,M,A,B,C,Y,TWOCOS,D,W,S)                    
      DIMENSION A(1),B(1),C(1),Y(1),TWOCOS(1),D(1),W(1) ,S(1)           
C     SUBROUTINE TO SOLVE PERIODIC TRIDIAGONAL SYSTEM                   
      MM1=M-1                                                           
C ********** THE NEXT 3 CARDS WERE CHANGED ON NOV. 9, 1983 ************ 
      FB=    KR+1                                                       
      FC=    LR+1                                                       
      L=  FB/FC                                                         
      LINT=1                                                            
      DO 300 K=1,KR                                                     
      X =-TWOCOS(K)                                                     
      IF(K.NE.L) GOTO 90                                                
      XX = -X-TWOCOS(KR+LINT)                                           
      DO 20 I=1,M                                                       
      S(I)=Y(I)                                                         
   20 Y(I)=XX*Y(I)                                                      
   90 CONTINUE                                                          
      DEN=1./(B(1)+X)                                                   
      D(1)=C(1)*DEN                                                     
      W(1)=A(1)*DEN                                                     
      Y(1)=Y(1)*DEN                                                     
      BM = B(M)                                                         
      Z = C(M)                                                          
      IF(MM1.EQ.1) GOTO 101                                             
      DO 100 I=2,MM1                                                    
      DEN=1./(B(I)+X-A(I)*D(I-1))                                       
      D(I)=C(I)*DEN                                                     
      W(I)=-A(I)*W(I-1)*DEN                                             
      Y(I)=(Y(I)-A(I)*Y(I-1))*DEN                                       
      Y(M) = Y(M) - Z*Y(I-1)                                            
      BM = BM - Z*W(I-1)                                                
  100 Z = -Z*D(I-1)                                                     
  101 D(MM1) = D(MM1) + W(MM1)                                          
      Z = A(M) + Z                                                      
      DEN = BM   + X - Z*D(MM1)                                         
      Y(M) = Y(M) - Z*Y(M-1)                                            
      Y(M) = Y(M)/DEN                                                   
      Y(MM1) = Y(MM1) - D(MM1)*Y(M)                                     
      IF(MM1.EQ.1) GOTO 201                                             
      DO 200 IP=2,MM1                                                   
      I = M - IP                                                        
  200 Y(I) = Y(I) - D(I)*Y(I+1) - W(I)*Y(M)                             
  201 IF(K.NE.L) GOTO 300                                               
      LINT=LINT+1                                                       
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L = (FLOAT(LINT)*FB)/FC                                           
      DO 250 I=1,M                                                      
  250 Y(I)=Y(I)+S(I)                                                    
  300 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRI3(M,A,B,C,K,Y1,Y2,Y3,TCOS,D,W1,W2,W3,S)             
      DIMENSION A(1),B(1),C(1),K(4),TCOS(1),Y1(1),Y2(1),Y3(1)           
     1         ,D(1),W1(1),W2(1),W3(1),S(1)                             
C                                                                       
C     SUBROUTINE TO SOLVE TRIDIAGONAL SYSTEMS                           
C                                                                       
      MM1 = M - 1                                                       
      K1=K(1)                                                           
      K2=K(2)                                                           
      K3=K(3)                                                           
      K4=K(4)                                                           
      K2K3K4=K2+K3+K4                                                   
      IF(K2K3K4.EQ.0) GOTO 299                                          
C ********** THE NEXT 7 CARDS WERE CHANGED ON NOV. 9, 1983 ************ 
      F1=K1+1                                                           
      F2=K2+1                                                           
      F3=K3+1                                                           
      F4=K4+1                                                           
      L1 = F1/F2                                                        
      L2 = F1/F3                                                        
      L3 = F1/F4                                                        
      LINT1=1                                                           
      LINT2=1                                                           
      LINT3=1                                                           
      KINT1=K1                                                          
      KINT2=KINT1+K2                                                    
      KINT3=KINT2+K3                                                    
  299 CONTINUE                                                          
      DO 300 N=1,K1                                                     
      X=TCOS(N)                                                         
      IF(K2K3K4.EQ.0) GOTO 90                                           
      IF (N .NE. L1) GO TO 30                                           
      DO 20 I=1,M                                                       
   20 W1(I)=Y1(I)                                                       
   30 IF (N .NE. L2)  GO TO 50                                          
      DO 40 I=1,M                                                       
   40 W2(I)=Y2(I)                                                       
   50 IF (N .NE. L3)  GO TO 90                                          
      DO 60 I=1,M                                                       
   60 W3(I)=Y3(I)                                                       
   90 CONTINUE                                                          
      Z=1./(B(1)-X)                                                     
      D(1)=C(1)*Z                                                       
      Y1(1)=Y1(1)*Z                                                     
      Y2(1)=Y2(1)*Z                                                     
      Y3(1)=Y3(1)*Z                                                     
      IM1=1                                                             
      DO 100 I=2,M                                                      
      Z=1./(B(I)-X-A(I)*D(IM1))                                         
      D(I)=C(I)*Z                                                       
      Y1(I)=(Y1(I)-A(I)*Y1(IM1))*Z                                      
      Y2(I)=(Y2(I)-A(I)*Y2(IM1))*Z                                      
      Y3(I)=(Y3(I)-A(I)*Y3(IM1))*Z                                      
  100 IM1=I                                                             
      IP1=M                                                             
      DO 200 IP=1,MM1                                                   
      I = M - IP                                                        
      Y1(I)=Y1(I)-D(I)*Y1(IP1)                                          
      Y2(I)=Y2(I)-D(I)*Y2(IP1)                                          
      Y3(I)=Y3(I)-D(I)*Y3(IP1)                                          
  200 IP1=I                                                             
      IF(K2K3K4.EQ.0) GOTO 300                                          
      IF (N .NE. L1)  GO TO 230                                         
      XX=X-TCOS(LINT1+KINT1)                                            
      DO 220 I=1,M                                                      
  220 Y1(I)=XX*Y1(I)+W1(I)                                              
      LINT1 = LINT1 + 1                                                 
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L1 = (FLOAT(LINT1)*F1)/F2                                         
  230 IF (N .NE. L2) GO TO 250                                          
      XX=X-TCOS(LINT2+KINT2)                                            
      DO 240 I=1,M                                                      
  240 Y2(I)=XX*Y2(I)+W2(I)                                              
      LINT2 = LINT2 + 1                                                 
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L2 = (FLOAT(LINT2)*F1)/F3                                         
  250 IF (N .NE. L3)  GO TO 300                                         
      XX=X-TCOS(LINT3+KINT3)                                            
      DO 260 I=1,M                                                      
  260 Y3(I)=XX*Y3(I)+W3(I)                                              
      LINT3 = LINT3 + 1                                                 
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L3 = (FLOAT(LINT3)*F1)/F4                                         
  300 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRI3P(M,A,B,C,K,Y1,Y2,Y3,TCOS,D,W1,W2,W3,S)            
      DIMENSION A(1),B(1),C(1),K(4),TCOS(1),Y1(1),Y2(1),Y3(1)           
     1         ,D(1),W1(1),W2(1),W3(1),S(1)                             
C                                                                       
C     SUBROUTINE TO SOLVE TRIDIAGONAL SYSTEMS                           
C                                                                       
      MM1 = M - 1                                                       
      K1=K(1)                                                           
      K2=K(2)                                                           
      K3=K(3)                                                           
      K4=K(4)                                                           
      K2K3K4=K2+K3+K4                                                   
      IF(K2K3K4.EQ.0) GOTO 299                                          
C ********** THE NEXT 7 CARDS WERE CHANGED ON NOV. 9, 1983 ************ 
      F1=K1+1                                                           
      F2=K2+1                                                           
      F3=K3+1                                                           
      F4=K4+1                                                           
      L1=F1/F2                                                          
      L2=F1/F3                                                          
      L3=F1/F4                                                          
      LINT1=1                                                           
      LINT2=1                                                           
      LINT3=1                                                           
      KINT1=K1                                                          
      KINT2=KINT1+K2                                                    
      KINT3=KINT2+K3                                                    
  299 CONTINUE                                                          
      DO 300 N=1,K1                                                     
      X=TCOS(N)                                                         
      IF(K2K3K4.EQ.0) GOTO 90                                           
      IF (N .NE. L1) GO TO 30                                           
      DO 20 I=1,M                                                       
   20 W1(I)=Y1(I)                                                       
   30 IF (N .NE. L2)  GO TO 50                                          
      DO 40 I=1,M                                                       
   40 W2(I)=Y2(I)                                                       
   50 IF (N .NE. L3)  GO TO 90                                          
      DO 60 I=1,M                                                       
   60 W3(I)=Y3(I)                                                       
   90 CONTINUE                                                          
      Z=1./(B(1)-X)                                                     
      D(1)=C(1)*Z                                                       
      S(1)=A(1)*Z                                                       
      Y1(1)=Y1(1)*Z                                                     
      Y2(1)=Y2(1)*Z                                                     
      Y3(1)=Y3(1)*Z                                                     
      BM=B(M)                                                           
      CM=C(M)                                                           
      IF(MM1.EQ.1) GOTO 101                                             
      IM1=1                                                             
      DO 100 I=2,MM1                                                    
      Z=1./(B(I)-X-A(I)*D(IM1))                                         
      D(I)=C(I)*Z                                                       
      S(I)=-A(I)*S(IM1)*Z                                               
      Y1(I)=(Y1(I)-A(I)*Y1(IM1))*Z                                      
      Y2(I)=(Y2(I)-A(I)*Y2(IM1))*Z                                      
      Y3(I)=(Y3(I)-A(I)*Y3(IM1))*Z                                      
      Y1(M)=Y1(M)-CM*Y1(IM1)                                            
      Y2(M)=Y2(M)-CM*Y2(IM1)                                            
      Y3(M)=Y3(M)-CM*Y3(IM1)                                            
      BM=BM-CM*S(IM1)                                                   
      CM=-CM*D(IM1)                                                     
  100 IM1=I                                                             
  101 D(MM1)=D(MM1)+S(MM1)                                              
      CM=A(M)+CM                                                        
      Z=1./(BM-X-CM* D(MM1))                                            
      Y1(M)=(Y1(M)-CM*Y1(MM1))*Z                                        
      Y2(M)=(Y2(M)-CM*Y2(MM1))*Z                                        
      Y3(M)=(Y3(M)-CM*Y3(MM1))*Z                                        
      Y1(MM1)=Y1(MM1)-D(MM1)*Y1(M)                                      
      Y2(MM1)=Y2(MM1)-D(MM1)*Y2(M)                                      
      Y3(MM1)=Y3(MM1)-D(MM1)*Y3(M)                                      
      IF(MM1.EQ.1) GOTO 201                                             
      IP1=MM1                                                           
      DO 200 IP=2,MM1                                                   
      I = M - IP                                                        
      Y1(I)=Y1(I)-D(I)*Y1(IP1)-S(I)*Y1(M)                               
      Y2(I)=Y2(I)-D(I)*Y2(IP1)-S(I)*Y2(M)                               
      Y3(I)=Y3(I)-D(I)*Y3(IP1)-S(I)*Y3(M)                               
  200 IP1=I                                                             
  201 IF(K2K3K4.EQ.0) GOTO 300                                          
      IF(N. NE. L1) GO TO 230                                           
      XX=X-TCOS(LINT1+KINT1)                                            
      DO 220 I=1,M                                                      
  220 Y1(I)=XX*Y1(I)+W1(I)                                              
      LINT1 = LINT1 + 1                                                 
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L1 = (FLOAT(LINT1)*F1)/F2                                         
  230 IF (N .NE. L2) GO TO 250                                          
      XX=X-TCOS(LINT2+KINT2)                                            
      DO 240 I=1,M                                                      
  240 Y2(I)=XX*Y2(I)+W2(I)                                              
      LINT2 = LINT2 + 1                                                 
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L2 = (FLOAT(LINT2)*F1)/F3                                         
  250 IF (N .NE. L3)  GO TO 300                                         
      XX=X-TCOS(LINT3+KINT3)                                            
      DO 260 I=1,M                                                      
  260 Y3(I)=XX*Y3(I)+W3(I)                                              
      LINT3 = LINT3 + 1                                                 
C ********** THE NEXT CARD WAS CHANGED ON NOV. 9, 1983 ************     
      L3 = (FLOAT(LINT3)*F1)/F4                                         
  300 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRIXS(IDEGBR,IDEGCR,M,A,B,C,Y,TCOS,D,W,SING,S)         
      DIMENSION A(1),B(1),C(1),Y(1),  TCOS(1),D(1),W(1),S(1)            
C                                                                       
C     SUBROUTINE TO SOLVE TRIDIAGONAL SYSTEMS                           
C                                                                       
C      SAME AS TRIX EXCEPT FOR MM WHICH IS REPLACED BY                  
C               MM.(A-2 I)                                              
C      THIS IS A SINGULAR MATRIX WHICH NEEDS SPECIAL HANDLING           
      CALL TRIX       (IDEGBR,IDEGCR,M,A,B,C,Y,TCOS,D,W,S)              
      MM1 = M - 1                                                       
      Z=1./(B(1)-2.)                                                    
      D(1)=C(1)*Z                                                       
      Y(1)=Y(1)*Z                                                       
      IF(MM1.LT.2) GOTO 501                                             
      IM1=1                                                             
      DO 500 I=2,MM1                                                    
      Z=1./(B(I)-2.-A(I)*D(IM1))                                        
      D(I)=C(I)*Z                                                       
      Y(I)=(Y(I)-A(I)*Y(IM1))*Z                                         
  500 IM1=I                                                             
  501 Z=B(M)-2.-A(M)*D(MM1)                                             
      X=Y(M)-A(M)*Y(MM1)                                                
      IF(SING*Z.NE.0.) GOTO 550                                         
      Y(M)=0.                                                           
      GOTO 575                                                          
  550 Y(M)=X/Z                                                          
  575 IP1=M                                                             
      DO 600 IP=1,MM1                                                   
      I=M-IP                                                            
      Y(I)=Y(I)-D(I)*Y(IP1)                                             
  600 IP1=I                                                             
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRIXPS(KR,LR        ,M,A,B,C,Y,TWOCOS,D,W,SING,S)      
      DIMENSION A(1),B(1),C(1),Y(1),TWOCOS(1),D(1),W(1),S(1)            
      CALL       TRIXP(KR,LR,M,A,B,C,Y,TWOCOS,D,W,S)                    
      MM1=M-1                                                           
      D(1) = C(1)/(B(1) - 2.)                                           
      W(1) = A(1)/(B(1) - 2.)                                           
      Y(1) = Y(1)/(B(1) - 2.)                                           
      BM = B(M)                                                         
      Z = C(M)                                                          
      IF(MM1.EQ.1) GOTO 501                                             
      DO 500 I=2,MM1                                                    
      DEN = B(I) - 2.- A(I)*D(I-1)                                      
      D(I) = C(I)/DEN                                                   
      W(I) = -A(I)*W(I-1)/DEN                                           
      Y(I) = (Y(I) - A(I)*Y(I-1))/DEN                                   
      Y(M) = Y(M) - Z*Y(I-1)                                            
      BM = BM - Z*W(I-1)                                                
  500 Z = -Z*D(I-1)                                                     
  501 D(MM1) = D(MM1) + W(MM1)                                          
      Z = A(M) + Z                                                      
      DEN = BM   - 2.- Z*D(MM1)                                         
      Y(M) = Y(M) - Z*Y(M-1)                                            
      IF(SING*DEN.NE.0.) GOTO 700                                       
  600 Y(M) = 0.                                                         
      GO TO 800                                                         
  700 Y(M) = Y(M)/DEN                                                   
  800 Y(MM1) = Y(MM1) - D(MM1)*Y(M)                                     
      IF(MM1.EQ.1) GOTO 901                                             
      DO 900 IP=2,MM1                                                   
      I = M - IP                                                        
  900 Y(I) = Y(I) - D(I)*Y(I+1) - W(I)*Y(M)                             
  901 RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRI3S (M,A,B,C,K,Y,Y2,Y3,TCOS,D,W1,W2,W3,SING,S)       
      DIMENSION A(1),B(1),C(1),K(4),Y(1),Y2(1),Y3(1),TCOS(1),D(1),S(1)  
     1, W1(1),W2(1),W3(1)                                               
      CALL TRI3(M,A,B,C,K,Y,Y2,Y3,TCOS,D,W1,W2,W3,S)                    
      MM1=M-1                                                           
      Z=1./(B(1)-2.)                                                    
      D(1)=C(1)*Z                                                       
      Y(1)= (Y(1)+Y2(1)+Y3(1))*Z                                        
      IF(MM1.LT.2)GO TO 501                                             
      IM1=1                                                             
      DO 500 I=2,MM1                                                    
      Z=1./(B(I)-2.-A(I)*D(IM1))                                        
      D(I)=C(I)*Z                                                       
      Y(I)=(Y(I)+Y2(I)+Y3(I)-A(I)*Y(IM1))*Z                             
  500 IM1=I                                                             
  501 Z=B(M)-2.-A(M)*D(MM1)                                             
      X=(Y(M)+Y2(M)+Y3(M))-A(M)*Y(MM1)                                  
      IF(SING*Z.NE.0.) GOTO 550                                         
      Y(M)=0.                                                           
      GO TO 575                                                         
  550 Y(M)=X/Z                                                          
  575 X=Y(M)                                                            
      DO 600 IP=1,MM1                                                   
      I=M-IP                                                            
      X=Y(I)-D(I)*X                                                     
  600 Y(I)=X                                                            
      RETURN                                                            
      END                                                               
C                                                                       
C ********************************************************************* 
C                                                                       
      SUBROUTINE TRI3PS(M,A,B,C,K,Y,Y2,Y3,TCOS,D,W1,W2,W3,SING,S)       
      DIMENSION A(1),B(1),C(1),K(4),Y(1),Y2(1),Y3(1),TCOS(1),D(1)       
     1, W1(1),W2(1),W3(1),S(1)                                          
      CALL TRI3P(M,A,B,C,K,Y,Y2,Y3,TCOS,D,W1,W2,W3,S)                   
      MM1=M-1                                                           
      Z=1./(B(1)-2.)                                                    
      D(1)=C(1)*Z                                                       
      S(1)=A(1)*Z                                                       
      Y(1)= (Y(1)+Y2(1)+Y3(1))*Z                                        
      Y(M)=Y(M)+Y2(M)+Y3(M)                                             
      BM=B(M)                                                           
      CM=C(M)                                                           
      IF(MM1.LT.2)GO TO 501                                             
      IM1=1                                                             
      DO 500 I=2,MM1                                                    
      Z=1./(B(I)-2.-A(I)*D(IM1))                                        
      D(I)=C(I)*Z                                                       
      S(I)=-A(I)*S(IM1)*Z                                               
      Y(I)=(Y(I)+Y2(I)+Y3(I)-A(I)*Y(IM1))*Z                             
      Y(M)=Y(M)-CM*Y(IM1)                                               
      BM=BM-CM*S(IM1)                                                   
      CM=-CM*D(IM1)                                                     
  500 IM1=I                                                             
  501 D(MM1)=D(MM1)+S(MM1)                                              
      CM=A(M)+CM                                                        
      Z=BM-2.-CM*D(MM1)                                                 
      X=Y(M)-CM*Y(MM1)                                                  
      IF(SING*Z.NE.0.) GOTO 550                                         
      Y(M)=0.                                                           
      GO TO 575                                                         
  550 Y(M)=X/Z                                                          
  575 Y(MM1)=Y(MM1)-D(MM1)*Y(M)                                         
      IF (MM1 .EQ. 1) GO TO 601                                         
      X=Y(MM1)                                                          
      BM=Y(M)                                                           
      DO 600 IP=2,MM1                                                   
      I=M-IP                                                            
      X=Y(I)-D(I)*X-S(I)*BM                                             
  600 Y(I)=X                                                            
  601 RETURN                                                            
      END                                                               
