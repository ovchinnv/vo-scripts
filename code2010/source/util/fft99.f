
CcPA54.JOB.CNTL(FFT99)                                                  
      SUBROUTINE FFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)            
C***********************************************************************
C                                                                       
C**   *SUBROUTINE* *FFT99 * - ROUTINE FOR REAL/COMPLEX MULTIPLE FFT     
C                                                                       
C     CLIVE TEMPERTON         ECMWF NOVEMBER, 1978                      
C     CHANGED BY RUSS REW     NCAR  SEPTEMBER, 1980                     
C     MINOR ADJUSTMENTS BY ULRICH SCHUMANN, DFVLR JUNE 1984             
C       I.E.: N=4 ADMITTED, ELIMINATION OF Q8QST4 AND ULIBER            
C                                                                       
C     PURPOSE                                                           
C     -------                                                           
C                                                                       
C              PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS ROUTINE 
C              WILL PERFORM A NUMBER OF SIMULTANEOUS REAL/HALF-COMPLEX  
C              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE     
C              TRANSFORMS, I.E.  GIVEN A SET OF REAL DATA VECTORS, THE  
C              PACKAGE RETURNS A SET OF 'HALF-COMPLEX' FOURIER          
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE   
C              TRANSFORMS MUST BE AN EVEN NUMBER GREATER THAN 3 THAT HAS
C              NO OTHER FACTORS EXCEPT POSSIBLY POWERS OF 2, 3, AND 5.  
C              THIS IS AN ALL FORTRAN VERSION OF THE CRAYLIB PACKAGE    
C              THAT IS MOSTLY WRITTEN IN CAL.                           
C                                                                       
C   FFT99 COMPUTES FOR                                                  
C   ISIGN =  1: TRANSFORM FROM FOURIER COEFFICIENTS TO DATA VALUES      
C         = -1: TRANSFORM FROM DATA VALUES TO FOURIER COEFFICIENTS      
C                                                                       
C     FOR ISIGN = -1 FFT99 COMPUTES                                     
C     C(K,L) = (1/N)*SUM(J=0,...,N-1)(X(J,L)*EXP(-2*I*J*K*PI/N))        
C               K=0,1,2,...,(N/2)                                       
C               L=1,2,...,M                                             
C                                                                       
C     FOR ISIGN = 1  FFT99 COMPUTES                                     
C     X(J,L) =       SUM(K=0,...,N-1)(C(K,L)*EXP(+2*I*J*K*PI/N))        
C               J=1,2,3,...,N-1                                         
C               L=1,2,...,M                                             
C                                                                       
C     WHERE C(K,L) = AA(K,L) + I*BB(K,L) AND I=SQRT(-1.)                
C     AND   C(N-K,L)=AA(K,L) - I*BB(K,L)                                
C                                                                       
C            AS A CONSEQUENCE OF X BEING REAL, WE HAVE                  
C            BB(0)=BB(N/2)=0. FOR A CALL WITH ISIGN=-1 IT IS NOT        
C                             ACTUALLY NECESSARY TO SUPPLY THESE ZEROS. 
C     WITH THE FOLLOWING MAPPING                                        
C           AA(K,L) -> A(K*INC*2+(L-1)*JUMP+1)                          
C           BB(K,L) -> A(K*INC*2+(L-1)*JUMP+2)                          
C           C(N-K,L) = CONJG(C(K),L)                                    
C               K=0,1,2,...,(N/2)                                       
C               L=1,2,...,M                                             
C     AND    X(J,L) -> A(J*INC+(L-1)*JUMP+1)                            
C            X(0,L) = X(N,L)                                            
C            X(N+1,L)=X(1,L)                                            
C               J=1,2,...,N                                             
C               L=1,2,...,M                                             
C                                                                       
C NOTE: A CALL WITH ISIGN=1, FOLLOWED BY A CALL WITH ISIGN=-1           
C (OR VICE VERSA) RETURNS THE ORIGINAL DATA                             
C (THIS IS NOT THE CASE FOR CFFT99)                                     
C                                                                       
C     INTERFACE                                                         
C     ---------                                                         
C                                                                       
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 1, 
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF      
C              CALLS TO TRANSFORM A GIVEN SET OF REAL VECTORS OF LENGTH 
C              N TO A SET OF 'HALF-COMPLEX' FOURIER COEFFICIENT VECTORS 
C              OF LENGTH N IS                                           
C                                                                       
C                   DIMENSION IFAX(13),TRIGS(3*N/2+1),A(M*(N+2)),       
C                  +          WORK(M*(N+1))                             
C                                                                       
C                   CALL FFTFAX (N, IFAX, TRIGS)                        
C                   CALL FFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)   
C                                                                       
C              SEE THE INDIVIDUAL WRITE-UPS                             
C              BELOW, FOR A DETAILED DESCRIPTION OF THE                 
C              ARGUMENTS.                                               
C                                                                       
C                                                                       
C ARGUMENT     A(LA), WORK(M*(N+1)), TRIGS(3*N/2+1), IFAX(13)           
C DIMENSIONS     LA=(N+1)*INC+(M-1)*JUMP+1 >= M*(N+2)                   
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT     A                                                        
C               AN ARRAY OF LENGTH LA>=M*(N+2) CONTAINING THE INPUT DATA
C               OR COEFFICIENT VECTORS.  THIS ARRAY IS OVERWRITTEN BY   
C               THE RESULTS.                                            
C                                                                       
C              WORK                                                     
C               A WORK ARRAY OF DIMENSION M*(N+1)                       
C                                                                       
C              TRIGS                                                    
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.  
C                                                                       
C              IFAX                                                     
C               AN ARRAY SET UP BY FFTFAX, WHICH MUST BE CALLED FIRST.  
C                                                                       
C              INC                                                      
C               THE INCREMENT (IN WORDS) BETWEEN SUCCESSIVE ELEMENTS OF 
C               EACH DATA OR COEFFICIENT VECTOR (E.G. INC=1 FOR         
C               CONSECUTIVELY STORED DATA).                             
C                                                                       
C              JUMP                                                     
C               THE INCREMENT (IN WORDS) BETWEEN THE FIRST ELEMENTS OF  
C               SUCCESSIVE DATA OR COEFFICIENT VECTORS. ON THE CRAY-1,  
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8 
C               (TO AVOID MEMORY BANK CONFLICTS). FOR CLARIFICATION OF  
C               INC AND JUMP, SEE THE MAPPING EXPLAINED UNDER PURPOSE.  
C                                                                       
C              N                                                        
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF         
C               TRANSFORMS, BELOW).                                     
C                                                                       
C              M                                                        
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.     
C                                                                       
C              ISIGN                                                    
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO       
C                    GRIDPOINT VALUES.                                  
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER   
C                    COEFFICIENTS.                                      
C                                                                       
C                                                                       
C     RESULTS                                                           
C     -------                                                           
C                                                                       
C ON OUTPUT    A                                                        
C               IF ISIGN = +1, AND M COEFFICIENT VECTORS ARE SUPPLIED   
C               EACH CONTAINING THE SEQUENCE:                           
C                                                                       
C               A(0),B(0),A(1),B(1),...,A(N/2),B(N/2)  (N+2 VALUES)     
C                                                                       
C               THEN THE RESULT CONSISTS OF M DATA VECTORS EACH         
C               CONTAINING THE CORRESPONDING N+2 GRIDPOINT VALUES:      
C                                                                       
C               X(N-1),X(0),X(1),X(2),...,X(N-1),X(0).                  
C                   (EXPLICIT CYCLIC CONTINUITY)                        
C                                                                       
C               IF ISIGN = -1, AND M DATA VECTORS ARE SUPPLIED EACH     
C               CONTAINING A SEQUENCE OF GRIDPOINT VALUES X(J) AS       
C               DEFINED ABOVE, THEN THE RESULT CONSISTS OF M VECTORS    
C               EACH CONTAINING THE CORRESPONDING FOURIER COFFICIENTS   
C               A(K), B(K), 0 .LE. K .LE N/2.                           
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     METHOD                                                            
C     ------                                                            
C                                                                       
C     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM        
C     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12       
C     (1970), 315-337)                                                  
C                                                                       
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN      
C     PARALLEL                                                          
C                                                                       
C     EXTERNALS                                                         
C     ---------                                                         
C     SUBROUTINE TO BE CALLED BEFORE FFT99: FFTFAX                      
C     SUBROUTINES CALLED BY FFTFAX: FAX, FFTRIG                         
C     SUBROUTINES CALLED BY FFT99: FFT99A, FFT99B, VPASSM               
C                                                                       
C     REFERENCE                                                         
C     ---------                                                         
C                                                                       
C     C.TEMPERTON: FAST FOURIER TRANSFORMS ON CRAY-1.                   
C                  EUROPEAN CENTRE FOR MEDIUM RANGE WEATHER FORECASTS,  
C                  SHINFIELD PARK, INTERNAL REPORT NO. 21, JAN. 1979    
C                                                                       
C***********************************************************************
      DIMENSION A(*),WORK(*),TRIGS(*),IFAX(*)                           
      NFAX=IFAX(1)                                                      
      NX=N+1                                                            
      NH=N/2                                                            
      INK=INC+INC                                                       
      IF (ISIGN.EQ.+1) GO TO 30                                         
C                                                                       
C     IF NECESSARY, TRANSFER DATA TO WORK AREA                          
      IGO=50                                                            
      IF (MOD(NFAX,2).EQ.1) GOTO 40                                     
      IBASE=INC+1                                                       
      JBASE=1                                                           
      DO 20 L=1,M                                                       
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 10 MM=1,N                                                      
      WORK(J)=A(I)                                                      
      I=I+INC                                                           
      J=J+1                                                             
   10 CONTINUE                                                          
      IBASE=IBASE+JUMP                                                  
      JBASE=JBASE+NX                                                    
   20 CONTINUE                                                          
C                                                                       
      IGO=60                                                            
      GO TO 40                                                          
C                                                                       
C     PREPROCESSING (ISIGN=+1)                                          
C     ------------------------                                          
C                                                                       
   30 CONTINUE                                                          
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,M)                            
      IGO=60                                                            
C                                                                       
C     COMPLEX TRANSFORM                                                 
C     -----------------                                                 
C                                                                       
   40 CONTINUE                                                          
      IA=INC+1                                                          
      LA=1                                                              
      DO 80 K=1,NFAX                                                    
      IF (IGO.EQ.60) GO TO 60                                           
   50 CONTINUE                                                          
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,                
     *   INK,2,JUMP,NX,M,NH,IFAX(K+1),LA)                               
      IGO=60                                                            
      GO TO 70                                                          
   60 CONTINUE                                                          
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,                
     *    2,INK,NX,JUMP,M,NH,IFAX(K+1),LA)                              
      IGO=50                                                            
   70 CONTINUE                                                          
      LA=LA*IFAX(K+1)                                                   
   80 CONTINUE                                                          
C                                                                       
      IF (ISIGN.EQ.-1) GO TO 130                                        
C                                                                       
C     IF NECESSARY, TRANSFER DATA FROM WORK AREA                        
      IF (MOD(NFAX,2).EQ.1) GO TO 110                                   
      IBASE=1                                                           
      JBASE=IA                                                          
      DO 100 L=1,M                                                      
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 90 MM=1,N                                                      
      A(J)=WORK(I)                                                      
      I=I+1                                                             
      J=J+INC                                                           
   90 CONTINUE                                                          
      IBASE=IBASE+NX                                                    
      JBASE=JBASE+JUMP                                                  
  100 CONTINUE                                                          
C                                                                       
C     FILL IN CYCLIC BOUNDARY POINTS                                    
  110 CONTINUE                                                          
      IA=1                                                              
      IB=N*INC+1                                                        
CVD$  NODEPCHK                                                          
      DO 120 L=1,M                                                      
      A(IA)=A(IB)                                                       
      A(IB+INC)=A(IA+INC)                                               
      IA=IA+JUMP                                                        
      IB=IB+JUMP                                                        
  120 CONTINUE                                                          
      GO TO 140                                                         
C                                                                       
C     POSTPROCESSING (ISIGN=-1):                                        
C     --------------------------                                        
C                                                                       
  130 CONTINUE                                                          
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,M)                            
C                                                                       
  140 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,M)                      
C     *SUBROUTINE* *FFT99A* - ROUTINE USED BY FFT99                     
C                                                                       
C     SUBROUTINE FFT99A - PREPROCESSING STEP FOR FFT99, ISIGN=+1        
C     (SPECTRAL TO GRIDPOINT TRANSFORM)                                 
C                                                                       
      DIMENSION A(*),WORK(*),TRIGS(*)                                   
      NH=N/2                                                            
      NX=N+1                                                            
      INK=INC+INC                                                       
C     A(0) AND A(N/2)                                                   
      IA=1                                                              
      IB=N*INC+1                                                        
      JA=1                                                              
      JB=2                                                              
CVD$  NODEPCHK                                                          
      DO 10 L=1,M                                                       
      WORK(JA)=A(IA)+A(IB)                                              
      WORK(JB)=A(IA)-A(IB)                                              
      IA=IA+JUMP                                                        
      IB=IB+JUMP                                                        
      JA=JA+NX                                                          
      JB=JB+NX                                                          
   10 CONTINUE                                                          
C     REMAINING WAVENUMBERS                                             
      IABASE=2*INC+1                                                    
      IBBASE=(N-2)*INC+1                                                
      JABASE=3                                                          
      JBBASE=N-1                                                        
C                                                                       
      DO 30 K=3,NH,2                                                    
      IA=IABASE                                                         
      IB=IBBASE                                                         
      JA=JABASE                                                         
      JB=JBBASE                                                         
      C=TRIGS(N+K)                                                      
      S=TRIGS(N+K+1)                                                    
CVD$  NODEPCHK                                                          
      DO 20 L=1,M                                                       
      WORK(JA)=(A(IA)+A(IB))-                                           
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))                     
      WORK(JB)=(A(IA)+A(IB))+                                           
     *    (S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))                     
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+             
     *    (A(IA+INC)-A(IB+INC))                                         
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-             
     *    (A(IA+INC)-A(IB+INC))                                         
      IA=IA+JUMP                                                        
      IB=IB+JUMP                                                        
      JA=JA+NX                                                          
      JB=JB+NX                                                          
   20 CONTINUE                                                          
      IABASE=IABASE+INK                                                 
      IBBASE=IBBASE-INK                                                 
      JABASE=JABASE+2                                                   
      JBBASE=JBBASE-2                                                   
   30 CONTINUE                                                          
C                                                                       
      IF (IABASE.NE.IBBASE) GO TO 50                                    
C     WAVENUMBER N/4 (IF IT EXISTS)                                     
      IA=IABASE                                                         
      JA=JABASE                                                         
CVD$  NODEPCHK                                                          
      DO 40 L=1,M                                                       
      WORK(JA)=2.0*A(IA)                                                
      WORK(JA+1)=-2.0*A(IA+INC)                                         
      IA=IA+JUMP                                                        
      JA=JA+NX                                                          
   40 CONTINUE                                                          
C                                                                       
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,M)                      
C     *SUBROUTINE* *FFT99B* - ROUTINE USED BY FFT99                     
C                                                                       
C     SUBROUTINE FFT99B - POSTPROCESSING STEP FOR FFT99, ISIGN=-1       
C     (GRIDPOINT TO SPECTRAL TRANSFORM)                                 
C                                                                       
      DIMENSION WORK(*),A(*),TRIGS(*)                                   
      NH=N/2                                                            
      NX=N+1                                                            
      INK=INC+INC                                                       
C     A(0) AND A(N/2)                                                   
      SCALE=1.0/FLOAT(N)                                                
      IA=1                                                              
      IB=2                                                              
      JA=1                                                              
      JB=N*INC+1                                                        
CVD$  NODEPCHK                                                          
      DO 10 L=1,M                                                       
      A(JA)=SCALE*(WORK(IA)+WORK(IB))                                   
      A(JB)=SCALE*(WORK(IA)-WORK(IB))                                   
      A(JA+INC)=0.0                                                     
      A(JB+INC)=0.0                                                     
      IA=IA+NX                                                          
      IB=IB+NX                                                          
      JA=JA+JUMP                                                        
      JB=JB+JUMP                                                        
   10 CONTINUE                                                          
C     REMAINING WAVENUMBERS                                             
      SCALE=0.5*SCALE                                                   
      IABASE=3                                                          
      IBBASE=N-1                                                        
      JABASE=2*INC+1                                                    
      JBBASE=(N-2)*INC+1                                                
C                                                                       
      DO 30 K=3,NH,2                                                    
      IA=IABASE                                                         
      IB=IBBASE                                                         
      JA=JABASE                                                         
      JB=JBBASE                                                         
      C=TRIGS(N+K)                                                      
      S=TRIGS(N+K+1)                                                    
CVD$  NODEPCHK                                                          
      DO 20 L=1,M                                                       
      A(JA)=SCALE*((WORK(IA)+WORK(IB))                                  
     *   +(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))            
      A(JB)=SCALE*((WORK(IA)+WORK(IB))                                  
     *   -(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))            
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    +(WORK(IB+1)-WORK(IA+1)))                                     
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))
     *    -(WORK(IB+1)-WORK(IA+1)))                                     
      IA=IA+NX                                                          
      IB=IB+NX                                                          
      JA=JA+JUMP                                                        
      JB=JB+JUMP                                                        
   20 CONTINUE                                                          
      IABASE=IABASE+2                                                   
      IBBASE=IBBASE-2                                                   
      JABASE=JABASE+INK                                                 
      JBBASE=JBBASE-INK                                                 
   30 CONTINUE                                                          
C                                                                       
      IF (IABASE.NE.IBBASE) GO TO 50                                    
C     WAVENUMBER N/4 (IF IT EXISTS)                                     
      IA=IABASE                                                         
      JA=JABASE                                                         
      SCALE=2.0*SCALE                                                   
CVD$  NODEPCHK                                                          
      DO 40 L=1,M                                                       
      A(JA)=SCALE*WORK(IA)                                              
      A(JA+INC)=-SCALE*WORK(IA+1)                                       
      IA=IA+NX                                                          
      JA=JA+JUMP                                                        
   40 CONTINUE                                                          
C                                                                       
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE FFTFAX(N,IFAX,TRIGS)                                   
C***********************************************************************
C                                                                       
C**   *SUBROUTINE* *FFTFAX* - ROUTINE FOR SET-UP OF IFAX/TRIGS FOR FFT99
C                                                                       
C     CLIVE TEMPERTON         ECMWF NOVEMBER, 1978                      
C     CHANGED BY RUSS REW     NCAR  SEPTEMBER, 1980                     
C     MINOR ADJUSTMENTS BY ULRICH SCHUMANN, DFVLR JUNE 1984             
C       I.E.: N=4 ADMITTED, ELIMINATION OF Q8QST4 AND ULIBER            
C                                                                       
C     PURPOSE                                                           
C     -------                                                           
C                                                                       
C PURPOSE      A SET-UP ROUTINE FOR FFT99 AND FFT991.  IT NEED ONLY BE  
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO THE FFT        
C              ROUTINES (PROVIDED THAT N IS NOT CHANGED).               
C                                                                       
C     INTERFACE                                                         
C     ---------                                                         
C                                                                       
C ARGUMENT     IFAX(13),TRIGS(3*N/2+1)                                  
C DIMENSIONS                                                            
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT     N                                                        
C               AN EVEN NUMBER GREATER THAN 3 THAT HAS NO PRIME FACTOR  
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE 
C               THE DOCUMENTATION FOR FFT99 FOR THE                     
C               DEFINITIONS OF THE TRANSFORMS).                         
C                                                                       
C              IFAX                                                     
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED 
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING    
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN A MILLION.     
C                                                                       
C              TRIGS                                                    
C               A FLOATING POINT ARRAY OF DIMENSION 3*N/2 IF N/2 IS     
C               EVEN, OR 3*N/2+1 IF N/2 IS ODD.                         
C                                                                       
C ON OUTPUT    IFAX                                                     
C               CONTAINS THE FACTORIZATION OF N/2.  IFAX(1) IS THE      
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF FFTFAX IS CALLED WITH N ODD, 
C               OR IF N HAS ANY PRIME FACTORS GREATER THAN 5, IFAX(1)   
C               IS SET TO -99.                                          
C                                                                       
C              TRIGS                                                    
C               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY  
C               USED BY THE FFT ROUTINES.                               
C                                                                       
C     EXTERNALS                                                         
C     ---------                                                         
C     FAX, FFTRIG                                                       
C                                                                       
C***********************************************************************
      DIMENSION IFAX(13),TRIGS(*)                                       
C                                                                       
C MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE      
C TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT       
C DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE     
C WAS WRITTEN.                                                          
C                                                                       
      DATA MODE /3/                                                     
      CALL FAX (IFAX, N, MODE)                                          
      I = IFAX(1)                                                       
C     IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99                 
      IF (IFAX(I+1) .GT. 5 .OR. N .LT. 4) IFAX(1) = -99                 
      IF (IFAX(1) .GT. 0 )                                              
     +CALL FFTRIG (TRIGS, N, MODE)                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE FAX(IFAX,N,MODE)                                       
C     *SUBROUTINE* *FAX   * - ROUTINE CALLED BY FFTFAX                  
      DIMENSION IFAX(10)                                                
      NN=N                                                              
      IF (IABS(MODE).EQ.1) GO TO 10                                     
      IF (IABS(MODE).EQ.8) GO TO 10                                     
      NN=N/2                                                            
      IF ((NN+NN).EQ.N) GO TO 10                                        
      IFAX(1)=-99                                                       
      RETURN                                                            
   10 K=1                                                               
C     TEST FOR FACTORS OF 4                                             
   20 IF (MOD(NN,4).NE.0) GO TO 30                                      
      K=K+1                                                             
      IFAX(K)=4                                                         
      NN=NN/4                                                           
      IF (NN.EQ.1) GO TO 80                                             
      GO TO 20                                                          
C     TEST FOR EXTRA FACTOR OF 2                                        
   30 IF (MOD(NN,2).NE.0) GO TO 40                                      
      K=K+1                                                             
      IFAX(K)=2                                                         
      NN=NN/2                                                           
      IF (NN.EQ.1) GO TO 80                                             
C     TEST FOR FACTORS OF 3                                             
   40 IF (MOD(NN,3).NE.0) GO TO 50                                      
      K=K+1                                                             
      IFAX(K)=3                                                         
      NN=NN/3                                                           
      IF (NN.EQ.1) GO TO 80                                             
      GO TO 40                                                          
C     NOW FIND REMAINING FACTORS                                        
   50 L=5                                                               
      INC=2                                                             
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4                           
   60 IF (MOD(NN,L).NE.0) GO TO 70                                      
      K=K+1                                                             
      IFAX(K)=L                                                         
      NN=NN/L                                                           
      IF (NN.EQ.1) GO TO 80                                             
      GO TO 60                                                          
   70 L=L+INC                                                           
      INC=6-INC                                                         
      GO TO 60                                                          
   80 IFAX(1)=K-1                                                       
C     IFAX(1) CONTAINS NUMBER OF FACTORS                                
      NFAX=IFAX(1)                                                      
C     SORT FACTORS INTO ASCENDING ORDER                                 
      IF (NFAX.EQ.1) GO TO 110                                          
      DO 100 II=2,NFAX                                                  
      ISTOP=NFAX+2-II                                                   
      DO 90 I=2,ISTOP                                                   
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90                                
      ITEM=IFAX(I)                                                      
      IFAX(I)=IFAX(I+1)                                                 
      IFAX(I+1)=ITEM                                                    
   90 CONTINUE                                                          
  100 CONTINUE                                                          
  110 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE FFTRIG(TRIGS,N,MODE)                                   
C     *SUBROUTINE* *FFTRIG* - ROUTINE CALLED BY FFTFAX                  
      DIMENSION TRIGS(*)                                                
      PI=2.0*ASIN(1.0)                                                  
      IMODE=IABS(MODE)                                                  
      NN=N                                                              
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2                             
      DEL=(PI+PI)/FLOAT(NN)                                             
      L=NN+NN                                                           
      DO 10 I=1,L,2                                                     
      ANGLE=0.5*FLOAT(I-1)*DEL                                          
      TRIGS(I)=COS(ANGLE)                                               
      TRIGS(I+1)=SIN(ANGLE)                                             
   10 CONTINUE                                                          
      IF (IMODE.EQ.1) RETURN                                            
      IF (IMODE.EQ.8) RETURN                                            
      DEL=0.5*DEL                                                       
      NH=(NN+1)/2                                                       
      L=NH+NH                                                           
      LA=NN+NN                                                          
      DO 20 I=1,L,2                                                     
      ANGLE=0.5*FLOAT(I-1)*DEL                                          
      TRIGS(LA+I)=COS(ANGLE)                                            
      TRIGS(LA+I+1)=SIN(ANGLE)                                          
   20 CONTINUE                                                          
      IF (IMODE.LE.3) RETURN                                            
      DEL=0.5*DEL                                                       
      LA=LA+NN                                                          
      IF (MODE.EQ.5) GO TO 40                                           
      DO 30 I=2,NN                                                      
      ANGLE=FLOAT(I-1)*DEL                                              
      TRIGS(LA+I)=2.0*SIN(ANGLE)                                        
   30 CONTINUE                                                          
      RETURN                                                            
   40 CONTINUE                                                          
      DEL=0.5*DEL                                                       
      DO 50 I=2,N                                                       
      ANGLE=FLOAT(I-1)*DEL                                              
      TRIGS(LA+I)=SIN(ANGLE)                                            
   50 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE CFFT99(A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)           
C***********************************************************************
C                                                                       
C**   *SUBROUTINE* *CFFT99* - ROUTINE FOR COMPEX/COMPLEX MULTIPLE FFT   
C                                                                       
C     CLIVE TEMPERTON         ECMWF NOVEMBER, 1978                      
C     CHANGED BY RUSS REW     NCAR  SEPTEMBER, 1980                     
C            AND DAVE FULKER  NCAR  NOVEMBER, 1980                      
C     MINOR ADJUSTMENTS BY ULRICH SCHUMANN, DFVLR JUNE 1984             
C       I.E.: ELIMINATION OF Q8QST4 AND ULIBER                          
C                                                                       
C     PURPOSE                                                           
C     -------                                                           
C                                                                       
C                                                                       
C PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS. THIS ROUTINE  
C              WILL PERFORM A NUMBER OF SIMULTANEOUS COMPLEX PERIODIC   
C              FOURIER TRANSFORMS OR CORRESPONDING INVERSE TRANSFORMS.  
C              THAT IS, GIVEN A SET OF COMPLEX GRIDPOINT VECTORS, THE   
C              PACKAGE RETURNS A SET OF COMPLEX FOURIER                 
C              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE   
C              TRANSFORMS MUST BE A NUMBER GREATER THAN 1 THAT HAS      
C              NO PRIME FACTORS OTHER THAN 2, 3, AND 5.                 
C                                                                       
C                                                                       
C USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 0, 
C              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF      
C              CALLS TO TRANSFORM A GIVEN SET OF COMPLEX VECTORS OF     
C              LENGTH N TO A SET OF (UNSCALED) COMPLEX FOURIER          
C              COEFFICIENT VECTORS OF LENGTH N IS                       
C                                                                       
C                   DIMENSION IFAX(13),TRIGS(2*N)                       
C                   COMPLEX A(...), WORK(...)                           
C                                                                       
C                   CALL CFTFAX (N, IFAX, TRIGS)                        
C                   CALL CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,M,ISIGN)CF
C                                                                       
C              THE OUTPUT VECTORS OVERWRITE THE INPUT VECTORS, AND      
C              THESE ARE STORED IN A.  WITH APPROPRIATE CHOICES FOR     
C              THE OTHER ARGUMENTS, THESE VECTORS MAY BE CONSIDERED     
C              EITHER THE ROWS OR THE COLUMNS OF THE ARRAY A.           
C              SEE BELOW.                                               
C                                                                       
C     INTERFACE                                                         
C     ---------                                                         
C                                                                       
C ARGUMENT     COMPLEX A(N*INC+(M-1)*JUMP), WORK(N*M)                   
C DIMENSIONS   REAL TRIGS(2*N), INTEGER IFAX(13)                        
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT     A                                                        
C               A COMPLEX ARRAY OF LENGTH N*INC+(M-1)*JUMP CONTAINING   
C               THE INPUT GRIDPOINT OR COEFFICIENT VECTORS.  THIS ARRAY 
C               IS OVERWRITTEN BY THE RESULTS.                          
C                                                                       
C               N.B. ALTHOUGH THE ARRAY A IS USUALLY CONSIDERED TO BE OF
C               TYPE COMPLEX IN THE CALLING PROGRAM, IT IS TREATED AS   
C               REAL WITHIN THE TRANSFORM PACKAGE.  THIS REQUIRES THAT  
C               SUCH TYPE CONFLICTS ARE PERMITTED IN THE USER"S         
C               ENVIRONMENT, AND THAT THE STORAGE OF COMPLEX NUMBERS    
C               MATCHES THE ASSUMPTIONS OF THIS ROUTINE.  THIS ROUTINE  
C               ASSUMES THAT THE REAL AND IMAGINARY PORTIONS OF A       
C               COMPLEX NUMBER OCCUPY ADJACENT ELEMENTS OF MEMORY.  IF  
C               THESE CONDITIONS ARE NOT MET, THE USER MUST TREAT THE   
C               ARRAY A AS REAL (AND OF TWICE THE ABOVE LENGTH), AND    
C               WRITE THE CALLING PROGRAM TO TREAT THE REAL AND         
C               IMAGINARY PORTIONS EXPLICITLY.                          
C                                                                       
C              WORK                                                     
C               A COMPLEX WORK ARRAY OF LENGTH N*M OR A REAL ARRAY      
C               OF LENGTH 2*N*M.    SEE N.B. ABOVE.                     
C                                                                       
C              TRIGS                                                    
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.  
C                                                                       
C              IFAX                                                     
C               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.  
C                                                                       
C                                                                       
C               N.B. IN THE FOLLOWING ARGUMENTS, INCREMENTS ARE MEASURED
C               IN WORD PAIRS, BECAUSE EACH COMPLEX ELEMENT IS ASSUMED  
C               TO OCCUPY AN ADJACENT PAIR OF WORDS IN MEMORY.          
C                                                                       
C              INC                                                      
C               THE INCREMENT (IN WORD PAIRS) BETWEEN SUCCESSIVE ELEMENT
C               OF EACH (COMPLEX) GRIDPOINT OR COEFFICIENT VECTOR       
C               (E.G.  INC=1 FOR CONSECUTIVELY STORED DATA).            
C                                                                       
C              JUMP                                                     
C               THE INCREMENT (IN WORD PAIRS) BETWEEN THE FIRST ELEMENTS
C               OF SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-
C               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8 
C               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF 
C               INC AND JUMP, SEE THE EXAMPLES BELOW.                   
C                                                                       
C              N                                                        
C               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF         
C               TRANSFORMS, BELOW).                                     
C                                                                       
C              M                                                        
C               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.     
C                                                                       
C              ISIGN                                                    
C               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER   
C                    COEFFICIENTS.                                      
C               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO       
C                    GRIDPOINT VALUES.                                  
C                                                                       
C ON OUTPUT    A                                                        
C                                                                       
C               FOR ISIGN = -1 CFFT99 COMPUTES                          
C            C(K,L) =       SUM(J=0,...,N-1)(X(J,L)*EXP(-2*I*J*K*PI/N)) 
C                K=0,1,2,...,N-1                                        
C                L=1,2,...,M                                            
C                                                                       
C               FOR ISIGN = 1  CFFT99 COMPUTES                          
C            X(J,L) =       SUM(K=0,...,N-1)(C(K,L)*EXP(+2*I*J*K*PI/N)) 
C                J=0,1,2,...,N-1                                        
C                L=1,2,...,M                                            
C                                                                       
C               WHERE C(K) = AA(K) + I*BB(K) AND I=SQRT(-1.)            
C                  X(J) = YY(J) + I*ZZ(J)                               
C                                                                       
C               WITH  C(K,L) -> A(K*INC+(L-1)*JUMP+1)                   
C                K=0,1,2,...,N-1                                        
C                L=1,2,...,M                                            
C               WHERE  X(J,L) -> A(J*INC+(L-1)*JUMP+1)                  
C                J=0,1,...,N-1                                          
C                L=1,2,...,M                                            
C                                                                       
C                                                                       
C               A CALL WITH ISIGN=-1 FOLLOWED BY A CALL WITH ISIGN=+1   
C               (OR VICE VERSA) RETURNS THE ORIGINAL DATA, MULTIPLIED   
C               BY THE FACTOR N. (DIFFERS FROM FFT99 IN THIS RESPECT)   
C                                                                       
C                                                                       
C EXAMPLE       GIVEN A 64 BY 9 GRID OF COMPLEX VALUES, STORED IN       
C               A 66 BY 9 COMPLEX ARRAY, A, COMPUTE THE TWO DIMENSIONAL 
C               FOURIER TRANSFORM OF THE GRID.  FROM TRANSFORM THEORY,  
C               IT IS KNOWN THAT A TWO DIMENSIONAL TRANSFORM CAN BE     
C               OBTAINED BY FIRST TRANSFORMING THE GRID ALONG ONE       
C               DIRECTION, THEN TRANSFORMING THESE RESULTS ALONG THE    
C               ORTHOGONAL DIRECTION.                                   
C                                                                       
C               COMPLEX A(66,9), WORK(64,9)                             
C               REAL TRIGS1(128), TRIGS2(18)                            
C               INTEGER IFAX1(13), IFAX2(13)                            
C                                                                       
C               SET UP THE IFAX AND TRIGS ARRAYS FOR EACH DIRECTION:    
C                                                                       
C               CALL CFTFAX(64, IFAX1, TRIGS1)                          
C               CALL CFTFAX( 9, IFAX2, TRIGS2)                          
C                                                                       
C               IN THIS CASE, THE COMPLEX VALUES OF THE GRID ARE        
C               STORED IN MEMORY AS FOLLOWS (USING U AND V TO           
C               DENOTE THE REAL AND IMAGINARY COMPONENTS, AND           
C               ASSUMING CONVENTIONAL FORTRAN STORAGE):                 
C                                                                       
C   U(1,1), V(1,1), U(2,1), V(2,1),  ...  U(64,1), V(64,1), 4 NULLS,    
C                                                                       
C   U(1,2), V(1,2), U(2,2), V(2,2),  ...  U(64,2), V(64,2), 4 NULLS,    
C                                                                       
C   .       .       .       .         .   .        .        .           
C   .       .       .       .         .   .        .        .           
C   .       .       .       .         .   .        .        .           
C                                                                       
C   U(1,9), V(1,9), U(2,9), V(2,9),  ...  U(64,9), V(64,9), 4 NULLS.    
C                                                                       
C               WE CHOOSE (ARBITRARILY) TO TRANSORM FIRST ALONG THE     
C               DIRECTION OF THE FIRST SUBSCRIPT.  THUS WE DEFINE       
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 64, THE          
C               NUMBER OF TRANSFORMS, M, TO BE 9, THE INCREMENT         
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 1,       
C               AND THE INCREMENT BETWEEN THE STARTING POINTS           
C               FOR EACH TRANSFORM, JUMP, TO BE 66 (THE FIRST           
C               DIMENSION OF A).                                        
C                                                                       
C               CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, 66, 64, 9, -1)  
C                                                                       
C               TO TRANSFORM ALONG THE DIRECTION OF THE SECOND SUBSCRIPT
C               THE ROLES OF THE INCREMENTS ARE REVERSED.  THUS WE DEFIN
C               THE LENGTH OF THE TRANSFORMS, N, TO BE 9, THE           
C               NUMBER OF TRANSFORMS, M, TO BE 64, THE INCREMENT        
C               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 66,      
C               AND THE INCREMENT BETWEEN THE STARTING POINTS           
C               FOR EACH TRANSFORM, JUMP, TO BE 1                       
C                                                                       
C               CALL CFFT99( A, WORK, TRIGS2, IFAX2, 66, 1, 9, 64, -1)  
C                                                                       
C               THESE TWO SEQUENTIAL STEPS RESULTS IN THE TWO-DIMENSIONA
C               FOURIER COEFFICIENT ARRAY OVERWRITING THE INPUT         
C               GRIDPOINT ARRAY, A.  THE SAME TWO STEPS APPLIED AGAIN   
C               WITH ISIGN = +1 WOULD RESULT IN THE RECONSTRUCTION OF   
C               THE GRIDPOINT ARRAY (MULTIPLIED BY A FACTOR OF 64*9).   
C                                                                       
C                                                                       
C-----------------------------------------------------------------------
C                                                                       
C     METHOD                                                            
C     ------                                                            
C                                                                       
C     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN      
C     PARALLEL                                                          
C                                                                       
C     EXTERNALS                                                         
C     ---------                                                         
C     SUBROUTINE TO BE CALLED BEFORE CFFT99: CFTFAX                     
C     SUBROUTINES CALLED BY CFTFAX: FACT, CFTRIG                        
C     SUBROUTINES CALLED BY CFFT99:  VPASSM                             
C                                                                       
C     REFERENCE                                                         
C     ---------                                                         
C                                                                       
C     C.TEMPERTON: FAST FOURIER TRANSFORMS ON CRAY-1.                   
C                  EUROPEAN CENTRE FOR MEDIUM RANGE WEATHER FORECASTS,  
C                  SHINFIELD PARK, INTERNAL REPORT NO. 21, JAN. 1979    
C                                                                       
C***********************************************************************
      DIMENSION A(*),WORK(*),TRIGS(*),IFAX(*)                           
C                                                                       
      NN = N+N                                                          
      INK=INC+INC                                                       
      JUM = JUMP+JUMP                                                   
      NFAX=IFAX(1)                                                      
      JNK = 2                                                           
      JST = 2                                                           
      IF (ISIGN.GE.0) GO TO 30                                          
C                                                                       
C     THE INNERMOST TEMPERTON ROUTINES HAVE NO FACILITY FOR THE         
C     FORWARD (ISIGN = -1) TRANSFORM.  THEREFORE, THE INPUT MUST BE     
C     REARRANGED AS FOLLOWS:                                            
C                                                                       
C     THE ORDER OF EACH INPUT VECTOR,                                   
C                                                                       
C     G(0), G(1), G(2), ... , G(N-2), G(N-1)                            
C                                                                       
C     IS REVERSED (EXCLUDING G(0)) TO YIELD                             
C                                                                       
C     G(0), G(N-1), G(N-2), ... , G(2), G(1).                           
C                                                                       
C     WITHIN THE TRANSFORM, THE CORRESPONDING EXPONENTIAL MULTIPLIER    
C     IS THEN PRECISELY THE CONJUGATE OF THAT FOR THE NORMAL            
C     ORDERING.  THUS THE FORWARD (ISIGN = -1) TRANSFORM IS             
C     ACCOMPLISHED                                                      
C                                                                       
C     FOR NFAX ODD, THE INPUT MUST BE TRANSFERRED TO THE WORK ARRAY,    
C     AND THE REARRANGEMENT CAN BE DONE DURING THE MOVE.                
C                                                                       
      JNK = -2                                                          
      JST = NN-2                                                        
      IF (MOD(NFAX,2).EQ.1) GOTO 40                                     
C                                                                       
C     FOR NFAX EVEN, THE REARRANGEMENT MUST BE APPLIED DIRECTLY TO      
C     THE INPUT ARRAY.  THIS CAN BE DONE BY SWAPPING ELEMENTS.          
C                                                                       
      IBASE = 1                                                         
      ILAST = (N-1)*INK                                                 
      NH = N/2                                                          
      DO 20 L=1,M                                                       
      I1 = IBASE+INK                                                    
      I2 = IBASE+ILAST                                                  
CVD$  NODEPCHK								08600005
      DO 10 MM=1,NH                                                     
C     SWAP REAL AND IMAGINARY PORTIONS                                  
      HREAL = A(I1)                                                     
      HIMAG = A(I1+1)                                                   
      A(I1) = A(I2)                                                     
      A(I1+1) = A(I2+1)                                                 
      A(I2) = HREAL                                                     
      A(I2+1) = HIMAG                                                   
      I1 = I1+INK                                                       
      I2 = I2-INK                                                       
   10 CONTINUE                                                          
      IBASE = IBASE+JUM                                                 
   20 CONTINUE                                                          
      GOTO 100                                                          
C                                                                       
   30 CONTINUE                                                          
      IF (MOD(NFAX,2).EQ.0) GOTO 100                                    
C                                                                       
   40 CONTINUE                                                          
C                                                                       
C     DURING THE TRANSFORM PROCESS, NFAX STEPS ARE TAKEN, AND THE       
C     RESULTS ARE STORED ALTERNATELY IN WORK AND IN A.  IF NFAX IS      
C     ODD, THE INPUT DATA ARE FIRST MOVED TO WORK SO THAT THE FINAL     
C     RESULT (AFTER NFAX STEPS) IS STORED IN ARRAY A.                   
C                                                                       
      IBASE=1                                                           
      JBASE=1                                                           
      DO 60 L=1,M                                                       
C     MOVE REAL AND IMAGINARY PORTIONS OF ELEMENT ZERO                  
      WORK(JBASE) = A(IBASE)                                            
      WORK(JBASE+1) = A(IBASE+1)                                        
      I=IBASE+INK                                                       
      J=JBASE+JST                                                       
CVD$  NODEPCHK                                                          
      DO 50 MM=2,N                                                      
C     MOVE REAL AND IMAGINARY PORTIONS OF OTHER ELEMENTS (POSSIBLY IN   
C     REVERSE ORDER, DEPENDING ON JST AND JNK)                          
      WORK(J) = A(I)                                                    
      WORK(J+1) = A(I+1)                                                
      I=I+INK                                                           
      J=J+JNK                                                           
   50 CONTINUE                                                          
      IBASE=IBASE+JUM                                                   
      JBASE=JBASE+NN                                                    
   60 CONTINUE                                                          
C                                                                       
  100 CONTINUE                                                          
C                                                                       
C     PERFORM THE TRANSFORM PASSES, ONE PASS FOR EACH FACTOR.  DURING   
C     EACH PASS THE DATA ARE MOVED FROM A TO WORK OR FROM WORK TO A.    
C                                                                       
C     FOR NFAX EVEN, THE FIRST PASS MOVES FROM A TO WORK                
      IGO = 110                                                         
C     FOR NFAX ODD, THE FIRST PASS MOVES FROM WORK TO A                 
      IF (MOD(NFAX,2).EQ.1) IGO = 120                                   
      LA=1                                                              
      DO 140 K=1,NFAX                                                   
      IF (IGO.EQ.120) GO TO 120                                         
  110 CONTINUE                                                          
      CALL VPASSM(A(1),A(2),WORK(1),WORK(2),TRIGS,                      
     *   INK,2,JUM,NN,M,N,IFAX(K+1),LA)                                 
      IGO=120                                                           
      GO TO 130                                                         
  120 CONTINUE                                                          
      CALL VPASSM(WORK(1),WORK(2),A(1),A(2),TRIGS,                      
     *    2,INK,NN,JUM,M,N,IFAX(K+1),LA)                                
      IGO=110                                                           
  130 CONTINUE                                                          
      LA=LA*IFAX(K+1)                                                   
  140 CONTINUE                                                          
C                                                                       
C     AT THIS POINT THE FINAL TRANSFORM RESULT IS STORED IN A.          
C                                                                       
      RETURN                                                            
      END                                                               
      SUBROUTINE CFTFAX(N,IFAX,TRIGS)                                   
C***********************************************************************
C                                                                       
C**   *SUBROUTINE* *CFTFAX* -ROUTINE FOR SET-UP OF IFAX/TRIGS FOR CFFT99
C                                                                       
C     CLIVE TEMPERTON         ECMWF NOVEMBER, 1978                      
C     CHANGED BY RUSS REW     NCAR  SEPTEMBER, 1980                     
C            AND DAVE FULKER  NCAR  NOVEMBER, 1980                      
C     MINOR ADJUSTMENTS BY ULRICH SCHUMANN, DFVLR JUNE 1984             
C       I.E.: ELIMINATION OF Q8QST4 AND ULIBER                          
C                                                                       
C     PURPOSE                                                           
C     -------                                                           
C              A SET-UP ROUTINE FOR CFFT99.  IT NEED ONLY BE            
C              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO CFFT99,        
C              PROVIDED THAT N IS NOT CHANGED.                          
C                                                                       
C     INTERFACE                                                         
C     ---------                                                         
C                                                                       
C ARGUMENT     IFAX(13),TRIGS(2*N)                                      
C DIMENSIONS                                                            
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT     N                                                        
C               AN EVEN NUMBER GREATER THAN 1 THAT HAS NO PRIME FACTOR  
C               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE 
C               THE DOCUMENTATION FOR CFFT99 FOR THE DEFINITION OF      
C               THE TRANSFORMS).                                        
C                                                                       
C              IFAX                                                     
C               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED 
C               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING    
C               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN 1 MILLION.     
C                                                                       
C              TRIGS                                                    
C               A REAL ARRAY OF DIMENSION 2*N                           
C                                                                       
C ON OUTPUT    IFAX                                                     
C               CONTAINS THE FACTORIZATION OF N.  IFAX(1) IS THE        
C               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
C               IN IFAX(2),IFAX(3),...  IF N HAS ANY PRIME FACTORS      
C               GREATER THAN 5, IFAX(1) IS SET TO -99.                  
C                                                                       
C              TRIGS                                                    
C               AN ARRAY OF TRIGONOMETRIC FUNCTION VALUES SUBSEQUENTLY  
C               USED BY THE CFT ROUTINES.                               
C                                                                       
C     EXTERNALS                                                         
C     ---------                                                         
C     FACT, CFTRIG                                                      
C                                                                       
C***********************************************************************
      DIMENSION IFAX(13),TRIGS(1)                                       
C                                                                       
C     THIS ROUTINE WAS MODIFIED FROM TEMPERTON"S ORIGINAL               
C     BY DAVE FULKER.  IT NO LONGER PRODUCES FACTORS IN ASCENDING       
C     ORDER, AND THERE ARE NONE OF THE ORIGINAL 'MODE' OPTIONS.         
C                                                                       
      CALL FACT(N,IFAX)                                                 
      K = IFAX(1)                                                       
      IF (K .LT. 1 .OR. IFAX(K+1) .GT. 5) IFAX(1) = -99                 
      IF (IFAX(1) .GT. 0 )                                              
     +CALL CFTRIG (N, TRIGS)                                            
      RETURN                                                            
      END                                                               
      SUBROUTINE FACT(N,IFAX)                                           
C     *SUBROUTINE* *FACT  * -ROUTINE CALLED BY CFTFAX                   
C     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4        
      DIMENSION IFAX(13)                                                
      IF (N.GT.1) GO TO 10                                              
      IFAX(1) = 0                                                       
      IF (N.LT.1) IFAX(1) = -99                                         
      RETURN                                                            
   10 NN=N                                                              
      K=1                                                               
C     TEST FOR FACTORS OF 4                                             
   20 IF (MOD(NN,4).NE.0) GO TO 30                                      
      K=K+1                                                             
      IFAX(K)=4                                                         
      NN=NN/4                                                           
      IF (NN.EQ.1) GO TO 80                                             
      GO TO 20                                                          
C     TEST FOR EXTRA FACTOR OF 2                                        
   30 IF (MOD(NN,2).NE.0) GO TO 40                                      
      K=K+1                                                             
      IFAX(K)=2                                                         
      NN=NN/2                                                           
      IF (NN.EQ.1) GO TO 80                                             
C     TEST FOR FACTORS OF 3                                             
   40 IF (MOD(NN,3).NE.0) GO TO 50                                      
      K=K+1                                                             
      IFAX(K)=3                                                         
      NN=NN/3                                                           
      IF (NN.EQ.1) GO TO 80                                             
      GO TO 40                                                          
C     NOW FIND REMAINING FACTORS                                        
   50 L=5                                                               
      MAX = SQRT(FLOAT(NN))                                             
      INC=2                                                             
C     INC ALTERNATELY TAKES ON VALUES 2 AND 4                           
   60 IF (MOD(NN,L).NE.0) GO TO 70                                      
      K=K+1                                                             
      IFAX(K)=L                                                         
      NN=NN/L                                                           
      IF (NN.EQ.1) GO TO 80                                             
      GO TO 60                                                          
   70 IF (L.GT.MAX) GO TO 75                                            
      L=L+INC                                                           
      INC=6-INC                                                         
      GO TO 60                                                          
   75 K = K+1                                                           
      IFAX(K) = NN                                                      
   80 IFAX(1)=K-1                                                       
C     IFAX(1) NOW CONTAINS NUMBER OF FACTORS                            
      RETURN                                                            
      END                                                               
      SUBROUTINE CFTRIG(N,TRIGS)                                        
C     *SUBROUTINE* *CFTRIG* -ROUTINE CALLED BY CFTFAX                   
      DIMENSION TRIGS(*)                                                
      PI=2.0*ASIN(1.0)                                                  
      DEL=(PI+PI)/FLOAT(N)                                              
      L=N+N                                                             
      DO 10 I=1,L,2                                                     
      ANGLE=0.5*FLOAT(I-1)*DEL                                          
      TRIGS(I)=COS(ANGLE)                                               
      TRIGS(I+1)=SIN(ANGLE)                                             
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,M,N,IFAC,LA)  
C     *SUBROUTINE* *VPASSM* -ROUTINE CALLED BY FFT99 AND CFFT99         
      DIMENSION A(*),B(*),C(*),D(*),TRIGS(*)                            
C                                                                       
C     SUBROUTINE "VPASSM" - MULTIPLE VERSION OF "VPASSA"                
C     PERFORMS ONE PASS THROUGH DATA                                    
C     AS PART OF MULTIPLE COMPLEX (INVERSE) FFT ROUTINE                 
C     A IS FIRST REAL INPUT VECTOR                                      
C     B IS FIRST IMAGINARY INPUT VECTOR                                 
C     C IS FIRST REAL OUTPUT VECTOR                                     
C     D IS FIRST IMAGINARY OUTPUT VECTOR                                
C     TRIGS IS PRECALCULATED TABLE OF SINES & COSINES                   
C     INC1 IS ADDRESSING INCREMENT FOR A AND B                          
C     INC2 IS ADDRESSING INCREMENT FOR C AND D                          
C     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S                    
C     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S                    
C     M IS THE NUMBER OF VECTORS                                        
C     N IS LENGTH OF VECTORS                                            
C     IFAC IS CURRENT FACTOR OF N                                       
C     LA IS PRODUCT OF PREVIOUS FACTORS                                 
C                                                                       
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,           
     *     SIN72/0.951056516295154/,COS72/0.309016994374947/,           
     *     SIN60/0.866025403784437/                                     
C                                                                       
      MM=N/IFAC                                                         
      IINK=MM*INC1                                                      
      JINK=LA*INC2                                                      
      JUMP=(IFAC-1)*JINK                                                
      IBASE=0                                                           
      JBASE=0                                                           
      IGO=IFAC-1                                                        
      IF (IGO.GT.4) RETURN                                              
      GO TO (10,50,90,130),IGO                                          
C                                                                       
C     CODING FOR FACTOR 2                                               
C                                                                       
   10 IA=1                                                              
      JA=1                                                              
      IB=IA+IINK                                                        
      JB=JA+JINK                                                        
      DO 20 L=1,LA                                                      
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 15 IJK=1,M
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      D(JA+J)=B(IA+I)+B(IB+I)                                           
      C(JB+J)=A(IA+I)-A(IB+I)                                           
      D(JB+J)=B(IA+I)-B(IB+I)                                           
      I=I+INC3                                                          
      J=J+INC4                                                          
   15 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
   20 CONTINUE                                                          
      IF (LA.EQ.MM) RETURN                                              
      LA1=LA+1                                                          
      JBASE=JBASE+JUMP                                                  
      DO 40 K=LA1,MM,LA                                                 
      KB=K+K-2                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      DO 30 L=1,LA                                                      
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 25 IJK=1,M                                                     
      C(JA+J)=A(IA+I)+A(IB+I)                                           
      D(JA+J)=B(IA+I)+B(IB+I)                                           
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))                 
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))                 
      I=I+INC3                                                          
      J=J+INC4                                                          
   25 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
   30 CONTINUE                                                          
      JBASE=JBASE+JUMP                                                  
   40 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C     CODING FOR FACTOR 3                                               
C                                                                       
   50 IA=1                                                              
      JA=1                                                              
      IB=IA+IINK                                                        
      JB=JA+JINK                                                        
      IC=IB+IINK                                                        
      JC=JB+JINK                                                        
      DO 60 L=1,LA                                                      
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 55 IJK=1,M                                                     
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))                                 
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))                                 
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))) 
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))) 
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))) 
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))) 
      I=I+INC3                                                          
      J=J+INC4                                                          
   55 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
   60 CONTINUE                                                          
      IF (LA.EQ.MM) RETURN                                              
      LA1=LA+1                                                          
      JBASE=JBASE+JUMP                                                  
      DO 80 K=LA1,MM,LA                                                 
      KB=K+K-2                                                          
      KC=KB+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      DO 70 L=1,LA                                                      
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 65 IJK=1,M                                                     
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))                                 
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))                                 
      C(JB+J)=                                                          
     *    C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=                                                          
     *    S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))
     *   +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=                                                          
     *    C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=                                                          
     *    S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))
     *   +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3                                                          
      J=J+INC4                                                          
   65 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
   70 CONTINUE                                                          
      JBASE=JBASE+JUMP                                                  
   80 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C     CODING FOR FACTOR 4                                               
C                                                                       
   90 IA=1                                                              
      JA=1                                                              
      IB=IA+IINK                                                        
      JB=JA+JINK                                                        
      IC=IB+IINK                                                        
      JC=JB+JINK                                                        
      ID=IC+IINK                                                        
      JD=JC+JINK                                                        
      DO 100 L=1,LA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 95 IJK=1,M                                                     
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))                       
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))                       
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))                       
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))                       
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))                       
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))                       
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))                       
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))                       
      I=I+INC3                                                          
      J=J+INC4                                                          
   95 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
  100 CONTINUE                                                          
      IF (LA.EQ.MM) RETURN                                              
      LA1=LA+1                                                          
      JBASE=JBASE+JUMP                                                  
      DO 120 K=LA1,MM,LA                                                
      KB=K+K-2                                                          
      KC=KB+KB                                                          
      KD=KC+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      C3=TRIGS(KD+1)                                                    
      S3=TRIGS(KD+2)                                                    
      DO 110 L=1,LA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 105 IJK=1,M                                                    
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))                       
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))                       
      C(JC+J)=                                                          
     *    C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))                      
     *   -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))                      
      D(JC+J)=                                                          
     *    S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))                      
     *   +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))                      
      C(JB+J)=                                                          
     *    C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))                      
     *   -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))                      
      D(JB+J)=                                                          
     *    S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))                      
     *   +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))                      
      C(JD+J)=                                                          
     *    C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))                      
     *   -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))                      
      D(JD+J)=                                                          
     *    S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))                      
     *   +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))                      
      I=I+INC3                                                          
      J=J+INC4                                                          
  105 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
  110 CONTINUE                                                          
      JBASE=JBASE+JUMP                                                  
  120 CONTINUE                                                          
      RETURN                                                            
C                                                                       
C     CODING FOR FACTOR 5                                               
C                                                                       
  130 IA=1                                                              
      JA=1                                                              
      IB=IA+IINK                                                        
      JB=JA+JINK                                                        
      IC=IB+IINK                                                        
      JC=JB+JINK                                                        
      ID=IC+IINK                                                        
      JD=JC+JINK                                                        
      IE=ID+IINK                                                        
      JE=JD+JINK                                                        
      DO 140 L=1,LA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 135 IJK=1,M                                                    
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))               
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))               
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) 
     *  -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))              
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) 
     *  +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))              
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) 
     *  +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))              
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) 
     *  -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))              
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) 
     *  -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))              
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) 
     *  +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))              
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) 
     *  +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))              
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) 
     *  -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))              
      I=I+INC3                                                          
      J=J+INC4                                                          
  135 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
  140 CONTINUE                                                          
      IF (LA.EQ.MM) RETURN                                              
      LA1=LA+1                                                          
      JBASE=JBASE+JUMP                                                  
      DO 160 K=LA1,MM,LA                                                
      KB=K+K-2                                                          
      KC=KB+KB                                                          
      KD=KC+KB                                                          
      KE=KD+KB                                                          
      C1=TRIGS(KB+1)                                                    
      S1=TRIGS(KB+2)                                                    
      C2=TRIGS(KC+1)                                                    
      S2=TRIGS(KC+2)                                                    
      C3=TRIGS(KD+1)                                                    
      S3=TRIGS(KD+2)                                                    
      C4=TRIGS(KE+1)                                                    
      S4=TRIGS(KE+2)                                                    
      DO 150 L=1,LA                                                     
      I=IBASE                                                           
      J=JBASE                                                           
CVD$  NODEPCHK                                                          
      DO 145 IJK=1,M                                                    
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))               
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))               
      C(JB+J)=                                                          
     *    C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) 
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))         
     *   -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) 
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))         
      D(JB+J)=                                                          
     *    S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) 
     *      -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))         
     *   +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) 
     *      +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))         
      C(JE+J)=                                                          
     *    C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) 
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))         
     *   -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) 
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))         
      D(JE+J)=                                                          
     *    S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I))) 
     *      +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))         
     *   +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I))) 
     *      -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))         
      C(JC+J)=                                                          
     *    C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) 
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))         
     *   -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) 
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))         
      D(JC+J)=                                                          
     *    S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) 
     *      -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))         
     *   +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) 
     *      +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))         
      C(JD+J)=                                                          
     *    C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) 
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))         
     *   -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) 
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))         
      D(JD+J)=                                                          
     *    S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I))) 
     *      +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))         
     *   +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I))) 
     *      -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))         
      I=I+INC3                                                          
      J=J+INC4                                                          
  145 CONTINUE                                                          
      IBASE=IBASE+INC1                                                  
      JBASE=JBASE+INC2                                                  
  150 CONTINUE                                                          
      JBASE=JBASE+JUMP                                                  
  160 CONTINUE                                                          
      RETURN                                                            
      END                                                               
