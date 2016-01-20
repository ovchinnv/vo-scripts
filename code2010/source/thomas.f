      SUBROUTINE F04EAF(N,D,DU,DL,B,IFAIL)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14A REVISED. IER-688 (DEC 1989).
C
C     F04EAF SOLVES THE EQUATIONS
C
C     T*X = B ,
C
C     WHERE T IS AN N BY N TRIDIAGONAL MATRIX, BY GAUSSIAN ELIMINATION
C     WITH PARTIAL PIVOTING.
C
C     FOR A DESCRIPTION OF THE PARAMETERS AND USE OF THIS ROUTINE SEE
C     THE NAG LIBRARY MANUAL.
C
C     -- WRITTEN ON 14-JANUARY-1983.  S.J.HAMMARLING.
C
C     NAG FORTRAN 66 BLACK BOX ROUTINE.
C
C     .. Parameters ..
      CHARACTER*6       SRNAME
      PARAMETER         (SRNAME='F04EAF')
C     .. Scalar Arguments ..
      INTEGER           IFAIL, N
C     .. Array Arguments ..
      DOUBLE PRECISION  B(*), D(*), DL(*), DU(*)
C     .. Local Scalars ..
      DOUBLE PRECISION  MULT, TEMP, ZERO
      INTEGER           K, KK
C     .. Local Arrays ..
      CHARACTER*1       P01REC(1)
C     .. External Functions ..
      INTEGER           P01ABF
      EXTERNAL          P01ABF
C     .. Intrinsic Functions ..
      INTRINSIC         ABS
C     .. Data statements ..
      DATA              ZERO/0.0e+0/
C     .. Executable Statements ..
C
      IF (N.GT.0) GO TO 20
      IFAIL = P01ABF(IFAIL,1,SRNAME,0,P01REC)
      RETURN
   20 CONTINUE
C
      IF (N.EQ.1) GO TO 140
      DO 120 K = 2, N
         IF (DL(K).NE.ZERO) GO TO 40
         IF (D(K-1).EQ.ZERO) GO TO 220
         GO TO 100
   40    IF (ABS(D(K-1)).LT.ABS(DL(K))) GO TO 60
         MULT = DL(K)/D(K-1)
         D(K) = D(K) - MULT*DU(K)
         B(K) = B(K) - MULT*B(K-1)
         IF (K.LT.N) DL(K) = ZERO
         GO TO 100
   60    CONTINUE
         MULT = D(K-1)/DL(K)
         D(K-1) = DL(K)
         TEMP = D(K)
         D(K) = DU(K) - MULT*TEMP
         IF (K.EQ.N) GO TO 80
         DL(K) = DU(K+1)
         DU(K+1) = -MULT*DL(K)
   80    CONTINUE
         DU(K) = TEMP
         TEMP = B(K-1)
         B(K-1) = B(K)
         B(K) = TEMP - MULT*B(K)
  100    CONTINUE
  120 CONTINUE
  140 CONTINUE
C
      IF (D(N).EQ.ZERO) GO TO 200
      B(N) = B(N)/D(N)
      IF (N.GT.1) B(N-1) = (B(N-1)-DU(N)*B(N))/D(N-1)
      IF (N.LE.2) GO TO 180
      K = N - 2
      DO 160 KK = 3, N
         B(K) = (B(K)-DU(K+1)*B(K+1)-DL(K+1)*B(K+2))/D(K)
         K = K - 1
  160 CONTINUE
  180 CONTINUE
C
      IFAIL = 0
      RETURN
C
  200 K = N + 1
  220 IFAIL = P01ABF(IFAIL,K,SRNAME,0,P01REC)
      RETURN
C
C     END OF F04EAF.
C
      END
C
C
      INTEGER FUNCTION P01ABF(IFAIL,IERROR,SRNAME,NREC,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C     MARK 13 REVISED. IER-621 (APR 1988).
C     MARK 13B REVISED. IER-668 (AUG 1988).
C
C     P01ABF is the error-handling routine for the NAG Library.
C
C     P01ABF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are 
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ABF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ABF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 20 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   20       CONTINUE
            IF (IFAIL.NE.-13) THEN
               WRITE (MESS,FMT=99999) SRNAME, IERROR
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ABF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': IFAIL',
     *  ' =',I6)
      END
C
C
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
C
C
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/0/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
C
C
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
