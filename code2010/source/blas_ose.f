c     DO NOT MODIFY!     
      SUBROUTINE BLASIUS_ose(Y, U, V, DUDY, NY, NPP, NP, NU, UINF, X, XHAT, DSHAT)
C     
C///////////////////////////////////////////////////////////////////////
C     
C     Author        : Andrea Pascarelli
c     modification and specialization: Victor Ovchinnikov 
C     Date          : January 12, 1999
C     Last update   : March 12 2003
C     Info          :       ***Purpose***
C     
C     Boundary layer velocity profiles: Blasius solution.
C     Input Variables:
c     NY -- max y position at which velocities are desired
c     NPP -- dimension of grid
c     NP  -- dimension of velocity components   
c     XHAT -- reference x value at which dshat (reference delta star) is specified
C     X    -- streamwise location at which solution is desired
c     UINF -- U at INFINITY
C     Output variables:
C     
C     Local variables:
C     
C     External Subroutines:
c     
C     PLEASE NOTE: U AND V ARE COMPUTED AT THE CELL CORNER;  DUDY AT THE CELL CENTER
c     this version solves the original equation (with the 1/2 factor)
C///////////////////////////////////////////////////////////////////////
C     
C.......................................................................
C     Declaration of formal arguments and local variables
C.......................................................................
C     ...Parameters...
      INTEGER ITEMAX, NETA
      REAL EPS, ETAMAX
      PARAMETER (ITEMAX=50, EPS=1.0E-14, ETAMAX = 10.0, NETA = 20000)
C     ...Scalar Arguments...
      INTEGER  NY
      REAL DSHAT, NU, UINF, X, XHAT
C     ...Array Arguments...
      REAL   Y(1:NPP), U(1:NP),  V(1:NP), DUDY(1:NP)
C     ...Intrinsic Functions...
      INTRINSIC REAL, SQRT
C     ...Local Scalars...
      INTEGER COUNTER, I, J, J1, JCOUNT, JSTAT
      REAL ALPHA, DETA, DF1N, DSTAR, ETAMIN, HSHAP
      REAL FGUESS, F0, G0, H0, I0, L, DL
      REAL RATIO, THETA, UJ, UJP, X0, YJ, YJP, YSTAT, VJ, VJP,
     &     DUDYJ, DUDYJP 
C     ...Local Arrays...
      REAL ETA(1:NETA), F(1:NETA), F1(1:NETA), F2(1:NETA)
C.......................................................................
C     Executable Statements
C.......................................................................
      ETAMIN = 0.0
      DETA   = (ETAMAX - ETAMIN)/ REAL(NETA-1)
C.......................................................................
C     Boundary conditions @ eta = 0
C.......................................................................
      F (1)   = 0.0
      F1(1)   = 0.0
      ETA(1)  = ETAMIN
C.......................................................................
C     Guess the value of f''(0) (greather than 0)
C.......................................................................
      FGUESS = 0.33206
      COUNTER = 0
C.......................................................................
C     Repeat until construct: 
C     integrate Blasius equation from 0 to etamax 
C.......................................................................
 100  CONTINUE
c     WRITE(UNIT=*,FMT=1) ' Integration #', COUNTER
c     1   FORMAT(/,1X, A, I2)
C     
      F2(1) = FGUESS
      DO 10 I = 2, NETA
         ETA(I) = ETA(I-1)+ DETA
         F(I)   = F(I-1)  + F1(I-1)*DETA
         F1(I)  = F1(I-1) + F2(I-1)*DETA
         F2(I)  = F2(I-1) - 0.5*F(I-1)*F2(I-1)*DETA 
 10   CONTINUE
C.......................................................................
C     Check to see whether f'(infty) = 1 is fulfilled, and if not,
C     improve the guess for f''(0).
C.......................................................................
      ALPHA  = 1.0/SQRT(F1(NETA))
      DF1N   = ABS(1.0-F1(NETA))
      FGUESS = FGUESS*ALPHA
      COUNTER = COUNTER + 1
      IF (.NOT.((COUNTER.GT.ITEMAX).OR.(DF1N.LT.EPS))) GOTO 100
c      write(6,*) 'input:',fguess
C     Writing   eta, f , f' & f''
c     WRITE(UNIT=*,FMT=*) (I, ETA(I), F(I), F1(I), F2(I), I = 1, NETA) 
C.......................................................................
C     Compute delta*, theta and shape factor
C.......................................................................
      DSTAR = 0.0
      THETA = 0.0
      DO 20 I = 2, NETA
         F0 = 1.0 - F1(I-1)
         G0 = 1.0 - F1(I)
         H0 = (1.0 - F1(I-1))*F1(I-1)
         I0 = (1.0 - F1(I-1))*F1(I-1) 
         DSTAR = DSTAR  + 0.5*(ETA(I)-ETA(I-1))*(F0 + G0)
         THETA = THETA  + 0.5*(ETA(I)-ETA(I-1))*(H0 + I0)
 20   CONTINUE
c     print*,dstar,theta
      X0   = XHAT - UINF*((DSHAT/DSTAR)**2)/NU
c     write(6,*) 'x0=',x0
      L    = SQRT(NU*(X-X0)/UINF)
      DL   = 0.5*SQRT(NU/(UINF*(X-X0)))
      DSTAR = L * DSTAR 
      THETA = L * THETA 
      HSHAP = DSTAR/THETA 
C     
c     print*, ' '
c      print*, ' dstar    ', DSTAR 
c     print*, ' '
c      print*, ' theta    ', THETA 
c     print*, ' '
c     print*, ' Hshape   ', HSHAP
c     print*, ' '
c     print*, ' Re_dstar ', DSTAR*UINF/NU
c     print*, ' '
c     print*, ' Re_theta ', THETA*UINF/NU
c     print*, ' '
C.......................................................................
C     Interpolate selfsimilar profiles to mesh y(j)
C.......................................................................
C.......................................................................
C     Break iteration to locate the right wall-normal indeces 
C.......................................................................
      JCOUNT = 1
      J1     = 1
      DO 30 JSTAT = 2, NY
         DO 40 J = J1, NETA-1
C.......................................................................
C     Convert eta ---> y
C.......................................................................
            YJ  = ETA(J)
            YJP = ETA(J+1)
            YSTAT = 0.5*(Y(JSTAT)+Y(JSTAT-1))
            IF ( (YSTAT.GE.YJ) .AND. (YSTAT.LT.YJP) )  THEN
C.......................................................................
C     Convert normalized quantities  ---> dimensional quantities
C.......................................................................
               DUDYJ = UINF*F2(J)!/L
               DUDYJP= UINF*F2(J+1)!/L
               RATIO = (YSTAT - YJ)/(YJP - YJ)
C.......................................................................
C     Interpolate u velocities
C.......................................................................
               DUDY(JSTAT) = DUDYJ + RATIO*(DUDYJP -DUDYJ)
c     print*, JCOUNT, J, YJ, YSTAT, YJP
               JCOUNT = JCOUNT + 1
               J1     = J
               GOTO 30
            END IF
 40      CONTINUE
         J1 = NETA
         DUDY(JSTAT) = 0.
 30   CONTINUE
C     
c     finished with u velocity; note that it is known at cell centers
c     -- therefore the 1st derivative will be second order at the cell corner      
c     
      JCOUNT = 1
      J1     = 1
      DO 31 JSTAT = 1, NY
         DO 41 J = J1, NETA-1
C.......................................................................
C     Convert eta ---> y
C.......................................................................
            YJ  = ETA(J)
            YJP = ETA(J+1)
            YSTAT = Y(JSTAT)
            IF ( (YSTAT.GE.YJ) .AND. (YSTAT.LT.YJP) )  THEN
C.......................................................................
C     Convert normalized quantities  ---> dimensional quantities
C.......................................................................
               VJ  = UINF*DL* (ETA(J)*F1(J)-F(J))
               VJP = UINF*DL* (ETA(J+1)*F1(J+1)-F(J+1))
               UJ  = UINF*F1(J)
               UJP = UINF*F1(J+1)
               RATIO = (YSTAT - YJ)/(YJP - YJ)
C.......................................................................
C     Interpolate v velocities and dU/dY
C.......................................................................
               V(JSTAT) = VJ + RATIO*(VJP - VJ)
               U(JSTAT) = UJ + RATIO*(UJP - UJ)
c     print*, JCOUNT, J, YJ, YSTAT, YJP
               JCOUNT = JCOUNT + 1
               J1     = J
               GOTO 31
            END IF
 41      CONTINUE
         J1 = NETA
         V(JSTAT) = UINF*DL*(ETA(neta)*F1(neta)-F(neta))
         U(JSTAT)= UINF
 31   CONTINUE
C     
      RETURN
      END
C
C///////////////////////////////////////////////////////////////////////

