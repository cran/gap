      SUBROUTINE ALMIZE(FX,NALL,PAR,XALL,ITP,ITPSAV,SE,USCORE,BND,
     +   IBND,IHESS,IVAR,IHX,IQOB,IVERB,IQ,IP,H,TOL,TRUPB,
     +   NIT,NFE,PTG,IDG,TL,MAXPAR,MAXITP,WORK,
     +   INLC,NC,FC,NK,CC,IAC,WRKNLC,N,MSG,iop)
C---
C--- INTERFACE TO ALMINI FOR "PATH" PROGRAMS
C---
C---   - INITIALIZES FLAGS AND CONSTANTS
C---   - CALLS ALMINI SUBROUTINE
C---   - INDICATES PARAMETERS SET TO BOUNDS
C---   - RETURNS U-SCORES
C---   - RETURNS CHARACTER STRING CONTAINING TERMINATION MESSAGE
C---
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER   PAR(*)*(*),MSG*(*)
      DIMENSION   XALL(*),ITP(*),ITPSAV(*)
      DIMENSION   SE(*),USCORE(*),BND(MAXPAR,2)
      DIMENSION   TL(MAXITP,MAXITP),WORK(MAXITP,16)
      DIMENSION   CC(NC),IAC(NC),WRKNLC(NC,4)
C
      CHARACTER   IDGTXT(11)*62
C
      DATA IDGTXT /
     +'MAXIMUM POSSIBLE ACCURACY REACHED',
     +'SEARCH DIRECTION NOT DOWNHILL ANYMORE',
     +'ACCUMULATION OF ROUNDING ERRORS PREVENTS FURTHER PROGRESS',
     +'EXCESSIVE CANCELLATION IN OPTIMAL CONDITIONING',
     +'SPECIFIED TOLERANCE ON NORMALIZED GRADIENT WAS MET',
     +'EXCESSIVE CANCELLATION IN GRADIENT',
     +'MAXIMUM NUMBER OF ITERATIONS REACHED',
     +'EXCESSIVE CANCELLATION IN GRADIENT CALCULATION',
     +'NO ITERATED PARAMETERS (OR ALL PARAMETERS SET TO BOUNDS)',
     +'PARAMETER WAS SET TO A BOUND',
     +'UNABLE TO MEET SPECIFIED TOLERANCE ON CONSTRAINT VIOLATION'/
C
C     I N I T I A L I Z A T I O N
C
C START THE TIMER
C
      CALL TIMER(ELAPSE)
C
C INITIALIZATION FOR NON-LINEAR CONSTRAINTS
C
      RATE  = 0.333
      SCL   = 5.0
      EPSCC = .00001
      ISCAL = 0
C
      DO 110 I=1,NC
         IAC(I) = 1
110   CONTINUE
C
C INITIALIZE EQUAL PARAMETERS, SUCH THAT I-TH PARAMETER SET EQUAL
C TO SPECIFIED PARAMETER
C
      DO 120 I=1,NALL
         IF (ITP(I) .LT. 0) XALL(I) = XALL(-ITP(I))
120   CONTINUE
C
C ALLOW FOR SPECIFICATION OF ITERATED PARAMETERS, SUCH THAT I-TH
C PARAMETER IS ITERATED IF ITP(I)=1
C
      N = 0
      DO 130 I=1,NALL
         ITPSAV(I) = ITP(I)
         IF (ITP(I) .EQ. 1) N = N + 1
130   CONTINUE
C
C DEFINE MAXIMUM NUMBER OF ITERATIONS FOR GEMINI
C
      MAXIT = 20 * N
C
C CHECK STARTING VALUES OF ITERATED PARAMETERS AGAINST BOUNDARIES
C
      IBND  = 1
      DO 140 I=1,NALL
         IF (ITP(I) .EQ. 1) THEN
            B1 = BND(I,1) + (H+H)
            B2 = BND(I,2) - (H+H)
            IF (XALL(I) .LE. B1) THEN
               XALL(I) = B1 + H
               WRITE (IQ,1140) PAR(I), XALL(I)
            ELSE IF (XALL(I) .GE. B2) THEN
               XALL(I) = B2 - H
               WRITE (IQ,1140) PAR(I), XALL(I)
            END IF
         END IF
140   CONTINUE
 1140 FORMAT (' PARAMETER ',A,' VIOLATES BOUND -- RESET TO',G17.9)
C
C CALL ALM TO MINIMIZE F(X) SUBJECT TO CONSTRAINTS
C
      CALL ALMINI(FX,NALL,XALL,ITP,SE,BND,
     +   IBND,IHESS,IVAR,IHX,IQOB,IVERB,MAXIT,IP,H,TOL,TRUPB,
     +   NIT,NFE,PTG,IDG,TL,MAXPAR,MAXITP,WORK,
     +   INLC,NC,ISCAL,RATE,SCL,EPSCC,FC,FXC,NK,CC,IAC,WRKNLC,iop)
C
C     O U T P U T
C
200   WRITE (IQ,1200)
 1200 FORMAT (//' ',59('='),' FINAL OUTPUT ',58('='),//)
C
C INDICATE PARAMETERS AT BOUND
C
      DO 210 I=1,NALL
         IF (ITPSAV(I) .NE. ITP(I)) WRITE (IQ,1210) PAR(I)
1210     FORMAT (' PARAMETER ',A,' SET TO A BOUND')
210   CONTINUE
C
C UNPACK U SCORES FROM WORK ARRAY
C
      N = 0
      DO 220 I=1,NALL
         USCORE(I) = 0.0
         IF (ITP(I) .EQ. 1) THEN
            N = N + 1
            USCORE(I) = WORK(N,1)
         END IF
220   CONTINUE
C
C INDICATE REASON FOR TERMINATION
C
      IF (IDG .GT. 10) THEN
         call intpr('UNRECOGNIZED IDG',-1,0,0)
         STOP
      END IF
      MSG = IDGTXT(IDG+1)
C
      RETURN
      END
