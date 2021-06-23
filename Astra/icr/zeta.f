      SUBROUTINE ZETA(X,ZF,ZD,ZW)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  Plasma Dispersion function and first derivative for real
C  argument, adapted for use in RAYIC.
C
C         ZF        = -x*Z(x)
C         ZD        = x*x*Z (x)
C         ZW        = ZF - ZD
C (with the present compiler an apostrophe in a comment is an error!
C  where Z(x) is the real part of the Fried & Conte function.
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DATA  ZERR /1.E-9/
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
         X2     = X*X
C
      IF(ABS(X).LE.5.7)  THEN
C
C  Power series
C
         SUM    = 1.
         TERME  = 1.
         FN     = 0.
C
   10    FN     = FN + 1.
         TERME  = TERME*X2/FN
         TERM   = TERME/(FN+FN + 1.)
         SUM    = SUM + TERM
C
      IF(ABS(TERM/SUM).GT.ZERR)  GO TO 10
C
         ZF     = 2.*X2*EXP(-X2)*SUM
         ZD     = -2.*X2*(1. - ZF)
         ZW     = ZF - ZD
C
      ELSE
C
C  Asymptotic series
C
         K      = 0
         TERMZ  = 1.E0
         FMULT  = 0.5E0/X2
         ZD     = 0.
         ZW     = 0.
C
   20    K      = K+1
         FACZ   = FLOAT(2*K-1)*FMULT
C
      IF(FACZ.GT.1.)  GOTO 30
C
         TERMZ  = TERMZ*FACZ
         TERMD  = FLOAT(2*K+1)*TERMZ
         ZD     = ZD + TERMD
         ZW     = ZW + TERMZ - TERMD
C
      IF(K.LE.10 .AND. ABS(TERMD).GT.ZERR)  GO TO 20
C
   30    ZD     = ZD + 1.
         ZF     = ZD*FMULT + 1.
C
      END IF
C
      RETURN
      END

      SUBROUTINE IBESSN(Z,FIBN,FIBNM1,FIBNP1,IORDER)
C
C ------------------------------------------------------------------
C
C         MODIFIED BESSEL FUNCTIONS   EXP(-Z)*IB(N,Z)
C
C         Z       (REAL, INPUT) - ARGUMENT
C         FIBN    (REAL, OUTPUT) - FUNCTION OF ORDER IORDER
C         FIBNM1  (REAL, OUTPUT) - FUNCTION OF ORDER IORDER-1
C         FIBNP1  (REAL, OUTPUT) - FUNCTION OF ORDER IORDER+1
C
C ------------------------------------------------------------------
C
      IF(IORDER.LE.0)  IORDER=1
C
         ZHF = Z/2.
         ZTN = 1.
C
      DO 20  N=1,IORDER
         ZTNM1 = ZTN
   20    ZTN = ZTN*ZHF/FLOAT(N)
         ZTNP1 = ZTN*ZHF/FLOAT(IORDER+1)
C
         Z2 = ZHF*ZHF
C
         ZEXP = EXP(-Z)
         FIBN = ZTN*ZEXP*ZIBSUM(Z2,IORDER)
         FIBNM1 = ZTNM1*ZEXP*ZIBSUM(Z2,IORDER-1)
         FIBNP1 = ZTNP1*ZEXP*ZIBSUM(Z2,IORDER+1)
C
      RETURN
      END
C
      FUNCTION ZIBSUM(Z2,N)
C
         ZIBSUM = 1.
         ZTERM = 1.
         M = 0
C
   10    M = M+1
         ZTERM = ZTERM*Z2/FLOAT(M*(M+N))
         ZIBSUM = ZIBSUM + ZTERM
      IF(ZTERM/ZIBSUM.GT.1.E-9)  GOTO 10
C
      RETURN
      END
