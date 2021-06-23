      SUBROUTINE CUBSPL(H,C,NPROF)
C ********** ********** ********** ********** ********** **********
C                                                                 *
C  THIS SUBROUTINE EVALUATES THE COEFFICIENTS FOR                 *
C        THE CUBIC SPLINE INTERPOLATION OF THE PROFILES.          *
C        DERIVATIVES IN R = 0 AND R = RPLASM                      *
C        EVALUATED FROM DIVIDED DIFFERENCES                       *
C                                                                 *
C ********** ********** ********** ********** ********** **********
      DIMENSION C(NPROF,4)
C ---------- ---------- ---------- ---------- ---------- ----------
         N1 = NPROF-1
C ---------- ---------- ---------- ---------- ---------- ----------
C  EVALUATE FINITE DIFFERENCES AND STORE IN C(J,4)
      DO 10  N=1,N1
   10    C(N,4) = (C(N+1,1)-C(N,1))/H
C ---------- ---------- ---------- ---------- ---------- ----------
C  FORWARD EVALUATION OF E(J) (=C(J,2)) AND F(J) (=C(J,3))
         C(1,2) = 0.
         D2 = C(2,4) - C(1,4)
         C(1,3) = C(1,4) - D2/2.E0
      DO 20  N=2,N1
         C(N,2) = -1./(4.+C(N-1,2))
   20    C(N,3) = (C(N-1,3) - 3.*(C(N-1,4)+C(N,4)))*C(N,2)
C ---------- ---------- ---------- ---------- ---------- ----------
C   BACKWARD EVALUATION OF THE DERIVATIVES S(J) (=C(2,J))
C      AND OF HIGHER COEFFICIENTS
         C(NPROF,2) = C(N1,4)
         C(NPROF,3) = 0.
         C(NPROF,4) = 0.
      DO 30 K=1,N1
         N = NPROF-K
         C(N,2) = C(N,2)*C(N+1,2) + C(N,3)
         CH4 = (C(N,2) + C(N+1,2) - 2.*C(N,4))/H
         C(N,3) = 2.*((C(N,4) - C(N,2))/H - CH4)
   30    C(N,4) = 6.*CH4/H
      RETURN
      END

      FUNCTION VALSPL(H,C,I,JDERIV,NPROF)
C ********** ********** ********** ********** ********** **********
C                                                                 *
C  THIS SUBROUTINE EVALUATES THE CUBIC SPLINE VALUE               *
C        OF THE JDERIV-DERIVATIVE OF THE FUNCTION                 *
C        STORED IN THE TABLE C(J,1)                               *
C        (JDERIV.LE.3)                                            *
C                                                                 *
C ********** ********** ********** ********** ********** **********
      DIMENSION C(NPROF,4)
C ---------- ---------- ---------- ---------- ---------- ----------
      IF(I.LE.0 .OR. I.GT.NPROF)  THEN  
         WRITE(6,9000)  I
 9000 FORMAT(' Index for the spline interpolation out of range',I5)
C ====== test begs
         WRITE(6,9010)  JDERIV,H,NPROF
 9010 FORMAT(' deriv. =',I2,' h =',1P,E13.4,0P,' Nprof =',I5)
         WRITE(6,9020) (C(k,1),k=1,NPROF)
 9020 FORMAT(1P,5E13.4,0P)
C ====== test ends
         STOP  
      END IF
C
         V      = 0.
         MMJDR  = 4-JDERIV
         MM = MMJDR
      DO 10  M=1,MMJDR
         N = 5-M
         V = H*V/FLOAT(MM) + C(I,N)
   10    MM = MM - 1 
         VALSPL = V
C
      RETURN
      END
