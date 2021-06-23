      SUBROUTINE HPCSD(FUNCT,OUTPUT)
C
C *********************************************************************
C
C     DESCRIPTION OF PARAMETERS
C           PRMT(1)    - LOWER BOUND OF THE INTERVAL (INPUT).
C           PRMT(2)    - UPPER BOUND OF THE INTERVAL (INPUT).
C           PRMT(3)    - INITIAL INCREMENT OF THE INDEPENDENT VARIABLE
C                        (INPUT).
C           PRMT(4)    - UPPER ERROR BOUND (INPUT). IF ABSOLUTE ERROR
C                        IS GREATER THAN PRMT(4), INCREMENT GETS HALVED
C                        IF INCREMENT IS LESS THAN PRMT(3) AND ABSOLUTE
C                        ERROR LESS THAN PRMT(4)/50, INCREMENT GETS
C                        DOUBLED. THE USER MAY CHANGE PRMT(4) BY MEANS
C                        OF HIS OUTPUT SUBROUTINE
C           PRMT(5)    - NO INPUT PARAMETER. SUBROUTINE HPCSD INITIALI-
C                        ZES PRMT(5) TO ZERO. IF THE USER WANTS TO TER-
C                        MINATE SUBROUTINE HPCSD AT ANY OUTPUT POINT HE
C                        HAS TO CHANGE PRMT(5) TO NON-ZERO BY MEANS OF
C                        SUBROUTINE OUTPUT.
C           Y          - INPUT VECTOR OF INITIAL VALUES (DESTROYED).
C                        LATER ON Y IS THE RESULTING VECTOR OF DEPEN-
C                        DENT VARIABLES Y COMPUTED AT THE INTERMEDIATE
C                        POINTS X.
C           DY         - INPUT VECTOR OF ERROR WEIGHTS (DESTROYED),
C                        THE SUM OF ITS COMPONENTS MUST BE EQUAL TO 1.
C                        LATER ON DY IS THE VECTOR OF THE DERIVATIVES
C                        WHICH BELONG TO FUNCTION VALUES Y AT POINT X.
C           HLAST      - OUTPUT CONTAINING THE CURRENT VALUE OF THE
C                        INDEP. VARIABLE STEP
C           NDIM       - AN INPUT VALUE WHICH SPECIFIES THE NUMBER OF
C                        EQUATIONS IN THE SYSTEM.
C           JHLF       - AN OUTPUT VALUE WHICH SPECIFIES THE NUMBER OF
C                        BISECTIONS OF THE INITIAL INCREMENT. IF JHLF
C                        GETS GREATER THAN 10, SUBROUTINE HPCSD RETURNS
C                        WITH ERROR MESSAGE JHLF=11 TO MAIN PROGRAM.
C                        ERROR MESSAGE JHLF=12 OR JHLF=13 APPEARS IN
C                        CASE PRMT(3)=0 OR IN CASE SIGN(PRMT(3)).NE.
C                        SIGN(PRMT(2)-PRMT(1)) RESPECTIVELY.
C           IND        - IF IND >= 0 TEST IS MADE ON ABSOLUTE ERROR
C                        IF IND < 0 TEST IS MADE ON RELATIVE ERROR
C           FUNCT      - THE NAME OF AN EXTERNAL SUBROUTINE USED. IT
C                        COMPUTES THE RIGHT HAND SIDE DY OF THE SYSTEM
C                        TO GIVEN VALUES OF X AND Y. ITS PARAMETER LIST
C                        IS X. FURTHERMORE IT MUST CONTAIN THE LABELED
C                        COMMON/CHPCSD/PRMT(5),Y(*),DY(*),W(16,*),NDIM
C                        ,JHLF. THIS SUBROUTINE SHOULD NOT DESTROY:
C                        X,Y,W,NDIM,JHLF.
C           OUTPUT     - THE NAME OF AN EXTERNAL OUTPUT SUBROUTINE USED
C                        ITS PARAMETER LIST IS X. FURTHERMORE IT MUST
C                        CONTAIN THE LABELED COMMON/CHPCSD/PRMT(5),Y(*)
C                        ,DY(*),W(16,*),NDIM,JHLF. NONE OF THESE PARA-
C                        METERS (EXCEPT PRMT(4), PRMT(5), IF NECESSA-
C                        RY) SHOULD BE CHANGED BY SUBROUTINE OUTPUT. IF
C                        PRMT(5) IS CHANGED TO NON-ZERO, SUBROUTINE
C                        HPCSD IS TERMINATED
C
C
C     REMARKS
C           THE PROCEDURE TERMINATES AND RETURNS TO CALLING PROGRAM IF
C           1) MORE THAN 10 BISECTIONS OF THE INITIAL INCREMENT ARE
C              NECESSARY TO GET SATISFACTORY ACCURACY (ERROR MESSAGE
C              JHLF=11).
C           2) INITIAL INCREMENT IS EQUAL TO ZERO OR HAS WRONG SIGN
C              (ERROR MESSAGES JHLF=12 OR 13 RESPECTIVELY).
C           3) THE WHOLE INTEGRATION INTERVAL IS WORKED THROUGH.
C           4) SUBROUTINE OUTPUT HAS CHANGED PRMT(5) TO NON-ZERO.
C
C
C *********************************************************************
C
      COMMON /CHPCSD/    PRMT(5),        Y(10),          DY(10),
     +                   HLAST,          NDIM,           JHLF,
     +                   IND
C
C *********************************************************************
C
      DIMENSION          CK(40),         W(16,10)
C
C *********************************************************************
C
      DATA  CK / .0                  E0  , .5                  E0  ,
     +           .66666666666667     E0  , .375                E0  ,
     +           .79166666666667     E0  , .20833333333333     E0  ,
     +           .41666666666667     E-1 , .33333333333333     E0  ,
     +           .4                  E0  , .29697760924775     E0  ,
     +           .15875964497103     E0  , .45573725421879     E0  ,
     +           .21810038822592     E0  , .30509651486929     E1  ,
     +           .38328647604670     E1  , .17476028226269     E0  ,
     +           .55148066287873     E0  , .12055355993965     E1  ,
     +           .17118478121952     E0  , .13333333333333     E1  ,
     +           .92561983471074     E0  , .125                E0  ,
     +           .9                  E1  , .3                  E1  ,
     +           .74380165289256     E-1 , .01                 E0  ,
     +           .2                  E-1 , .89629629629630     E1  ,
     +           .33611111111111     E1  , .390625             E-2 ,
     +           .8                  E2  , .135                E3  ,
     +           .4                  E2  , .1171875            E0  ,
     +           .6                  E1  , .2                  E1  ,
     +           .12                 E2  , .108                E3  ,
     +           .234375             E-1 , .18                 E2  /
C
C *********************************************************************
C
            N          = 1
            JHLF       = 0
            X          = PRMT(1)
            H          = PRMT(3)
            PRMT(5)    = CK(1)
      DO 100        K0 = 1,NDIM
            W(16,K0)   = CK(1)
            W(15,K0)   = DY(K0)
            W(1,K0)    = Y(K0)
 100  CONTINUE
                         IF(H*(PRMT(2)-X)) 2,1,3
C
C           ERROR RETURNS
C           -------------
 1          JHLF       = 12
                         GO TO 3
 2          JHLF       = 13
C
C           COMPUTATION OF DY FOR STARTING VALUES
C           -------------------------------------
 3                                       CALL FUNCT(X)
C
C           RECORDING OF STARTING VALUES
C           ----------------------------
                                         CALL OUTPUT(X)
                         IF(PRMT(5)) 5,4,5
 4                       IF(JHLF)    6,6,5
 5    RETURN
C     <><><>
 6    DO 110        K1 = 1,NDIM
            W(8,K1)    = DY(K1)
 110  CONTINUE
C
C           COMPUTATION OF W(2,K2)
C           ----------------------
            JSW        = 1
                         GO TO 21
 7          X          = X + H
      DO 120        K2 = 1,NDIM
            W(2,K2)    = Y(K2)
 120  CONTINUE
C
C           INCREMENT IS TESTED BY MEANS OF BISECTION
C           -----------------------------------------
 8          JHLF       = JHLF + 1
            X          = X - H
      DO 130        K3 = 1,NDIM
           W(4,K3)     = W(2,K3)
 130  CONTINUE
            H          = H*CK(2)
            HLAST = H
            N          = 1
            JSW        = 2
                         GO TO 21
 9          X          = X + H
                                         CALL FUNCT(X)
            N          = 2
      DO 140        K4 = 1,NDIM
            W(2,K4)    = Y(K4)
            W(9,K4)    = DY(K4)
 140  CONTINUE
            JSW        = 3
                         GO TO 21
C
C           COMPUTATION OF TEST VALUE DELTA.
C           --------------------------------
C           Q IS AN AUXILIARY STORAGE LOCATION.
 10         DELTA      = CK(1)
      DO 150        K5 = 1,NDIM
            Q          = Y(K5) - W(4,K5)
C ****** modified 16/07/93
C *old*                  IF(IND.GT.0.AND.Y(K5).NE.CK(1)) Q=Q/Y(K5)
                         IF(IND.GT.0) Q=Q/(1. + ABS(Y(K5)))
C ****** end
                         IF(Q.GT.CK(1)) GO TO 11
            Q          = -Q
 11         DELTA      = DELTA + W(15,K5)*Q
 150  CONTINUE
            DELTA      = DELTA*CK(3)
                         IF(DELTA-PRMT(4)) 14,14,12
 12                      IF(JHLF-10)        8,13,13
C
C           NO SATISFACTORY ACCURACY AFTER 10 BISECTIONS:
C           ---------------------------------------------
C           ERROR MESSAGE
 13         JHLF       = 11
            X          = X + H
                         GO TO 3
C
C           THERE IS SATISFACTORY ACCURACY AFTER < 11 BISECTIONS
C           ----------------------------------------------------
 14         X          = X + H
                                         CALL FUNCT(X)
      DO 160        K6 = 1,NDIM
            W(3,K6)    = Y(K6)
            W(10,K6)   = DY(K6)
 160  CONTINUE
            N          = 3
            JSW        = 4
                         GO TO 21
 15         N          = 1
            X          = X + H
                                         CALL FUNCT(X)
            X          = PRMT(1)
      DO 170        K7 = 1,NDIM
            W(11,K7)   = DY(K7)
            Y(K7)      = W(1,K7) + H*(CK(4)*W(8,K7)+CK(5)*W(9,K7)-
     +                                CK(6)*W(10,K7)+CK(7)*DY(K7))
 170  CONTINUE
 16         X          = X + H
            N          = N + 1
                                         CALL FUNCT(X)
                                         CALL OUTPUT(X)
                         IF(PRMT(5)) 5,17,5
 17                      IF(N-4)    18,22,22
 18   DO 180        K8 = 1,NDIM
            W(N,K8)    = Y(K8)
            W(N+7,K8)  = DY(K8)
 180  CONTINUE
                         IF(N-3) 19,20,22
 19   DO 190        K9 = 1,NDIM
            DELTA      = W(9,K9) + W(9,K9)
            DELTA      = DELTA + DELTA
            Y(K9)      = W(1,K9) + CK(8)*H*(W(8,K9)+DELTA+W(10,K9))
 190  CONTINUE
                         GO TO 16
 20   DO 200        L0 = 1,NDIM
            DELTA      = W(9,L0) + W(10,L0)
            DELTA      = DELTA + DELTA + DELTA
            Y(L0)      = W(1,L0) + CK(4)*H*(W(8,L0)+DELTA+W(11,L0))
 200  CONTINUE
                         GO TO 16
C
C     THE FOLLOWING PART OF THE SUBROUTINE HPCSD COMPUTES BY MEANS OF
C     A RUNGE-KUTTA METHOD STARTING VALUES FOR THE NOT SELF-STARTING
C     PREDICTOR-CORRECTOR METHOD. T IS AN AUXILIARY STORAGE LOCATION.
 21   DO 210        L1 = 1,NDIM
            T          = H*W(N+7,L1)
            W(5,L1)    = T
            Y(L1)      = W(N,L1) + CK(9)*T
 210  CONTINUE
            T          = X + CK(9)*H
                                         CALL FUNCT(T)
      DO 220        L2 = 1,NDIM
            T          = H*DY(L2)
            W(6,L2)    = T
            Y(L2)      = W(N,L2) + CK(10)*W(5,L2) + CK(11)*T
 220  CONTINUE
            T          = X + CK(12)*H
                                         CALL FUNCT(T)
      DO 230        L3 = 1,NDIM
            T          = H*DY(L3)
            W(7,L3)    = T
            Y(L3)      = W(N,L3) + CK(13)*W(5,L3) - CK(14)*W(6,L3) +
     +                             CK(15)*T
 230  CONTINUE
            T          = X + H
                                         CALL FUNCT(T)
      DO 240        L4 = 1,NDIM
            Y(L4)      = W(N,L4) + CK(16)*W(5,L4) - CK(17)*W(6,L4) +
     +                             CK(18)*W(7,L4) + CK(19)*H*DY(L4)
 240  CONTINUE
C
C           POSSIBLE BREAK-POINT FOR LINKAGE
C           --------------------------------
                                         GO TO (7,9,10,15,41), JSW
C
C     STARTING VALUES ARE COMPUTED. NOW STARTS HAMMING-S MODIFIED
C     -----------------------------------------------------------
C     PREDICTOR-CORRECTOR METHOD.
C     ---------------------------
 22         JSTEP      = 3
 23                      IF(N-8) 25,24,25
C
C           N=8 CAUSES THE ROWS OF W TO CHANGE THEIR STORAGE LOCATIONS
C           ----------------------------------------------------------
 24   DO 260        L6 = 2,7
      DO 250        L5 = 1,NDIM
            W(L6-1,L5) = W(L6,L5)
            W(L6+6,L5) = W(L6+7,L5)
 250  CONTINUE
 260  CONTINUE
            N          = 7
C
C           N<8 CAUSES N+1 TO GET N
C           -----------------------
 25         N          = N + 1
C
C           COMPUTATION OF NEXT VECTOR Y
C           ----------------------------
      DO 270        L7 = 1,NDIM
            W(N-1,L7)  = Y(L7)
            W(N+6,L7)  = DY(L7)
 270  CONTINUE
            X          = X + H
 26         JSTEP      = JSTEP + 1
      DO 280        L8 = 1,NDIM
            DELTA      = W(N-4,L8) + CK(20)*H*(W(N+6,L8)+W(N+6,L8)-
     +                   W(N+5,L8)+W(N+4,L8)+W(N+4,L8))
            Y(L8)      = DELTA - CK(21)*W(16,L8)
            W(16,L8)   = DELTA
 280  CONTINUE
C
C           PREDICTOR IS NOW GENERATED IN ROW 16 OF W, MODIFIED PREDIC-
C           TOR IS GENERATED IN Y. DELTA MEANS AN AUXILIARY LOCATION.
C           THE DERIVATIVES OF THE MODIFIED PREDICTOR WILL NOW BE GENE-
C           RATED IN DY.
                                         CALL FUNCT(X)
      DO 290        L9 = 1,NDIM
            DELTA      = CK(22)*(CK(23)*W(N-1,L9)-W(N-3,L9)+CK(24)*H*
     +                   (DY(L9)+W(N+6,L9)+W(N+6,L9)-W(N+5,L9)))
            W(16,L9)   = W(16,L9) - DELTA
            Y(L9)      = DELTA + CK(25)*W(16,L9)
 290  CONTINUE
C
C           TEST WHETHER H MUST BE HALVED OR DOUBLED
C           ----------------------------------------
            DELTA      = CK(1)
      DO 300        M0 =1,NDIM
            Q          = W(16,M0)
C ****** modified 16/07/93
C *old*                  IF(IND.GT.0.AND.Y(M0).NE.CK(1)) Q=Q/Y(M0)
                         IF(IND.GT.0) Q=Q/(1. + ABS(Y(M0)))
C ****** end
                         IF(Q.LT.CK(1)) Q=-Q
            DELTA      = DELTA + W(15,M0)*Q
 300  CONTINUE
            IF(DELTA-PRMT(4)) 27,38,38
C
C           H MUST NOT BE HALVED. THAT MEANS Y(K) ARE GOOD
C           ------------------------------------------
 27                                      CALL FUNCT(X)
                                         CALL OUTPUT(X)
                         IF(PRMT(5)) 5,28,5
 28                      IF(JHLF-11) 30,5,5
 29         QDELTA     = (X-PRMT(2))/H
                         IF(ABS(QDELTA).GT.CK(26)) GOTO 40
      DO 350       L10 = 1,NDIM
            Y(L10)     = Y(L10)+QDELTA*(W(N-1,L10)-Y(L10))
 350  CONTINUE
            X          = X-H*QDELTA
                                         CALL FUNCT(X)
                                         CALL OUTPUT(X)
      RETURN
 30                      IF(H*(X-PRMT(2))) 31,29,29
 31         Q1         = X - PRMT(2)
                         IF(Q1.LT.CK(1)) Q1=-Q1
            Q2         = H
                         IF(Q2.LT.CK(1)) Q2=-Q2
                         IF(Q1-CK(26)*Q2) 29,32,32
 32                      IF(DELTA-CK(27)*PRMT(4)) 33,33,23
C
C           H COULD BE DOUBLED IF ALL NECESSARY PRECEEDING VALUES EXIST
C           -----------------------------------------------------------
 33                      IF(JHLF) 23,23,34
 34                      IF(N-7)  23,35,35
 35                      IF(JSTEP-4) 23,36,36
 36         JMOD       = JSTEP/2
                         IF(JSTEP-JMOD-JMOD) 23,37,23
 37         H          = H + H
            HLAST = H
            JHLF       = JHLF - 1
            JSTEP      = 0
      DO 310        M1 = 1,NDIM
            W(N-1,M1)  = W(N-2,M1)
            W(N-2,M1)  = W(N-4,M1)
            W(N-3,M1)  = W(N-6,M1)
            W(N+6,M1)  = W(N+5,M1)
            W(N+5,M1)  = W(N+3,M1)
            W(N+4,M1)  = W(N+1,M1)
            DELTA      = W(N+6,M1) + W(N+5,M1)
            DELTA      = DELTA + DELTA + DELTA
            W(16,M1)   = CK(28)*(Y(M1)-W(N-3,M1)) -
     +                   CK(29)*(DY(M1)+DELTA+W(N+4,M1))*H
 310  CONTINUE
                         GO TO 23
C
C           H MUST BE HALVED
C           ----------------
 38         JHLF       = JHLF + 1
                         IF(JHLF-10) 39,39,27
 39         H          = CK(2)*H
            HLAST = H
            JSTEP      = 0
      DO 320        M2 = 1,NDIM
            Y(M2)      = CK(30)*(CK(31)*W(N-1,M2)+CK(32)*W(N-2,M2)+
     +                           CK(33)*W(N-3,M2)+W(N-4,M2))-CK(34)*H*
     +                           (W(N+6,M2)-CK(35)*W(N+5,M2)-W(N+4,M2))
            W(N-4,M2)  = CK(30)*(CK(37)*W(N-1,M2)+CK(32)*W(N-2,M2)+
     +                           CK(38)*W(N-3,M2)+W(N-4,M2)) -
     +                   CK(39)*(W(N+6,M2)+CK(40)*W(N+5,M2)-CK(23)*
     +                                         W(N+4,M2))*H
            W(N-3,M2)  = W(N-2,M2)
            W(N+4,M2)  = W(N+5,M2)
 320  CONTINUE
            X          = X - H
            DELTA      = X - (H+H)
                                         CALL FUNCT(DELTA)
      DO 330        M3 = 1,NDIM
            W(N-2,M3)  = Y(M3)
            W(N+5,M3)  = DY(M3)
            Y(M3)      = W(N-4,M3)
 330  CONTINUE
            DELTA      = DELTA - (H+H)
                                         CALL FUNCT(DELTA)
      DO 340        M4 = 1,NDIM
            DELTA      = W(N+5,M4) + W(N+4,M4)
            DELTA      = DELTA + DELTA + DELTA
            W(16,M4)   = CK(28)*(W(N-1,M4)-Y(M4)) -
     +                   CK(29)*(W(N+6,M4)+DELTA+DY(M4))*H
            W(N+3,M4)  = DY(M4)
 340  CONTINUE
                         GO TO 26
 40         JSW        = 5
            H          = PRMT(2)-X
      DO 360       L11 = 1,NDIM
            W(1,L11)   = Y(L11)
            W(8,L11)   = DY(L11)
 360  CONTINUE
            N          = 1
                         GOTO 21
 41         X          = X+H
                                        CALL FUNCT(X)
                                        CALL OUTPUT(X)
      RETURN
C
C *********************************************************************
C
      END
