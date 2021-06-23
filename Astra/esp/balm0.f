C
C
C...    BALLOONING AND MERCIER CRITERION
C...    PERIODICITY ASSUMED
C...    2*NPER PERIODS IN POLOIDAL DIRECTION SYMMETRIC AROUND J0
C...    NB > 2*NPER*M
C
c  version  ==18.01.09
        SUBROUTINE BALM0(NPER,J0,INCR,IOPT,EPSOPT,
     &                  KPG,PS0,PS1,PI2,
     &                  N,AP,PF,PP,FR,
     &                  M,UM,VM,UK,VK,
     &                  NAZ,NTZ,RO,
     &                  DMERC,RMERC,RV,SV,
     &                  DINCR,PPOPT,
     &                  W,WT,WJ,WQ,
     &                  NBZ,WB)
C
        INCLUDE'implic.inc'
        DIMENSION AP(*),PF(*),PP(*),FR(*)
        DIMENSION UM(*),VM(*),UK(*),VK(*)
        DIMENSION RO(1:NAZ,1:NTZ)
c        DIMENSION RO(1:NAZ,0:NTZ-1)
        DIMENSION DMERC(*),RMERC(*),RV(*),SV(*)
        DIMENSION DINCR(*),PPOPT(*)
c----------------------------------------------
        DIMENSION W(1:N,1:3),WT(0:M,1:7)
        DIMENSION WJ(1:N,1:M),WQ(1:N,1:M)
        DIMENSION WB(1:NBZ,1:5)
C
        AMAX1(XARG,YARG)=DMAX1(XARG,YARG)
C
        PSM=PS0-PS1
C
        EIT=.125
        QUA=.25
        HAL=.5
C
C...    J*R AT PSI GRIDS TO 2D ARRAY
C...    J/R AT PSI GRIDS TO 2D ARRAY
        DO 601 J=1,M-1
C...        JACOBIAN AT THE CENTERS
          DO 611 I=1,N-1
            UKP=UK(J+1)-UM(J+1)
            UKM=UK(J  )-UM(J  )
            VKP=VK(J+1)-VM(J+1)
            VKM=VK(J  )-VM(J  )
            U1=UM(J  )+RO(I  ,J  )*UKM
            V1=VM(J  )+RO(I  ,J  )*VKM
            U2=UM(J+1)+RO(I  ,J+1)*UKP
            V2=VM(J+1)+RO(I  ,J+1)*VKP
            U3=UM(J  )+RO(I+1,J  )*UKM
            V3=VM(J  )+RO(I+1,J  )*VKM
            U4=UM(J+1)+RO(I+1,J+1)*UKP
            V4=VM(J+1)+RO(I+1,J+1)*VKP           
            U21=U2-U1
            U43=U4-U3
            U31=U3-U1
            U42=U4-U2
            V21=V2-V1
            V43=V4-V3
            V31=V3-V1
            V42=V4-V2
            SC=QUA*(-(U21+U43)*(V31+V42)+(U31+U42)*(V21+V43))
            HPS=-PSM*( AP(I+1)**2-AP(I)**2 )
            UC=QUA*(U1+U2+U3+U4)
            G3C=1.+KPG*(UC-1.)
            DJ=SC/HPS
            W(I,1)=DJ*G3C
            W(I,3)=DJ/G3C
C...        ARCLENGTH ALONG THE RAY
            UL=HAL*(U1+U2)
            VL=HAL*(V1+V2)
            UR=HAL*(U3+U4)
            VR=HAL*(V3+V4)
            W(I,2)=SQRT((UR-UL)**2+(VR-VL)**2)
 611      CONTINUE
C...        TO THE MAGNETIC SURFACES
          DO 612 I=2,N-1
            WJ(I,J)=(W(I-1,2)*W(I,1)+W(I,2)*W(I-1,1))
     &             /(W(I-1,2)+W(I,2))
            WQ(I,J)=(W(I-1,2)*W(I,3)+W(I,2)*W(I-1,3))
     &             /(W(I-1,2)+W(I,2))
 612      CONTINUE
C...        EDGES (QUADRATIC)
          ARC0=0.
          ARC1=ARC0+HAL*W(1,2)
          ARC2=ARC1+HAL*(W(1,2)+W(2,2))
          ARC3=ARC2+HAL*(W(2,2)+W(3,2))
          SPW1=(ARC0-ARC2)*(ARC0-ARC3)
     &        /(ARC1-ARC2)/(ARC1-ARC3)
          SPW2=(ARC0-ARC1)*(ARC0-ARC3)
     &        /(ARC2-ARC1)/(ARC2-ARC3)
          SPW3=(ARC0-ARC1)*(ARC0-ARC2)
     &        /(ARC3-ARC1)/(ARC3-ARC2)
          WJ(1,J)=SPW1*W(1,1)+SPW2*W(2,1)+SPW3*W(3,1)
          WQ(1,J)=SPW1*W(1,3)+SPW2*W(2,3)+SPW3*W(3,3)
C
          ARC0=0.
          ARC1=ARC0+HAL*W(N-1,2)
          ARC2=ARC1+HAL*(W(N-1,2)+W(N-2,2))
          ARC3=ARC2+HAL*(W(N-2,2)+W(N-3,2))
          SPW1=(ARC0-ARC2)*(ARC0-ARC3)
     &        /(ARC1-ARC2)/(ARC1-ARC3)
          SPW2=(ARC0-ARC1)*(ARC0-ARC3)
     &        /(ARC2-ARC1)/(ARC2-ARC3)
          SPW3=(ARC0-ARC1)*(ARC0-ARC2)
     &        /(ARC3-ARC1)/(ARC3-ARC2)
          WJ(N,J)=SPW1*W(N-1,1)+SPW2*W(N-2,1)+SPW3*W(N-3,1)
          WQ(N,J)=SPW1*W(N-1,3)+SPW2*W(N-2,3)+SPW3*W(N-3,3)
C
 601    CONTINUE
C
        DO 500 I=1,N
          WJ(I,M)=WJ(I,1)
          WQ(I,M)=WQ(I,1)
 500    CONTINUE
C
C...    THE LOOP OVER PSI INTERVALS
C
        DO 700 I=1,N-1
C
          HPS=-PSM*(AP(I+1)**2 - AP(I)**2)
          HPS2=HPS**2
          FRC=HAL*(FR(I+1)+FR(I))
          PFC=HAL*(PF(I+1)+PF(I))
          PPC=HAL*(PP(I+1)+PP(I))
C
C...      TETA GRID VALUES
          DO 711 J=1,M-1      
            UKP=UK(J+1)-UM(J+1)
            UKM=UK(J  )-UM(J  )
            VKP=VK(J+1)-VM(J+1)
            VKM=VK(J  )-VM(J  )
            U1=UM(J  )+RO(I  ,J  )*UKM
            V1=VM(J  )+RO(I  ,J  )*VKM
            U2=UM(J+1)+RO(I  ,J+1)*UKP
            V2=VM(J+1)+RO(I  ,J+1)*VKP
            U3=UM(J  )+RO(I+1,J  )*UKM
            V3=VM(J  )+RO(I+1,J  )*VKM
            U4=UM(J+1)+RO(I+1,J+1)*UKP
            V4=VM(J+1)+RO(I+1,J+1)*VKP           
            U21=U2-U1
            U43=U4-U3
            U31=U3-U1
            U42=U4-U2
            V21=V2-V1
            V43=V4-V3
            V31=V3-V1
            V42=V4-V2
            SC=QUA*(-(U21+U43)*(V31+V42)+(U31+U42)*(V21+V43))
            UC=QUA*(U1+U2+U3+U4)
            G3C=1.+KPG*(UC-1.)
C...        G11CONTR
            WT(J,1)=QUA*((U21+U43)**2+(V21+V43)**2)/SC**2*HPS2
C...        J*G12CONTR
            WT(J,3)=-QUA*((U21+U43)*(U31+U42)+(V21+V43)*(V31+V42))/SC
C...        ARCLENGTH ALONG THE SURFACE
            UB=HAL*(U3+U1)
            VB=HAL*(V3+V1)
            UO=HAL*(U4+U2)
            VO=HAL*(V4+V2)
            WT(J,2)=SQRT((UO-UB)**2+(VO-VB)**2)
711       CONTINUE
          WT(0,1)=WT(M-1,1)
          WT(M,1)=WT(1,1)
          WT(0,2)=WT(M-1,2)
          WT(M,2)=WT(1,2)
          WT(0,3)=WT(M-1,3)
          WT(M,3)=WT(1,3)
          DO 712 J=1,M-1
            WT(J,4)=(WT(J-1,2)*WT(J,1)+WT(J,2)*WT(J-1,1))
     &             /(WT(J-1,2)+WT(J,2))
            WT(J,5)=(WT(J-1,2)*WT(J,3)+WT(J,2)*WT(J-1,3))
     &             /(WT(J-1,2)+WT(J,2))
 712      CONTINUE
          WT(0,4)=WT(M-1,4)
          WT(M,4)=WT(1,4)
          WT(0,5)=WT(M-1,5)
          WT(M,5)=WT(1,5)
C
          A3R=0.
          B2A3R=0.
          B2MB1=0.
          B2SA3R=0.
C
          D1=0.
          D2=0.
          D3=0.
          D4=0.
          D5=0.
          D6=0.                                
          D7=0.
          D8=0.
C                                
          DO 701 J=1,M-1
            UKP=UK(J+1)-UM(J+1)
            UKM=UK(J  )-UM(J  )
            VKP=VK(J+1)-VM(J+1)
            VKM=VK(J  )-VM(J  )
            U1=UM(J  )+RO(I  ,J  )*UKM
            V1=VM(J  )+RO(I  ,J  )*VKM
            U2=UM(J+1)+RO(I  ,J+1)*UKP
            V2=VM(J+1)+RO(I  ,J+1)*VKP
            U3=UM(J  )+RO(I+1,J  )*UKM
            V3=VM(J  )+RO(I+1,J  )*VKM
            U4=UM(J+1)+RO(I+1,J+1)*UKP
            V4=VM(J+1)+RO(I+1,J+1)*VKP           
            U21=U2-U1
            U43=U4-U3
            U31=U3-U1
            U42=U4-U2
            V21=V2-V1
            V43=V4-V3
            V31=V3-V1
            V42=V4-V2
            SC=QUA*(-(U21+U43)*(V31+V42)+(U31+U42)*(V21+V43))
            UC=QUA*(U1+U2+U3+U4)
            G3C=1.+KPG*(UC-1.)
            UL=HAL*(U1+U2)
            UR=HAL*(U3+U4)
            G3L=1.+KPG*(UL-1.)
            G3R=1.+KPG*(UR-1.)
            UB=HAL*(U1+U3)
            UO=HAL*(U2+U4)
            G3B=1.+KPG*(UB-1.)
            G3O=1.+KPG*(UO-1.)
C
            RJ=SC*G3C/HPS
            G11K=QUA*((U21+U43)**2+(V21+V43)**2)/SC**2*HPS2
            BS=(G11K+FRC**2)/G3C**2
            RNP=-(FR(I+1)*WQ(I+1,J)-FR(I)*WQ(I,J))/HPS
            DJG12K=-QUA*((U21+U43)*(U31+U42)+(V21+V43)*(V31+V42))/SC
            RNG=-FRC*DJG12K/G3C/G11K
            RNGT=-FRC*(WT(J+1,5)/WT(J+1,4)/G3O-WT(J,5)/WT(J,4)/G3B)
            SH=-(RNP+RNGT)/RJ
            WELL= FRC*SH/BS
     &           -(WJ(I+1,J)-WJ(I,J))/HPS/RJ
     &           -(WT(J+1,5)/WT(J+1,4)*G3O-WT(J,5)/WT(J,4)*G3B)/RJ
     &           +PPC/BS
            BJ=FRC*PPC+BS*PFC/FRC
            DD=BJ/G11K-SH
            BSPS=BS/G11K
            DK=DD*BJ/BS+WELL*PPC
C
            BSRT=G3O**2/(WT(J+1,4)+FRC**2)-G3B**2/(WT(J,4)+FRC**2)
            EINT=-RNP
            EIAD=-RNG
C
C...        PERIODIC COEFFICIENTS FOR BALLOONING EQUATION
C...        WT(J,3) AND WT(J,4) TO BE MULTIPLIED BY PPC
C
            WT(J,5)=EINT
            WT(J,6)=EIAD
            WT(J,7)=RJ**2
C
C...        NB: MINUS TO AVOID NEGATIVE JACOBIAN
            WT(J,1)=-(1./(RJ*G11K))
            WT(J,2)=-(G11K/(RJ*BS))
            WT(J,3)=-(-RJ*WELL)
            WT(J,4)=-(-FRC*BSRT)
C
C...        COEFFICIENTS FOR MERCIER CONSISTENT WITH BALLOONING
            A3R=A3R+BSPS*RJ
            B2A3R=B2A3R-RJ*FRC/G11K
            B2MB1=B2MB1-RJ*WELL+FRC*RJ*SH/BS
            B2SA3R=B2SA3R+FRC**2*RJ/(G11K*BS)
C
            D1=D1+DD*RJ
            D2=D2+SH*RJ
            D3=D3+BSPS*RJ
            D4=D4+DK*RJ
            D5=D5+RJ/G11K
            D6=D6+BS*RJ
            D7=D7+RJ
            D8=D8-RJ*FRC/G3C**2
C
 701      CONTINUE
C
C...      VOLUME
          RV(I+1)=D7*HPS
C...      DQDP/Q
          QPC=-D2
          SV(I+1)=QPC/D8
C
          A3R=A3R/QPC**2
          B2A3R=B2A3R/QPC
C
          DMERC(I+1)=(D1/D2+.5)**2-D3*D4/D2**2
          HVAL=FRC*PPC/D2*(D5-D3*D7/D6)
          RMERC(I+1)=HVAL
C
C...      SHIFT THE COEFFICIENTS TO MAKE J0=1
          DO 7001 K=1,J0-1
            DO 7002 L=1,7
              W1=WT(1,L)
              DO 7003 J=2,M-1
                WT(J-1,L)=WT(J,L)
 7003         CONTINUE
              WT(M-1,L)=W1
 7002       CONTINUE
 7001     CONTINUE
          DO 7000 L=1,7
            WT(0,L)=WT(M-1,L)
            WT(M,L)=WT(1,L)
 7000     CONTINUE
C
C...      MERCIER FROM BALLOONING COEFFICIENTS
          SA3=0.
          SB2MB1=0.
          SB2=0.
          SB1=0.
          SB2A3=0.
          SB2SA3=0.
          CB2L=0.
          CEIL=0.
          RNRM=1./(M-1)
          QPCM=QPC*RNRM
          QPCM2=QPCM**2
          DO 7300 J=1,M-1
C
            CB2R=CB2L-(-WT(J,4))*QPCM
            CB2C=HAL*(CB2L+CB2R)
            CEIR=CEIL+QPCM+WT(J,5)
            CEIC=HAL*(CEIL+CEIR)
C
            CA3=1./WT(J,2)/QPCM2
            CB1=-WT(J,3)-WT(J,4)*(CEIC+WT(J,6))
            CB2A3=CB2C*CA3
            CB2SA3=CB2C**2*CA3
C
            SA3=SA3+CA3
            SB2MB1=SB2MB1+CB2C-CB1
            SB2=SB2+CB2C
            SB1=SB1+CB1
            SB2A3=SB2A3+CB2A3
            SB2SA3=SB2SA3+CB2SA3
C
            CB2L=CB2R  
            CEIL=CEIR
C
 7300     CONTINUE
C
          SA3=SA3*RNRM**2
          SB2A3=SB2A3*RNRM
          DMBALS=(.5-PPC*SB2A3)**2
     &          +PPC*SA3*(SB2MB1-PPC*SB2SA3)
          RMERC(I+1)=DMBALS
C  
C...      INITIAL STABILITY AND STORE NONPERIODIC FOR OPT
C
          PPT=PPC
C
          JBF=NPER*(M-1)+1
          JBB=JBF-1
          NBP=2*NPER*(M-1)+1
C
          DO 7111 JB=1,NBP
            WB(JB,1)=0.
            WB(JB,2)=0.
            WB(JB,3)=0.
 7111     CONTINUE
C
          EIKF0=0.
          EIKB0=0.
          DO 7100 JP=1,NPER
C...        COEFFICIENTS FORWARD IN TETA
            DO 7101 J=1,M-1
C...          INTEGRATE EIKONAL
              EIKF1=EIKF0+WT(J,5)
              EIKC=HAL*(EIKF1+EIKF0)+WT(J,6)
              EIKF0=EIKF1
              C1=WT(J,1)+EIKC**2*WT(J,2)
              C2=HAL*(WT(J,3)+EIKC*WT(J,4))
              C3=HAL*C1*WT(J,7)
              CC=C1+PPT*C2
              WB(JBF,1)=WB(JBF,1)+CC
              WB(JBF+1,1)=WB(JBF+1,1)+CC
              WB(JBF+1,2)=WB(JBF+1,2)-C1
C...          NORM FOR INCREMENT
              WB(JBF,3)=WB(JBF,3)+C3
              WB(JBF+1,3)=WB(JBF+1,3)+C3
C...          STORE FOR OPTIMIZATION
              WB(JBF,4)=C1
              WB(JBF,5)=C2
              JBF=JBF+1
 7101       CONTINUE
C...        COEFFICIENTS BACK IN TETA
            EIK0=0.
            DO 7102 J=M-1,1,-1
C...          INTEGRATE EIKONAL
              EIKB1=EIKB0-WT(J,5)
              EIKC=HAL*(EIKB1+EIKB0)+WT(J,6)
              EIKB0=EIKB1
              C1=WT(J,1)+EIKC**2*WT(J,2)
              C2=HAL*(WT(J,3)+EIKC*WT(J,4))
              C3=HAL*C1*WT(J,7)
              CC=C1+PPT*C2
              WB(JBB,1)=WB(JBB,1)+CC
              WB(JBB+1,1)=WB(JBB+1,1)+CC
              WB(JBB+1,2)=WB(JBB+1,2)-C1
C...          NORM FOR INCREMENT
              WB(JBB,3)=WB(JBB,3)+C3
              WB(JBB+1,3)=WB(JBB+1,3)+C3
C...          STORE FOR OPTIMIZATION
              WB(JBB,4)=C1
              WB(JBB,5)=C2
              JBB=JBB-1
 7102       CONTINUE
 7100     CONTINUE
C...      BOUNDARY DUMMIES
          WB(1,1)=1.
          WB(NBP,1)=1.
          WB(2,2)=0.
          WB(1,3)=1.E-10
          WB(NBP,2)=0.
          WB(NBP,3)=1.E-10
          DO 7200 JB=1,NBP
            WB(JB,2)=WB(JB,2)**2
 7200     CONTINUE
C
C...      MERCIER INDEX
C
          DMBALM=(.5-PPT*B2A3R)**2
     &          +PPT*A3R*(B2MB1-PPT*B2SA3R)
C
C...      BALLOONING INDEX
C
          NI=0
          SI=WB(1,1)
          DO 7410 JB=2,NBP
            SI=WB(JB,1)-WB(JB,2)/SI
            IF(SI.LE.0) NI=NI+1
 7410     CONTINUE
          DINCR(I+1)=-NI
C
C...      INCREMENT BY RATQR
          IF(INCR.EQ.1) THEN
C...        RATQR MATRIX FORMAT
            DO 7401 JB=1,NBP
              WB(JB,1)=WB(JB,1)/WB(JB,3)
 7401       CONTINUE
            DO 7402 JB=1,NBP-1
              WB(JB+1,2)=WB(JB+1,2)/WB(JB,3)/WB(JB+1,3)
 7402       CONTINUE
C
            EPM=1.E-15
            DLA=0.
C...        NO BOUNDARY INVOLVED
            CALL RATQR(NBP-2,DLA,EPM,WB(2,1),WB(2,2))
            DINCR(I+1)=WB(2,1)
C
C...        RESTORE WB(JB,2)
            DO 7413 JB=1,NBP
              WB(JB,2)=0.
 7413       CONTINUE
            DO 7403 JB=1,NBP-1
              C1=WB(JB,4)
              WB(JB+1,2)=WB(JB+1,2)-C1
 7403       CONTINUE
            WB(2,2)=0.
            WB(NBP,2)=0.
            DO 7400 JB=1,NBP
              WB(JB,2)=WB(JB,2)**2
 7400       CONTINUE
          ENDIF
C
C...      TOTAL STABILITY INDEX
          IF(NI.NE.0 .OR. DMBALM.LT.0.) THEN
            ISTAB=-1
          ELSE
            ISTAB=1
          ENDIF
C...      INITIAL STABILITY INDEX
          ISTAB0=ISTAB
          W(I+1,1)=ISTAB0
          W(I+1,2)=ISTAB
C
C...    OPTIMIZATION LOOP BEGINS HERE
C
        PPL=PPC
        HPP=.5*PPC*EPSOPT
        HPP=.5*DABS(PPC*EPSOPT)
        HZE=0.
        KSTAB=ISTAB0+ISTAB0
C
        DO 9999 KOPT=1,IOPT
C
          IF(KSTAB.EQ.0) THEN
C...      INSTABILITY INTERVAL FOUND 
            IF(ISTAB.LT.0) THEN
               PPR=PPT
            ELSE
               PPL=PPT
            ENDIF
            PPT=.5*(PPR+PPL)
          ELSE IF(KSTAB.EQ.-2) THEN
C...      STILL UNSTABLE
            HPP=2.*HPP    
            PPL=AMAX1(PPT-HPP,HZE)
            PPR=PPT
            PPT=.5*(PPR+PPL)
          ELSE IF(KSTAB.EQ. 2) THEN
C...      STILL STABLE    
            HPP=2.*HPP    
            PPL=PPT
            PPR=PPT+HPP
            PPT=.5*(PPR+PPL)
          ENDIF
C
C...      MERCIER INDEX
C
          DMBALM=(.5-PPT*B2A3R)**2
     &          +PPT*A3R*(B2MB1-PPT*B2SA3R)
C
C...      BALLOONING INDEX
C
          DO 7113 JB=1,NBP
            WB(JB,1)=0.
 7113     CONTINUE
C
          DO 7103 JB=1,NBP-1
            CC=WB(JB,4)+PPT*WB(JB,5)
            WB(JB,1)=WB(JB,1)+CC
            WB(JB+1,1)=WB(JB+1,1)+CC
 7103     CONTINUE
C              
C...      BOUNDARY DUMMIES
          WB(1,1)=1.
          WB(NBP,1)=1.
C
          NI=0
          SI=WB(1,1)
          DO 7210 JB=2,NBP
            SI=WB(JB,1)-WB(JB,2)/SI
            IF(SI.LE.0) NI=NI+1
 7210     CONTINUE
C
C...      TOTAL STABILITY INDEX
          IF(NI.NE.0 .OR. DMBALM.LT.0.) THEN
            ISTAB=-1
          ELSE
            ISTAB=1
          ENDIF
          W(I+1,2)=ISTAB
C
C...      REMEMBER ABOUT ISTAB0 UNTIL INTERVAL FOUND
          IF(KSTAB.NE.0) KSTAB=ISTAB+ISTAB0
C
C...      ZERO PRESSURE LEAVE
          IF(PPT.EQ.0.) GO TO 9998
C
C...      CHECK ACCURACY FOR FOUND INTERVAL AND LEAVE
          IF(KSTAB.EQ.0) THEN
            IF(.5*(PPR-PPL)/PPT.LE.EPSOPT) GO TO 9998
          ENDIF            
C
 9999   CONTINUE
C
 9998   CONTINUE
C
        PPOPT(I+1)=PPL
C              
 700    CONTINUE
C
        RV(1)=0.
        DO 801 I=1,N-1
          RV(I+1)=RV(I)+RV(I+1)
 801    CONTINUE
        RVN=RV(N)
        DO 802 I=1,N
          RV(I)=SQRT(RV(I)/RVN)
 802    CONTINUE
C          
        DO 803 I=1,N-1
          HPS=-PSM*(AP(I+1)**2 - AP(I)**2)
          DPRV=HPS/(RV(I+1)-RV(I))
          SV(I+1)=SV(I+1)*HAL*(RV(I+1)+RV(I))*DPRV
 803    CONTINUE
C
        RETURN
        END
C
C
      SUBROUTINE RATQR(N,DLA,EPM,D ,B2)
      INCLUDE'implic.inc'
C      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION D(N),B2(N)
      B2(1)=0.
      S=0.
      Q=0.
      TOT=D(1)
      DO 1 I=1,N
      P=Q
      Q=DSQRT(B2(N-I+1))
      E = D(N-I+1)-P-Q
    1 IF(E.LT.TOT) TOT=E
      DO 2 I=1,N
    2 D(I)=D(I)-TOT
    3 TOT=TOT+S
      EPT=DABS(EPM*TOT)
      IF(DLA.LT.EPT) DLA=EPT
      DE = D(N)-S
      E = B2(N)/DE
      QP=DE+E
      P=1.
      DO 5 I=2,N
      Q=D(N-I+1)-S-E
      R=Q/QP
      P=P*R+1
      EP=E*R
      D(N-I+2)=QP+EP
      DE=Q-EP
      IF(DE.LE.DLA) GOTO6
      E=B2(N-I+1)/Q
      QP=DE+E
      B2(N-I+2)=QP*EP
    5 CONTINUE
      D(1)=QP
       S=QP/P
      IF(TOT+S.GE.TOT) GO TO 3
    6 CONTINUE
      D(1)=TOT
    7 FORMAT(E12.6)
      B2(1)=DABS(DE)
      RETURN
      END
        


                


   
        
