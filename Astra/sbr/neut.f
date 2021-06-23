	subroutine	NEUT
C-----------------------------------------------------------22.01.97---|
C	Input:	ABC,NA1,NA,AMJ,NAB,NNCX
C		AMAIN(j),ZMAIN(j),TE(j),TI(j),NE(j),NI(j),ZEF(j),SNNBM(j)
C		ENCL,ENWM or wall neutral distribution
C		NNCL,NNWM
C	Warning:	ENCL > 0.5;	ENWM =/= 0 if NNWM =/= 0
C	Output:	NN,	TN,    ALBPL
C---------------------------------------------CHANGED BY POLEVOY-------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
C	double precision ABC,AMJ,ENCL,ENWM,NNCL,NNWM,ALBPL
C	integer NA1,NA,NAB,NNCX
	include 'for/status.inc'
	double precision
     1		SCXNI(NRD),TEN(NRD),SRCNN(NRD),YVI(NRD),NIN(NRD),
     2		NN0(NRD),THICKN(NRD),YKERN(NRD,NRD),YKERNF(NRD,NRD),
     3		AKDLT,Y,SVIE,SVREC,SVCX,YNN0,YXI1,YXI2,
     4		YV1,YV2,YHA,YTI,YN0,YZZN,YZ1,YZ2,YZZ,YZZT
	integer	J,JJ,JN
	equivalence	(WORK1(1,1),SCXNI(1)),	(WORK1(1,2),TEN(1)),
     1			(WORK1(1,3),SRCNN(1)),	(WORK1(1,4),YVI(1)),
     2			(WORK1(1,5),NIN(1)),	(WORK1(1,6),NN0(1)),
     3			(WORK1(1,7),THICKN(1)), (WORK1(1,8),YKERN(1,1)),
     4			(WORK1(1,7+NRD),YKERNF(1,1))
C----------------------------------------------------------------------|
	if (ENCL .lt. 0.0001)	then
	   write(*,'(A,F6.3,A)')
     1	      ">>> NEUT >>> Too low energy of incoming neutrals ENCL ="
     2	      ,1000.*ENCL,"eV"
	   write(*,*) "             Setting ENCL = 2 eV"
	   ENCL = 0.002
	endif
	if (ENWM .lt. 0.002)	then
	   write(*,'(A,F6.3,A)')
     1	      ">>> NEUT >>> Too low energy of incoming neutrals ENWM ="
     2	      ,1000.*ENWM,"eV"
	   write(*,*) "             Setting ENWM = 2 eV"
	   ENWM = 0.002
	endif
C Neutral generations
	YHA=ABC/(NA1-.5)
	do	1	J=1,NA1
	include 'fml/svcx'
		SCXNI(J)=SVCX*NI(J)
	include 'fml/svie'
		YKERN(J,4)=SVIE*NE(J)
	include 'fml/svrec'
		SRCNN(J)=SVREC*NE(J)*NI(J)+SNNBM(J)
C V_i is divided by SQRT(3)
		YTI = max(TI(j),1.d-3)
		YVI(J) = 2.52E5*SQRT(YTI/AMAIN(J))
 1	continue
	YVI(NA1) = YVI(NA)
C...	THICKN(1)=0.5*(SCXNI(1)+YKERN(1,4))*YHA
						THICKN(1)=0.
	do	2	J=1,NA
 2	THICKN(J+1)=THICKN(J)+(SCXNI(J)+YKERN(J,4)+
     +	SCXNI(J+1)+YKERN(J+1,4))*YHA*0.5
	do	3	JJ=1,NA1
 	SCXNI(JJ)	=0.5*SCXNI(JJ)/YVI(JJ)
	do	3	J=1,NA1
		YZ1	=(THICKN(J)+THICKN(JJ))/YVI(JJ)
		if (J.GE.JJ)	then
		YZ2	=(THICKN(J)-THICKN(JJ))/YVI(JJ)
				else
		YZ2	=(THICKN(JJ)-THICKN(J))/YVI(JJ)
				endif
	 	YKERN(J,JJ)	=AKDLT(YZ1)+AKDLT(YZ2)
 3		YKERNF(J,JJ)	=YKERN(J,JJ)*SCXNI(JJ)
C-------------------------------	Zero generation:
	YNN0	= NNCL+NNWM+1.E-11
	YXI1	= (NNCL+1.E-11)/YNN0
	YXI2	=NNWM/YNN0
	YV1	=4.37E5*SQRT(ENCL/AMJ)
	YV2	=4.37E5*SQRT(ENWM/AMJ)
	do	J=1,NA1
	   NIN(J) = 0.5*SRCNN(J)/YNN0/YVI(J)
	   YVI(J) = NIN(J)*TI(J)
	enddo
	do	5	J=1,NA1
	   YZZ	=0.5*(NIN(1)*YKERN(J,1)+NIN(NA1)*YKERN(J,NA1))
	   YZZT	=0.5*(YVI(1)*YKERN(J,1)+YVI(NA1)*YKERN(J,NA1))
	   do	JJ=2,NA
	      YZZ	=YZZ+NIN(JJ)*YKERN(J,JJ)
	      YZZT	=YZZT+YVI(JJ)*YKERN(J,JJ)
	   enddo
	   NN0(J)	=YHA*YZZ
	   TEN(J)	=YHA*YZZT
	   YZ1	=(THICKN(NA1)+THICKN(J))/YV1
	   YZ2	=(THICKN(NA1)-THICKN(J))/YV1
	   YN0	=YXI1*(AKDLT(YZ1)+AKDLT(YZ2))
	   NN0(J)	=NN0(J)+YN0
	   TEN(J)	=TEN(J)+ENCL*YN0
	   if (NNWM.LE.0.)	goto	51
	   YZ1	=(THICKN(NA1)+THICKN(J))/YV2
	   YZ2	=(THICKN(NA1)-THICKN(J))/YV2
	   YN0	=YXI2*(AKDLT(YZ1)+AKDLT(YZ2))
	   NN0(J)	=NN0(J)+YN0
	   TEN(J)	=TEN(J)+ENWM*YN0
 51	   continue
 5	continue
C-----------------------	Zero generation density is ready
C-------------------------------	Iterations:
	do	J=1,NA1
	   NN(J)	=NN0(J)
	enddo
	do	10	JN=1,NNCX
	   do	J=1,NA1
	      NIN(J)	=NN0(J)
	   enddo
C	YZ1 = 0.
C	YZ2 = 1.
	   do	10	J=1,NA1
	      YZZ = 0.5*(NIN(1)*YKERNF(J,1)+NIN(NA1)*YKERNF(J,NA1))
	      do       JJ=2,NA
		 YZZ = YZZ+NIN(JJ)*YKERNF(J,JJ)
	      enddo
	      NN0(J) = YHA*YZZ
	      NN(J)  = NN(J)+NN0(J)
C	      YZ1 = max(YZ1,abs(NN0(j)))
C	      YZ2 = min(YZ2,NN(j))
 10	continue
C	CNEUT3 = YZ1
C	CNEUT4 = YZ2
C------------------- End of iterations
	do	J=1,NA1
	   YVI(J)	=TI(J)*NN(J)
	enddo
	do	J=1,NA1
	   YZZ	= 0.5*(YVI(1)*YKERNF(J,1)+YVI(NA1)*YKERNF(J,NA1))
	   YZZN	= 0.5*(NN(1)*YKERNF(J,1)+NN(NA1)*YKERNF(J,NA1))
	   do	JJ=2,NA
	      YZZN = YZZN+NN(JJ)*YKERNF(J,JJ)
	      YZZ  = YZZ+YVI(JJ)*YKERNF(J,JJ)
	   enddo
	   TN(J) = (YZZ*YHA+TEN(J))/NN(J)
	enddo
C albpl: [d/l]	Plasma albedo
C		Pereverzev	15-05-95
C (Neutral_outflux)/(Neutral_influx) =
C	= sqrt{(NN-N1-N2)(NN*TN-N1*E1-N2*E2)}/(N1*sqrt(E1)+N2*sqrt(E2))
	ALBPL=sqrt((NNCL+NNWM)*(NN(NA1)-1)*
     .		((NNCL+NNWM)*NN(NA1)*TN(NA1)-NNCL*ENCL-NNWM*ENWM))/
     .		(NNCL*sqrt(ENCL)+NNWM*sqrt(ENWM))
	if (NA1 .ge. NAB)	return
	do	J=NA1+1,NAB
	   NN(J) = NN(NA1)
	   TN(J) = TN(NA1)
	enddo
	end
C======================================================================|
	double precision	function	AKDLT(X)
	implicit none
	double precision	X
	if (X.GT.30.)	then
		AKDLT	=0.
		RETURN
			else
		AKDLT=EXP(-X)
			endif
	end
C======================================================================|
