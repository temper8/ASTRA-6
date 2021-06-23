C=======================================================================
C  NOTE:  All functions in this file should be decribed as
C         double precision and have a name no longer than 6 characters
C-----------------------------------------------------------------------
C    Modules:
C	RFA	RFAN	XFA	XFAN	AFR	AFX	FRMAX	FRMIN
C	RFMIN	RFMAX	RFVAL	AFVAL	RFVEX	AFVEX	RFVIN	AFVIN
C	RECR	GAUSS	ASTEP	RSTEP	XSTEP	GRAD	STEP	VINT
C	IINT	RADIAL	RADINT	ATR	ATX	CUT	FRAMP	FJUMP
C	TIMDER	TIMINT	TIMAVG	FIXVAL	FTAV	FTMIN	FTMAX
C=======================================================================
	double precision function RFA(YA)
C-----------------------------------------------------------------------
C rho=f(a) "rho" at a given radius "a"
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	YA,DRHODA,QUADIN
	integer	j
	RFA = QUADIN(NA1,AMETR,RHO,YA,DRHODA,j)
	end
C=======================================================================
	double precision function RFAN(YAN)
C-----------------------------------------------------------------------
C rho=f(a) "rho" at a given radius "a"
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YAN,RFA
	integer	j
	RFAN = RFA(YAN*ABC)
	end
C=======================================================================
	double precision function XFA(YA)
C-----------------------------------------------------------------------
C x=f(a) "x"=rho/roc at a given radius "a"
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YA,RFA
	integer	j
	XFA = RFA(YA)/ROC
	end
C=======================================================================
	double precision function XFAN(YAN)
C-----------------------------------------------------------------------
C rho=f(a) "rho" at a given radius "a"
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YAN,XFA
	integer	j
	XFAN = XFA(YAN*ABC)
	end
C=======================================================================
	double precision function AFR(YR)
C-----------------------------------------------------------------------
C a=f(rho) "a" at a given radius "rho"
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	YR,DADRHO,QUADIN
	integer	j
	AFR = QUADIN(NA1,RHO,AMETR,YR,DADRHO,j)
	end
C=======================================================================
	double precision function AFX(YX)
C-----------------------------------------------------------------------
C a=f(rho/roc) "a" at a given radius "x"=rho/roc
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YX,AFR
	AFX = AFR(YX*ROC)
	end
C=======================================================================
C FRMAX []: max_r(astra_array). Maximal value of argument array on [1,NA1]
C			(Pereverzev 17-JUL-97)
C	Usage:	...=FRMAX(TE);	qmin_1/FRMAX(MU);
C
	double precision function FRMAX(YA)
	implicit none
	double precision	YA(*)
	integer j
	include	'for/parameter.inc'
	include	'for/const.inc'
	FRMAX=-1.E37
	do 1 J=1,NA1
 1	FRMAX=max(FRMAX,YA(j))
	end
C=======================================================================
C FRMIN []: min_r(astra_array). Minimal value of argument array on [1,NA1]
C			(Pereverzev 17-JUL-97)
C	Usage:	...=FRMIN(TE);	qmax_1/FRMIN(MU);
C
	double precision function FRMIN(YA)
	implicit none
	double precision	YA(*)
	integer j
	include	'for/parameter.inc'
	include	'for/const.inc'
	FRMIN=1.E37
	do 1 J=1,NA1
 1	FRMIN=min(FRMIN,YA(j))
	end
C=======================================================================
C RFMIN [m]: Position of a minimal value of array on 0<=rho/ROC<=1
C			(Pereverzev 17-JUL-97)
C	Usage:	CAR1=NUES;	...=RFMIN(CAR1);
C		rmnH_RFMIN(HE);
C		Further exmples see in RFMAX
C
	double precision function RFMIN(YA)
	implicit none
	integer j
	include  'for/parameter.inc'
	include  'for/const.inc'
	include	'for/status.inc'
	double precision	YMIN,YA(*)
	YMIN = 1.E37
	RFMIN = 0.
	do 1 J=1,NA1
	if (YA(j).gt.YMIN)	goto	1
	YMIN = YA(j)
	if (j.ne.1)	RFMIN = RHO(j)
 1	continue
	end
C=======================================================================
C RFMAX [m]: Position of a minimal value of array on 0<=rho/ROC<=1
C			(Pereverzev 17-JUL-97)
C	Usage:	CAR1=NUES;	...=RFMAX(CAR1);
C		rmnH_RFMAX(MU);		returns position on "rho" in [m]
C		rmnH_RFMAX(MU)/ROC;	returns position on "rho" normalized
C		rmnH_AMETR(RFMAX(MU));	returns position on "a" in [m]
C
C		Two expressions:	MU(RFMAX(MU))
C				and	FRMAX(MU)
C				return the same result.
C
C	Using in FORTRAN:	RFMAX(TE)   or   RFMAX(TE(1))
C
	double precision function RFMAX(YA)
	implicit none
	integer j
	include  'for/parameter.inc'
	include  'for/const.inc'
	include	'for/status.inc'
	double precision	YMAX,YA(*)
	YMAX = -1.E37
	RFMAX = 0.
	do 1 J=1,NA1
	if (YA(j).lt.YMAX)	goto	1
	YMAX = YA(j)
	if (j.ne.1)	RFMAX = RHO(j)
 1	continue
	end
C=======================================================================
C RFVAL [m]: Outermost radial (0<=rho<=ROC) position where
C		 array YA(*) takes a given value (YVAL) 
C			(Pereverzev 17-MAR-98)
C	Usage:	(1) CAR1=NUES;	       ! Coordinate of a point where NUES=0.1
C		    ...=RFVAL(CAR1,.1);
C	     (2) rq2_RFVAL(MU,.5)/ROC; ! ------------------------ MU=0.5
C	     (3) RFVAL(AMETR,a0)  ! Recalculates "a0" in "rho"
	double precision function RFVAL(YA,YVAL)
	implicit none
	integer j
	include  'for/parameter.inc'
	include  'for/const.inc'
	include	'for/status.inc'
	double precision	YVAL,YA(*),YA1,YA2
	RFVAL = 0.
	YA1 = YA(1)-YVAL
	do	j=2,NA1
	   YA2 = YA(j)-YVAL
	   if (YA1*YA2 .gt. 0.)				goto	1
	   YA1 = abs(YA(j-1))+abs(YA(j))
	   if (YA1 .eq. 0.)				goto	1
	   if (1000.*abs(YA(j-1)-YA(j)) .lt. YA1)	then
	      RFVAL = RHO(j)
	      goto	1
	   endif
	   RFVAL = (RHO(j-1)*YA(j)-RHO(j)*YA(j-1)+
     >		   YVAL*(RHO(j)-RHO(j-1)))/(YA(j)-YA(j-1))
 1	   continue
	   YA1 = YA2
	enddo
	end
C=======================================================================
C AFVAL [m]: Outermost radial (0<=a<=ABC) position where
C		 array YA(*) takes a given value (YVAL) 
C			(Pereverzev 17-MAR-98)
C	Usage:	CAR1=NUES;	...=AFVAL(CAR1,.1);
C		aq2_AFVAL(MU,.5)/ROC;
C
	double precision function AFVAL(YA,YVAL)
	implicit none
	integer j
	include  'for/parameter.inc'
	include  'for/const.inc'
	include	'for/status.inc'
	double precision	YVAL,YA(*),YA1,YA2
	AFVAL = 0.
	YA1 = YA(1)-YVAL
	do	j=2,NA1
	   YA2 = YA(j)-YVAL
	   if (YA1*YA2 .gt. 0.)				goto	1
	   YA1 = abs(YA(j-1))+abs(YA(j))
	   if (YA1 .eq. 0.)				goto	1
	   if (1000.*abs(YA(j-1)-YA(j)) .lt. YA1)	then
	      AFVAL = AMETR(j)
	      goto	1
	   endif
	   AFVAL = (AMETR(j-1)*YA(j)-AMETR(j)*YA(j-1)+
     >		   YVAL*(AMETR(j)-AMETR(j-1)))/(YA(j)-YA(j-1))
 1	   continue
	   YA1 = YA2
	enddo
	end
C=======================================================================
C RFVEX [m]: Outermost radial (0<=rho<=ROC) position where
C		 array YA(*) takes a given value (YVAL) 
C			(Pereverzev 17-MAR-98)
C	Usage:	CAR1=NUES;	...=RFVEX(CAR1,.1);
C		rq2_RFVEX(MU,.5)/ROC;
C
	double precision function RFVEX(YA,YVAL)
	implicit none
	double precision	YVAL,YA(*),RFVAL
	RFVEX = RFVAL(YA,YVAL)
	end
C=======================================================================
C AFVEX [m]: Outermost radial (0<=rho<=ROC) position where
C		 array YA(*) takes a given value (YVAL) 
C			(Pereverzev 17-MAR-98)
C	Usage:	CAR1=NUES;	...=AFVEX(CAR1,.1);
C		rq2_AFVEX(MU,.5)/ROC;
C
	double precision function AFVEX(YA,YVAL)
	implicit none
	double precision	YVAL,YA(*),AFVAL
	AFVEX = AFVAL(YA,YVAL)
	end
C=======================================================================
C RFVIN [m]: Innermost radial (0<=rho<=ROC) position where
C		 array YA(*) takes a given value (YVAL) 
C			(Pereverzev 17-MAR-98)
C	Usage:	CAR1=NUES;	...=RFVIN(CAR1,.1);
C		rq2_RFVIN(MU,.5)/ROC;
C
	double precision function RFVIN(YA,YVAL)
	implicit none
	integer j
	include  'for/parameter.inc'
	include  'for/const.inc'
	include	'for/status.inc'
	double precision	YVAL,YA(*),YA1,YA2
	RFVIN = 0.
	YA1 = YA(NA1)-YVAL
	do 1 j=NA1,2,-1
	   YA2 = YA(j-1)-YVAL
	   if (YA1*YA2 .gt. 0.)				goto	1
	   YA1 = abs(YA(j-1))+abs(YA(j))
	   if (YA1 .eq. 0.)				goto	1
	   if (1000.*abs(YA(j-1)-YA(j)) .lt. YA1)	then
	      RFVIN = RHO(j)
	      goto	1
	   endif
	   RFVIN = (RHO(j-1)*YA(j)-RHO(j)*YA(j-1)+
     >		YVAL*(RHO(j)-RHO(j-1)))/(YA(j)-YA(j-1))
 1	   continue
	   YA1 = YA2
	end
C=======================================================================
C AFVIN [m]: Innermost radial (0<=a<=ABC) position where
C		 array YA(*) takes a given value (YVAL) 
C			(Pereverzev 17-MAR-98)
C	Usage:	CAR1=NUES;	...=AFVIN(CAR1,.1);
C		rq2_AFVIN(MU,.5)/ROC;
C
	double precision function AFVIN(YA,YVAL)
	implicit none
	integer j
	include  'for/parameter.inc'
	include  'for/const.inc'
	include	'for/status.inc'
	double precision	YVAL,YA(*),YA1,YA2
	AFVIN = 0.
	YA1 = YA(NA1)-YVAL
	do 1 j=NA1,2,-1
	   YA2 = YA(j-1)-YVAL
	   if (YA1*YA2 .gt. 0.)				goto	1
	   YA1 = abs(YA(j-1))+abs(YA(j))
	   if (YA1 .eq. 0.)				goto	1
	   if (1000.*abs(YA(j-1)-YA(j)) .lt. YA1)	then
	      AFVIN = AMETR(j)
	      goto	1
	   endif
	   AFVIN = (AMETR(j-1)*YA(j)-AMETR(j)*YA(j-1)+
     >		YVAL*(AMETR(j)-AMETR(j-1)))/(YA(j)-YA(j-1))
 1	   continue
	   YA1 = YA2
	end
C=======================================================================
	double precision function RECR(YZ,N)
C-----------------------------------------------------------------------
C The function translates ECR frequency into a radius of resonance surface
C			for N-th ECR harmonic at z-plane z=YZ
C Input:
C	YZ    [m]	vertical pozition with respect to the chamber centre
C	FECR  [GHz]	ECR frequency 
C	N     [ ]	ECR harmonic number
C	BTOR  [T]	toridal field at the chamber centre
C	RTOR  [m]	major radius of the chamber centre
C Output:
C	RECR  [m]	rho_res
C    Optional output:	a_res, shift, elongation, triangularity
C Example: CF3=RECR(0.,2)/ROC;	PEECR=QECR*GAUSS(CF3,0.1);
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	YA,YR,YZ,YY,RZ2A,QUADIN
	integer	N,j
	FECR = 140.
	YR = 28.*N*BTOR*RTOR/FECR
	YA = RZ2A(YR,YZ,NA1)
C	YD = QUADIN(NA1,AMETR,SHIF,YA,YY,j)
C	YE = QUADIN(NA1,AMETR,ELON,YA,YY,j)
C	YT = QUADIN(NA1,AMETR,TRIA,YA,YY,j)
	RECR = QUADIN(NA1,AMETR,RHO,YA,YY,j)
	end
C=======================================================================
C GAUSS:	Gauss distribution function 
C	Note: Both arguments are dimensionless and the distribution are
C	      given with respect to the 
C	      normalized effective minor radius "rho".
C	      The 3rd argument is added by the ASTRA compiler
C Examples:
C	PE=GAUSS(.1,0.2);	Gaus\GAUSS(0.,.1);
C						(Pereverzev 01-AUG-96)
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    Does not work for multiple calls
	double precision function GAUSS(YX,YW,j)
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer jj,j
	double precision	YPOW,YX,YW,YR
	save	YPOW
	data	YPOW/0./
	if(j.eq.1)	then
	YPOW	= 0.
	do	1	jj = 1,NA1-1
	   YR	= (RHO(jj)/ROC-YX)/YW 
	   YPOW	= YPOW + exp(-YR*YR)*VR(jj)
 1	continue
	YPOW	= YPOW*HRO
	endif

	if(YPOW.le.0.)	then
Can happen when the function is called from "tmp/detvar.inc" before METRIC
	   GAUSS = 0.
	   return
	endif
	GAUSS=exp(-((RHO(j)/ROC-YX)/YW)**2)/YPOW
	end
C=======================================================================
C ASTEP:	Step function with dimensional [m] radial argument "a"
C	Returns 0. for a < YA
C	    and 1. otherwise
C Examples:
C	PE=ASTEP(.1)-ASTEP(0.2);	Step\ASTEP(.1);
C	CV1=XQMINB;	HE=HE*ASTEP(CV1);
C						(Pereverzev 01-AUG-96)
	double precision function ASTEP(YA,j)
	implicit none
	include	'for/parameter.inc'
C	include	'for/const.inc'
	include	'for/status.inc'
	integer j
	double precision    YA
	if(AMETR(j) .lt. YA)	then
	   ASTEP=0.
	else
	   ASTEP=1.
	endif
	end
C=======================================================================
C RSTEP:	Step function with dimensional [m] radial argument "rho"
C	Returns 0. for rho < YR
C	    and 1. otherwise
C Examples:
C	PE=RSTEP(.1)-RSTEP(0.2);	Step\RSTEP(.1);
C						(Pereverzev 01-AUG-96)
	double precision function RSTEP(YR,j)
	implicit none
	include	'for/parameter.inc'
	include	'for/status.inc'
	integer j
	double precision    YR
	if(RHO(j) .lt. YR)	then
	   RSTEP=0.
	else
	   RSTEP=1.
	endif
	end
C=======================================================================
C XSTEP:	Step function with a dimensionless radial argument
C	Returns 0. for rho < YX*ROC
C	    and 1. otherwise
C Examples:
C	PE=XSTEP(.1)-XSTEP(0.2);	Step\XSTEP(.1);
C	CV1=XQMINB;	HE=HE*XSTEP(CV1);
C						(Pereverzev 01-AUG-96)
	double precision function XSTEP(YX,j)
	implicit none
	integer j
	double precision    YX
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	if(RHO(j)/ROC .lt. YX)	then
	   XSTEP=0.
	else
	   XSTEP=1.
	endif
	end
C=======================================================================
C GRAD :	Gradient
C Only a radially dependent array may be the 1st parameter of the function
C Examples:
C    out\GRAD(CAR3)	!Radial profile of gradient CAR3 dCAR3/dRo
C    out_GRAD(CAR3B)	!Gradient CAR3 dCAR3/dRo at the boundary
C    out_GRAD(CAR3C)	!Gradient CAR3 dCAR3/dRo at the center
	double precision function GRAD(Y,j)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	Y(*)
	integer j
	if(j.lt.NA)	then
	   GRAD	=(Y(j+1)-Y(j))/HRO
	else
	   GRAD	=(Y(NA1)-Y(NA))/HROA
	endif
C	write(*,*)j,Y(j),GRAD
	end
C=======================================================================
C GRADS :	Gradient
C Only a radially dependent array may be the 1st parameter of the function
C Examples:
C    out\GRADS(CAR3)	!Radial profile of gradient CAR3 dCAR3/dRo
C    out_GRADS(CAR3B)	!Gradient CAR3 dCAR3/dRo at the boundary
C    out_GRADS(CAR3C)	!Gradient CAR3 dCAR3/dRo at the center
	double precision function GRADS(Y,j)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	Y(*)
	integer j
	if(j.lt.NA)	then
	   GRADS	=(Y(j+1)-Y(j))/HRO
	else
	   GRADS	=(Y(NA1)-Y(NA))/(ROC-NA*HRO)
	endif
C	write(*,*)j,Y(j),GRADS
	end
C======================================================================|
C STEP:	Step function with a dimensionless radial argument
C	Returns 0. for YX < 0
C	    and 1. otherwise
C Examples:
C	PE=STEP(CAR1);	Step\STEP(CAR1);
C						(Pereverzev 17-FEB-99)
	double precision function STEP(YX)
	implicit none
	double precision    YX
	if (YX .lt. 0.)	then
	   STEP=0.
	else
	   STEP=1.
	endif
	end
C=======================================================================
C VINT:	Volume integral {0,R} of any array
C Only a radially dependent array may be the 1st parameter of the function
C Examples:
C    out\Vint(CAR3)	!Radial profile of CAR3 volume integral
C    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
C    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
C			(Yushmanov 26-DEC-90)
	double precision function VINT(ARR,YR)
	implicit none
	double precision	ARR(*)
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	VINT=0.
	do 1 J=1,JK
 1	VINT=VINT+ARR(J)*VR(J)
	VINT=HRO*(VINT-ARR(JK)*YDR)
	end
C======================================================================|
C IINT:	Integral {0,R} of current density
C 	Iint=integral {0,R} (ARR/IPOL**2)dV*IPOL/(GP2*Ro)
C 	Only radial dependent array may be a parameter of the function
C 	Examples:
C    out\Iint(CU)	!Radial profile of toroidal current
C    out_Iint(CD,Ro)    !Toroidal driven current inside {0,Ro}
C    out_Iint(CUB)      !Total toroidal current =Iint(CU,ROC); (=IPL)
C			(Pereverzev 23-OCT-99)
	double precision function IINT(ARR,YR)
	implicit none
	integer J,JK
	double precision	ARR(*),YR,YDR,YA
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	IINT = 0.
	if (YR .le. 0.)	return
	JK = YR/HRO+1.-1.E-4
	if (JK .gt. NA)	JK = NA
	YA = 0.
	do   1	J=1,JK
	IINT = IINT+YA
	YA   = ARR(J)*RHO(J)/(G33(J)*IPOL(J)**3)
	YDR = YR-JK*HRO+HRO
 1	continue
	if (JK .ge. NA)	then
	   YDR = min(YDR,0.5*(HRO+HROA))
	endif
	IINT = GP2*IPOL(JK)*(HRO*IINT+YDR*YA)
	end
C=======================================================================
C NODE	Radial node number nearest to  YR
	integer function NODE(YR)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YR
	if(YR.le.0.)	NODE=1
	NODE=YR/HRO+1
	if (YR .gt. ROC-0.5*HROA)	NODE = NA1
	NODE=min(NA1,NODE)
	end
C=======================================================================
C RADIAL	Interpolation of ARR to the radial position R
C 	 It is assumed that the array ARR(j) is given on the grid RHO(j)
	double precision function RADIAL(ARR,YR)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	ARR(*),YR,YDR
	integer JK,node
	JK = node(YR)
	if (JK.ge.NA1)	RADIAL=ARR(NA1)
	if (JK.le.1)	RADIAL=ARR(1)
	if (JK.ge.NA1 .or. JK.le.1)	return
	if (YR.gt.RHO(NA))	then
	RADIAL=((ARR(NA1)-ARR(NA))*YR-RHO(NA)*ARR(NA1)+ROC*ARR(NA))/HROA
	return
	endif
	YDR=YR-RHO(JK)
	RADIAL=(YDR*ARR(JK+1)+(HRO-YDR)*ARR(JK))/HRO
	end
C=======================================================================
C RADINT	Interpolation of ARR to the radial position R
C 	 It is assumed that the array ARR(j) is given at points j*HRO
	double precision function RADINT(ARR,YR)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	ARR(*),YR,YDR
	integer JK
	JK=YR/HRO
	if(JK.ge.NA1)	RADINT=ARR(NA1)
	if(JK.eq.0)	RADINT=(4.*ARR(1)-ARR(2))/3.*(1.-YR*YR)
	if(JK.lt.0)	RADINT=(4.*ARR(1)-ARR(2))/3.
	if(JK.ge.NA1.or.JK.le.0)	return
	YDR=YR-JK*HRO
	RADINT=(YDR*ARR(JK+1)+(HRO-YDR)*ARR(JK))/HRO
	end
C=======================================================================
C ATR	synonym for RADIAL
	double precision function ATR(ARR,YR)
	implicit none
	double precision	ARR(*),YR,RADIAL
	ATR=RADIAL(ARR,YR)
	end
C=======================================================================
C ATX		Interpolation of ARR to the radial position R
C 	 It is assumed that the array ARR(j) is given on the grid RHO(j)
	double precision function ATX(ARR,YR)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	ARR(*),YR,RADIAL
	integer JK
	ATX=RADIAL(ARR,YR*ROC)
	end
C=======================================================================
C CUT	The function cuts off values of Y out of the range:
C				-X < Y < X
C	Usage in a model:	Name\cut(.005,CAR1)\.005;
C		or		nues\cut(100.,HEXP)\100
	double precision function CUT(X,Y)
	implicit none
	double precision X,Y
	CUT=max(-X,min(X,Y))
	end
C=======================================================================
C FRAMP:	linear ramp from 0 at Time<T1 to 1 at Time>T2
C Example:
C	IPL=.1+.2*framp(0.21,0.27)
C			(Yushmanov 26-DEC-90)
	double precision function FRAMP(T1,T2)
	implicit none
	double precision    T1,T2
	include	'for/parameter.inc'
	include	'for/const.inc'
	if(T1.ge.T2)	then
		write(*,*)' Function FRAMP(t1,t2) error: t1 > t2'
		FRAMP=0.
		return
			endif
	if(TIME.le.T1)	then
		FRAMP=0.
		return
			endif
	if(TIME.ge.T2)	then
		FRAMP=1.
		return
			endif
	FRAMP	=(TIME-T1)/(T2-T1)
	end
C=======================================================================
C FJUMP:	jump from 0 to 1 at Time=T1
C Example:
C	IPL=.1+.2*fjump(.21d0)
C			(Yushmanov 26-DEC-90)
	double precision function FJUMP(T1)
	implicit none
	double precision    T1
	include	'for/parameter.inc'
	include	'for/const.inc'
	if(TIME.le.T1)	then
		FJUMP=0.
			else
		FJUMP=1.
			endif
	end
C=======================================================================
C FTBOX:		Box-like time dependence:
C       FTBOX = | 0   if  time <= T1  or time >= T2
C		| 1   otherwise
C Example:
C	CV3=.2*ftbox(0.21d0,5.d-1)
C
C       CAR3=FBOX(
C			(Pereverzev 26-JAN-06)
	double precision function FTBOX(T1,T2)
	implicit none
	double precision    T1,T2
	include	'for/parameter.inc'
	include	'for/const.inc'
	if(TIME.le.T1 .or. TIME.ge.T2)	then
		FTBOX=0.
			else
		FTBOX=1.
			endif
	end
C=======================================================================
C FXBOX:	Box-like radial (argument x=rho_norm)  dependence:
C       FXBOX = | 0   if  x <= X1  or x >= X2
C		| 1   otherwise
C Example:
C	       CAR3=FXBOX(5.d-1,7.d-1)
C					(Pereverzev 26-JAN-06)
	double precision function FXBOX(X1,X2,j)
	implicit none
	double precision    X1,X2
	integer		    j
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	if(RHO(j).le.X1*ROC .or. RHO(j).ge.X2*ROC)	then
		FXBOX=0.
			else
		FXBOX=1.
			endif
	end
C======================================================================|
C TimDer:	time derivative
C 	Examples:
C	    cv1=TimDer(IPL);		cv2=TIMDER(WTOTB);
C	    dIdt_cv1;			dWdt_cv2;
C Note! the function cannot be used in the time output directly !!!!!!!!
C			(Yushmanov 13-FEB-91)
Changed by Pereverzev 15.10.98
C-----------------------------------------------------------------------
	double precision function TIMDER(Y)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer  NLOC,IY,ITSME,ICALL
	parameter	(NLOC=2200)
	double precision	Y,YO(NLOC),YT(NLOC)
	save	ITSME,ICALL,YO,YT
	data	ITSME/0/	ICALL/0/	YO/NLOC*0./	IY/0/
	call	getid(Y,itsme,IY)
	if (IY .lt. 0)		goto	10		! Error !
	if (IY .gt. NLOC)	goto	11

	if (ICALL .eq. 0.)	then
	   ICALL  = 1
	   TIMDER = 0.
	   YO(IY) = Y
	   YT(IY) = TIME
	   return
	endif
	if (TIME .le. YT(IY))	return		! no time step in between

	TIMDER = (Y-YO(IY))/(TIME-YT(IY))
	YO(IY) = Y
	YT(IY) = TIME

	return
 10	if (IY.eq.-2)	write(*,*)"            Calling from TIMDER"
	if (IY.eq.-1)	goto	11
	TIMDER=0.
	return
 11	write(*,*)' too many time derivatives >',NLOC
	TIMDER=0.
	end
C======================================================================|
C TimInt:	time integral
C 	Examples:
C		cv3=TimInt(QITOTB)
C    		A_cv3	! Total ion energy input
C
C Note! the function cannot be used in the time output directly !!!!!!!!
C			(Yushmanov 13-FEB-91)
Changed by Pereverzev 15.10.98
C-----------------------------------------------------------------------
	double precision function TIMINT(Y)
	implicit none
	integer  NLOC,IY,ITSME,ICALL
	parameter	(NLOC=2200)
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	Y,YO(NLOC),YT(NLOC)
	save	ITSME,YO,YT,ICALL
	data	ITSME/0/	ICALL/0/	YO/NLOC*0./	IY/0/
	call	getid(Y,itsme,IY)
	if (IY .lt. 0)		goto	10		! Error !
	if (IY .gt. NLOC)	goto	11

	if (ICALL .eq. 0.)	then
	   ICALL  = 1
C	   TIMINT = Y*TAU
	   TIMINT = 0.
C	   YO(IY) = TIMINT
	   YT(IY) = TIME
	   return
	endif
C	if (TIME .le. YT(IY))	return		! no time step in between

	TIMINT = YO(IY)+Y*(TIME-YT(IY))
	YO(IY) = TIMINT
	YT(IY) = TIME
	return

 10	TIMINT=0.
	if (IY.eq.-1)	goto	11
	if (IY.eq.-2)	write(*,*)"            Calling from TIMINT"
	return
 11	write(*,*)' too many time integrals >',NLOC
	end
C======================================================================|
	double precision function TIMAVG(Y,YTINT)
C----------------------------------------------------------------------|
C Timavg:	time average
C 	Examples:
C		cv3=Timavg(NNCL,1.)
C    		N_cv3	! Neutral gas density averaged over 1 sec
C					Pereverzev 01.06.02
C----------------------------------------------------------------------|
C   Y      - quantity to be averaged over YTINT sec
C   YTINT  - sampling time 
C----------------------------------------------------------------------|
C   IY     - ID of the input quantity Y
C   YI(IY) - integral of Y
C   YT(IY) - time of the previous calling
C----------------------------------------------------------------------|
	implicit none
	integer  NLOC,IY,ITSME,j
	parameter	(NLOC=2200)
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	Y,YTINT,YI(NLOC),YT(NLOC),YDT,YST
	save	ITSME,YI,YT
	data	ITSME/0/	YI/NLOC*-1.E9/		IY/0/
	call	getid(Y,itsme,IY)
	if (IY .lt. 0)		goto	10	! Error !
	if (IY .gt. NLOC)	goto	11

	if (YI(IY) .lt. -0.9E9)	then		! 1st call for IY
	   YI(IY) = 0.
	   YST = TIME
	   TIMAVG = Y
	   YT(IY) = TIME
	   return
	endif					! all subsequent calls
	j = TIME/YTINT
	YST = j*YTINT
	if (YST .lt. YT(IY))	then		! intermediate call
	   YI(IY) = YI(IY)+Y*(TIME-YT(IY))
	   TIMAVG = YI(IY)/(TIME-YST)
	else					! YTINT starting call
	   YI(IY) = YI(IY)+Y*(YST-YT(IY))
	   TIMAVG = YI(IY)/YTINT
	   YI(IY) = Y*(TIME-YST)
	endif
	YT(IY) = TIME
	return

 10	TIMAVG=0.
	if (IY.eq.-1)	goto	11
	if (IY.eq.-2)	write(*,*)"            Calling from TIMAVG"
	return
 11	write(*,*)' too many time integrals >',NLOC
	end
C======================================================================|
C The function returns the value Yvar(time<=Ytime)
C 	Example:
C	    dli_LINTB-FIXVAL(LINTB,1.)
C	This output produces an increment of l_i(t) with respect to t=1sec
C			(Pereverzev 15-OCT-98)
C
	double precision function FIXVAL(Y,YTIME)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer  NLOC
	parameter	(NLOC=2200)
	integer	IY,ITSME
	double precision	Y,YO(NLOC),YTIME
	save	ITSME,YO
	data	ITSME/0/	IY/0/
C IY is the ID (ordinal number) of "Y"
	call	getid(Y,itsme,IY)
	if (IY .lt. 0)		goto	10		! Error !
	if (IY .gt. NLOC)	goto	11

	if (TIME .le. YTIME)	YO(IY)	= Y
	FIXVAL	= YO(IY)
	return
 10	if (IY.eq.-2)	write(*,*)"            Calling from FIXVAL"
	if (IY.eq.-1)	goto	11
	FIXVAL=0.
	return
 11	write(*,*)' >>> FIXVAL >>> too many calls: >',NLOC
	FIXVAL=0.
	end
C======================================================================|
	double precision function FTAV(Y,YTAV)
C----------------------------------------------------------------------|
C The function returns g_i=g_{i-1}exp(-tau/YTAV)+f_i*[1-exp(-tau/YTAV)]
C 	i.e.
C		g_i=f_i		if	 tau << YTAV
C		g_i=g_{i-1}	if	 tau >> YTAV
C----------------------------------------------------------------------|
C Floating Time AVerage, exponential decay with YTAV
C 
C 	Examples:
C		CV2=FTAV(cv1,CF1)
C		CAR2=FTAV(UPL,3*tau);
C		CAR3=HARL
C		HE=FTAV(CAR3,.01);
C	Note:
C	     Do not use	...=FTAV(UPLB,.1)    (UPLB - ASTRA abbreviation)
C	     Use CV1=UPLB;	CV2=FTAV(CV1,.1)
C
C	     Do not use	...=FTAV(HARL,.1)	(HARL - ASTRA formula)
C	     Use CAR3=HARL;	HE=...+FTAV(CAR3,.01)
C
C	G.W. Pacher (18/01/1994)
Changed by Pereverzev 15.10.98
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer  NLOC
	parameter	(NLOC=2200)
	integer	IY,ITSME,ICALL
	double precision	Y,YO(NLOC),YTAV
	save	ITSME,ICALL,YO
	data	ITSME/0/	ICALL/0/	YO/NLOC*0./	IY/0/
	call	getid(Y,itsme,IY)
C IY is the ID (ordinal number) of "Y"
C itsme is the ID (ordinal number) of the calling function (FTAV)
	if (IY .lt. 0)		goto	10		! Error !
	if (IY .gt. NLOC)	goto	11

	if (ICALL .eq. 0)	then
	    FTAV=Y
	    YO(IY)=FTAV
	    ICALL = 1
	    return
	endif

	if (YTAV .le. .1*TAU) then
	    FTAV = Y				! Return the input value
	 else
	    FTAV=Y+(YO(IY)-Y)*EXP(-TAU/YTAV)	! Return a weighted value
	endif
	YO(IY) = FTAV				! Save the previous value
	return
 10	if (IY.eq.-2)	write(*,*)"            Calling from FTAV"
	if (IY.eq.-1)	goto	11
	FTAV = 0.
	return
 11	write(*,*)' >>> FTAV >>> buffer overflow: >',NLOC
	FTAV = 0.
	end
C======================================================================|
	double precision function FTAV3(Y)
C----------------------------------------------------------------------|
C The function returns an average of the quantity Y over JBASE time steps
C The argument to FTAV3 can be either an ASTRA variable or ASTRA array
C 
C 	Examples:
C		CV2=FTAV3(cv1)
C	Note:
C	     Do not use	...=FTAV3(UPLB)      (UPLB - ASTRA abbreviation)
C	     Use CV1=UPLB;	CV2=FTAV3(CV1)
Created by Pereverzev 8.08.2007
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer  JBASE
	parameter	(JBASE=1000)
	integer	j,JY
	double precision	Y,YY(JBASE),YTIME,YAV
	save	YY, YTIME, YAV, JY
	data	JY/0/ YTIME/-1.d37/ YAV/0.d0/

	if (JY .eq. 0)	 goto	1
	FTAV3 = YAV
	if (YTIME .eq. TIME)	return
 1	continue
	if (jy .lt. JBASE) then
	   JY = JY+1
	   YY(jy) = Y
	   YAV = 0.
	   do j=1,JY
	      YAV = YAV+YY(j)
	   enddo
	   YAV = YAV/JY
	else
	   YAV = YAV+(Y-YY(1))/JY
	   do j=2,JY
	      YY(j-1) = YY(j)
	   enddo
	   YY(jy) = Y
	endif
	FTAV3 = YAV
	YTIME = TIME
C	write(*,'(6F10.6)')(1.d3*YY(j),j=1,5),1.d3*FTAV3
	end
C======================================================================|
	double precision function FTAV2(Y)
C----------------------------------------------------------------------|
C The function returns an average of the quantity Y over JBASE time steps
C The argument to FTAV2 can be either an ASTRA variable or ASTRA array
C 
C 	Examples:
C		CV2=FTAV2(cv1)
C	Note:
C	     Using <dt>_FTAV2(TAU) is allowed although can give a wrong
C		result: It averages over time output times instead of
C		all calculation times.
C	     Do not use	...=FTAV2(UPLB)      (UPLB - ASTRA abbreviation)
C	     Use CV1=UPLB;	CV2=FTAV2(CV1)
Created by Pereverzev 8.08.2007
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer  JBASE
	parameter	(JBASE=100)
	integer	j,JY
	double precision	Y,YY(JBASE),YTIME,YAV
	save	YY, YTIME, YAV, JY
	data	JY/0/ YTIME/-1.d37/ YAV/0.d0/

	if (JY .eq. 0)	 goto	1
	FTAV2 = YAV
	if (YTIME .eq. TIME)	return
 1	continue
	if (jy .lt. JBASE) then
	   JY = JY+1
	   YY(jy) = Y
	   YAV = 0.
	   do j=1,JY
	      YAV = YAV+YY(j)
	   enddo
	   YAV = YAV/JY
	else
	   YAV = YAV+(Y-YY(1))/JY
	   do j=2,JY
	      YY(j-1) = YY(j)
	   enddo
	   YY(jy) = Y
	endif
	FTAV2 = YAV
	YTIME = TIME
C	write(*,'(6F10.6)')(1.d3*YY(j),j=1,5),1.d3*FTAV2
	end
C======================================================================|
	double precision function FTMIN(Y)
C----------------------------------------------------------------------|
C The function returns min_t[Y]
C The argument to FTMIN can be either an ASTRA variable or ASTRA array
C 
C 	Examples:
C		CV2=FTMIN(cv1)
C		CAR2=FTMIN(UPL);
C		CAR3=HARL
C		HE=FTMIN(CAR3);
C	Note:
C	     Do not use	...=FTMIN(UPLB)      (UPLB - ASTRA abbreviation)
C	     Use CV1=UPLB;	CV2=FTMIN(CV1)
C
C	     Do not use	...=FTMIN(HARL)	    (HARL - ASTRA formula)
C	     Use CAR3=HARL;	HE=...+FTMIN(CAR3)
C
Changed by Pereverzev 15.10.98
C----------------------------------------------------------------------|
	implicit none
	integer  NLOC
	parameter	(NLOC=2200)
	integer	IY,ITSME
	double precision	Y,YO(NLOC)
	save	ITSME,YO
	data	ITSME/0/	YO/NLOC*1.E37/		IY/0/

C IY is the ID (ordinal number) of "Y"
	call	getid(Y,itsme,IY)
	if (IY .lt. 0)		goto	10		! Error !
	if (IY .gt. NLOC)	goto	11

	FTMIN=1.E37
	if (YO(IY).le.Y)	then
	   FTMIN = YO(IY)
	   return
	endif
	FTMIN=Y
	YO(IY)=Y
	return

 10	if (IY.eq.-2)	write(*,*)"            Calling from FTMIN"
	if (IY.eq.-1)	goto	11
	FTMIN=0.
	return
 11	write(*,*)' >>> FTMIN >>> too many calls: >',NLOC
	FTMIN=0.
	end
C======================================================================|
	double precision function FTMAX(Y)
C----------------------------------------------------------------------|
C The function returns max_t[Y]
C The argument to FTMAX can be either an ASTRA variable or ASTRA array
C 
C 	Examples:
C		CV2=FTMAX(cv1)
C		CAR2=FTMAX(UPL);
C		CAR3=HARL
C		HE=FTMAX(CAR3);
C	Note:
C	     Do not use	...=FTMAX(UPLB)      (UPLB - ASTRA abbreviation)
C	     Use CV1=UPLB;	CV2=FTMAX(CV1)
C
C	     Do not use	...=FTMAX(HARL)	    (HARL - ASTRA formula)
C	     Use CAR3=HARL;	HE=...+FTMAX(CAR3)
C
C	G.V.Pereverzev	15.10.98
C----------------------------------------------------------------------|
	implicit none
	integer  NLOC
	parameter	(NLOC=2200)
	integer	IY,ITSME
	double precision	Y,YO(NLOC)
	save	ITSME,YO
	data	ITSME/0/	YO/NLOC*-1.E37/		IY/0/
C IY is the ID (ordinal number) of "Y"
	call	getid(Y,itsme,IY)
	if (IY .lt. 0)		goto	10		! Error !
	if (IY .gt. NLOC)	goto	11

	FTMAX=-1.E37
	if (YO(IY).ge.Y)	then
	   FTMAX = YO(IY)
	   return
	endif
	FTMAX=Y
	YO(IY)=Y
	return
 10	if (IY.eq.-2)	write(*,*)"            Calling from FTMAX"
	if (IY.eq.-1)	goto	11
	FTMAX=0.
	return
 11	write(*,*)' >>> FTMAX >>> too many calls: >',NLOC
	FTMAX=0.
	end
C======================================================================|
