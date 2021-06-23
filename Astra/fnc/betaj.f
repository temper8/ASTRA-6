C======================================================================|
C BETAJR []:	Beta poloidal (r)		(Pereverzev 12-JAN-90)
C
C In this subroutine, the dimensionless \beta_J is defined as
C
C   2*c*c		       2*c*c
C  ------*{Int(p*dS) - pS} = - ------ * Int[S*dp]
C    I*I		        I*I
C
C Equivalent definition:
C
C    c*c           p*dV              c*c         V*dp
C  -------- * {Int[----] - pV} = - ------- * Int[----]
C  \pi*I*I          r              \pi*I*I        r
C
C----------------------------------------------------------------------|
C The same in SI units:
C
C     8\pi                           8\pi
C  ---------*{Int(p*dS) - pS} = - --------- * Int[S*dp]
C  \mu_0*I*I                      \mu_0*I*I
C
C or
C
C      4            p*dV                4           V*dp
C  --------- * {Int[----] - pV} = - --------- * Int[----]
C  \mu_0*I*I         r              \mu_0*I*I        r
C
C----------------------------------------------------------------------|
C Note!
C   The simplified representation S ~= \pi\lambda*a^2 is employed here.
C   It is valid for 3-moment equilibrium only.
C======================================================================|
	double precision function BETAJR(YR)
	implicit  none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	integer   J,JR
	double precision	  YR,YB,YP,YP1
C----------------------------------------------------------------------|
	JR = YR/HRO+0.5
	if (JR .lt. 2)	JR=2
	if (JR .gt. NA)	JR=NA
	YB = 0.
	YP = NE(1)*TE(1)+NI(1)*TI(1)
	do	1	J=1,JR
	   YP1 = NE(J+1)*TE(J+1)+NI(J+1)*TI(J+1)
	   YB = YB+(YP-YP1)*ELON(J)*AMETR(J)**2
	   YP = YP1
 1	continue
	BETAJR = 6.4E-4*GP2*YB*
     .		 (RTOR/(G22(JR)*IPOL(JR)*BTOR*JR*HRO*MU(JR)))**2
	end
C======================================================================|
C----------------------------------------------------------------------|
C Another definition:
C
C     2*c*c		               2*c*c      1 
C  ------------ * {Int(p*dV) - pV} = - ----- * -------- {Int(V*dp)}
C  2\pi*r_c*I*I		                I*I    2\pi*r_c
C
C where "r_c" can be defined in different ways, e.g. r_c = RTOR
C More reasonable definition is r_c = RTOR+SHIF(1)
C
C Similarly, in SI units:
C
C        4                                 4
C  -------------*{Int(p*dV) - pV} = - ------------- * Int[V*dp]
C  \mu_0*r_c*I*I                      \mu_0*r_c*I*I
C
C----------------------------------------------------------------------|
C======================================================================|
