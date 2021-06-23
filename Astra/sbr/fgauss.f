C======================================================================|
C Gaussian deposition function P=P0*exp(-((r-r0*a)/a/width)^2)
C  r    [m]	- coordinate in the equatorial cross-section 
C  r0   [d/l]	- relative position of the centre of distribution 
C width [d/l]	- relative width, normalized to the plasma radius 
C  a	[m]	- minor radius
C Normalization:	INT(PdV)=1
C
C Usage example:	FGAUSS(CF1,CF2,CAR1):;	PE=...+CF3*CAR1;
C
	subroutine FGAUSS(YCENTR,YWIDTH,YPROF)
	implicit none
	double precision	YCENTR,YWIDTH,YPROF(*),YR,YRO,YWH,YPOW
	integer	j
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	YRO	= YCENTR*ABC
	YWH	= YWIDTH*ABC
	YPOW	= 0.
	do	1	j = 1,NA1
		YR	= (AMETR(j)-YRO)/YWH 
		YPROF(j)= exp(-YR*YR)
		if (j .eq. NA1)	goto 1
		YPOW	= YPOW + YPROF(j)*VR(j)
 1	continue
	YPOW	= YPOW*HRO
	do	2	j = 1,NA1
		YPROF(j)= YPROF(j)/YPOW
 2	continue
	end
