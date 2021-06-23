C======================================================================|
	subroutine	MYCONV(ARR,DELMAX)
C----------------------------------------------------------------------|
C The subroutine enables user's timestep control:
C    if a relative variation of ARR per time step is larger than DELMAX
C    then the time step TAU is reduced to the minimal value TAUMIN
C
C Input:
C	ARR(NA1),TAUMIN,TAUMAX,TAU,TAUINC,DELMAX,DTOUT,DPOUT
C Output:
C	TAU
C Example of usage for a stiff diffusivity model.
C (1) strict control:
C	CAR1=grad(TI)-grad(TI)_critical;	CF1=1;
C	MYCONV(CAR1,CF1):;
C (2) soft control:
C	MYCONV(XI,CF1):;
C
C Warning:
C	Only one call from a model is allowed!
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	double precision	CTAU,ARR(NA1),ARRO(NRD),DELMAX
	integer j,ICALL
	save	ARRO,ICALL
	data	ICALL/0/
	if (ICALL .eq. 0)	goto	1
	CTAU = 1./TAUINC
	do	j = 1,NA1
		CTAU = MAX(CTAU,ABS(ARRO(j)/ARR(j)-1.)/DELMAX)
	enddo
	TAU = MIN(TAU/CTAU,TAU)
	TAU = MAX(TAUMIN/TAUINC,TAU)
 1	continue
	do	j=1,NA1
		ARRO(j) = ARR(j)
	enddo
	ICALL = 1
	end
C======================================================================|
