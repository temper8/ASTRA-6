C======================================================================|
	subroutine	TSCTRL(ARR1,ARR2,ARR3,DELMAX)
C----------------------------------------------------------------------|
C The subroutine enables user's timestep control:
C    if a relative variation of ARR per time step is larger than DELMAX
C    then the time step TAU is reduced to the minimal value TAUMIN
C
C Input:
C	ARR1(NA1),ARR2(NA1),ARR3(NA1),  
C	TAUMIN,TAUMAX,TAU,TAUINC,DELMAX,DTOUT,DPOUT
C Output:
C	TAU
C Example of usage for a stiff diffusivity model.
C (1) strict control:
C	CAR1=grad(TI)-grad(TI)_critical;	CF1=1;
C	TSCTRL(CAR1,CAR1,CAR1,CF1):;	! One array only is checked
C (2) soft control:
C	TSCTRL(XI,HE,DN,CF1):;
C
C Warning:
C	Only one call from a model is allowed!
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	double precision CTAU,DELMAX,ARR1(NA1),ARR2(NA1),ARR3(NA1)
	double precision ARR1O(NRD),ARR2O(NRD),ARR3O(NRD),YDA
	integer j,ICALL,jj(3)
	save	ARR1O,ARR2O,ARR3O,ICALL
	data	ICALL/0/jj/3*0/
	if (DELMAX .le. 0.)	return
	if (ICALL  .eq. 0)	goto	1
	CTAU = 1./TAUINC
	do	j = 1,NA
           YDA  = max(1.d-6,abs(ARR1O(j)+ARR1(j)))
           CTAU = MAX(CTAU,ABS(ARR1O(j)-ARR1(j))/YDA/DELMAX)
           YDA  = max(1.d-6,abs(ARR2O(j)+ARR2(j)))
           CTAU = MAX(CTAU,ABS(ARR2O(j)-ARR2(j))/YDA/DELMAX)
           YDA  = max(1.d-6,abs(ARR3O(j)+ARR3(j)))
           CTAU = MAX(CTAU,ABS(ARR3O(j)-ARR3(j))/YDA/DELMAX)
	enddo
	TAU = MIN(TAU/CTAU,TAU)
	TAU = MAX(TAUMIN,TAU)
 1	ICALL = 1
	do	j=1,NA
		ARR1O(j) = ARR1(j)
		ARR2O(j) = ARR2(j)
		ARR3O(j) = ARR3(j)
	enddo
	end
C======================================================================|
