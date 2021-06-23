C QMIN []: min(q). Minimal value of safety factor "q" on [0,a]
C			(Pereverzev 13-JAN-97)
C	Usage:	qmin_QMINB;
C
	double precision FUNCTION QMINR(YR)
	implicit none
	double precision YR
	integer  J
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	QMINR=1.E5
	do 1 J=1,NA1
 1	QMINR=min(QMINR,1./MU(j))
	end
