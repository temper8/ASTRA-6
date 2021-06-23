C XQMIN []: Position of minimal "q" value on 0<=a/ABC<=1
C			(Pereverzev 13-JAN-97)
C	Usage:	xmin_XQMINB;
C		CV1=XQMINB;	HE=HE*XSTEP(CV1);
C
	double precision FUNCTION XQMINR(YR)
	implicit none
	double precision YR,YMIN
	integer  J
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	YMIN = 0.
	XQMINR = 0.
	do 1 J=1,NA1
	if (MU(j).lt.YMIN)	goto	1
	YMIN = MU(j)
	if (j.ne.1)	XQMINR = AMETR(j)/ABC
 1	continue
	end
