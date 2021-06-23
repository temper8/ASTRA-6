C BETR  Troyon factor:	beta_toroidal[%]/(I[MA]/a[m]B[T])
C		
C Example:	beTr_BETRB;
C					(Pereverzev 02-MAY-2006)
	double precision function BETRR(YR)
	implicit  none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
        double precision YR,AFR,BETAR,ITOTR
        external AFR,BETAR,ITOTR
	BETRR = BETAR(YR)*AFR(YR)*BTOR/ITOTR(YR)
	return
	end
