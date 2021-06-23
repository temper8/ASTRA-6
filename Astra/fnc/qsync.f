C QSYNC [MW]:	 Integral {0,R} ( PSYNC ) dV
C			(Pereverzev 9-AUG-02)
	double precision function QSYNCR(YR)
	implicit none
	double precision YR,PSYNC,RADIAL,TEAVR,NEAVR
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	PSYNC=1.32E-7*(TEAVR(ROC)*BTOR)**2.5
     .	*SQRT(NEAVR(ROC)/AB*(1.+18.*AB/(RTOR*SQRT(TEAVR(ROC)))))
	QSYNCR=PSYNC*RADIAL(VOLUM,YR)
	end
