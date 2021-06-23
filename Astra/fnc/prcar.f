C PRCARR [MW*m^3]:	Carbon Radiation power
C       POST D.E. JENSEN R.V. e.a., Atomic Data and Nuclear Tables,
C				    Vol.20 (1977) 397.
C						Pereverzev 29-07-96
C	Usage:		PE=...+PRCAR*NE*NIZ1
	double precision FUNCTION PRCARR(YR)
	implicit none
	double precision YR,YT,YZ,Y
	integer  JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	JK=YR/HRO+.5
	if(JK.GE.NA1)	JK = NA1
	if(JK.le.0)	JK = 1
        YT=TE(JK)*1000.
        YZ=LOG10(TE(JK))
        if(YT.le.13.)                    Y=-14.4*(YZ+2.2)**2
        if(YT.le.25..and.YT.gt.13.)       Y=-1.3-3.67*(YZ+1.9)
        if(YT.le.63..and.YT.gt.25.)       Y=-2.8+10.*(YZ+1.4)**2
        if(YT.le.159..and.YT.gt.63.)      Y=-2.12-7.2*(YZ+1.)**2
        if(YT.le.500..and.YT.gt.159.)     Y=-2.4-1.2*(YZ+.8)
        if(YT.le.4.E3.and.YT.gt.500.)     Y=-3.3+.957*(YZ-.26)**2
        if(YT.gt.4.E3)                   Y=-3.2+0.357*(YZ-0.6)
        PRCARR	=10.**(Y+1.)
	end
