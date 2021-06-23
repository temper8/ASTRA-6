Changed 12.05.2009 (Pereverzev)
	double precision	YR,YDR,YR1
	integer JK,J
	if(YR.le.0.)	return
	if (YR .le. RHO(NA)) 	then
	   JK = YR/HRO+1
	   YR1 =  YR
	else
	   JK = NA
	   YR1 = min(YR,RHO(NA)+.5d0*HROA)
	endif
	YDR=(JK-YR1/HRO)*VR(JK)
