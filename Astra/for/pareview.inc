C		NCURVS max No. of radial profiles stored
C		NRMAX max length of rad. arrays, max No. of rad. slices 
C		NTMAX  max length of time arrays, max No. of time slices 
	integer	  NRMAX, NTMAX, NCURVS
	parameter(NRMAX=512, NCURVS=NRMAX*(NRW+7))
C	parameter(NTMAX=NRD*(4*NRD+7)/(NRW+20))
	parameter(NTMAX=NRD*(4*NRD+7)/(NRW+2))
