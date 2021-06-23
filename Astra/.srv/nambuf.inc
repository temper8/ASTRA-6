C----------------------------------------------------------------------|
	include	'for/parameter.inc'
	integer	NCVA,KEYOUT,NOUT,NFINMX
	parameter(NCVA=512, NFINMX=100)
C NCVA - maximum number of Constants, Variables, Internal_variables, Arrays
C        and Funcions (file for/stdfun.f) recognized by the Astra compiler
C	 i.e. for   NCONS, NMAIVA, NSRV, NARR
C 	similar to NCONST in the ASTRA core
	integer		NFIN,NARR,NFML,NFNC,NMAIVA,NSRV,NCONS,
     1			NSTDF,ARXNAM,JARX,Tflag
	character*6	FINNAM,ARRNAM,FMLNAM,CONNAM,FNCNAM,MAIVAR,STDF,
     1			SRVNAM
	common /NAMLEN/ NFIN,NARR,NFML,NFNC,NMAIVA,NSRV,NCONS,NSTDF
	common /NAMBUF/ FINNAM(NFINMX),ARRNAM(NCVA),FMLNAM(500),
     1			CONNAM(NCVA),FNCNAM(200),MAIVAR(NCVA),
     2			SRVNAM(NCVA),STDF(NCVA)
	common /KEYOUT/	KEYOUT,NOUT,ARXNAM(NARRX),JARX,Tflag
C----------------------------------------------------------------------|
