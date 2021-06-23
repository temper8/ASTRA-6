	integer	NP,       ITER
	parameter  (NP=101)
	double precision	WSACC
	double precision	
     1		WGBJ,      WGB,       WGB0,      WGBST,     WGBD,
     2	        WGMJ,      WGMJEX,    WSQC,      WBBS0,     WBJ0,
     3		WBR00,     WBR0,      WSA(NP),   WSAA(NP),
     4		WSD1(NP),  WGL(NP),   WDGL(NP),
     5		WDSD1(NP), WD2SD1(NP),WBS1(NP),  WSAC1,
     6		WD2GL(NP), WBS2(NP),  WSAC2,
     7		WSD3(NP),  WDSD3(NP), WD2SD3(NP),WBS3(NP),  WSAC3,
     8		WBK02(NP), WBK0(NP),  WDBK00(NP),WBD02(NP),
     9		WBD0(NP),  WBG02(NP), WBG0(NP),  WBD12(NP),
     &		WBG332(NP),WBG33(NP), WBG222(NP),WBG22(NP)
	double precision
     1		WSL0(NP),  WSV0(NP),  WSU0(NP),  WSW0(NP),
     2		WGMC(NP),  WDGMC(NP), WBA(NP),   WDBA(NP),
     3		WBB(NP),   WDBB(NP),
     4		WSCJ1(NP), WSCI1(NP), WSCJ3(NP), WSCI3(NP),
     5		WSCJ5(NP), WSCI5(NP), WSCJ7(NP), WSCI7(NP),
     6		WBK10(NP), WBK11(NP), WBK20(NP), WBK22(NP),
     7		WBK13(NP), WBK30(NP), WBK31(NP), WBK33(NP),
     8		WSL1(NP),  WSV1(NP),  WSU1(NP),
     9		WSL22(NP), WSL2(NP),  WSV2(NP),  WSU2(NP),
     &		WSL3(NP),  WSV3(NP),  WSU3(NP)
	double precision	
     1		WSP(NP),   WDSP(NP),  WSJP(NP),  WDSJP(NP),
     2		WSJ(NP),   WDSJ(NP),  WSQ(NP),   WDSQ(NP),
     3		WDSQRQ(NP),WBF(NP),   WBFF(NP),  WSJSL(NP),
     4		WSJSR(NP), WGP(NP)

C Total length 76*NP (+21*NP)
	common/EQACC/ WSACC, ITER
	common/EQDCH/ WGBJ,  WGB,   WGB0,  WGBST, WGBD,
     *	              WGMJ,  WGMJEX,WSQC,  WBBS0, WBJ0
	common/EQG0/  WBR00, WBR0,  WSA,   WSAA,  WSD1,  WGL,   WDGL
	common/EQG1/  WDSD1, WD2SD1,WBS1,  WSAC1
	common/EQG2/  WD2GL, WBS2,  WSAC2
	common/EQG3/  WSD3,  WDSD3, WD2SD3,WBS3,  WSAC3
	common/EQGA/  WBK02, WBK0,  WDBK00,WBD02, WBD0,  WBG02, WBG0,
     *		      WBD12, WBG332,WBG33, WBG222,WBG22,
     *		      WSL0,  WSV0,  WSU0,  WSW0
	common/EQEQP/ WGMC,  WDGMC, WBA,   WDBA,  WBB,   WDBB
	common/EQIC0/ WSCJ1, WSCI1, WSCJ3, WSCI3, WSCJ5, WSCI5
	common/EQLIC1/WSCJ7, WSCI7
	common/EQLK1/ WBK10, WBK11
	common/EQLK2/ WBK20, WBK22
	common/EQLK3/ WBK13, WBK30, WBK31, WBK33
	common/EQLL1/ WSL1,  WSV1,  WSU1
	common/EQLL2/ WSL22, WSL2,  WSV2,  WSU2
	common/EQLL3/ WSL3,  WSV3,  WSU3
	common/EQPP/  WSP,   WDSP,  WSJP,  WDSJP, WSJ,   WDSJ,  WSQ,
     *		      WDSQ,  WDSQRQ,WBF,   WBFF,  WSJSL, WSJSR, WGP
