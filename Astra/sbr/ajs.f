C======================================================================|
	subroutine	AJSP(NCH,ICALL)
C----------------------------------------------------------------------|
C   The subroutine is called from IFKEY and writes the binary file
C   ~/astra/runs/runID/a4jsp (unit 1)
C   with account of the control file
C   ~/astra/runs/a2jsp (unit NCH)
C----------------------------------------------------------------------|
	implicit  none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
      double precision ALMHD,BETE,BETI,BETPL,BETR,CCMHD,CCNEU,CCSP
      double precision CCSPX,CERL,CHOTF,CNHH,CNHR,CNSA,COULG,CS,CUOHM
      double precision D2TI,DBOHM,DCHA,DCHH,DCHR,DCHR2,DCKIM,DCKM1
      double precision DCSA,DHKIM,DIDT,DNDIF,DNEXP,EDR,EFLHN,EFLHW
      double precision ENHH0,ENHH1,ENHH2,EPAR,EPL,ETAE,ETAI,ETAN,FA
      double precision FLIN,FOWC,FPA,FPR,FR,FRS,FTE,FX,GAITG,GITG
      double precision GITG0,HAALC,HABM,HABMS,HABOM,HACTE,HAED,HAETI
      double precision HAETS,HAGB,HAGBS,HAGBS1,HAITF,HAITR,HAMM,HANAG
      double precision HANAL,HAPA,HAPUE,HAPYU,HAQ1,HAQ1C,HARL,HARLS
      double precision HARNQ,HARPL,HASCL,HASDL,HATL,HATLI,HATLS,HBJET
      double precision HCHA,HCHGP,HCHH,HCHII,HCHR,HCHR2,HCKIM,HCKM1
      double precision HCSA,HEEFF,HEGN,HETAI,HETIS,HEXP,HGBEJ,HGBIJ
      double precision HIBOM,HIT89,HMHD1,HMHD2,HNASI,HNCHI,HNGSB,HNGSE
      double precision HNGSI,HNGSP,HNPSI,IBM,IBS,ICD,IOHM,ITOT,LNE,LNI
      double precision LNZ1,LTE,LTI,NUE,NUEE,NUES,NUFE,NUI,NUIS,NUPP
      double precision PAION,PAION1,PBEIE,PBICX,PBR1,PBR2,PBR3,PBRAD
      double precision PDT,PDTF,PEDT,PEDT1,PEDTF,PEGN,PEHCL,PEI,PEICL
      double precision PEIGN,PENEU,PENLI,PETSL,PHICL,PICX,PIDT,PIDT1
      double precision PIDTF,PIGN,PINEU,PIONZ,PIREC,PITCX,PITSL,PJOUL
      double precision POH,PRARG,PRBER,PRCAR,PRESE,PRESI,PREST,PRFER
      double precision PRNEO,PRNIT,PROXI,PROXY,PRWOL,PSYNC,PTOT,QBM
      double precision QBREM,QBS,QCD,QDT,QDTF,QEDT,QEDTF,QEGN,QEICL
      double precision QEIGN,QENEU,QETOT,QEX,QIDT,QIDTF,QIGN,QINEU
      double precision QITOT,QIX,QJOUL,QNIND,QNTOT,QNX,QOH,QRAD,QRADX
      double precision QSYNC,QTOT,RLI,RLS,RLTCR,RLTCZ,RLTKD,RLTWN
      double precision RNODE,ROTSH,SHAT,SNNEU,SNNIE,SNNII,SNNR
      double precision SQZ,SVCX,SVCXX,SVD1,SVD2,SVDBH,SVDHE,SVDT,SVIE
      double precision SVIEP,SVIEX,SVII,SVREC,SVRECX,TAU89,TAUNA,TECRJ
      double precision TECRL,TECRL1,TEF,THQ99,TITER,TLQ97,TPF,VDIA
      double precision VPSWW,VRAS,VRHH,VSI,VTE,VTI,WE,WEX,WI,WIX,WTOT
      double precision XCH86,XCHA,XCHH,XCHR,XCHR2,XCKIM,XCKM1,XCSA
      double precision XEXP,XIEFF,XIGN,XIRL,XMHD1,XMHD2,XTEXP,ZIAR
      double precision ZIBE,ZICAR,ZIFER,ZINE,ZINEO,ZINIT,ZIOXI,ZIWOL
      double precision ZZEF
      double precision VINT,IINT,GRAD,FRMAX,FRMIN,RFMIN,RFMAX
      double precision RFVAL,AFVAL,RFVEX,AFVEX,RFVIN,AFVIN,RFA
      double precision RFAN,XFA,XFAN,AFR,AFX,RECR,ATR,ATX
      double precision TIMINT,TIMDER,TIMAVG,GAUSS
      double precision RADIAL,RADINT,ASTEP,RSTEP,XSTEP,STEP,CUT
      double precision FIXVAL,FTAV,FTMIN,FTMAX,FRAMP,FJUMP
      external IINT
      double precision BETAJR,BETBMR,BETRR,FTLLMR,HNASCR,HNCHCR,IBMR
      double precision IBSR,ICDR,IECRR,IFIR,IFWR,IICRR,ILHR,IOHMR,ITOTR
      double precision IXR,LICDR,LINTR,NEAVR,NECHR,NELAR,NEXAVR,PRCARR
      double precision PTOTR,QBRADR,QBTOTR,QDTR,QDTFR,QEDTR,QEDTFR
      double precision QEDWTR,QEGNR,QEICLR,QEIGNR,QENEUR,QETOTR,QEXR
      double precision QICXR,QIDTR,QIDTFR,QIDWTR,QIGNR,QINEUR,QITOTR
      double precision QIXR,QJOULR,QMINR,QNDNTR,QNTOTR,QNXR,QOHR,QRADR
      double precision QRADXR,QSYNCR,QTOTR,SHINER,TAUER,TAUEER
      double precision TAUEIR,TAUGR,TAUIGR,TAUPR,TEAVR,TENDNR,TEXAVR
      double precision TIAVR,TINDNR,TIXAVR,WALFR,WBPOLR,WCR,WCER,WCEXR
      double precision WCIR,WCIXR,WCXR,WER,WEXR,WIR,WIXR,WTOTR,WTOTXR
      double precision XQMINR,ZNDNR,BETAR
      double precision Y,YR,SNNI,QNNEU
	character*4 STRI*132,DATATYPE(264),DATATABLE(264)
	integer	  j,j1,j2,NCH,ICALL,IERR,length,jdt,jn
	double precision YVC,YVE,YY6,YY7,YEPS,YVALP,YLLAME,YLLAMI,YLLAMA
	save	DATATYPE,DATATABLE,jdt
	data	DATATABLE(1)	/'XRHO'/
	data	DATATABLE(2)	/'R   '/
	data	DATATABLE(3)	/'PSI '/
	data	DATATABLE(4)	/'RHO '/
	data	DATATABLE(5)	/'XA  '/
	data	DATATABLE(6)	/'XPSI'/
	data	DATATABLE(7)	/'XPSQ'/
	data	DATATABLE(8)	/'RI  '/
	data	DATATABLE(9)	/'ZU  '/
	data	DATATABLE(10)	/'ZL  '/
	data	DATATABLE(11)	/'RZMA'/
	data	DATATABLE(12)	/'VOL '/
	data	DATATABLE(13)	/'DVEQ'/
	data	DATATABLE(14)	/'SURF'/
	data	DATATABLE(15)	/'AREA'/
	data	DATATABLE(16)	/'LP  '/
	data	DATATABLE(17)	/'RGE '/
	data	DATATABLE(18)	/'FA  '/
	data	DATATABLE(19)	/'ELO '/
	data	DATATABLE(20)	/'BP  '/
	data	DATATABLE(21)	/'POLB'/
	data	DATATABLE(22)	/'EPS '/
	data	DATATABLE(23)	/'A   '/
	data	DATATABLE(24)	/'K   '/
	data	DATATABLE(25)	/'F   '/
	data	DATATABLE(26)	/'B0B '/
	data	DATATABLE(27)	/'BB0 '/
	data	DATATABLE(28)	/'GRHQ'/
	data	DATATABLE(29)	/'NE  '/
	data	DATATABLE(30)	/'NI  '/
	data	DATATABLE(31)	/'NH  '/
	data	DATATABLE(32)	/'NI1 '/
	data	DATATABLE(33)	/'NI2 '/
	data	DATATABLE(34)	/'NIMP'/
	data	DATATABLE(35)	/'N01 '/
	data	DATATABLE(36)	/'N02 '/
	data	DATATABLE(37)	/'NALF'/
	data	DATATABLE(38)	/'TE  '/
	data	DATATABLE(39)	/'TI  '/
	data	DATATABLE(40)	/'T01 '/
	data	DATATABLE(41)	/'T02 '/
	data	DATATABLE(42)	/'PRE '/
	data	DATATABLE(43)	/'PRI '/
	data	DATATABLE(44)	/'PR  '/
	data	DATATABLE(45)	/'Q   '/
	data	DATATABLE(46)	/'ZEFF'/
	data	DATATABLE(47)	/'NUES'/
	data	DATATABLE(48)	/'NUIS'/
	data	DATATABLE(49)	/'BPOL'/
	data	DATATABLE(50)	/'BP2 '/
	data	DATATABLE(51)	/'ETA '/
	data	DATATABLE(52)	/'EZ  '/
	data	DATATABLE(53)	/'JZ  '/
	data	DATATABLE(54)	/'JZBS'/
	data	DATATABLE(55)	/'JZNB'/
	data	DATATABLE(56)	/'JZLH'/
	data	DATATABLE(57)	/'JZEC'/
	data	DATATABLE(58)	/'JZRF'/
	data	DATATABLE(59)	/'CUR '/
	data	DATATABLE(60)	/'CURA'/
	data	DATATABLE(61)	/'CUBS'/
	data	DATATABLE(62)	/'CUNB'/
	data	DATATABLE(63)	/'CULH'/
	data	DATATABLE(64)	/'CUEC'/
	data	DATATABLE(65)	/'CURF'/
	data	DATATABLE(66)	/'VOLT'/
	data	DATATABLE(67)	/'WALD'/
	data	DATATABLE(68)	/'S0  '/
	data	DATATABLE(69)	/'S01 '/
	data	DATATABLE(70)	/'S02 '/
	data	DATATABLE(71)	/'SRC1'/
	data	DATATABLE(72)	/'SRC2'/
	data	DATATABLE(73)	/'SCX1'/
	data	DATATABLE(74)	/'SCX2'/
	data	DATATABLE(75)	/'SII1'/
	data	DATATABLE(76)	/'SII2'/
	data	DATATABLE(77)	/'S0D '/
	data	DATATABLE(78)	/'S01D'/
	data	DATATABLE(79)	/'S02D'/
	data	DATATABLE(80)	/'SDRC'/
	data	DATATABLE(81)	/'SDII'/
	data	DATATABLE(82)	/'SDCX'/
	data	DATATABLE(83)	/'SB  '/
	data	DATATABLE(84)	/'SB1 '/
	data	DATATABLE(85)	/'SB2 '/
	data	DATATABLE(86)	/'SBD1'/
	data	DATATABLE(87)	/'SBD2'/
	data	DATATABLE(88)	/'SPD '/
	data	DATATABLE(89)	/'RTHD'/
	data	DATATABLE(90)	/'RTH '/
	data	DATATABLE(91)	/'RRFD'/
	data	DATATABLE(92)	/'RRF '/
	data	DATATABLE(93)	/'RNBD'/
	data	DATATABLE(94)	/'RNB '/
	data	DATATABLE(95)	/'RRD '/
	data	DATATABLE(96)	/'TORQ'/
	data	DATATABLE(97)	/'VTOR'/
	data	DATATABLE(98)	/'VPOL'/
	data	DATATABLE(99)	/'VPOZ'/
	data	DATATABLE(100)	/'VPAR'/
	data	DATATABLE(101)	/'VPAZ'/
	data	DATATABLE(102)	/'TOPR'/
	data	DATATABLE(103)	/'TVPO'/
	data	DATATABLE(104)	/'TVTO'/
	data	DATATABLE(105)	/'ERAD'/
	data	DATATABLE(106)	/'ETOR'/
	data	DATATABLE(107)	/'EPOL'/
	data	DATATABLE(108)	/'EPRE'/
	data	DATATABLE(109)	/'RR  '/
	data	DATATABLE(110)	/'NT  '/
	data	DATATABLE(111)	/'WNBD'/
	data	DATATABLE(112)	/'WNB '/
	data	DATATABLE(113)	/'DNBD'/
	data	DATATABLE(114)	/'DNB '/
	data	DATATABLE(115)	/'WRD '/
	data	DATATABLE(116)	/'WR  '/
	data	DATATABLE(117)	/'PERF'/
	data	DATATABLE(118)	/'PARN'/
	data	DATATABLE(119)	/'WDE '/
	data	DATATABLE(120)	/'WDI '/
	data	DATATABLE(121)	/'WE  '/
	data	DATATABLE(122)	/'WI  '/
	data	DATATABLE(123)	/'WBP '/
	data	DATATABLE(124)	/'QRFE'/
	data	DATATABLE(125)	/'PRFE'/
	data	DATATABLE(126)	/'QRFI'/
	data	DATATABLE(127)	/'PRFI'/
	data	DATATABLE(128)	/'QLHE'/
	data	DATATABLE(129)	/'PLHE'/
	data	DATATABLE(130)	/'QLHI'/
	data	DATATABLE(131)	/'PLHI'/
	data	DATATABLE(132)	/'QECE'/
	data	DATATABLE(133)	/'PECE'/
	data	DATATABLE(134)	/'QNBE'/
	data	DATATABLE(135)	/'PNBE'/
	data	DATATABLE(136)	/'QNBI'/
	data	DATATABLE(137)	/'PNBI'/
	data	DATATABLE(138)	/'QRAD'/
	data	DATATABLE(139)	/'PRAD'/
	data	DATATABLE(140)	/'QBRE'/
	data	DATATABLE(141)	/'PBRE'/
	data	DATATABLE(142)	/'QTHX'/
	data	DATATABLE(143)	/'PTHX'/
	data	DATATABLE(144)	/'PION'/
	data	DATATABLE(145)	/'Q0RC'/
	data	DATATABLE(146)	/'Q0II'/
	data	DATATABLE(147)	/'Q0CX'/
	data	DATATABLE(148)	/'Q0E '/
	data	DATATABLE(149)	/'Q0I '/
	data	DATATABLE(150)	/'P0E '/
	data	DATATABLE(151)	/'P0I '/
	data	DATATABLE(152)	/'Q0  '/
	data	DATATABLE(153)	/'P0  '/
	data	DATATABLE(154)	/'QOH '/
	data	DATATABLE(155)	/'POH '/
	data	DATATABLE(156)	/'QALE'/
	data	DATATABLE(157)	/'PALE'/
	data	DATATABLE(158)	/'QALI'/
	data	DATATABLE(159)	/'PALI'/
	data	DATATABLE(160)	/'QSYR'/
	data	DATATABLE(161)	/'PSYR'/
	data	DATATABLE(162)	/'TGEQ'/
	data	DATATABLE(163)	/'TGCE'/
	data	DATATABLE(164)	/'TGCI'/
	data	DATATABLE(165)	/'TGC '/
	data	DATATABLE(166)	/'TGRE'/
	data	DATATABLE(167)	/'TGRI'/
	data	DATATABLE(168)	/'TGR '/
	data	DATATABLE(169)	/'TGLO'/
	data	DATATABLE(170)	/'GROQ'/
	data	DATATABLE(171)	/'GRTC'/
	data	DATATABLE(172)	/'GRTE'/
	data	DATATABLE(173)	/'GRTI'/
	data	DATATABLE(174)	/'GRNE'/
	data	DATATABLE(175)	/'GRPE'/
	data	DATATABLE(176)	/'GRPI'/
	data	DATATABLE(177)	/'GRQ '/
	data	DATATABLE(178)	/'SH  '/
	data	DATATABLE(179)	/'ETAE'/
	data	DATATABLE(180)	/'ETAI'/
	data	DATATABLE(181)	/'XE  '/
	data	DATATABLE(182)	/'XI  '/
	data	DATATABLE(183)	/'XE0 '/
	data	DATATABLE(184)	/'XE1 '/
	data	DATATABLE(185)	/'XE2 '/
	data	DATATABLE(186)	/'XE3 '/
	data	DATATABLE(187)	/'XE4 '/
	data	DATATABLE(188)	/'XE5 '/
	data	DATATABLE(189)	/'XI0 '/
	data	DATATABLE(190)	/'XI1 '/
	data	DATATABLE(191)	/'XI2 '/
	data	DATATABLE(192)	/'XI3 '/
	data	DATATABLE(193)	/'XI4 '/
	data	DATATABLE(194)	/'XI5 '/
	data	DATATABLE(195)	/'XEF '/
	data	DATATABLE(196)	/'XIF '/
	data	DATATABLE(197)	/'XF1 '/
	data	DATATABLE(198)	/'XF2 '/
	data	DATATABLE(199)	/'DFI '/
	data	DATATABLE(200)	/'QEC '/
	data	DATATABLE(201)	/'QIC '/
	data	DATATABLE(202)	/'QE  '/
	data	DATATABLE(203)	/'QE2 '/
	data	DATATABLE(204)	/'QEW '/
	data	DATATABLE(205)	/'QTE '/
	data	DATATABLE(206)	/'QI  '/
	data	DATATABLE(207)	/'QI2 '/
	data	DATATABLE(208)	/'QTI '/
	data	DATATABLE(209)	/'DEDT'/
	data	DATATABLE(210)	/'DIDT'/
	data	DATATABLE(211)	/'GAME'/
	data	DATATABLE(212)	/'GAMI'/
	data	DATATABLE(213)	/'QEAV'/
	data	DATATABLE(214)	/'QIAV'/
	data	DATATABLE(215)	/'D1  '/
	data	DATATABLE(216)	/'D2  '/
	data	DATATABLE(217)	/'GANE'/
	data	DATATABLE(218)	/'GINE'/
	data	DATATABLE(219)	/'GWE '/
	data	DATATABLE(220)	/'GTE '/
	data	DATATABLE(221)	/'GANI'/
	data	DATATABLE(222)	/'GINI'/
	data	DATATABLE(223)	/'GWI '/
	data	DATATABLE(224)	/'GTI '/
	data	DATATABLE(225)	/'VW  '/
	data	DATATABLE(226)	/'VIN1'/
	data	DATATABLE(227)	/'VIN2'/
	data	DATATABLE(228)	/'TEO '/
	data	DATATABLE(229)	/'TSHE'/
	data	DATATABLE(230)	/'TBET'/
	data	DATATABLE(231)	/'TRHO'/
	data	DATATABLE(232)	/'OMEB'/
	data	DATATABLE(233)	/'ALFV'/
	data	DATATABLE(234)	/'OMOG'/
	data	DATATABLE(235)	/'OMOS'/
	data	DATATABLE(236)	/'SHOM'/
	data	DATATABLE(237)	/'SHOI'/
	data	DATATABLE(238)	/'GA01'/
	data	DATATABLE(239)	/'GA02'/
	data	DATATABLE(240)	/'GA03'/
	data	DATATABLE(241)	/'GA04'/
	data	DATATABLE(242)	/'GA05'/
	data	DATATABLE(243)	/'GA06'/
	data	DATATABLE(244)	/'GA07'/
	data	DATATABLE(245)	/'GA08'/
	data	DATATABLE(246)	/'GA09'/
	data	DATATABLE(247)	/'GA10'/
	data	DATATABLE(248)	/'GA11'/
	data	DATATABLE(249)	/'FR01'/
	data	DATATABLE(250)	/'FR02'/
	data	DATATABLE(251)	/'FR03'/
	data	DATATABLE(252)	/'FR04'/
	data	DATATABLE(253)	/'FR05'/
	data	DATATABLE(254)	/'FR06'/
	data	DATATABLE(255)	/'FR07'/
	data	DATATABLE(256)	/'FR08'/
	data	DATATABLE(257)	/'FR09'/
	data	DATATABLE(258)	/'FR10'/
	data	DATATABLE(259)	/'FR11'/
	data	DATATABLE(260)	/'CHEW'/
	data	DATATABLE(261)	/'CHIW'/
	data	DATATABLE(262)	/'DPAW'/
	data	DATATABLE(263)	/'RTR '/
	data	DATATABLE(264)	/'RTRS'/
	if (ICALL .eq. 1)	goto	902
	if (ICALL .ne. 0)	then
	   write(*,*)'>>> ERROR calling AJSP. Call ignored'
	endif
	jdt = 0
C Read file '$HOME/astra/runs/a2jsp'
C----------------------------------------------------------------------|
	if (NCH .eq. 1)	then
	   write(*,*)'>>> IFKEY >>> Calling AJSP ignored.'
	   return
	endif
	call	OPENRD(1,FILEX(1:index(FILEX,'runs')+4)//'a2jsp',0,IERR)
	if (IERR .gt. 0)	stop
	do	j=1,8
	   read(1,fmt='(1A132)',ERR=908,END=908)STRI
	enddo
 900	read(1,fmt='(1A132)',ERR=907,END=901)STRI
	read(STRI(1:4),*)j
	if (j .eq. 0)	goto	900
	if (j .ne. 1)	goto	908
	if (jdt .eq. 264)	then
	   write(*,*)'>>> ERROR >>> JSP control table is too long'
	   return
	endif
	jdt = jdt+1
	DATATYPE(jdt) = STRI(7:10)
	goto	900
 901	close(1)
	write(NCH)jdt
	write(NCH)(DATATYPE(j),j=1,jdt)
 902	continue
	write(NCH)TIME,NA1
	do 905	j2=1,jdt
	   jn = 0
	   do	j1=1,264
	      if (DATATYPE(j2) .eq. DATATABLE(j1))	then
		 jn = j1
		 goto	903
	      endif
	   enddo
 903	   continue
C 	   write(*,*)jn,char(9),DATATABLE(jn)
	   if (jn .le. 0)	then
	      write(*,*)'Unknown datatype requested: "',DATATYPE(j2),'"'
	      goto	905
	   endif
	   goto( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     >		11,12,13,14,15,16,17,18,19,20,
     >		21,22,23,24,25,26,27,28,29,30,
     >		31,32,33,34,35,36,37,38,39,40,
     >		41,42,43,44,45,46,47,48,49,50,
     >		51,52,53,54,55,56,57,58,59,60,
     >		61,62,63,64,65,66,67,68,69,70,
     >		71,72,73,74,75,76,77,78,79,80,
     >		81,82,83,84,85,86,87,88,89,90,
     >		91,92,93,94,95,96,97,98,99,100,
     >		101,102,103,104,105,106,107,108,109,110,
     >		111,112,113,114,115,116,117,118,119,120,
     >		121,122,123,124,125,126,127,128,129,130,
     >		131,132,133,134,135,136,137,138,139,140,
     >		141,142,143,144,145,146,147,148,149,150,
     >		151,152,153,154,155,156,157,158,159,160,
     >		161,162,163,164,165,166,167,168,169,170,
     >		171,172,173,174,175,176,177,178,179,180,
     >		181,182,183,184,185,186,187,188,189,190,
     >		191,192,193,194,195,196,197,198,199,200,
     >		201,202,203,204,205,206,207,208,209,210,
     >		211,212,213,214,215,216,217,218,219,220,
     >		221,222,223,224,225,226,227,228,229,230,
     >		231,232,233,234,235,236,237,238,239,240,
     >		241,242,243,244,245,246,247,248,249,250,
     >		251,252,253,254,255,256,257,258,259,260,
     >		261,262,263,264),jn
	   write(*,*)jn,char(9),DATATABLE(jn)
	   write(*,*)'Unrecognized datatype: "',DATATYPE(j2),'"'
	   goto	905
 1	   write(NCH)(RHO(j)/ROC,j=1,NA1)
	   goto	905
 2	   write(NCH)(RTOR+SHIFT+SHIF(j)+AMETR(j),j=1,NA1)
	   goto	905
 3	   write(NCH)(FP(j),j=1,NA1)
	   goto	905
 4	   write(NCH)(RHO(j),j=1,NA1)
	   goto	905
 5	   write(NCH)(AMETR(j)/ABC,j=1,NA1)
	   goto	905
 6	   write(NCH)(FP(j)/FP(NA1),j=1,NA1)
	   goto	905
 7	   write(NCH)(sqrt(FP(j)/FP(NA1)),j=1,NA1)
	   goto	905
 8	   write(NCH)(RTOR+SHIFT+SHIF(j)-AMETR(j),j=1,NA1)
	   goto	905
 9	   goto	904
 10	   goto	904
 11	   goto	904
 12	   write(NCH)(VOLUM(j),j=1,NA1)
	   goto	905
 13	   write(NCH)(VR(j),j=1,NA1)
	   goto	905
 14	   goto	904
 15	   goto	904
 16	   goto	904
 17	   goto	904
 18	   goto	904
 19	   write(NCH)(ELON(j),j=1,NA1)	   
	   goto	905
 20	   goto	904
 21	   goto	904
 22	   goto	904
 23	   goto	904
 24	   goto	904
 25	   goto	904
 26	   goto	904
 27	   goto	904
 28	   goto	904
 29	   write(NCH)(1.d19*NE(j),j=1,NA1)
	   goto	905
 30	   write(NCH)(1.d19*NI(j),j=1,NA1)
	   goto	905
 31	   write(NCH)(1.d19*NHYDR(j),j=1,NA1)
	   goto	905
 32	   write(NCH)(1.d19*NDEUT(j),j=1,NA1)
	   goto	905
 33	   write(NCH)(1.d19*NTRIT(j),j=1,NA1)
	   goto	905
 34	   write(NCH)(1.d19*NIZ1(j),j=1,NA1)
	   goto	905
 35	   write(NCH)(1.d19*(NNCL+NNWM)*NN(j),j=1,NA1)
	   goto	905
 36	   goto	904
 37	   write(NCH)(1.d19*NALF(j),j=1,NA1)
	   goto	905
 38	   write(NCH)(1.d3*TE(j),j=1,NA1)
	   goto	905
 39	   write(NCH)(1.d3*TI(j),j=1,NA1)
	   goto	905
 40	   write(NCH)(1.d3*TN(j),j=1,NA1)
	   goto	905
 41	   goto	904
 42	   write(NCH)(1.6d3*NE(j)*TE(j),j=1,NA1)
	   goto	905
 43	   write(NCH)(1.6d3*NI(j)*TI(j),j=1,NA1)
	   goto	905
 44	   write(NCH)(1.6d3*(NE(j)*TE(j)+NI(j)*TI(j)),j=1,NA1)
	   goto	905
 45	   write(NCH)(1./MU(j),j=1,NA1)
	   goto	905
 46	   write(NCH)(ZEF(j),j=1,NA1)
	   goto	905
 47	   do	j=1,na1
	      include 'fml/nues'
	      work1(j,1)=NUES
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 48	   do	j=1,na1
	      include 'fml/nuis'
	      work1(j,1)=NUIS
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 49	   goto	904
 50	   goto	904
 51	   write(NCH)(1.d-6/max(CC(j),0.1d0),j=1,NA1)
	   goto	905
 52	   write(NCH)(ULON(j)/gp2/RTOR,j=1,NA1)
	   goto	905
 53	   write(NCH)(1.d6*CU(j),j=1,NA1)
	   goto	905
 54	   write(NCH)(1.d6*CUBS(j),j=1,NA1)
	   goto	905
 55	   write(NCH)(1.d6*CUBM(j),j=1,NA1)
	   goto	905
 56	   write(NCH)(1.d6*CULH(j),j=1,NA1)
	   goto	905
 57	   write(NCH)(1.d6*CUECR(j),j=1,NA1)
	   goto	905
 58	   write(NCH)(1.d6*CUICR(j),j=1,NA1)
	   goto	905
 59	   write(NCH)(1.d6*ITOTR(RHO(j)),j=1,NA1)
	   goto	905
 60	   goto	904
 61	   write(NCH)(1.d6*IBSR(RHO(j)),j=1,NA1)
	   goto	905
 62	   write(NCH)(1.d6*IBMR(RHO(j)),j=1,NA1)
	   goto	905
 63	   write(NCH)(1.d6*ILHR(RHO(j)),j=1,NA1)
	   goto	905
 64	   write(NCH)(1.d6*IECRR(RHO(j)),j=1,NA1)
	   goto	905
 65	   write(NCH)(1.d6*IICRR(RHO(j)),j=1,NA1)
	   goto	905
 66	   write(NCH)(UPL(j),j=1,NA1)
	   goto	905
 67	   goto	904
 68	   QNNEU=0.
	   do	j=1,na1
	      include 'fml/snneu'
	      QNNEU=QNNEU+SNNEU*NE(j)*VR(j)*HRO
	      work1(j,1)=1.D19*QNNEU
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 69	   goto	904
 70	   goto	904
 71	   goto	904
 72	   goto	904
 73	   goto	904
 74	   goto	904
 75	   goto	904
 76	   goto	904
 77	   do	j=1,na1
	      include 'fml/snneu'
	      work1(j,1)=1.D19*SNNEU*NE(j)
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 78	   goto	904
 79	   goto	904
 80	   goto	904
 81	   goto	904
 82	   goto	904
 83	   write(NCH)(1.d19*VINT(SNIBM1,RHO(j)),j=1,NA1)
	   goto	905
 84	   goto	904
 85	   goto	904
 86	   write(NCH)(1.d19*SNIBM1(j),j=1,NA1)
	   goto	905
 87	   goto	904
 88	   goto	904
 89	   goto	904
 90	   goto	904
 91	   goto	904
 92	   goto	904
 93	   goto	904
 94	   goto	904
 95	   goto	904
 96	   goto	904
 97	   write(NCH)(VTOR(j),j=1,NA1)
	   goto	905
 98	   write(NCH)(VPOL(j),j=1,NA1)
	   goto	905
 99	   goto	904
 100	   goto	904
 101	   goto	904
 102	   goto	904
 103	   goto	904
 104	   goto	904
 105	   write(NCH)(ER(j),j=1,NA1)
	   goto	905
 106	   write(NCH)(-BTOR*VPOL(j),j=1,NA1)
	   goto	905
 107	   write(NCH)(BTOR*j*HRO*MU(j)*VTOR(j)/RTOR,j=1,NA1)
	   goto	905
 108	   do	j=1,na1
	      include 'fml/vdia'
	      work1(j,1)=VDIA
	   enddo
	   write(NCH)(BTOR*work1(j,1),j=1,NA1)
	   goto	905
 109	   goto	904
 110	   goto	904
 111	   goto	904
 112	   goto	904
 113	   goto	904
 114	   goto	904
 115	   goto	904
 116	   goto	904
 117	   goto	904
 118	   goto	904
 119	   goto	904
 120	   goto	904
 121	   write(NCH)(1.d6*WER(RHO(j)),j=1,NA1)
	   goto	905
 122	   write(NCH)(1.d6*WIR(RHO(j)),j=1,NA1)
	   goto	905
 123	   write(NCH)(1.d6*WBPOLR(RHO(j)),j=1,NA1)
	   goto	905
 124	   write(NCH)(1.d6*PEICR(j),j=1,NA1)
	   goto	905
 125	   write(NCH)(1.d6*VINT(PEICR,RHO(j)),j=1,NA1)
	   goto	905
 126	   write(NCH)(1.d6*PIICR(j),j=1,NA1)
	   goto	905
 127	   write(NCH)(1.d6*VINT(PIICR,RHO(j)),j=1,NA1)
	   goto	905
 128	   goto	904
 129	   goto	904
 130	   goto	904
 131	   goto	904
 132	   write(NCH)(1.d6*PEECR(j),j=1,NA1)
	   goto	905
 133	   write(NCH)(1.d6*VINT(PEECR,RHO(j)),j=1,NA1)
	   goto	905
 134	   write(NCH)(1.d6*PEBM(j),j=1,NA1)
	   goto	905
 135	   write(NCH)(1.d6*VINT(PEBM,RHO(j)),j=1,NA1)
	   goto	905
 136	   write(NCH)(1.d6*PIBM(j),j=1,NA1)
	   goto	905
 137	   write(NCH)(1.d6*VINT(PIBM,RHO(j)),j=1,NA1)
	   goto	905
 138	   write(NCH)(1.d6*PRAD(j),j=1,NA1)
	   goto	905
 139	   write(NCH)(1.d6*VINT(PRAD,RHO(j)),j=1,NA1)
	   goto	905
 140	   goto	904
 141	   goto	904
 142	   do	j=1,na1
	      include 'fml/peicl'
	      work1(j,1)=1.d6*PEICL
	   enddo
 	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 143	   do	j=1,na1
	      include 'fml/qeicl'
	      work1(j,1)=1.d6*QEICL
	   enddo
 	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 144	   goto	904
 145	   goto	904
 146	   goto	904
 147	   goto	904
 148	   goto	904
 149	   goto	904
 150	   goto	904
 151	   goto	904
 152	   goto	904
 153	   goto	904
 154	   do	j=1,na1
	      include 'fml/poh'
	      work1(j,1)=1.d6*POH
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 155	   goto	904
 156	   do	j=1,na1
	      include 'fml/pedt'
	      work1(j,1)=1.d6*PEDT
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 157	   goto	904
 158	   do	j=1,na1
	      include 'fml/pidt'
	      work1(j,1)=1.d6*PIDT
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 159	   goto	904
 160	   goto	904
 161	   goto	904
 162	   goto	904
 163	   goto	904
 164	   goto	904
 165	   goto	904
 166	   goto	904
 167	   goto	904
 168	   goto	904
 169	   goto	904
 170	   goto	904
 171	   goto	904
 172	   goto	904
 173	   goto	904
 174	   goto	904
 175	   goto	904
 176	   goto	904
 177	   goto	904
 178	   write(NCH)(SHEAR(j),j=1,NA1)
	   goto	905
 179	   goto	904
 180	   goto	904
 181	   write(NCH)(HE(j),j=1,NA1)
	   goto	905
 182	   write(NCH)(XI(j),j=1,NA1)
	   goto	905
 183	   goto	904
 184	   write(NCH)(work(j,305),j=1,NA1)
	   goto	905
 185	   write(NCH)(HE(j)-work(j,102),j=1,NA1)! XE2 all except for GLF
	   goto	905
 186	   write(NCH)(work(j,102),j=1,NA1)	! XE3 used for GLF
	   goto	905
 187	   goto	904
 188	   goto	904
 189	   goto	904
 190	   write(NCH)(work(j,365),j=1,NA1)
	   goto	905
 191	   write(NCH)(XI(j)-work(j,101),j=1,NA1)! XI2 all except for GLF
	   goto	905
 192	   write(NCH)(work(j,101),j=1,NA1)	! XI3 used for GLF
	   goto	905
 193	   goto	904
 194	   goto	904
 195	   do	j=1,na1
	      include 'fml/heeff'
	      work1(j,1)=HEEFF
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 196	   do	j=1,na1
	      include 'fml/xieff'
	      work1(j,1)=XIEFF
	   enddo
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 197	   goto	904
 198	   goto	904
 199	   goto	904
 200	   goto	904
 201	   goto	904
 202	   goto	904
 203	   goto	904
 204	   goto	904
 205	   goto	904
 206	   goto	904
 207	   goto	904
 208	   goto	904
 209	   goto	904
 210	   goto	904
 211	   goto	904
 212	   goto	904
 213	   goto	904
 214	   goto	904
 215	   write(NCH)(DN(j),j=1,NA1)	   
	   goto	905
 216	   write(NCH)(work(j,362),j=1,NA1)
	   goto	905
 217	   goto	904
 218	   goto	904
 219	   goto	904
 220	   goto	904
 221	   goto	904
 222	   goto	904
 223	   goto	904
 224	   goto	904
 225	   goto	904
 226	   write(NCH)(VF2(j),j=1,NA1)
	   goto	905
 227	   goto	904
 228	   goto	904
 229	   goto	904
 230	   goto	904
 231	   goto	904
 232	   do	j=1,NA
	      include 'fml/rotsh'
	      work1(j,1)=ROTSH
	   enddo
	   work1(NA1,1)=work1(na,1)
	   write(NCH)(work1(j,1),j=1,NA1)
	   goto	905
 233	   goto	904
 234	   goto	904
 235	   goto	904
 236	   goto	904
 237	   goto	904
 238	   goto	904
 239	   goto	904
 240	   goto	904
 241	   goto	904
 242	   goto	904
 243	   goto	904
 244	   goto	904
 245	   goto	904
 246	   goto	904
 247	   goto	904
 248	   goto	904
 249	   goto	904
 250	   goto	904
 251	   goto	904
 252	   goto	904
 253	   goto	904
 254	   goto	904
 255	   goto	904
 256	   goto	904
 257	   goto	904
 258	   goto	904
 259	   goto	904
 260	   goto	904
 261	   goto	904
 262	   goto	904
 263	   goto	904
 264	   goto	904
 904	   write(*,*)'Sorry, the requested datatype "',
     >			DATATYPE(j2),'" is not implemented'
 905	continue
	return
 907	write(*,*)'>>> File "a2jsp" reading error'
	stop
 908	write(*,*)'>>> File "a2jsp" reading error'
	write(*,*)'    Line "',STRI(1:length(STRI)),'"'
	stop
 910	format(1p,5e13.5)
	end
C======================================================================|
	subroutine	AJST(NCH,ICALL)
C----------------------------------------------------------------------|
C   The subroutine is called from IFKEY and writes the binary file
C   ~/astra/runs/runID/a4jst (unit 1)
C   with account of the control file
C   ~/astra/runs/a2jst (unit NCH)
C----------------------------------------------------------------------|
	implicit  none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
      double precision ALMHD,BETE,BETI,BETPL,BETR,CCMHD,CCNEU,CCSP
      double precision CCSPX,CERL,CHOTF,CNHH,CNHR,CNSA,COULG,CS,CUOHM
      double precision D2TI,DBOHM,DCHA,DCHH,DCHR,DCHR2,DCKIM,DCKM1
      double precision DCSA,DHKIM,DIDT,DNDIF,DNEXP,EDR,EFLHN,EFLHW
      double precision ENHH0,ENHH1,ENHH2,EPAR,EPL,ETAE,ETAI,ETAN,FA
      double precision FLIN,FOWC,FPA,FPR,FR,FRS,FTE,FX,GAITG,GITG
      double precision GITG0,HAALC,HABM,HABMS,HABOM,HACTE,HAED,HAETI
      double precision HAETS,HAGB,HAGBS,HAGBS1,HAITF,HAITR,HAMM,HANAG
      double precision HANAL,HAPA,HAPUE,HAPYU,HAQ1,HAQ1C,HARL,HARLS
      double precision HARNQ,HARPL,HASCL,HASDL,HATL,HATLI,HATLS,HBJET
      double precision HCHA,HCHGP,HCHH,HCHII,HCHR,HCHR2,HCKIM,HCKM1
      double precision HCSA,HEEFF,HEGN,HETAI,HETIS,HEXP,HGBEJ,HGBIJ
      double precision HIBOM,HIT89,HMHD1,HMHD2,HNASI,HNCHI,HNGSB,HNGSE
      double precision HNGSI,HNGSP,HNPSI,IBM,IBS,ICD,IOHM,ITOT,LNE,LNI
      double precision LNZ1,LTE,LTI,NUE,NUEE,NUES,NUFE,NUI,NUIS,NUPP
      double precision PAION,PAION1,PBEIE,PBICX,PBR1,PBR2,PBR3,PBRAD
      double precision PDT,PDTF,PEDT,PEDT1,PEDTF,PEGN,PEHCL,PEI,PEICL
      double precision PEIGN,PENEU,PENLI,PETSL,PHICL,PICX,PIDT,PIDT1
      double precision PIDTF,PIGN,PINEU,PIONZ,PIREC,PITCX,PITSL,PJOUL
      double precision POH,PRARG,PRBER,PRCAR,PRESE,PRESI,PREST,PRFER
      double precision PRNEO,PRNIT,PROXI,PROXY,PRWOL,PSYNC,PTOT,QBM
      double precision QBREM,QBS,QCD,QDT,QDTF,QEDT,QEDTF,QEGN,QEICL
      double precision QEIGN,QENEU,QETOT,QEX,QIDT,QIDTF,QIGN,QINEU
      double precision QITOT,QIX,QJOUL,QNIND,QNTOT,QNX,QOH,QRAD,QRADX
      double precision QSYNC,QTOT,RLI,RLS,RLTCR,RLTCZ,RLTKD,RLTWN
      double precision RNODE,ROTSH,SHAT,SNNEU,SNNIE,SNNII,SNNR
      double precision SQZ,SVCX,SVCXX,SVD1,SVD2,SVDBH,SVDHE,SVDT,SVIE
      double precision SVIEP,SVIEX,SVII,SVREC,SVRECX,TAU89,TAUNA,TECRJ
      double precision TECRL,TECRL1,TEF,THQ99,TITER,TLQ97,TPF,VDIA
      double precision VPSWW,VRAS,VRHH,VSI,VTE,VTI,WE,WEX,WI,WIX,WTOT
      double precision XCH86,XCHA,XCHH,XCHR,XCHR2,XCKIM,XCKM1,XCSA
      double precision XEXP,XIEFF,XIGN,XIRL,XMHD1,XMHD2,XTEXP,ZIAR
      double precision ZIBE,ZICAR,ZIFER,ZINE,ZINEO,ZINIT,ZIOXI,ZIWOL
      double precision ZZEF
      double precision VINT,IINT,GRAD,FRMAX,FRMIN,RFMIN,RFMAX
      double precision RFVAL,AFVAL,RFVEX,AFVEX,RFVIN,AFVIN,RFA
      double precision RFAN,XFA,XFAN,AFR,AFX,RECR,ATR,ATX
      double precision TIMINT,TIMDER,TIMAVG,GAUSS
      double precision RADIAL,RADINT,ASTEP,RSTEP,XSTEP,STEP,CUT
      double precision FIXVAL,FTAV,FTMIN,FTMAX,FRAMP,FJUMP
      external IINT
      double precision BETAJR,BETBMR,BETRR,FTLLMR,HNASCR,HNCHCR,IBMR
      double precision IBSR,ICDR,IECRR,IFIR,IFWR,IICRR,ILHR,IOHMR,ITOTR
      double precision IXR,LICDR,LINTR,NEAVR,NECHR,NELAR,NEXAVR,PRCARR
      double precision PTOTR,QBRADR,QBTOTR,QDTR,QDTFR,QEDTR,QEDTFR
      double precision QEDWTR,QEGNR,QEICLR,QEIGNR,QENEUR,QETOTR,QEXR
      double precision QICXR,QIDTR,QIDTFR,QIDWTR,QIGNR,QINEUR,QITOTR
      double precision QIXR,QJOULR,QMINR,QNDNTR,QNTOTR,QNXR,QOHR,QRADR
      double precision QRADXR,QSYNCR,QTOTR,SHINER,TAUER,TAUEER
      double precision TAUEIR,TAUGR,TAUIGR,TAUPR,TEAVR,TENDNR,TEXAVR
      double precision TIAVR,TINDNR,TIXAVR,WALFR,WBPOLR,WCR,WCER,WCEXR
      double precision WCIR,WCIXR,WCXR,WER,WEXR,WIR,WIXR,WTOTR,WTOTXR
      double precision XQMINR,ZNDNR,BETAR
	character*4 STRI*132,DATATYPE(189),DATATABLE(189)
	integer	  j,j1,j2,NCH,ICALL,IERR,length,jdt,jn
	double precision YVC,YVE,YEPS,YVALP,YLLAME,YLLAMI,YLLAMA
	save	DATATYPE,DATATABLE,jdt
C The data list below is made making use of 
C /u/sim/docs/source/default/config/jst
	data	DATATABLE(1)	/'TEAX'/
	data	DATATABLE(2)	/'TIAX'/
	data	DATATABLE(3)	/'TEAV'/
	data	DATATABLE(4)	/'TIAV'/
	data	DATATABLE(5)	/'TE08'/
	data	DATATABLE(6)	/'TEBO'/
	data	DATATABLE(7)	/'TIBO'/
	data	DATATABLE(8)	/'NEAX'/
	data	DATATABLE(9)	/'NHAX'/
	data	DATATABLE(10)	/'NEAV'/
	data	DATATABLE(11)	/'NIAV'/
	data	DATATABLE(12)	/'NI1 '/
	data	DATATABLE(13)	/'NI2 '/
	data	DATATABLE(14)	/'NIMP'/
	data	DATATABLE(15)	/'NNB1'/
	data	DATATABLE(16)	/'NNB2'/
	data	DATATABLE(17)	/'NIN1'/
	data	DATATABLE(18)	/'NIN2'/
	data	DATATABLE(19)	/'NPEL'/
	data	DATATABLE(20)	/'NEBO'/
	data	DATATABLE(21)	/'NHBO'/
	data	DATATABLE(22)	/'NEL '/
	data	DATATABLE(23)	/'ZEFF'/
	data	DATATABLE(24)	/'S0  '/
	data	DATATABLE(25)	/'S01 '/
	data	DATATABLE(26)	/'S02 '/
	data	DATATABLE(27)	/'S0II'/
	data	DATATABLE(28)	/'S0RC'/
	data	DATATABLE(29)	/'S0CX'/
	data	DATATABLE(30)	/'SII1'/
	data	DATATABLE(31)	/'SII2'/
	data	DATATABLE(32)	/'SRC1'/
	data	DATATABLE(33)	/'SRC2'/
	data	DATATABLE(34)	/'SCX1'/
	data	DATATABLE(35)	/'SCX2'/
	data	DATATABLE(36)	/'SNB '/
	data	DATATABLE(37)	/'SNB1'/
	data	DATATABLE(38)	/'SNB2'/
	data	DATATABLE(39)	/'SIN1'/
	data	DATATABLE(40)	/'SIN2'/
	data	DATATABLE(41)	/'PIN1'/
	data	DATATABLE(42)	/'PIN2'/
	data	DATATABLE(43)	/'SOUA'/
	data	DATATABLE(44)	/'FLXE'/
	data	DATATABLE(45)	/'FLX1'/
	data	DATATABLE(46)	/'FLX2'/
	data	DATATABLE(47)	/'FLXI'/
	data	DATATABLE(48)	/'NT  '/
	data	DATATABLE(49)	/'RR  '/
	data	DATATABLE(50)	/'RRTH'/
	data	DATATABLE(51)	/'RRRF'/
	data	DATATABLE(52)	/'RRNB'/
	data	DATATABLE(53)	/'LVOL'/
	data	DATATABLE(54)	/'VLP '/
	data	DATATABLE(55)	/'RFLU'/
	data	DATATABLE(56)	/'WTHE'/
	data	DATATABLE(57)	/'WTHI'/
	data	DATATABLE(58)	/'WTH '/
	data	DATATABLE(59)	/'WTOT'/
	data	DATATABLE(60)	/'WERL'/
	data	DATATABLE(61)	/'POH '/
	data	DATATABLE(62)	/'PLHE'/
	data	DATATABLE(63)	/'PLHI'/
	data	DATATABLE(64)	/'PLH '/
	data	DATATABLE(65)	/'PRFE'/
	data	DATATABLE(66)	/'PRFI'/
	data	DATATABLE(67)	/'PRF '/
	data	DATATABLE(68)	/'PNBE'/
	data	DATATABLE(69)	/'PNBI'/
	data	DATATABLE(70)	/'PNB '/
	data	DATATABLE(71)	/'PALE'/
	data	DATATABLE(72)	/'PALI'/
	data	DATATABLE(73)	/'PALF'/
	data	DATATABLE(74)	/'PTOT'/
	data	DATATABLE(75)	/'PLHC'/
	data	DATATABLE(76)	/'PRFC'/
	data	DATATABLE(77)	/'PNBC'/
	data	DATATABLE(78)	/'PFEL'/
	data	DATATABLE(79)	/'PFIO'/
	data	DATATABLE(80)	/'PTHX'/
	data	DATATABLE(81)	/'PION'/
	data	DATATABLE(82)	/'P0E '/
	data	DATATABLE(83)	/'P0I '/
	data	DATATABLE(84)	/'P0  '/
	data	DATATABLE(85)	/'P0RC'/
	data	DATATABLE(86)	/'P0II'/
	data	DATATABLE(87)	/'P0CX'/
	data	DATATABLE(88)	/'PRAD'/
	data	DATATABLE(89)	/'PBRE'/
	data	DATATABLE(90)	/'PSYR'/
	data	DATATABLE(91)	/'PFUS'/
	data	DATATABLE(92)	/'TAUP'/
	data	DATATABLE(93)	/'TACE'/
	data	DATATABLE(94)	/'TACI'/
	data	DATATABLE(95)	/'TACT'/
	data	DATATABLE(96)	/'TART'/
	data	DATATABLE(97)	/'TAUG'/
	data	DATATABLE(98)	/'TA98'/
	data	DATATABLE(99)	/'TAIT'/
	data	DATATABLE(100)	/'QDT '/
	data	DATATABLE(101)	/'Q5  '/
	data	DATATABLE(102)	/'WRRF'/
	data	DATATABLE(103)	/'WRNB'/
	data	DATATABLE(104)	/'DRNB'/
	data	DATATABLE(105)	/'WR  '/
	data	DATATABLE(106)	/'PERF'/
	data	DATATABLE(107)	/'CUR '/
	data	DATATABLE(108)	/'CURA'/
	data	DATATABLE(109)	/'CUBS'/
	data	DATATABLE(110)	/'CUNB'/
	data	DATATABLE(111)	/'CULH'/
	data	DATATABLE(112)	/'CURF'/
	data	DATATABLE(113)	/'WBP '/
	data	DATATABLE(114)	/'BETP'/
	data	DATATABLE(115)	/'BETT'/
	data	DATATABLE(116)	/'BETI'/
	data	DATATABLE(117)	/'LI  '/
	data	DATATABLE(118)	/'Q0  '/
	data	DATATABLE(119)	/'Q95 '/
	data	DATATABLE(120)	/'QBO '/
	data	DATATABLE(121)	/'ROQ1'/
	data	DATATABLE(122)	/'BTOR'/
	data	DATATABLE(123)	/'DRBA'/
	data	DATATABLE(124)	/'WDE '/
	data	DATATABLE(125)	/'WDI '/
	data	DATATABLE(126)	/'WDOT'/
	data	DATATABLE(127)	/'RONC'/
	data	DATATABLE(128)	/'XI1N'/
	data	DATATABLE(129)	/'QEB '/
	data	DATATABLE(130)	/'QIB '/
	data	DATATABLE(131)	/'NUB '/
	data	DATATABLE(132)	/'NUBA'/
	data	DATATABLE(133)	/'PBAL'/
	data	DATATABLE(134)	/'PBAF'/
	data	DATATABLE(135)	/'PLOS'/
	data	DATATABLE(136)	/'QBCI'/
	data	DATATABLE(137)	/'QBCE'/
	data	DATATABLE(138)	/'QBPI'/
	data	DATATABLE(139)	/'QBPE'/
	data	DATATABLE(140)	/'QBNI'/
	data	DATATABLE(141)	/'QBNE'/
	data	DATATABLE(142)	/'QBVI'/
	data	DATATABLE(143)	/'QBVE'/
	data	DATATABLE(144)	/'QBAI'/
	data	DATATABLE(145)	/'QBAE'/
	data	DATATABLE(146)	/'QBMI'/
	data	DATATABLE(147)	/'QBME'/
	data	DATATABLE(148)	/'QBTI'/
	data	DATATABLE(149)	/'QBTE'/
	data	DATATABLE(150)	/'QNEO'/
	data	DATATABLE(151)	/'PFPI'/
	data	DATATABLE(152)	/'PFP '/
	data	DATATABLE(153)	/'PFN '/
	data	DATATABLE(154)	/'PFV '/
	data	DATATABLE(155)	/'PFA '/
	data	DATATABLE(156)	/'PFM '/
	data	DATATABLE(157)	/'ATM '/
	data	DATATABLE(158)	/'DTMX'/
	data	DATATABLE(159)	/'QMIN'/
	data	DATATABLE(160)	/'ROQM'/
	data	DATATABLE(161)	/'NHDO'/
	data	DATATABLE(162)	/'XNBA'/
	data	DATATABLE(163)	/'XABA'/
	data	DATATABLE(164)	/'DBAR'/
	data	DATATABLE(165)	/'NEBA'/
	data	DATATABLE(166)	/'TEBA'/
	data	DATATABLE(167)	/'TIBA'/
	data	DATATABLE(168)	/'XEBA'/
	data	DATATABLE(169)	/'XIBA'/
	data	DATATABLE(170)	/'RLBA'/
	data	DATATABLE(171)	/'XI1B'/
	data	DATATABLE(172)	/'CBAR'/
	data	DATATABLE(173)	/'RECY'/
	data	DATATABLE(174)	/'GI1 '/
	data	DATATABLE(175)	/'GI2 '/
	data	DATATABLE(176)	/'N0AX'/
	data	DATATABLE(177)	/'N0BO'/
	data	DATATABLE(178)	/'T0AX'/
	data	DATATABLE(179)	/'T0BO'/
	data	DATATABLE(180)	/'CN01'/
	data	DATATABLE(181)	/'CN02'/
	data	DATATABLE(182)	/'VELE'/
	data	DATATABLE(183)	/'VELI'/
	data	DATATABLE(184)	/'VELP'/
	data	DATATABLE(185)	/'KSIE'/
	data	DATATABLE(186)	/'KSII'/
	data	DATATABLE(187)	/'KSIP'/
	data	DATATABLE(188)	/'ROBA'/
	data	DATATABLE(189)	/'JTOB'/
	if (ICALL .eq. 1)	goto	902
	if (ICALL .ne. 0)	then
	   write(*,*)'>>> ERROR calling AJST. Call ignored'
	endif
	jdt = 0
C Read file '$HOME/astra/runs/a2jst'
C----------------------------------------------------------------------|
	if (NCH .eq. 1)	then
	   write(*,*)'>>> IFKEY >>> Calling AJST ignored.'
	   return
	endif
	call	OPENRD(1,FILEX(1:index(FILEX,'runs')+4)//'a2jst',0,IERR)
	if (IERR .gt. 0)	stop
	do	j=1,8
	   read(1,fmt='(1A132)',ERR=908,END=908)STRI
	enddo
 900	read(1,fmt='(1A132)',ERR=907,END=901)STRI
	read(STRI(1:4),*)j
	if (j .eq. 0)	goto	900
	if (j .ne. 1)	goto	908
	if (jdt .eq. 189)	then
	   write(*,*)'>>> ERROR >>> JST control table is too long'
	   return
	endif
	jdt = jdt+1
	DATATYPE(jdt) = STRI(7:10)
	goto	900
 901	close(1)
	write(NCH)jdt
	write(NCH)(DATATYPE(j),j=1,jdt)
 902	continue
	write(NCH)TIME
	do 905	j2=1,jdt
	   jn = 0
	   do	j1=1,189
	      if (DATATYPE(j2) .eq. DATATABLE(j1))	then
		 jn = j1
		 goto	903
	      endif
	   enddo
 903	   continue
C	   write(*,*)jn,char(9),DATATABLE(jn)
	   if (jn .le. 0)	then
	      write(*,*)'Unknown datatype requested: "',DATATYPE(j2),'"'
	      goto	905
	   endif
	   goto( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     >		11,12,13,14,15,16,17,18,19,20,
     >		21,22,23,24,25,26,27,28,29,30,
     >		31,32,33,34,35,36,37,38,39,40,
     >		41,42,43,44,45,46,47,48,49,50,
     >		51,52,53,54,55,56,57,58,59,60,
     >		61,62,63,64,65,66,67,68,69,70,
     >		71,72,73,74,75,76,77,78,79,80,
     >		81,82,83,84,85,86,87,88,89,90,
     >		91,92,93,94,95,96,97,98,99,100,
     >		101,102,103,104,105,106,107,108,109,110,
     >		111,112,113,114,115,116,117,118,119,120,
     >		121,122,123,124,125,126,127,128,129,130,
     >		131,132,133,134,135,136,137,138,139,140,
     >		141,142,143,144,145,146,147,148,149,150,
     >		151,152,153,154,155,156,157,158,159,160,
     >		161,162,163,164,165,166,167,168,169,170,
     >		171,172,173,174,175,176,177,178,179,180,
     >		181,182,183,184,185,186,187,188,189),jn
	   write(*,*)jn,char(9),DATATABLE(jn)
	   write(*,*)'Unrecognized datatype: "',DATATYPE(j2),'"'
	   goto	905
 1	   write(NCH)1.d3*TE(1)
	   goto	905
 2	   write(NCH)1.d3*TI(1)
	   goto	905
 3	   write(NCH)1.d3*TEAVR(ROC)
	   goto	905
 4	   write(NCH)1.d3*TIAVR(ROC)
	   goto	905
 5	   write(NCH)1.d3*RADIAL(TE,RFA(AFX(8.d-1)))
	   goto	905
 6	   write(NCH)1.d3*TE(NA1)
	   goto	905
 7	   write(NCH)1.d3*TI(NA1)
	   goto	905
 8	   write(NCH)1.d19*NE(1)
	   goto	905
 9	   write(NCH)1.d19*NHYDR(1)
	   goto	905
 10	   write(NCH)1.d19*NEAVR(ROC)
	   goto	905
 11	   write(NCH)1.d19*VINT(NI,ROC)/VOLUME
	   goto	905
 12	   write(NCH)1.d19*VINT(NIZ1,ROC)/VOLUME
	   goto	905
 13	   write(NCH)1.d19*VINT(NIZ2,ROC)/VOLUME
	   goto	905
 14	   goto	904
 15	   goto	904
 16	   goto	904
 17	   goto	904
 18	   goto	904
 19	   goto	904
 20	   write(NCH)1.d19*NE(NA1)
	   goto	905
 21	   write(NCH)1.d19*NHYDR(NA1)
	   goto	905
 22	   write(NCH)1.d19*NECHR(ROC)
	   goto	905
 23	   write(NCH)1.d19*VINT(ZEF,ROC)/VOLUME
	   goto	905
 24	   write(NCH)1.d19*(VINT(SNTOT,ROC)-VINT(SNIBM1,ROC))
	   goto	905
 25	   goto	904
 26	   goto	904
 27	   goto	904
 28	   goto	904
 29	   goto	904
 30	   goto	904
 31	   goto	904
 32	   goto	904
 33	   goto	904
 34	   goto	904
 35	   goto	904
 36	   write(NCH)1.d19*VINT(SNIBM1,ROC)
	   goto	905
 37	   goto	904
 38	   goto	904
 39	   write(NCH)1.d19*((VINT(SNTOT,ROC)-VINT(SNIBM1,ROC))*
     .                    (1-CBND1*ALBPL)-CBND1*QN(NA1))
	   goto	905
 40	   goto	904
 41	   goto	904
 42	   goto	904
 43	   goto	904
 44	   write(NCH)1.d19*QN(NA1)
	   goto	905
 45	   goto	904
 46	   goto	904
 47	   goto	904
 48	   goto	904
 49	   goto	904
 50	   goto	904
 51	   goto	904
 52	   goto	904
 53	   write(NCH)UPL(NA)
	   goto	905
 54	   write(NCH)ULON(NA)
	   goto	905
 55	   goto	904
 56	   write(NCH)1.d6*WER(ROC)
	   goto	905
 57	   write(NCH)1.d6*WIR(ROC)
	   goto	905
 58	   write(NCH)1.d6*WTOTR(ROC)
	   goto	905
 59	   goto	904
 60	   goto	904
 61	   write(NCH)1.d6*QOHR(ROC)
	   goto	905
 62	   write(NCH)1.d6*VINT(PELH,ROC)
	   goto	905
 63	   goto	904
 64	   write(NCH)QLH
	   goto	905
 65	   write(NCH)1.d6*VINT(PEICR,ROC)
	   goto	905
 66	   write(NCH)1.d6*VINT(PIICR,ROC)
	   goto	905
 67	   write(NCH)1.d6*QICR
	   goto	905
 68	   write(NCH)1.d6*VINT(PEBM,ROC)
	   goto	905
 69	   write(NCH)1.d6*VINT(PIBM,ROC)
	   goto	905
 70	   write(NCH)1.d6*QNBI
	   goto	905
 71	   write(NCH)1.d6*QEDTR(ROC)
	   goto	905
 72	   write(NCH)1.d6*QIDTR(ROC)
	   goto	905
 73	   write(NCH)1.d6*QDTR(ROC)
	   goto	905
 74	   write(NCH)1.d6*QTOTR(ROC)
	   goto	905
 75	   goto	904
 76	   goto	904
 77	   goto	904
 78	   write(NCH)1.d6*QE(NA1)
	   goto	905
 79	   write(NCH)1.d6*QI(NA1)
	   goto	904
 80	   goto	904
 81	   goto	904
 82	   goto	904
 83	   goto	904
 84	   goto	904
 85	   goto	904
 86	   goto	904
 87	   goto	904
 88	   write(NCH)1.d6*QRADR(ROC)
	   goto	905
 89	   write(NCH)1.d6*QBRADR(ROC)
	   goto	905
 90	   write(NCH)1.d6*QSYNCR(ROC)
	   goto	905
 91	   write(NCH)4.997d6*QDTR(ROC)
	   goto	905
 92	   write(NCH)TAUPR(ROC)
	   goto	905
 93	   write(NCH)TAUEER(ROC)
	   goto	905
 94	   write(NCH)TAUEIR(ROC)
	   goto	905
 95	   write(NCH)TAUER(ROC)
	   goto	905
 96	   write(NCH)TAUER(ROC)
	   goto	905
 97	   goto	904
 98	   continue
	   include 'fml/thq99'
	   write(NCH)THQ99
	   goto	905
 99	   continue
	   include 'fml/tau89'
	   write(NCH)TAU89
	   goto	905
 100	   goto	904
 101	   goto	904
 102	   goto	904
 103	   goto	904
 104	   goto	904
 105	   goto	904
 106	   goto	904
 107	   write(NCH)1.d6*ITOTR(ROC)
	   goto	905
 108	   write(NCH)1.d6*ICDR(ROC)
	   goto	905
 109	   write(NCH)1.d6*IBSR(ROC)
	   goto	905
 110	   write(NCH)1.d6*IBMR(ROC)
	   goto	905
 111	   write(NCH)1.d6*ILHR(ROC)
	   goto	905
 112	   write(NCH)1.d6*IICRR(ROC)
	   goto	905
 113	   write(NCH)1.d6*WBPOLR(ROC)
	   goto	905
 114	   write(NCH)BETAJR(ROC)
	   goto	905
 115	   write(NCH)BETAR(ROC)
	   goto	905
 116	   write(NCH)BETAJR(ROC)
	   goto	905
 117	   write(NCH)LINTR(ROC)
	   goto	905
 118	   write(NCH)1./MU(1)
	   goto	905
 119	   write(NCH)1./RADIAL(MU,RFA(AFX(5.d-1)))
	   goto	905
 120	   write(NCH)1./MU(NA1)
	   goto	905
 121	   goto	904
 122	   write(NCH)BTOR
	   goto	905
 123	   goto	904
 124	   goto	904
 125	   goto	904
 126	   goto	904
 127	   goto	904
 128	   goto	904
 129	   goto	904
 130	   goto	904
 131	   goto	904
 132	   goto	904
 133	   goto	904
 134	   goto	904
 135	   goto	904
 136	   goto	904
 137	   goto	904
 138	   goto	904
 139	   goto	904
 140	   goto	904
 141	   goto	904
 142	   goto	904
 143	   goto	904
 144	   goto	904
 145	   goto	904
 146	   goto	904
 147	   goto	904
 148	   goto	904
 149	   goto	904
 150	   goto	904
 151	   goto	904
 152	   goto	904
 153	   goto	904
 154	   goto	904
 155	   goto	904
 156	   goto	904
 157	   write(NCH)AMJ
	   goto	905
 158	   goto	904
 159	   goto	904
 160	   goto	904
 161	   goto	904
 162	   goto	904
 163	   goto	904
 164	   goto	904
 165	   goto	904
 166	   goto	904
 167	   goto	904
 168	   goto	904
 169	   goto	904
 170	   goto	904
 171	   goto	904
 172	   goto	904
 173	   goto	904
 174	   goto	904
 175	   goto	904
 176	   goto	904
 177	   goto	904
 178	   goto	904
 179	   goto	904
 180	   goto	904
 181	   goto	904
 182	   goto	904
 183	   goto	904
 184	   goto	904
 185	   goto	904
 186	   goto	904
 187	   goto	904
 188	   goto	904
 189	   goto	904
 904	   write(*,*)'Sorry, the requested datatype "',
     >			DATATYPE(j2),'" is not implemented'
 905	continue
	return
 907	write(*,*)'>>> File "a2jst" reading error'
	stop
 908	write(*,*)'>>> File "a2jst" reading error'
	write(*,*)'    Line "',STRI(1:length(STRI)),'"'
	stop
 910	format(1p,5e13.5)
	end
C======================================================================|
	subroutine	AJSE(NCH,ICALL)
C----------------------------------------------------------------------|
	implicit  none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer	  j,j1,j2,jdt,NCH,ICALL,IERR,length,jn,jnb
	character*4 STRI*132,DATATYPE(2),DATATABLE(2)
	double	precision	YRB(250),YZB(250),Y,YS
	save	DATATYPE,DATATABLE,jdt
	data	DATATABLE(1)	/'RBND'/
	data	DATATABLE(2)	/'ZBND'/
	if (ICALL .eq. 1)	goto	902
	if (ICALL .ne. 0)	then
	   write(*,*)'>>> ERROR calling AJSE. Call ignored'
	endif
	jdt = 0
C Read file '$HOME/astra/runs/a2jsp'
C----------------------------------------------------------------------|
	if (NCH .eq. 1)	then
	   write(*,*)'>>> IFKEY >>> Calling AJSE ignored.'
	   return
	endif
	call	OPENRD(1,FILEX(1:index(FILEX,'runs')+4)//'a2jse',0,IERR)
	if (IERR .gt. 0)	stop
	do	j=1,8					! Skip 8 lines
	   read(1,fmt='(1A132)',ERR=908,END=908)STRI
	enddo
 900	read(1,fmt='(1A132)',ERR=907,END=901)STRI
	read(STRI(1:4),*)j
	if (j .eq. 0)	goto	900
	if (j .ne. 1)	goto	908
	if (jdt .eq. 2)	then
	   write(*,*)'>>> ERROR >>> Unknown signal in JSE control table'
	   return
	endif
	jdt = jdt+1
	DATATYPE(jdt) = STRI(7:10)
	goto	900
 901	close(1)
C Information from .../runs/a2jse is not used for the records
cc	write(NCH)jdt			! Ordinal No. of the record
cc	write(NCH)(DATATYPE(j),j=1,jdt)
C	write(*,*)jdt
C	write(*,'(a,x,a)')(DATATYPE(j),j=1,jdt)
 902	continue

	goto(10,20,30),NBNT
 10	continue			! No boundary points given
	jnb = 18
	do	j=1,jnb
	   y  = GP2*(j-1.)/jnb
	   ys = sin(y)
	   YRB(j) = RTOR+SHIFT+ABC*(cos(y)-TRIAN*ys**2)
	   YZB(j) = UPDWN+ABC*ELONG*ys
	enddo
	goto	40
 20	continue			! Fixed boundary NBNT==1
	jnb = NBND
	goto	40
 30	continue			! NBNT>1
	jnb = NBND
 40	continue

C	write(*,*)
C	write(*,*)TIME,jnb
	write(NCH)TIME,jnb
	do 905	j2=1,jdt
	   jn = 0
	   do	j1=1,2
	      if (DATATYPE(j2) .eq. DATATABLE(j1))	then
		 jn = j1
		 goto	903
	      endif
	   enddo
 903	   continue
	   if (jn .le. 0)	then
	      write(*,*)'Unknown datatype requested: "',DATATYPE(j2),'"'
	      goto	905
	   endif
	   goto( 1, 2 ),jn
	   write(*,*)jn,char(9),DATATABLE(jn)
	   write(*,*)'Unrecognized datatype: "',DATATYPE(j2),'"'
	   goto	905
 1	   continue
C	   write(*,*)jn,char(9),DATATABLE(jn)
C	   write(*,*)(YRB(j),j=1,jnb)
	   write(NCH)(YRB(j),j=1,jnb)
	   goto	905
 2	   continue
C	   write(*,*)jn,char(9),DATATABLE(jn)
C	   write(*,*)(YZB(j),j=1,jnb)
	   write(NCH)(YZB(j),j=1,jnb)
	   goto	905
 905	continue

	return
 907	write(*,*)'>>> File "a2jse" reading error'
	stop
 908	write(*,*)'>>> File "a2jse" reading error'
	write(*,*)'    Line "',STRI(1:length(STRI)),'"'
	stop
 910	format(1p,5e13.5)
	end
C======================================================================|
