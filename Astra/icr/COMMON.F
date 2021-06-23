C ----------------------------------------------------------------------
C |            INPUT-OUTPUT COMMON BLOCKS FOR RAYIC                    |
C ----------------------------------------------------------------------
C  Parameters
C -----------------------------
C  Maximum number of ion species
      PARAMETER (NSPION =   3)
C  Maximum number of radial mesh points
      PARAMETER (NPSIPT =  51)
C  Maximum number of rays in the poloidal plane
      PARAMETER (NTHRAY =  21)
C  Maximum number of toroidal modes
      PARAMETER (NPHRAY =  21)
C ----------------------------------------------------------------------
C  MHD Cconfiguration
C -----------------------------
      COMMON /CMHDPS/ RTORUS, RPLASM, BZERO,  TCURR,
     +                RMAJ(NPSIPT),   TRIANG(NPSIPT), ELONG(NPSIPT),
     +                NMHD,   NEQUIL
C ----------------------------------------------------------------------
C  Plasma parameters
C -----------------------------
      COMMON /CPLPAR/ DENC,   TEMPEC, ATM(NSPION),    AZI(NSPION),
     +                DENIC(NSPION),  TEMPIC(NSPION), ACONC(NSPION),
     +                NSPEC,  MAINSP
C ----------------------------------------------------------------------
C  Profile tables
C -----------------------------
      COMMON /CPROIC/ PPSI(NPSIPT),   PNE(NPSIPT),    PTE(NPSIPT),
     +                PNI(NPSIPT,NSPION),     PTI(NPSIPT,NSPION),
     +                PJAVG(NPSIPT),  NPROF,  IPROEQ
C ----------------------------------------------------------------------
C  Frequency and spectrum
C -----------------------------
      COMMON /CPOWIC/ FREQCY, POWER,  POYNT(NPHRAY,NTHRAY),   
     +                THSTRT(NTHRAY), PWCPL(NPHRAY),  PWTRA(NPHRAY),
     +                NPHI(NPHRAY),   NUMPHI, NPHI1,  JUMPHI, NTHIN,
     +                MAXHRM, MAXREF, ISMOOT
C ----------------------------------------------------------------------
C  Antenna
C -----------------------------
      COMMON /CANTIC/ WIDTH,  HEIGTH, DISTAP, DISTAW, WGAP,   THANTN,
     +                EPHASE, DPHASE, ANTKY,  NANTS,  JPOLE,  JALIM
C ----------------------------------------------------------------------
C  Control of the output
C -----------------------------
      COMMON /CPRTIC/ DEYKPR, IPRINT, JOUTA,  IGRAPH, IPLOPR,
     +                IPLOXY(2,7)
C ----------------------------------------------------------------------
C  Output Tables
C -----------------------------
      COMMON /OUTPIC/ PWE(NPSIPT),            PWI(NPSIPT,NSPION),
     +                PWPARI(NPHRAY,NSPION),  PWPARE(NPHRAY),
     +                PWPARB(NPHRAY),         PWTOTI(NSPION),
     +                PWIBWS, PWTOTE, PWEDGE
C
C ----------------------------------------------------------------------
C |                  PRIVATE COMMONS FOR RAYIC                         |
C ----------------------------------------------------------------------
C
      COMMON /CCONST/ PI
      COMMON /PROTAB/ TQR(NPSIPT,4),  TNE(NPSIPT,4),  TTE(NPSIPT,4),
     +                TNI(NPSIPT,4,NSPION),
     +                TTIX(NPSIPT,4,NSPION),  TTIZ(NPSIPT,4,NSPION)
      COMMON /COMEQL/ BX0(NPSIPT,4),  BX2(NPSIPT,4),  BZ0(NPSIPT,4),
     +                VOLUME(NPSIPT), PQR(NPSIPT),    PSIMHD(NPSIPT),
     +                POLASP(2),      RSHIFT, DTRIAN, DPSMHD
      COMMON /CPARMX/ DENS,   DERDEN, TEMPE, DERTEM,  DENSI(NSPION),
     +                DERDNI(NSPION), TEMPIX(NSPION), TEMPIZ(NSPION),
     +                DERTIX(NSPION), DERTIZ(NSPION),
     +                BTOT,   BTOR,   BPOL
      COMMON /CADIMX/ OPI2(NSPION),   OHI(NSPION),    OPCYI2(NSPION),
     +                OMCYI(NSPION),  BETAI(NSPION),  OXI(NSPION),
     +                XI(5,NSPION),   ZVTHI2(NSPION),
     +                ZVTHE2, OPE2,   OHE,    OPCYE2, XE,     BETAE
      COMMON /COINIT/ ZPSINI(NTHRAY), ZTHINI(NTHRAY), ZKPSIN(NTHRAY),
     +                ZKTHIN(NTHRAY), INZWFR(NTHRAY)
      COMMON /CCOORD/ EIKON,  PSI,    THETA,  UX,     UZ,     URHS,
     +                XU,     ZU,     ACOSTH, ASINTH, ACOSTQ, ASINTQ
      COMMON /COMTRD/ AJ,     ANTAU,  AG,     ANTAU2, DXPSI,  DXTH,
     +                DZPSI,  DZTH,   DJP,    DJT,    DGP,    DGT,
     +                DNTAUP, DNTAUT, DOMPSI, DOMTH,  DQPSI,  DQTH
      COMMON /CWAVEC/ AKPSI,  AKTH,   ANPHI,  PKPSI,  PKTAU,  PKPHI,
     +                PKHETA, PKX,    PKZ,    PKPAR,  PKPAR2, PKP,
     +                PKP2
      COMMON /CWVDER/ DKPSIP, DKPSIT, DKHEP,  DKHET,  DKPZP,  DKPZT,
     +                DKP2P,  DKP2T
      COMMON /CDIREL/ AR,     AL,     AS,     ARWARM, ALWARM, ATWARM,
     +                ANR,    ANL,    ANS,    HA,     HB,     HC
      COMMON /CHDISP/ DISPH,  DHKP2,  DHKZ,   DHOM,   DHNE,   DHTE,
     +                DHNI(NSPION),   DHTX(NSPION),   DHTZ(NSPION)
      COMMON /CPOWER/ ZEXPI(2,NSPION),        GAMI(NSPION),   GAME(4),
     +                EPOL(3),        ZEXPE,  ZEXPTT, POYNTN, PWX,
     +                GAMMA,  EPOLX,  EPOLY
      COMMON /ICONTR/ ISTOP,  ISTART, INDGR,   ISXZ,   LEFT,   IPWDP,
     +                IPWDPE, IPWDPI, KRES,    IASYM,  NHLAST, NUREFL,
     +                ITH,    INPHI,  IXIEQ(NPSIPT)
      COMMON /AUXIL1/ PZL(NSPION),    ACCUR,  DSTEPH, DPSI,   AKRINI,
     +                AKXINI, DKDRIN, EIKREF, PWINIT, POWOLD, PWCONF,
     +                XIJUMP, HH,     PKPMIN, PSIMAX, ZDIVY,  ZNUMY,
     +                ZERR
      COMMON /AUXIL2/ UKZERO, ULENGT, URTOR,  URPLAS, XMIN,   XMAX,
     +                URSHFT, BAXIS,  RAXMAG, DXIEQ,  PSASYM, ZDEYKP
      COMMON /AUXIL3/ ATZI(NSPION),   ATOPI(NSPION),  PRECDY(5),
     +                UXRES(20),      UXCOF(20),      UXCYCL(2,NSPION),
     +                NRES,           NCOF,           NCYCLR(2,NSPION)
      COMMON /AUXIL4/ HTOPT1, HTOPT2, REFTWH, ABSTWH, TRATWH, UPWABE,
     +                PWMODC, UPWABI(NSPION),
     +                MINOR(NSPION),  IRES(NSPION),   IHARM(NSPION),
     +                IDOPPL, IMIN,   ITWIHR, MHARM
C
C ----------------------------------------------------------------------
C |         COMMON BLOCKS FOR THE GRAPHICAL OUTPUT                     |
C ----------------------------------------------------------------------
C
      PARAMETER (NSEQMX = 200)
C
      COMMON /STABL1/ SEIKON(NSEQMX,NTHRAY),  SPSI(NSEQMX,NTHRAY),
     +                SXU(NSEQMX,NTHRAY),     SZU(NSEQMX,NTHRAY),
     +                SKP(NSEQMX,NTHRAY),     SKPAR(NSEQMX,NTHRAY),
     +                SKR(NSEQMX,NTHRAY),     SKTH(NSEQMX,NTHRAY),
     +                SKX(NSEQMX,NTHRAY),     SKZ(NSEQMX,NTHRAY),
     +                SPOWER(NSEQMX,NTHRAY)
      COMMON /SINTGS/ SXMIN,  SXMAX,  SYMIN,  SYMAX,  AXVERT, AXHORZ,
     +                ISIMPL(NTHRAY), NGR,    NPS,    NXVARB, NYVARB,
     +                KGRAPH
