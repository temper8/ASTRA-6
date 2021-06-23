C ----------------------------------------------------------------
C   PRIVATE COMMON OF THE ANTENNA EVALUATION PACKAGE
C      (Include COMINOUT or COMPARAM)
C ----------------------------------------------------------------
C
      PARAMETER (NYDIM = 26,  NZDIM = NPHRAY)
C
      COMMON /CMUNIT/ CZERO,          CRE,            CIM
      COMMON /CSPECT/ CY(NYDIM),              CYASIM(NYDIM),
     +                CJY(NYDIM),             CJZ(51),
     +                CEXNYY(NYDIM,NTHRAY),   CZEDG(NYDIM),
     +                ANY(NYDIM),             ANZ(NZDIM),
     +                ADFLUX(NYDIM),          PWNEAR(NZDIM),
     +                YY(NTHRAY)
      COMMON /CAUXAN/ PSIINT,         ULENGR,
     +                DUAW,           DUPW,           ANUCOL,
     +                ARASYM,         ALASYM,         ASASYM,
     +                ADASYM,         BZEDGE,         BZASYM,
     +                DENEDG,         DNASYM,         ACOPE2,
     +                ACOPH2,         ACOHE,          ACHIH,
     +                IOUTA,          NUMNY,          IREANT,
     +                IBREAK
      COMMON /CANTXV/ VNZ2,           VANY,           VANY2,
     +                VR,             VL,             VS,
     +                VNR,            VNL,            VNS,
     +                VD,             DENSX,          BTOTX,
     +                NZ
