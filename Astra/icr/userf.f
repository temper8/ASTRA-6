      FUNCTION PRONE(PSI)
C--------------------------------------------------------
C ANALYTIC ELECTRON DENSITY PROFILE
C--------------------------------------------------------
      DATA       PINT  /3.0/,       PEXT /1.0/
      DATA       EDGE  /0.25/
C--------------------------------------------------------
         PRONE = EDGE + (1.-EDGE)*(1.-PSI**PINT)**PEXT
      RETURN
      END

      FUNCTION PROTE(PSI)
C--------------------------------------------------------
C ANALYTIC ELECTRON TEMPERATURE PROFILE
C--------------------------------------------------------
      DATA       PINT  /2.0/,       PEXT /1.0/
      DATA       EDGE  /0.20/
C--------------------------------------------------------
         PROTE = EDGE + (1.-EDGE)*(1.-PSI**PINT)**PEXT
      RETURN
      END

      FUNCTION PROTI(PSI)
C--------------------------------------------------------
C ANALYTIC ION TEMPERATURE PROFILE
C--------------------------------------------------------
      DATA       PINT  /2.0/,       PEXT /1.0/
      DATA       EDGE  /0.20/
C--------------------------------------------------------
         PROTI = EDGE + (1.-EDGE)*(1.-PSI**PINT)**PEXT
      RETURN
      END

      FUNCTION PSIWFR(PSASYM,THETA,DPSDTH)
C
C =================================================================
C  This user-defined function must provide the shape of the
C  initial wavefront psi(theta) and its derivative d.psi/d.theta
C           PSIWFR = PSI(THETA)
C           DPSDTH = D(PSI)/D(THETA)
C =================================================================
      DATA     PSI0 / 0.8/,   DPSI0 /0./
C =================================================================
C
      IF(PSASYM.LE.0. .OR. PSASYM.GT.1.)  THEN
         PSIWFR = PSI0*(1.+ DPSI0*THETA*THETA)
         DPSDTH = 2.*PSI0*DPSI0*THETA
      ELSE
         PSIWFR = PSASYM
         DPSDTH = 0.
      END IF
C
      RETURN
      END

      FUNCTION CFUNJZ(NPH)
C
C =================================================================
C  This user-defined function must provide the complex Fourier
C     transform of the toroidal current distribution in the antenna
C      Input         :  NPH = Toroidal wavenumber (integer)
C      Output        :  CJFUNZ = Complex Fourier component
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
C =================================================================
      DATA ZPI /3.14159265359/
C =================================================================
C
         ZRANT = RTORUS + (RPLASM+DISTAP)*COS(ZPI*THANTN/180.)
         ZNZ   = FLOAT(NPH)/ZRANT
C
C   Single element factor
C
         ZELEFC = 1.
         ZARG = 0.5*WIDTH*ZNZ
         IF(ZARG.NE.0.)
     >      ZELEFC = SIN(ZARG)/ZARG
C
C   multipole antenna form factor
C
            ZWITH = 1.
      IF(JPOLE.GT.1)  THEN
         IF(NPH.EQ.0 .AND. ABS(EPHASE).EQ.180.)  THEN
            ZWITH = 0.
         ELSE
            ZARG = 0.5*((WIDTH + WGAP)*ZNZ - ZPI*EPHASE/180.)
            IF(ZARG.NE.0.)
     >         ZWITH = SIN(FLOAT(JPOLE)*ZARG)/(FLOAT(JPOLE)*SIN(ZARG))
         END IF
      END IF
C
         CFUNJZ = CMPLX(ZWITH*ZELEFC,0.)
C
      RETURN
      END

      FUNCTION CFUNJY(ZNY)
C
C =================================================================
C  This user-defined function must provide the complex Fourier
C     transform of the poloidal current distribution in the antenna
C      Input         :  ZNY = Poloidal wavevector c k_y/w
C      Output        :  CJFUNY = Complex Fourier component
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
C =================================================================
      DATA     ZPI    /3.14159265359/
C =================================================================
C
         ZHFH = 2.*ZPI*FREQCY*HEIGTH/2.9979E8
C
         ZARGS = (ZNY + ANTKY)*ZHFH
      IF(ZARGS.NE.0.)  THEN
         ZSS = SIN(ZARGS)/ZARGS
         ZCS = (COS(ZARGS)-1.)/ZARGS
      ELSE
         ZSS = 1.
         ZCS = 0.
      END IF
C
         ZARGD = (ZNY - ANTKY)*ZHFH
      IF(ZARGD.NE.0.)  THEN
         ZSD = SIN(ZARGD)/ZARGD
         ZCD = (COS(ZARGD)-1.)/ZARGD
      ELSE
         ZSD = 1.
         ZCD = 0.
      END IF
C
      IF(JALIM.EQ.1)  THEN
C
C  Feeders up and down, equatorial short
C
         ZEPH = (180.-EPHASE)*ZPI/360.
         ZFJ = (ZCS+ZCD)*SIN(ZEPH) - (ZSS+ZSD)*COS(ZEPH)
C
      ELSE
C
C  Shorts up and down, equatorial feeder
C
         ZKYH = ANTKY*ZHFH
         ZFJ = (ZSS + ZSD)*COS(ZKYH) - (ZCS - ZCD)*SIN(ZKYH)
C
      END IF
C
         CFUNJY = CMPLX(ZFJ,0.)
C
      RETURN
      END
