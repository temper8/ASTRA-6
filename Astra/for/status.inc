C	The total number of arrays below cannot exceed the dimension 
C	of the array ARRNAM from the common block /NAMBUF/ in the
C	Astra model compiler (presently 512)
C
C   Here follows the list of arrays recognizable by the ASTRA model compiler
C Array name     [   ] - Comment
C AMAIN            [ ] - Mass/Mp of main ions
C AMETR            [m] - Equatorial magnetic surface size
C B0DB2            []  - <B0**2/B**2> (for Neoclassical)
C BDB0             []  - <B/B0> (for Neoclassical)
C BDB02            []  - <B**2/B0**2> (for Neoclassical)
C BMAXT            [T] - Maximum B field on the surface
C BMINT            [T] - Minimum B field on the surface
C CAR1                 - Through channel
C CAR10                - Through channel
C CAR10X               - Through channel
C CAR11                - Through channel
C CAR11X               - Through channel
C CAR12                - Through channel
C CAR12X               - Through channel
C CAR13                - Through channel
C CAR13X               - Through channel
C CAR14                - Through channel
C CAR14X               - Through channel
C CAR15                - Through channel
C CAR15X               - Through channel
C CAR16                - Through channel
C CAR16X               - Through channel
C CAR17                - Through channel
C CAR17X               - Through channel
C CAR18                - Through channel
C CAR18X               - Through channel
C CAR19                - Through channel
C CAR19X               - Through channel
C CAR1X                - Through channel
C CAR2                 - Through channel
C CAR20                - Through channel
C CAR20X               - Through channel
C CAR21                - Through channel
C CAR21X               - Through channel
C CAR22                - Through channel
C CAR22X               - Through channel
C CAR23                - Through channel
C CAR23X               - Through channel
C CAR24                - Through channel
C CAR24X               - Through channel
C CAR25                - Through channel
C CAR25X               - Through channel
C CAR26                - Through channel
C CAR26X               - Through channel
C CAR27                - Through channel
C CAR27X               - Through channel
C CAR28                - Through channel
C CAR28X               - Through channel
C CAR29                - Through channel
C CAR29X               - Through channel
C CAR2X                - Through channel
C CAR3                 - Through channel
C CAR30                - Through channel
C CAR30X               - Through channel
C CAR31                - Through channel
C CAR31X               - Through channel
C CAR32                - Through channel
C CAR32X               - Through channel
C CAR3X                - Through channel
C CAR4                 - Through channel
C CAR4X                - Through channel
C CAR5                 - Through channel
C CAR5X                - Through channel
C CAR6                 - Through channel
C CAR6X                - Through channel
C CAR7                 - Through channel
C CAR7X                - Through channel
C CAR8                 - Through channel
C CAR8X                - Through channel
C CAR9                 - Through channel
C CAR9X                - Through channel
C CD          [MA/m^2] - Driven current density			!all
C CU          [MA/m^2] - Longitudinal current density		!currents
C CUBM        [MA/m^2] - Beam driven current			!here
C CUBS        [MA/m^2] - Bootstrep current			
C CUECR       [MA/m^2] - EC driven current			!except
C CUFI        [MA/m^2] - Fast ion current			!for CUTOR
C CUFW        [MA/m^2] - FW driven current			
C CUICR       [MA/m^2] - IC driven current			!are
C CULH        [MA/m^2] - LH driven current			!assumed
C CUTOR       [MA/m^2] - Toroidal current density
C CUX         [MA/m^2] - External current density		!to be
C CV          [MA/m^2] - Current density creating MV		!parallel to B
C DIMP1        [m^2/s] - 1st impurity diffusivity (STRAHL)
C DIMP2        [m^2/s] - 2nd impurity diffusivity (STRAHL)
C DIMP3        [m^2/s] - 3rd impurity diffusivity (STRAHL)
C DRODA             [] - d\rho/da
C DRODAX            [] - d\rho/da as read from a data file
C ELON              [] - Magnetic surface elongation
C ELX               [] - ELON as taken from a data file
C EQFF        [MA/m^2] - Toroidal current density j(r,z) is:
C EQPF        [MA/m^2] - j(r,z) = R0/r*EQFF(rho) + r/R0*EQPF(rho)
C ER             [V/m] - Radial electric field
C F0                   -
C F0O                  - 
C F0X                  - 
C F1                   -
C F1O                  - 
C F1X                  - F0,  F1, etc. dummy unknowns at the current time
C F2                   - F0O, F1O, etc. -"-   --"--   at the previous time step
C F2O                  - F0X, F1X, etc. -"-   --"--   as taken from a data file
C F2X                  -                dimensionality depends on usage
C F3                   - 
C F3O                  - 
C F3X                  -
C F4                   -
C F4O                  - 
C F4X                  - 
C F5                   -
C F5O                  - 
C F5X                  - 
C F6                   -
C F6O                  - 
C F6X                  - 
C F7                   -
C F7O                  - 
C F7X                  - 
C F8                   -
C F8O                  - 
C F8X                  - 
C F9                   -
C F9O                  - 
C F9X                  - 
C FOFB              [] - <(B0/B)**2*(1.-sqrt(1-B/Bmax)*(1+.5B/Bmax)> (f_t-LLM)
C FP             [V*s] - Poloidal flux
C FPO            [V*s] - Poloidal flux from the previous time step
C FV             [V*s] - Vacuum poloidal flux
C G11            [m^2] - <(nabla(Rho))**2>*dV/dRho
C G11X           [m^2] - G11 as taken from a data file
C G22              [m] - V'/(4pi^2)*<(nabla(Rho)/r)**2>*RTOR/IPOL
C G22X             [m] - G22 as taken from a data file
C G33               [] - <(RTOR/r)^2>
C G33X              [] - G33 as taken from a data file
C GN     [10^19/s*m^2] - Electron flux density
C GNX    [10^19/m^2/s] - External particle flux density
C GRADRO            [] - <|grad(rho)|> (lateral area -> VR*GRADRO)
C IPOL              [] - Normalized poloidal current (=1+O(betator))
C IPOLX             [] - IPOL as taken from a data file
C MU                [] - Inversed safety factor (1/q)
C MUX               [] - Inversed safety factor (1/q) 
C MV                [] - Vacuum rotational transform for a stellarator (iota)
C MVX               [] - Vacuum rotational transform as taken from a data file
C NALF     [10^19/m^3] - Alpha ion density
C NDEUT    [10^19/m^3] - Deuterium ion density
C NE       [10^19/m^3] - Electron density
C NEO      [10^19/m^3] - Old electron density
C NEX      [10^19/m^3] - External electron density
C NHE3     [10^19/m^3] - He3 ion density
C NHYDR    [10^19/m^3] - Hydrogen ion density
C NI       [10^19/m^3] - Main ion density
C NIBM     [10^19/m^3] - Beam ion density
C NIO      [10^19/m^3] - Old main ion density
C NIX      [10^19/m^3] - External ion density
C NIZ1     [10^19/m^3] - Impurity ions 1 density
C NIZ2     [10^19/m^3] - Impurity ions 2 density
C NIZ3     [10^19/m^3] - Impurity ions 3 density
C NN               [ ] - Neutral relative density (NN/NN_B)
C NNBM1    [10^19/m^3] - Hot neutrals densities
C NNBM2    [10^19/m^3] - Hot neutrals densities
C NNBM3    [10^19/m^3] - Hot neutrals densities
C NTRIT    [10^19/m^3] - Tritium ion density
C PBEAM       [MW/m^3] - Total beam power source
C PBLON [10^19keV/m^3] - Longitudinal beam pressure
C PBOL1       [MW/m^3] - Bolometric power, impurity No.1 (STRAHL output)
C PBOL2       [MW/m^3] - Bolometric power, impurity No.2 (STRAHL output)
C PBOL3       [MW/m^3] - Bolometric power, impurity No.3 (STRAHL output)
C PBPER [10^19keV/m^3] - Perpendicular beam pressure
C PDE         [MW/m^3] - Artificial energy source introduced by DVE
C PDI         [MW/m^3] - Artificial energy source introduced by DVI
C PEBM        [MW/m^3] - Beam energy source into electrons
C PEECR       [MW/m^3] - Electron heating due to ECR
C PEFW        [MW/m^3] - Electron heating due to FW
C PEICR       [MW/m^3] - Electron heating due to ICR
C PELH        [MW/m^3] - Electron heating due to LHCD
C PELON [10^19keV/m^3] - Longitudinal electron pressure due to H & CD
C PEPER [10^19keV/m^3] - Perpendicular electron pressure due to H & CD
C PETOT       [MW/m^3] - Total electron power input
C PEX         [MW/m^3] - External electron power source
C PFAST [10^19keV/m^3] - Fast particle pressure
C PIBM        [MW/m^3] - Beam energy source into ions
C PIFW        [MW/m^3] - Ion heating due to FW
C PIICR       [MW/m^3] - Ion heating due to ICR
C PITOT       [MW/m^3] - Total ion power input
C PIX         [MW/m^3] - External ion power source
C PRAD        [MW/m^3] - Impurity radiation power
C PRADX       [MW/m^3] - External radiation power
C PRES  [10^19keV/m^3] - Plasma pressure
C PSXR1       [MW/m^3] - Soft X-ray power impurity No.1 (STRAHL output)
C PSXR2       [MW/m^3] - Soft X-ray power impurity No.2 (STRAHL output)
C PSXR3       [MW/m^3] - Soft X-ray power impurity No.3 (STRAHL output)
C QE              [MW] - Total electron energy flux 
C QF0         [[F0]/s] - Total flux for the equation F0
C QF1         [[F1]/s] - Total flux for the equation F1
C QF2         [[F2]/s] - Total flux for the equation F2
C QF3         [[F3]/s] - Total flux for the equation F3
C QF4         [[F4]/s] - Total flux for the equation F4
C QF5         [[F5]/s] - Total flux for the equation F5
C QF6         [[F6]/s] - Total flux for the equation F6
C QF7         [[F7]/s] - Total flux for the equation F7
C QF8         [[F8]/s] - Total flux for the equation F8
C QF9         [[F9]/s] - Total flux for the equation F9
C QI              [MW] - Total ion energy flux
C QN         [10^19/s] - Total electron flux
C RHO              [m] - Current magnetic surface equivalent radius (Rho)
C SCUBM        [N/m^3] - Momentum source density from NBI [kg m/s^2/m^3]
C SD0     [[F0]/m^3/s] - Artificial source introduced by DVF0
C SD1     [[F0]/m^3/s] - Artificial source introduced by DVF1
C SD2     [[F0]/m^3/s] - Artificial source introduced by DVF2
C SD3     [[F0]/m^3/s]  -Artificial source introduced by DVF3
C SD4     [[F0]/m^3/s] - Artificial source introduced by DVF4
C SD5     [[F0]/m^3/s] - Artificial source introduced by DVF5
C SD6     [[F0]/m^3/s] - Artificial source introduced by DVF6
C SD7     [[F0]/m^3/s] - Artificial source introduced by DVF7
C SD8     [[F0]/m^3/s] - Artificial source introduced by DVF8
C SD9     [[F0]/m^3/s] - Artificial source introduced by DVF9
C SDN    [10^19/m^3/s] - Artificial electron source introduced by DVN
C SF0TOT
C SF1TOT
C SF2TOT
C SF3TOT
C SF4TOT
C SF5TOT
C SF6TOT
C SF7TOT
C SF8TOT
C SF9TOT
C SHEAR            [ ] - Magnetic shear (replaces fnc/shearr.f and fml/shear)
C SHIF             [m] - Magnetic surface Shafranov shift (alias SHIH)
C SHIV             [m] - Magnetic surface vertical (up-down) shift
C SHX              [m] - SHIF as taken from a data file
C SLAT           [m^2] - Lateral area of a magnetic surface
C SLATX          [m^2] - Xarray for SLAT
C SNEBM  [10^19/m^3/s] - Electron source from NB
C SNIBM1 [10^19/m^3/s] - Source of ions from NB (energy EBEAM)
C SNIBM2 [10^19/m^3/s] - Source of ions from NB (energy EBEAM/2)
C SNIBM3 [10^19/m^3/s] - Source of ions from NB (energy EBEAM/3)
C SNNBM  [10^19/m^3/s] - Warm neutrals source
C SNTOT  [10^19/m^3/s] - Total particle source
C SNX    [10^19/m^3/s] - External particle source
C SQEPS            [ ] - SQRT(eps)=sqrt(a/R)=sqrt(AMETR(j)/RTOR)
C TE             [keV] - Electron temperature
C TEO            [keV] - Old electron temperature
C TEX            [keV] - External electron temperature
C TI             [keV] - Ion temperature
C TIO            [keV] - Old ion temperature
C TIX            [keV] - External ion temperature
C TN             [keV] - Neutral temperature
C TRIA             [ ] - Current magnetic surface triangularity
C TRX              [ ] - TRIA as taken from a data file
C ULON             [V] - Parallel loop voltage (<UB>/Bo=ULON)
C UPL              [V] - Toroidal loop voltage (d(FP)/dt=ULON/IPOL/G33)
C VIMP1          [m/s] - 1st impurity pinch velocity (STRAHL)
C VIMP2          [m/s] - 2nd impurity pinch velocity (STRAHL)
C VIMP3          [m/s] - 3rd impurity pinch velocity (STRAHL)
C VOLUM          [m^3] - Volume within flux surface
C VP             [m/s] - ULON/(dFP/dRho) (= c*E||/Bp) pinch velocity
C VPFP           [m/s] - Flux surface velocity (for adiabatic compression)
C VPOL           [m/s] - Poloidal rotation speed
C VPOLX          [m/s] - Poloidal rotation speed as taken from data
C VR             [m^2] - dV/dRho, main grid, current time step
C VRO            [m^2] - dV/dRho, main grid, previous time step
C VRS            [m^2] - dV/dRho, intermediate grid
C VRX            [m^2] - VR (dV/dRho) as taken from a data file
C VTOR           [m/s] - Toroidal rotation speed
C VTORX          [m/s] - Toroidal rotation speed as taken from data
C ZEF              [ ] - Plasma Z effective
C ZEF1             [ ] - 1st impurity Zeff contribution
C ZEF2             [ ] - 2nd impurity Zeff contribution
C ZEF3             [ ] - 3rd impurity Zeff contribution
C ZEFX             [ ] - External Z effective
C ZIM1             [ ] - 1st impurity Z (used in NBI package)
C ZIM2             [ ] - 2nd impurity Z (used in NBI package)
C ZIM3             [ ] - 3rd impurity Z (used in NBI package)
C ZMAIN            [ ] - Z of main ions
C End array			[Total number 274 (must be < 512)]
C
Changing the list of arrays recognizable in "equ" file
C    change also: (i) dimension of the array "STAARR" in this file,
C		 (ii) "DO" loops in the subroutine DEFARR
C
	double precision
     &		TE(NRD),   TEO(NRD),  TI(NRD),   TIO(NRD),  NE(NRD),
     1		NEO(NRD),  FP(NRD),   FPO(NRD),  MU(NRD),   ULON(NRD),
     2		CU(NRD),   VP(NRD),   GN(NRD),   SQEPS(NRD),
     3		PRES(NRD), PELH(NRD), PETOT(NRD),PELON(NRD),CULH(NRD),
     4		PITOT(NRD),SNTOT(NRD),ZEF(NRD),  PEECR(NRD),PEFW(NRD),
     5		PEICR(NRD),PIICR(NRD),PIFW(NRD), PEPER(NRD),QE(NRD),
     6		QI(NRD),   QN(NRD),   QF0(NRD),  QF1(NRD),  QF2(NRD),
     7		QF3(NRD),  QF4(NRD),  QF5(NRD),  QF6(NRD),  QF7(NRD),
     8		QF8(NRD),  QF9(NRD),
     9		CUECR(NRD),CUFW(NRD), CUICR(NRD),CUTOR(NRD),UPL(NRD),
     &		NI(NRD),   NIO(NRD),  ZMAIN(NRD),AMAIN(NRD),VRO(NRD)
	double precision
     1		VR(NRD),   SHIF(NRD), SHIV(NRD), ELON(NRD), TRIA(NRD),
     2		G11(NRD),  G22(NRD),  G33(NRD),  VPFP(NRD), IPOL(NRD),
     3		CUBS(NRD), VOLUM(NRD),CV(NRD),   DRODA(NRD),
     4		ER(NRD),   VRS(NRD),  PDE(NRD),  PDI(NRD),  SDN(NRD),
     5	        SD0(NRD),  SD1(NRD),  SD2(NRD),  SD3(NRD),  SD4(NRD),
     6		SD5(NRD),  SD6(NRD),  SD7(NRD),  SD8(NRD),  SD9(NRD),
     7		STAARR(NRD,80)
	common /A_VECTORS/
     1	    TE,   TEO,  TI,   TIO,  NE,   NEO,  FP,   FPO,  MU,   ULON,
     2	    CU,   VP,   GN,   SQEPS,PELH, PETOT,PELON,CULH, PITOT,SNTOT,
     3	    ZEF,  PEECR,PEFW, PEICR,PIICR,PIFW, PEPER,QE,   QI,   QN,
     4	    QF0,  QF1,  QF2,  QF3,  QF4,  QF5,  QF6,  QF7,  QF8,  QF9,
     5	    CUECR,CUFW, CUICR,CUTOR,UPL,  ER,	VRS,  NI,   NIO,  ZMAIN,
     6	    AMAIN,VRO,  VR,   SHIF, SHIV, ELON, TRIA, G11,  G22,  G33,
     7	    VPFP, IPOL, CUBS, VOLUM,CV,   DRODA,PRES, PDE,  PDI,  SDN,
     8	    SD0,  SD1,  SD2,  SD3,  SD4,  SD5,  SD6,  SD7,  SD8,  SD9
	equivalence     (TE(1),STAARR(1,1))
	double precision
     1		B0DB2(NRD),BDB02(NRD),BDB0(NRD), BMAXT(NRD),BMINT(NRD),
     2		FOFB(NRD),GRADRO(NRD),EQFF(NRD), EQPF(NRD), SLAT(NRD),
     3		FV(NRD),   MV(NRD),   AMETR(NRD),RHO(NRD),  SHEAR(NRD)
	common /A_EQUIL/   EQFF,  EQPF,  B0DB2, BDB02, BDB0,
     1			   BMAXT, BMINT, FOFB,  GRADRO,SLAT,
     2			   FV,    MV,    AMETR, RHO,   SHEAR
	double precision
     1		NN,   TN,   NIZ1, NIZ2, NIZ3,
     2		NALF, PRAD, VTOR, NHYDR,NDEUT,  VPOL,
     3		NTRIT,NHE3, ZIM1, ZIM2, ZIM3,
     4		PBOL1,PBOL2,PBOL3,PSXR1,PSXR2,
     5		PSXR3,ZEF1, ZEF2, ZEF3, DIMP1,
     6		DIMP2,DIMP3,VIMP1,VIMP2,VIMP3,  EXTARR(NRD,31)
	common /A_IONS/
     1		NN(NRD),   TN(NRD),   NIZ1(NRD), NIZ2(NRD), NIZ3(NRD),
     2		NALF(NRD), PRAD(NRD), VTOR(NRD),
     3		NHYDR(NRD),NDEUT(NRD),NTRIT(NRD),NHE3(NRD), 
     4		ZIM1(NRD), ZIM2(NRD), ZIM3(NRD), 
     5		PBOL1(NRD),PBOL2(NRD),PBOL3(NRD),
     6		PSXR1(NRD),PSXR2(NRD),PSXR3(NRD),
     7		ZEF1(NRD), ZEF2(NRD), ZEF3(NRD), 
     8		DIMP1(NRD),DIMP2(NRD),DIMP3(NRD),
     9		VIMP1(NRD),VIMP2(NRD),VIMP3(NRD), VPOL(NRD)
! Note! Array EXTARR is used in NBI routine
	equivalence     (NN(1),EXTARR(1,1))

Changing dimensionality of arrays STAARR, EXTARR, CAR 
change also the DO loop limits in DEFARR (file for/defarr.f)
	double precision	CAR(NRD,32)
	double precision
     1		CAR1, CAR2, CAR3, CAR4, CAR5,
     2		CAR6, CAR7, CAR8, CAR9, CAR10,
     3		CAR11,CAR12,CAR13,CAR14,CAR15,
     4		CAR16,CAR17,CAR18,CAR19,CAR20,
     5		CAR21,CAR22,CAR23,CAR24,CAR25,
     6		CAR26,CAR27,CAR28,CAR29,CAR30,
     7		CAR31,CAR32
	equivalence     (CAR1(1),CAR(1,1))
	common /A_CARS/
     1		CAR1(NRD), CAR2(NRD), CAR3(NRD), CAR4(NRD), CAR5(NRD),
     2		CAR6(NRD), CAR7(NRD), CAR8(NRD), CAR9(NRD), CAR10(NRD),
     3		CAR11(NRD),CAR12(NRD),CAR13(NRD),CAR14(NRD),CAR15(NRD),
     4		CAR16(NRD),CAR17(NRD),CAR18(NRD),CAR19(NRD),CAR20(NRD),
     5		CAR21(NRD),CAR22(NRD),CAR23(NRD),CAR24(NRD),CAR25(NRD),
     6		CAR26(NRD),CAR27(NRD),CAR28(NRD),CAR29(NRD),CAR30(NRD),
     7		CAR31(NRD),CAR32(NRD)
C
	double precision
     1		CAR1X,  CAR2X,  CAR3X,  CAR4X,
     2		CAR5X,  CAR6X,  CAR7X,  CAR8X,
     3		CAR9X,  CAR10X, CAR11X, CAR12X,
     4		CAR13X, CAR14X, CAR15X, CAR16X,
     5		CAR17X, CAR18X, CAR19X, CAR20X,
     6		CAR21X, CAR22X, CAR23X, CAR24X,
     7		CAR25X, CAR26X, CAR27X, CAR28X,
     8		CAR29X, CAR30X, CAR31X, CAR32X,
     9		F0X,  F1X,  F2X,  F3X,  F4X,
     &		F5X,  F6X,  F7X,  F8X,  F9X,
     1		MUX,  MVX,  GNX,  SNX,  PEX,
     2		PIX,  PRADX,TEX,  TIX,  NEX,
     3		CUX,  ZEFX, VRX,  SHX,  ELX,
     4		TRX,  G11X, G22X, G33X, DRODAX,
     5		IPOLX,NIX,  VPOLX,VTORX,SLATX, EXT(NRD,NARRX)
Changing the list of the arrays readable from an "exp" file "EXARNM"
C    change also: (i) "NARRX" in "parameter.inc",
C Note: The order of arrays in this COMMON block A_CARX must be the same
C	as in the BLOCKDATA statement file for/init.f
	common /A_CARX/
     1		CAR1X(NRD),  CAR2X(NRD),  CAR3X(NRD),  CAR4X(NRD),
     2		CAR5X(NRD),  CAR6X(NRD),  CAR7X(NRD),  CAR8X(NRD),
     3		CAR9X(NRD),  CAR10X(NRD), CAR11X(NRD), CAR12X(NRD),
     4		CAR13X(NRD), CAR14X(NRD), CAR15X(NRD), CAR16X(NRD),
     5		CAR17X(NRD), CAR18X(NRD), CAR19X(NRD), CAR20X(NRD),
     6		CAR21X(NRD), CAR22X(NRD), CAR23X(NRD), CAR24X(NRD),
     7		CAR25X(NRD), CAR26X(NRD), CAR27X(NRD), CAR28X(NRD),
     8		CAR29X(NRD), CAR30X(NRD), CAR31X(NRD), CAR32X(NRD),
     9		F0X(NRD),  F1X(NRD),  F2X(NRD),  F3X(NRD),  F4X(NRD),
     &		F5X(NRD),  F6X(NRD),  F7X(NRD),  F8X(NRD),  F9X(NRD),
     1		MUX(NRD),  MVX(NRD),  GNX(NRD),  SNX(NRD),  PEX(NRD),
     2		PIX(NRD),  PRADX(NRD),TEX(NRD),  TIX(NRD),  NEX(NRD),
     3		CUX(NRD),  ZEFX(NRD), VRX(NRD),  SHX(NRD),  ELX(NRD),
     4		TRX(NRD),  G11X(NRD), G22X(NRD), G33X(NRD), DRODAX(NRD),
     5		IPOLX(NRD),NIX(NRD),  VPOLX(NRD),VTORX(NRD),SLATX(NRD)
	equivalence     (CAR1X(1),EXT(1,1))
	double precision
     1		PBEAM,SNNBM,SNEBM,NNBM1,NNBM2,
     2		NNBM3,SCUBM,CUBM, PEBM, PIBM,
     3		NIBM, PBLON,PBPER,CUFI, PFAST, PBFUS,
     4		SNIBM1, SNIBM2, SNIBM3,NNBM(NRD,3),IOBMN(NRD,17)
	common /CFAST/ PFAST(NRD),PBFUS(NRD)
	common /IOBEAM/
     1		PBEAM(NRD),SNNBM(NRD),SNEBM(NRD),NNBM1(NRD),NNBM2(NRD),
     2		NNBM3(NRD),SCUBM(NRD),CUBM(NRD), PEBM(NRD), PIBM(NRD),
     3		NIBM(NRD), PBLON(NRD),PBPER(NRD),CUFI(NRD),
     4		SNIBM1(NRD), SNIBM2(NRD), SNIBM3(NRD)	
	equivalence  (NNBM(1,1),NNBM1(1)),(PBEAM(1),IOBMN(1,1))
	double precision
     1		CC,   CD,   CE,   CI,   CN,
     2		DC,   DE,   DI,   DN,   HC,
     3		HE,   HI,   HN,   PE,   PET,
     4		PI,   PIT,  SN,   SNN,  XC,
     5		XE,   XI,   XN,   DVN,  DSN,
     6		DVE,  DSE,  DVI,  DSI
	common /A_TRCOEF/
     1		CC(NRD),   CD(NRD),   CE(NRD),   CI(NRD),   CN(NRD),
     2		DC(NRD),   DE(NRD),   DI(NRD),   DN(NRD),
     3		HC(NRD),   HE(NRD),   HI(NRD),   HN(NRD),
     4		PE(NRD),   PET(NRD),  PI(NRD),   PIT(NRD),
     5		SN(NRD),   SNN(NRD),  DVN(NRD),  DSN(NRD),
     6		DVE(NRD),  DSE(NRD),  DVI(NRD),  DSI(NRD),
     7		XC(NRD),   XE(NRD),   XI(NRD),   XN(NRD)
C The length of WORK1 is set accordingly to NEUT
	double precision	WORK1(NRD,2*NRD+7),WORK(NRD,2*NRD),
     1				WORK1D(NRD*(4*NRD+7))
	common /A_WORKAR/ WORK1D
	equivalence (WORK1D(1),WORK(1,1)),
     1		    (WORK1D(NRD*(2*NRD+7)+1),WORK1(1,1))
	double precision
     1		F0(NRD),DF0(NRD),VF0(NRD),SF0(NRD),SFF0(NRD),F0O(NRD),
     2		F1(NRD),DF1(NRD),VF1(NRD),SF1(NRD),SFF1(NRD),F1O(NRD),
     3		F2(NRD),DF2(NRD),VF2(NRD),SF2(NRD),SFF2(NRD),F2O(NRD),
     4		F3(NRD),DF3(NRD),VF3(NRD),SF3(NRD),SFF3(NRD),F3O(NRD),
     5		F4(NRD),DF4(NRD),VF4(NRD),SF4(NRD),SFF4(NRD),F4O(NRD),
     6		F5(NRD),DF5(NRD),VF5(NRD),SF5(NRD),SFF5(NRD),F5O(NRD),
     7		F6(NRD),DF6(NRD),VF6(NRD),SF6(NRD),SFF6(NRD),F6O(NRD),
     8		F7(NRD),DF7(NRD),VF7(NRD),SF7(NRD),SFF7(NRD),F7O(NRD),
     9		F8(NRD),DF8(NRD),VF8(NRD),SF8(NRD),SFF8(NRD),F8O(NRD),
     &		F9(NRD),DF9(NRD),VF9(NRD),SF9(NRD),SFF9(NRD),F9O(NRD)
	double precision
     1		GF0(NRD),GF1(NRD),GF2(NRD),GF3(NRD),GF4(NRD),
     2		GF5(NRD),GF6(NRD),GF7(NRD),GF8(NRD),GF9(NRD),
     3		SF0TOT(NRD),SF1TOT(NRD),SF2TOT(NRD),SF3TOT(NRD),
     4		SF4TOT(NRD),SF5TOT(NRD),SF6TOT(NRD),SF7TOT(NRD),
     5		SF8TOT(NRD),SF9TOT(NRD),
     6		DVF0(NRD),DVF1(NRD),DVF2(NRD),DVF3(NRD),DVF4(NRD),
     7		DVF5(NRD),DVF6(NRD),DVF7(NRD),DVF8(NRD),DVF9(NRD),
     6		DSF0(NRD),DSF1(NRD),DSF2(NRD),DSF3(NRD),DSF4(NRD),
     7		DSF5(NRD),DSF6(NRD),DSF7(NRD),DSF8(NRD),DSF9(NRD),
     8		FJ(NRD,0:9) ,DFJ(NRD,0:9),VFJ(NRD,0:9),GFJ(NRD,0:9),
     9		SFJ(NRD,0:9),SFFJ(NRD,0:9),FJO(NRD,0:9),SFJTOT(NRD,0:9),
     &		DVF(NRD,0:9),DSF(NRD,0:9)
	common /A_EQS/
     1		F0,  F1,  F2,  F3,  F4,  F5,  F6,  F7,  F8,  F9, 
     2		F0O, F1O, F2O, F3O, F4O, F5O, F6O, F7O, F8O, F9O,
     3		DF0, DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9,
     4		VF0, VF1, VF2, VF3, VF4, VF5, VF6, VF7, VF8, VF9,
     5		GF0, GF1, GF2, GF3, GF4, GF5, GF6, GF7, GF8, GF9,
     6		SF0, SF1, SF2, SF3, SF4, SF5, SF6, SF7, SF8, SF9,
     7		SFF0,SFF1,SFF2,SFF3,SFF4,SFF5,SFF6,SFF7,SFF8,SFF9,
     8		SF0TOT,   SF1TOT,   SF2TOT,   SF3TOT,   SF4TOT,
     9		SF5TOT,   SF6TOT,   SF7TOT,   SF8TOT,   SF9TOT,
     &		DVF0,DVF1,DVF2,DVF3,DVF4,DVF5,DVF6,DVF7,DVF8,DVF9,
     &		DSF0,DSF1,DSF2,DSF3,DSF4,DSF5,DSF6,DSF7,DSF8,DSF9
	equivalence     ( FJ(1,0), F0(1)), (SFJTOT(1,0),SF0TOT(1))
	equivalence     (DFJ(1,0),DF0(1)), (VFJ(1,0),VF0(1))
	equivalence     (GFJ(1,0),GF0(1)), (FJO(1,0),F0O(1))
	equivalence     (SFJ(1,0),SF0(1)), (SFFJ(1,0),SFF0(1))
	equivalence     (DVF(1,0),DVF0(1)),(DSF(1,0),DSF0(1))
	double precision VOLDR(NRD),VNEWR(NRD),CUBTS(NRD),R2MTR(NRD),
     2			 G2GR2(NRD),G22G(NRD), SHIH(NRD)
	equivalence	(VRO(1),VOLDR(1)),	(VR(1),VNEWR(1)),
     1			(CUBS(1),CUBTS(1)),	(G33(1),R2MTR(1)),
     2			(G11(1),G2GR2(1)),	(G22(1),G22G(1)),
     3			(SHIF(1),SHIH(1))
