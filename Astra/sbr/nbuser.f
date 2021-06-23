C======================================================================|
C User defined functions for the NBI package
C   The input parameters CVER1,CVER2,CHOR1,CHOR2 are given in
C   the NBI configuration file for each NBI source (last line)
C======================================================================|
	double precision function NBFHZ(Z,JSRC,CVER1,CVER2) 
C====== Vertical beam power distribution in the footprint =========
C 	Is precribed by User
C	Z is normalized height in respect of the footprint center
C	Z = (H - HBEAM)/H_half_width,		|Z|<=1
C	H_half_width = CBMS4*(RBMAX-RBMIN)/2
C	JSRC is the NBI source number
C	NBFHZ [a.u.] is renormalized in the calling routine
C==================================================================
	implicit none
	integer JSRC
	double precision	Z,CVER1,CVER2
	NBFHZ=EXP(-CVER1*abs(Z)**CVER2)
	return
	end
C======================================================================|
	double precision function NBFRY(Y,JSRC,CHOR1,CHOR2)
C====== Horizontal beam power distribution in the footprint =======
C 	Is precribed by User
C	Y is normalized width in respect of the footprint center
C	Y = (R - RBtang)/R_half_width,		|Y|<=1
C	R_half_width = (RBMAX-RBMIN)/2;	RBtang= (RBMAX+RBMIN)/2
C	JSRC is	the NBI source number
C	NBFHY [a.u.] is renormalized in the calling routine
C==================================================================
	implicit none
	integer JSRC
	double precision	Y,CHOR1,CHOR2
	NBFRY =EXP(-CHOR1*abs(Y)**CHOR2)
	return
	end
C======================================================================|
	double precision	function	RIPRAD(YUPDWN,J)
C
C Ripple loss boundary R[m] for each magnetic surface
C YUPDWN [m] 	plasma midplane shift in respect to the plane
C 		of the ripple losses simmetry
C J 	- magnetic surface number
	implicit none
	integer J
	double precision	YUPDWN
C Default (No ripple losses
	RIPRAD	= 99999.
	return
	end
C=======================================================================
	double precision	FUNCTION	RNB2(X2)
C-----------------------------------------22.07.89
C           X     
C RNB2=     S du/(1+u3)
C           0
	implicit none
	double precision	X2,X
	X	=SQRT(X2)
	RNB2	=(-0.166666667*LOG((1-X+X2)/(1.+2.*X+X2))
     .		+0.57735026*(ATAN(0.57735026*(2.*X-1.))+0.52359874))
	RETURN
	END
C=================================================================
