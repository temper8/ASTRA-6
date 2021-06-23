	subroutine bln
	implicit none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	integer	IsUnstable(10),jn,i,j
	double precision	ya(10),XFA
C   If IsUnStable[i] is 
C       0 - surface is ballooning Stable; 
C   Otherwise, first 4 bits are essential: 
C       1 - surface is ballooning Untable; 
C       2 - surface is Mercier Untable; 
C       4 - large numerical error in the code; 
C       8 - no convergence in the code. 
C       double *a - Array of normalized minor radii, rho_tor 
C       0 <= a[i] <= 1 */ 
C       int *jn  Number of test points in arrays IsUnStable[], a[] 
	jn = 10
	i = 0
	do j=NA1-10,NA1
	   i = i+1
	   if (i .le. jn) ya(i) = XFA(AMETR(j))
	enddo
	call pestbal(IsUnstable,ya,jn)
	do i=1,jn
	   if(IsUnstable(i) .gt.3) write(*,100)
     >		'a(',i,')=',ya(i),'  - No convergence'
	   if(IsUnstable(i) .gt.0) write(*,100)
     >		'a(',i,')=',ya(i),'  - Unstable'
	   if(IsUnstable(i) .eq.0) write(*,100)
     >		'a(',i,')=',ya(i),'  - Stable'
	enddo
 100	format(a,i2,a,f8.3,a)
	end
