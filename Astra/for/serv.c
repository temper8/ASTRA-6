#include	<unistd.h>
#include	<time.h>
#include	<string.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<unistd.h>
#include	<errno.h>
#include	<sys/times.h>
#include	<sys/types.h>
#include	<sys/ipc.h>
#include	<sys/sem.h>
#include	<sys/shm.h>
void	getid (double*,int*,int*);	void	getid_(double*,int*,int*);
void	getval (int*,double*,double*);	void	getval_(int*,double*,double*);
void	ftdrv(double*);			void	ftdrv_(double*);
void	getgeo();			void	getgeo_();
void	drawcnf (char *file, int *type);
void	drawcnf_(char *file, int *type);
void	linter (double*,double*,int*,int*);
void	linter_(double*,double*,int*,int*);
void	timdat (int*,int*,int*,int*,int*,int*);
void	timdat_(int*,int*,int*,int*,int*,int*);
void	stcopy (char*,char*,int);
void	iroundA(char*,int*,int*,int*);
void	colovm (int*);
/* void	num2str (double*,char*,int*);
   void	num2str_(double*,char*,int*);
   void	numstrA (double,char*,int);  */
/**************************************************************************/
#include "A_X_common.h"
/* Format #include <A_X_common.h> requires option "-Ifor/" in cc command line */
/**************************************************************************/
void drawcnf (file,type)	char *file;	int *type;
{    drawcnf_(file,type);	}
void drawcnf_(file,type)	char *file;	int *type;
{  int	i, ierr, nSHOT, nEDIT=0, Ndim=750, NdGC=40, NGC;
   int	ixbeg[40], lenix[40], valix[40], EraseColor=31;
   char	GCnam[40];
   float	xyGC[750][2];

   colovm(&EraseColor); nSHOT = 16982;
/* kkGCd0 ( &ierr, "AUGD", "YGC", &nSHOT, &nEDIT, &Ndim, &xyGC[0][0],
  		 &NdGC, &NGC, &ixbeg[0], &lenix[0], &valix[0], &GCnam[0] );
   drawvm(&id, is, is+1, is, is+1);
*/ return;
}
/***************** Get ID of calling function and its parent  *************/
void ftdrv (double* var)
{    ftdrv_(var);	}
void ftdrv_(double* var)
{	/* No more than 16 functions for no more than 2200 calling each */
  static double	*varid[16][2200];
	static	int	*callid[16];
	int		i, j;
	printf("var = %d,  *var = %g\n", var, *var);
	return;
}
/***************** Get ID of calling function and its parent  *************/
void getid (var, mediator, id)	double *var;	int *mediator, *id;
{    getid_(var, mediator, id);	}
void getid_(var, mediator, id)	double *var;	int *mediator, *id;
{	/* No more than 16 functions for no more than 2200 calling each */
	static	double	*varid[16][2200];
	static	int	*callid[16];
	int		i, j;
	if   ( *mediator == 0 )
	     { for ( i=0; i<=15; ) { if ( callid[i] == 0) break; i++; }
		if ( i==16 )	goto Message1;
		callid[i] = mediator;	/* identify & store calling function */
		*mediator=i+1;
	     }				/* loop over calling functions */
	if   ( *mediator < 0 || *mediator >= 17)	goto Message2;

	i = *mediator-1;
	for (j=0; j<=2109;) { if (varid[i][j] == var) { *id=++j; return; }
	if ( varid[i][j] == 0  ) break; j++;
			  }	if ( j==2200 )	goto Message3;
	varid[i][j] = var; *id = ++j;
	return;

Message1:  printf(">>> GetID:  > 16 calling functions\n");	  return;
Message2:  printf(">>> GetID:  Illegal 2nd parameter\n"); *id=-2; return;
Message3:  /* Too many varibles > 2200 */		  *id=-1; return;
}
/****************** Get ID of calling function and its parent  *************/
void getval (shift, arr, value)	double *value, *arr;	int *shift;
{    getval_(shift, arr, value);	}
void getval_(shift, arr, value)	double *value, *arr;	int *shift;
{	*value = *(arr+*shift);
	/* printf("Shift = %d,    -> = %lf\n",*shift,*value); */
	return;
}
/**************************************************************************/
void linter (fun, arg, nd1, na1) double *fun, *arg; int *nd1, *na1;
{	linter_(fun, arg, nd1, na1);	}
void linter_(fun, arg, nd1, na1) double *fun, *arg; int *nd1, *na1;
{ int	j;	double y,x,y1;
  if (*nd1 >= *na1)	return;
  j = *nd1;  y1 = *(arg+*na1-1)-*(arg+*nd1-1); 
  for ( j = *nd1-1; j <= *na1-1; j++ )
    {	x = *(arg+j);
        y = *(fun+*na1-1)*(x-*(arg+*nd1-1))-*(fun+*nd1-1)*(x-*(arg+*na1-1));
	*(fun+j) = y/y1;	/* printf("%d   %f   %f\n",j,*(fun+j),y); */
    }
}
/**************************************************************************/
void timdat (iy, im, id, ih, mi, is)    int *iy, *im, *id, *ih, *mi, *is;
{	timdat_(iy, im, id, ih, mi, is);	}
void timdat_(iy, im, id, ih, mi, is)	int *iy, *im, *id, *ih, *mi, *is;
{
  long		ltime;
  struct	tm *td;
  time_t	timer;
	ltime = time(&timer);		td = localtime(&timer);
	*iy = (*td).tm_year;		*im = ++(*td).tm_mon;
	*id = (*td).tm_mday;		*ih = (*td).tm_hour;
	*mi = (*td).tm_min;		*is = (*td).tm_sec;
/*	  printf("Longtime %d\n", ltime);
	  printf("Localtime %s\n",asctime(td));
	  printf("Year %d Month %d Day %d Hour %d Min %d Sec %d \n",
	    *iy,*im,*id,*ih,*mi,*is);		*/
}
/****************************************************************************/
void stcopy (s1,s2,n) int n; char *s1, *s2;
/* Copy n characters from string s2 to string s1 
   equivalent to	{ strncpy(s1, s2, n);	s1[n] = '\0'; }  */
{	int	i;
	for(i=0; i <= n-1 ; ++i)	*(s1+i) = *(s2+i);
	*(s1+i) = '\0';
}
/**************************************************************************/
numstrA(x, str, l)	char *str; double x; int l;
/* Inpit:
	x   - number
	l   - string length
  Output:
	str - string of the length "l"
*/
{	static int      n = 12;
	int             i, ne, is, ll, k, i5;
	char            ch[30];

	for (i = 0; i < l; i++) str[i] = ' ';
	str[l] = '\0';
	if (x == 0.) { str[l - 1] = '0'; return (0); }
	sprintf(ch, "%+.*e", n, x);
	i = sscanf((ch + n + 4), "%d", &ne);
	if (ch[0] == '-') is = 1;
	else		  is = 0;
	ll = l - is;
	if (ll <= 0) { str[l - 1] = '*'; return (1); }
	k = n + 2;
	while (ch[k] == '0') k--;
	ch[0] = ch[1];
	for (i = 1; i < k; ch[i] = ch[i + 2], i++); k--;
	ch[k] = '\0'; ne += (1 - k); i5 = 0;
	if (k > ll) { ne = ne + k - ll; k = ll; iroundA(ch, &k, &ne, &i5); }
	while (k) {
	    if (ne >= 0) {
		if (ne + k <= ll) {
		    for (i = l - 1; ne > 0; ne--, str[i--] = '0');
		    for (k--; k >= 0; str[i--] = ch[k--]);
		    if (is) str[i] = '-'; return (0); } 
		else {	i = 2; if (ne > 9) i++; if (ne > 99) i++;
		    if (k + i > ll) { k--;
			if (k==0) { str[l-1]='*'; return (1); }
			ne++; iroundA(ch, &k, &ne, &i5); } 
		    else { i = l - i; sprintf(str + i, "e%d", ne);
			for (--k; k >= 0; str[--i]=ch[k--]);
			if (is) str[--i] = '-';
			return (0); } } } 
	    else { if (ll > -ne) {
		    if (k - ll) { i = k + ne;
			if (i >= 0) { i = l; while (ne) 
			    { k--; i--; str[i] = ch[k]; ne++; }
			    i--; str[i] = '.';
			    while (k) { k--; i--; str[i] = ch[k]; } } 
			else { i = l; while (k) 
			    { k--; i--; str[i]=ch[k]; ne++; }
			    while (ne) { i--; str[i]='0'; ne++; }
			    i--; str[i] = '.'; }
			if (is) { i--; str[i] = '-'; }
			return (0); } 
		    else { k--; ne++; iroundA(ch, &k, &ne, &i5); } }
		else {	i = 3; if (ne < -9) i = 4; if (ne < -99) i = 5;
		    if (ll >= k + i) { i = l - i;
			sprintf(str + i, "e%d", ne);
			while (k) { i--; k--; str[i] = ch[k]; }
			if (is) { i--; str[i] = '-'; }
			return (0); }
		    k--; if (k == 0) { str[l-1] = '*'; return (1); }
		    ne++; iroundA(ch, &k, &ne, &i5);
	}   }	}
	return (0);
}
/**************************************************************************/
num2str(x1, str, l1)	char *str; double *x1; int *l1;
	{ num2str_(x1, str, l1); }
num2str_(x1, str, l1)	char *str; double *x1; int *l1;
/* the same as numstrA but with float argument and callable from FORTRAN */
{	static int      n = 12;
	int             l, i, ne, is, ll, k, i5;
	char            ch[30];
	double		x;

	l = *l1;	x = *x1;
	for (i = 0; i < l; i++) str[i] = ' ';
	str[l] = '\0';
	if (x == 0.) { str[l - 1] = '0'; return (0); }
	sprintf(ch, "%+.*e", n, x);
	i = sscanf((ch + n + 4), "%d", &ne);
	if (ch[0] == '-') is = 1;
	else		  is = 0;
	ll = l - is;
	if (ll <= 0) { str[l - 1] = '*'; return (1); }
	k = n + 2;
	while (ch[k] == '0') k--;
	ch[0] = ch[1];
	for (i = 1; i < k; ch[i] = ch[i + 2], i++); k--;
	ch[k] = '\0'; ne += (1 - k); i5 = 0;
	if (k > ll) { ne = ne + k - ll; k = ll; iroundA(ch, &k, &ne, &i5); }
	while (k) {
	    if (ne >= 0) {
		if (ne + k <= ll) {
		    for (i = l - 1; ne > 0; ne--, str[i--] = '0');
		    for (k--; k >= 0; str[i--] = ch[k--]);
		    if (is) str[i] = '-'; return (0); } 
		else {	i = 2; if (ne > 9) i++; if (ne > 99) i++;
		    if (k + i > ll) { k--;
			if (k==0) { str[l-1]='*'; return (1); }
			ne++; iroundA(ch, &k, &ne, &i5); } 
		    else { i = l - i; sprintf(str + i, "e%d", ne);
			for (--k; k >= 0; str[--i]=ch[k--]);
			if (is) str[--i] = '-';
			return (0); } } } 
	    else { if (ll > -ne) {
		    if (k - ll) { i = k + ne;
			if (i >= 0) { i = l; while (ne) 
			    { k--; i--; str[i] = ch[k]; ne++; }
			    i--; str[i] = '.';
			    while (k) { k--; i--; str[i] = ch[k]; } } 
			else { i = l; while (k) 
			    { k--; i--; str[i]=ch[k]; ne++; }
			    while (ne) { i--; str[i]='0'; ne++; }
			    i--; str[i] = '.'; }
			if (is) { i--; str[i] = '-'; }
			return (0); } 
		    else { k--; ne++; iroundA(ch, &k, &ne, &i5); } }
		else {	i = 3; if (ne < -9) i = 4; if (ne < -99) i = 5;
		    if (ll >= k + i) { i = l - i;
			sprintf(str + i, "e%d", ne);
			while (k) { i--; k--; str[i] = ch[k]; }
			if (is) { i--; str[i] = '-'; }
			return (0); }
		    k--; if (k == 0) { str[l-1] = '*'; return (1); }
		    ne++; iroundA(ch, &k, &ne, &i5);
	}   }	}
	return (0);
}
/**************************************************************************/
void iroundA(ch, k, ne, i5)	char *ch;	int *k, *ne, *i5;
{	int             i;
	i = *k - 1;
	if (ch[*k] >= 53 + *i5) { ch[i]++; *i5 = 0; if (ch[i] == 53) *i5 = 1; }
	while (ch[i] == 58) {
		if (i > 0) { ch[i] = 48; i--; ch[i]++; *i5 = 0;
			if (ch[i] == 53) *i5 = 1;
			(*k)--; (*ne)++; } 
		else { ch[i] = 49; (*ne)++; } }
	while (ch[(*k) - 1] == 48) { if ((*k) > 1) { (*k)--; (*ne)++; } }
}
/**********************  not used  ********************************/
void waitas (sec) long int *sec;
{	sleep(*sec);	}
void waitas_ (sec) long int *sec;
{	sleep(*sec);	}
/**************************************************************************/
int Data2Spline(double*g0,double*g2,double*f,double*x,
		int I1,double ga0,double ga1,double ga2,double ga3)
{
  int i,ii;
  double h,r,p,q,H,ff;
  double h1,r1,p1,q1;
  double w0,w1,w2,w3,w4,w5;
  double gl[I1],V[9*I1],*v;
  double t0,t2;
  int I;
  static  double ccr2=0.5,ccr3=0.3333333333333,ccr6=0.16666666666666;

  I	=I1-1;
  h	=x[I]-x[0];
  h	*=h;
  ga1	*=h;
  ga2	*=h*h;
  ga3	*=h*h*h;

  ii	=0;
  i	=1;
  h	=x[i]-x[ii];
  H	=ccr2*h;
  r	=1./h;
  p	=ccr6*h;
  q	=ccr3*h;

  t0	=ga1*r;
  t2	=ga3*r;
  w0	=1./(H+ga0*H+t0);
  w2	=1./(ga2*H+t2);
  v	=V;
  v[0]	=w0*t0;
  v[1]	=0.;
  v[2]	=w0*r;
  v[3]	=0.;
  v[4]	=w2*t2;
  v[5]	=-w2*p;
  v[6]	=0.;
  v[7]	=0.;
  v[8]	=0.;
  g0[0]	=w0*H*f[0];
  g2[0]	=0.;
  gl[0]	=0.;
  while(i < I){
    ff	=f[i];
    g0[i]=t0*g0[ii]+r*gl[ii];
    g2[i]=t2*g2[ii]-p*gl[ii];
    gl[i]=r*g0[ii]-p*g2[ii];
    w0	=-(t0*v[0]+r*v[6]);
    w1	=p*v[6]-t2*v[3];
    w2	=p*v[7]-t2*v[4];
    w3	=p*v[3]-r*v[0];
    w4	=p*v[4]-r*v[1];
    w5	=p*v[5]-r*v[2];
    i++;
    ii++;
    v	+=9;

    h1=h;
    r1=r;
    p1=p;
    q1=q;
    h	=x[i]-x[ii];
    r	=1./h;
    p	=ccr6*h;
    q	=ccr3*h;
    H	=ccr2*(h+h1);
    w0=1./(w0+H+ga0*H+ga1*(r1+r));
    w2=1./(w2+ga2*H+ga3*(r1+r)-w1*w1*w0);
    w3+=r1+r;
    w4+=q1+q-w3*w1*w0;
    w5	=1./(w5-w3*w3*w0-w4*w4*w2);
    g0[ii]	+=H*ff;
    g0[ii]	*=w0;
    g2[ii]	=(g2[ii]-w1*g0[ii])*w2;
    gl[ii]	=(gl[ii]-w3*g0[ii]-w4*g2[ii])*w5;
    g2[ii]	-=w4*gl[ii]*w2;
    g0[ii]	-=(w1*g2[ii]+w3*gl[ii])*w0;

    t0	=ga1*r;
    t2	=ga3*r;
    v[0]=t0*w0;
    v[1]=0.;
    v[2]=r*w0;
    v[3]=0.;
    v[4]=t2;
    v[5]=-p;
    v[6]=r;
    v[7]=-p;
    v[8]=0.;
    if(i == I){
      v[2]=0.;
      v[5]=0.;
    }
    /*
      w0 w1 w3 | t0  0  r | v0 v1 v2 | g0
      w1 w2 w4 |  0 t2 -p | v3 v4 v5 | g2
      w3 w4 w5 |  r -p  0 | v6 v7 v8 | gl
    */
    v[3]	=(v[3]-w1*v[0])*w2;
    v[6]	=(v[6]-w3*v[0]-w4*v[3])*w5;
    v[3]	-=w4*v[6]*w2;
    v[0]	-=(w1*v[3]+w3*v[6])*w0;

    v[4]	=(v[4]-w1*v[1])*w2;
    v[7]	=(v[7]-w3*v[1]-w4*v[4])*w5;
    v[4]	-=w4*v[7]*w2;
    v[1]	-=(w1*v[4]+w3*v[7])*w0;

    v[5]	=(v[5]-w1*v[2])*w2;
    v[8]	=(v[8]-w3*v[2]-w4*v[5])*w5;
    v[5]	-=w4*v[8]*w2;
    v[2]	-=(w1*v[5]+w3*v[8])*w0;

    /*
      w0	=1./w0;
      w1	=w1;
      w3	=w3;
      w2	=1./(w2-w1*w1*w0);
      w4	=w4-w3*w1*w0;
      w5	=1./(w5-w3*w3*w0-w4*w4*w2);
      r0	*=w0;
      r1	=(r1-w1*r0)*w2;
      r2	=(r2-w3*r0-w4*r1)*w5;
      r1	=r1-w4*r2*w2;
      r0	=r0-(w1*r1+w3*r2)*w0;
    */
  }
  ff	=f[i];
  g0[i]	=t0*g0[ii]+r*gl[ii];
  g2[i]	=t2*g2[ii]-p*gl[ii];
  gl[i]	=0.;
  w0	=-(t0*v[0]+r*v[6]);
  w1	=p*v[6]-t2*v[3];
  w2	=p*v[7]-t2*v[4];
  w3	=0.;
  w4	=0.;
  w5	=0.;
  /*
    w0 w1 w3 | t0  0  r | v0 v1 v2 | g0
    w1 w2 w4 |  0 t2 -p | v3 v4 v5 | g2
    w3 w4 w5 |  0  0  0 | v6 v7 v8 | gl
  */
  i++;
  ii++;
  H	=ccr2*h;
  w0	=1./(w0+H+ga0*H+ga1*r);
  w2	=1./(w2+ga2*H+ga3*r-w1*w1*w0);
  w5	=1.;
  g0[ii]	+=H*ff;
  g0[ii]	*=w0;
  g2[ii]	=(g2[ii]-w1*g0[ii])*w2;
  g0[ii]	-=w1*g2[ii]*w0;
  while(ii){
    ii--;
    i--;
    g0[ii]	+=v[0]*g0[i]+v[1]*g2[i]+v[2]*gl[i];
    g2[ii]	+=v[3]*g0[i]+v[4]*g2[i]+v[5]*gl[i];
    gl[ii]	+=v[6]*g0[i]+v[7]*g2[i]+v[8]*gl[i];
    v	-=9;
  }

  return(0);
}
/**************************************************************************/
int spfit_(double*g0,double*g2,double*f,double*x,
		int*j1,double*fa0,double*fa1,double*fa2,double*fa3)
{  int I1;
   double ga0, ga1, ga2, ga3;
   I1 = *j1;
   ga0 = *fa0;
   ga1 = *fa1;
   ga2 = *fa2;
   ga3 = *fa3;
   Data2Spline(g0,g2,f,x,I1,ga0,ga1,ga2,ga3);
}
/**************************************************************************/
int spfit(double*g0,double*g2,double*f,double*x,
		int*j1,double*fa0,double*fa1,double*fa2,double*fa3)
{
  spfit_(g0,g2,f,x,j1,fa0,fa1,fa2,fa3);
}
/**************************************************************************/
/*	printf("%d,<%.20s><%.20s>\n",spos,stri,array+ind);
	printf("<Tab>,<%.20s><%.20s>\n",stri,array+ind);
	printf("%d,<%.20s><%.20s>\n",spos,stri,array+ind);
	printf("theEvent.type = %d\n",theEvent.type);
	for (j=0; j < nclmn; j++) printf("%d, ",nwid[j]);	printf("  ");
	for (j=0; j < nclmn; j++) printf("%d, ",nsta[j]);	printf("\n");
for (j=0; j < nclmn; j++) printf("%d + %d,  ",nsta[j],nwid[j]);	printf("\n");
	for (j=0; j < nclmn; j++) printf("\n[%10.10s] ",stri+nwid[j]+1);
	printf("\n%d\n",nclmn);
	printf("Width = %d, Height = %d;  %d x %d\n",Width,Height,nlines,nclmn);
*/
/*************************************************************************/
/* void markloc(char *s);			void markloc_(char *s); */
static char whereami[120] = "Start";
static int  l2where = 0, l3where = 0;

/* For all primary subroutines (those that are called directly from "stepon")
   use "markloc". For other subroutines (called by primary) use "add2loc" */
void markloc_(s)	char s[]; 
{
  strcpy(whereami,s);
  l2where = 0;	l3where = 0;
			/* printf("markloc:  %s\n",whereami); */
  return;
}
void markloc(char *s) /*Equivalent description: void markloc_(s) char s[];*/
{
  markloc_(s);
}

/****** Use "add2loc" for 2nd level subroutines (called by primary) ******/
void add2loc_(s)	char s[]; /* This function is called under Linux */
{
  if ( l2where )	/* repeated call -> overwrite */
    {
      strcpy(whereami+l2where,s);
    }
  else			/* 1st call after markloc -> append */
    {
      l2where = 6+strlen(whereami);
      strcat(whereami,"\n     ");
      strcat(whereami,s);
    }			/* printf("add2loc:  %s\n",whereami); */
  l3where = l2where;
  return;
}
void add2loc(char *s)
{
  add2loc_(s);
}
/************ Equivalent version *************
void add2loc(s)	 char s[];  { add2loc_(s); } */

void add3loc(char *s)
{
  if ( l3where )	/* call after add2loc -> append */
    {
    strcpy(whereami+l3where,s);
    }
  else			/* 1st call after markloc -> append */
    {
      l3where = strlen(whereami);
      strcat(whereami,s);
    }
  /*   printf("add3loc:  %s\n",whereami); */
  return;
}
/* Use "add3loc" for 3rd level subroutines (called by 2nd level) */
void add3loc_(s)	char s[]; 
{
  if ( l3where )
    {
    strcpy(whereami+l3where,s);
    }
  else
    {
      l3where = strlen(whereami);
      strcat(whereami,s);
    }
  /*  printf("%s\n",whereami); */
  return;
}

void locus(void)
{  printf("After %s\n",whereami);  return;  }
void locus_(void)
{  printf("After %s\n",whereami);  return;  }

void fenvex(void);			void fenvex_(void);
void fenvdx(void);			void fenvdx_(void);
#ifdef AFENV
#include <signal.h>
#include <fenv.h>
void float_error(int);
void sigint_handler(int sig); /* prototype */
void fenvex_(void)
{
//  feenableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO );		  /*
//  feenableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); */
  /* @ Linux FE_INVALID=1,  FE_OVERFLOW=8,  FE_DIVBYZERO=4,  FE_UNDERFLOW=16
  printf("%d %d %d %d\n",FE_INVALID,FE_OVERFLOW,FE_DIVBYZERO,FE_UNDERFLOW);
  */
  signal(SIGFPE, float_error);
  signal(SIGINT, sigint_handler);
  return;
}
void fenvex(void)
{ void fenvex_(); }

void fenvdx_(void)
{
 // fedisableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO );		  /*
//  fedisableexcept( FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO | FE_UNDERFLOW); */
 // signal(SIGFPE, float_error);
  return;
}
void fenvdx(void)
{ 
  //void fenvdx_();
   }

void sigint_handler(int sig)			/* this is the handler */
{
  const int j = 47;  
  /*  printf("Not this time!\n"); */
  /* system("if ( -e ${HOME}/bin/cln ) ${HOME}/bin/cln"); */
  /* if (system("/../ESC/cln") != -1)  ifkey_(&j); */
  ifkey_(&j);
}

void float_error(int i)
{
  int	j;
  fprintf(stderr,"\n >>> Floating point exception in %s\n",whereami);
  j = 257;	ifkey_(&j);
  exit(EXIT_FAILURE);
}
#else
void fenvex_(void)
{  return;	  }
void fenvex(void)
{ void fenvex_(); }
void fenvdx_(void)
{  return;	  }
void fenvdx(void)
{ void fenvdx_(); }
#endif
/*
void fatal(char *s)
{
  fprintf(stderr,"Error %s\n",s);
  exit(EXIT_FAILURE);
}
*/
