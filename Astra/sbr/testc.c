#include	<stdio.h>
void testc();			void testc_();
// void statushmid();		void statushmid_();
/*========================================================================
void testc(a)	double *a;	
{ testc_(a); }
void testc();	void testc_();
void testc_(a)	double *a;
{ system("pwd");
  printf("Parameter = %.2f\n",*a);
}
==========================================================================*/
void testc();	void testc_();
/*------------------------------------------------------------------------*/
void testc()
   { testc_(); }
/*------------------------------------------------------------------------*/
void testc_()
{  
/* system("echo CPU time: `top -d1 | grep kdbh4iter.exe | cut -c46-53`");*/
 system("echo CPU time: `top -d1 | grep *.exe`");
}
