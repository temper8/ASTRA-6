/* Example of a C subprogram for ASTRA */
/*========================================================================
void testc(double*);
void testc_(double*);
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
void testc()	{ testc_(); }
/*------------------------------------------------------------------------*/
void testc_()
{  system("echo CPU time: `top -d1 | grep kdbh4iter.exe | cut -c46-53`");
}
/*========================================================================*/
