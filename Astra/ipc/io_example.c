#include	<time.h>
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
extern int  A_NB1, A_ShmNum;
extern int  A_Nsems, A_SemID, A_ShmL[], A_ShmID[];	/* Optional */
extern void *A_ShmAdr[];
#define NC1 A_NB1
/*---------------------------------------------------------------------*/
int to_example (int* N)
{   to_example_(N); }
int to_example_(int* N)
{ int j;
  if (A_ShmNum < 0)	return(0);
}
/*---------------------------------------------------------------------*/
int ot_example (int* N, double* cpuse)
{   ot_example_(N, cpuse); }
int ot_example_(int* N, double* cpuse)
{ int j;
#include "A_vars.h"					/* Optional */
#include "A_arrs.h"					/* Optional */
  struct A_example
  {
    struct A_proc_info My;
  } *IO;
  if (A_ShmNum < 0)	return(0);
  AVARS = (struct A_vars *)A_ShmAdr[0];			/* Optional */
  AARRS = (struct A_arrs *)A_ShmAdr[1];			/* Optional */
  IO = (struct A_example *)A_ShmAdr[*N+1];		/* Example */
  *cpuse = (*IO).My.CPUse;
}
/*---------------------------------------------------------------------*/
