#include	<stdio.h>
#include	<stdlib.h>
#include	<sys/types.h>
#include	<sys/ipc.h>
#include	<sys/shm.h>
extern char  A_ChNa[][16];	/* Child process ID */
extern int   A_NB1, A_ShmNum, A_ShmL[], A_ShmID[];
extern void *A_ShmAdr[];
#define NC1 A_NB1
#include "A_vars.h"
/*---------------------------------------------------------------------*/
/* This function is called from Astra -> toneutex -> to_neutex
   "to_neutex" fills NEUTEX specific data in the dedicated memory segment   
   Similarly, Astra -> otneutex -> ot_neutex returns the data calculated
   by the separate process "neutex" and stored in shared memory.
*/
int to_neutex__(int* N, int* NNCX, double* SNNBM)
{   to_neutex_(N, NNCX, SNNBM); }
int to_neutex (int* N, int* NNCX, double* SNNBM)
{   to_neutex_(N, NNCX, SNNBM); }
int to_neutex_(int* N, int* NNCX, double* SNNBM)
{ int i;
  struct shmid_ds Myshmid_ds;
  struct A_neutIO
  {
    struct A_proc_info My;		/* General IO information */
    int	 Size;			/* Control: Size of the Shmem */
    int	 nncx;			/* Input:  Number of iterations */
    double snnbm[NC1];		/* Input:  Implicit source  */
    double albpl;		/* Output: Plasma albedo    */
    double nn[NC1];		/* Output: Neutral density  */
    double tn[NC1];		/* Output: Neutral temperature  */
  } *IONEUT;
  if (A_ShmNum < 0)	return(0);
  i = *N+1;
  if (shmctl(A_ShmID[i], IPC_STAT, &Myshmid_ds) < 0)
  {
      printf(">>> Process # %d: shmctl error >>>\n",i);
      exit(1);
  }
  else if (Myshmid_ds.shm_segsz != A_ShmL[i])
  {
      printf(">>> Process No.%d, \"%s\" >>>",*N,&A_ChNa[i][0]);
      // printf(">>> Process No.%d, \"%s\" >>>",*N,&(*IONEUT).My.Path);
      printf(" Size of shared memory mismatch\n");
      printf("    Allocated %d != %d(found)\n",
	     A_ShmL[i], Myshmid_ds.shm_segsz);
      exit(1);
  }
  /* printf("\n\nSizeof A_proc_info   %d bytes,   ",sizeof(My));
  printf("Sizeof ShMem[%d]   %d bytes, \n", i, Myshmid_ds.shm_segsz);  */
  IONEUT = (struct A_neutIO *)A_ShmAdr[i];
  (*IONEUT).nncx = *NNCX;
  for (i=0; i < NC1; i++)  IONEUT->snnbm[i]  = *(SNNBM+i);
}
/*---------------------------------------------------------------------*/
int ot_neutex__(int* N, double* nn, double* tn, double* albpl, double* cpuse)
{   ot_neutex_(N, nn, tn, albpl, cpuse); }
int ot_neutex (int* N, double* nn, double* tn, double* albpl, double* cpuse)
{   ot_neutex_(N, nn, tn, albpl, cpuse); }
int ot_neutex_(int* N, double* nn, double* tn, double* albpl, double* cpuse)
{ int j;
  struct A_neutIO
  {
    struct A_proc_info My;		/* General IO information */
    int	 Size;			/* Control: Size of the Shmem */
    int	 nncx;			/* Input:  Number of iterations */
    double snnbm[NC1];		/* Input:  Implicit source  */
    double albpl;		/* Output: Plasma albedo    */
    double nn[NC1];		/* Output: Neutral density  */
    double tn[NC1];		/* Output: Neutral temperature  */
  } *IONEUT;
  if (A_ShmNum < 0)	return(0);
  AVARS = (struct A_vars *)A_ShmAdr[0];
  IONEUT = (struct A_neutIO *)A_ShmAdr[*N+1];
  *albpl = IONEUT->albpl;
  *cpuse = (*IONEUT).My.CPUse;
  for (j=0; j < NC1; j++)  nn[j] = IONEUT->nn[j];
  for (j=0; j < NC1; j++)  tn[j] = IONEUT->tn[j];
/*
  printf("ot_neutex IO = %d,  %d\n", IONEUT, &(AVARS->ab));
  j = AVARS->na1;
  printf("ot_neutex IO = %d,  %d  %g\n", IONEUT, j, nn[j-1]);
*/
}
/*---------------------------------------------------------------------*/
