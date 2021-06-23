/* #include	<time.h>
   #include	<math.h>
   #include	<unistd.h>
   #include	<errno.h>
   #include	<string.h>
   #include	<sys/times.h>
   #include	<sys/sem.h>   */
#include	<stdio.h>
#include	<stdlib.h>
#include	<sys/types.h>
#include	<sys/ipc.h>
#include	<sys/shm.h>
extern char  A_ChNa[][16];			/* Child process ID */
extern int   A_NB1, A_ShmNum, A_ShmL[], A_ShmID[];
extern void *A_ShmAdr[];
#define NC1  A_NB1
#include "A_vars.h"					/* Optional */
/*---------------------------------------------------------------------*/
int to_template (int* N)
{   to_template_(N); }
int to_template_(int* N)
{
  struct A_template					/* IO data  */
  {
    struct A_proc_info My;
  } *IO;
  struct shmid_ds Myshmid_ds;				/* Optional */
  if (A_ShmNum < 0)	return(0);
/* Optional check
   (1) shmem segment is accessible? 
   (2) compare expected size_of_segment with existing one 
*/
  if (shmctl(A_ShmID[*N+1], IPC_STAT, &Myshmid_ds) < 0)
  {
      printf(">>> Process # %d: shmctl error >>>\n",*N+1);
      exit(1);
  }
  else if (Myshmid_ds.shm_segsz != A_ShmL[*N+1])
  {
      printf(">>> Process No.%d, \"%s\" >>>",*N,&A_ChNa[*N+1][0]);
      // printf(">>> Process No.%d, \"%s\" >>>",*N,&(*IONEUT).My.Path);
      printf(" Size of shared memory mismatch\n");
      printf("    Allocated %d != %d(found)\n",
	     A_ShmL[*N+1], Myshmid_ds.shm_segsz);
      exit(1);
  }						/* End of shmem check */
  /* Add IO data setting here */
}
/*---------------------------------------------------------------------*/
int ot_template (int* N, double* cpuse)
{   ot_template_(N, cpuse); }
int ot_template_(int* N, double* cpuse)
{
#include "A_arrs.h"					/* Optional */
  struct A_template					/* IO data  */
  {
    struct A_proc_info My;
  } *IO;
  if (A_ShmNum < 0)	return(0);
  AVARS = (struct A_vars *)A_ShmAdr[0];			/* Optional */
  AARRS = (struct A_arrs *)A_ShmAdr[1];			/* Optional */
  IO = (struct A_template *)A_ShmAdr[*N+1];		/* Example */
  *cpuse = (*IO).My.CPUse;
  /* Add IO data setting here */
}
/*---------------------------------------------------------------------*/
