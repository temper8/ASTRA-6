#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<sys/ipc.h>
#include	<sys/shm.h>
extern int  A_NB1, A_ShmNum;
extern int  A_Nsems, A_SemID, A_ShmL[], A_ShmID[];	/* Optional */
extern void *A_ShmAdr[];
#define NC1 A_NB1
#include "A_vars.h"
//double wa[1][1];
/*---------------------------------------------------------------------*/
int to_scope (
	      int* Nord, 		/* SCoPE ordinal number */
	      int* Ngrd, 		/* NA1 */
	      int* N1, 			/* NRD */
	      double* atime, 		/* TIME */
	      double* wa		/* work array */
	     )
{   to_scope_(Nord, Ngrd, N1, atime, wa); }
int to_scope_(int* Nord, int* Ngrd, int* N1, double* atime, double* wa)
{ int j, i;
#include "A_SCoPE.h"
  struct shmid_ds Myshmid_ds;
  if (A_ShmNum < 0)	return(0);
  j = *Nord+1;
  if (shmctl(A_ShmID[j], IPC_STAT, &Myshmid_ds) < 0)
  {
      printf(">>> Process # %d: shmctl error >>>\n",j);
      exit(1);
  }
  else if (Myshmid_ds.shm_segsz != A_ShmL[j])
  {
      printf(">>> Process No.%d, \"%s\" >>>",*Nord,&(*ASCOPE).My.Path);
      printf(" Size of shared memory mismatch\n");
      printf("    Allocated %d != %d(found)\n",
	     A_ShmL[j], Myshmid_ds.shm_segsz);
      exit(1);
  }
  /* printf("\n\nSizeof A_proc_info   %d bytes,   ",sizeof(My));
  printf("Sizeof ShMem[%d]   %d bytes, \n", i, Myshmid_ds.shm_segsz);  */
  ASCOPE = (struct A_scope *)A_ShmAdr[j];

// printf(">>> Process No.%d, \"%s\"\n",*Nord,(*ASCOPE).My.Path);
//  wa[0][0] = wa;
  (*ASCOPE).SCoPE_in.Nrho = *Ngrd;
  (*ASCOPE).SCoPE_in.time = *atime;
  for (j=0; j < NC1; j++)
  {
     ASCOPE->SCoPE_in.rho[j] = *(wa+j);
     ASCOPE->SCoPE_in.T_e[j] = *(wa+j+*N1);
     ASCOPE->SCoPE_in.T_i[j] = *(wa+j+(*N1)*2);
     ASCOPE->SCoPE_in.n_e[j] = *(wa+j+(*N1)*3);
     ASCOPE->SCoPE_in.n_i[j] = *(wa+j+(*N1)*4);
     ASCOPE->SCoPE_in.cc[j]  = *(wa+j+(*N1)*5);
     ASCOPE->SCoPE_in.jni[j] = *(wa+j+(*N1)*6);
     ASCOPE->SCoPE_in.pr[j]  = *(wa+j+(*N1)*7);
  }
/* printf("   %d %e %e  After setting shmem\n",ASCOPE->SCoPE_in.Nrho,
   ASCOPE->SCoPE_in.n_e[9],ASCOPE->SCoPE_in.n_e[0]); */
//printf(">>> Process No.%d\n >>>",(*ASCOPE).SCoPE_in.Nrho);
//printf(" %g  %g\n >>>",(*ASCOPE).SCoPE_in.T_e[0],(*ASCOPE).SCoPE_in.n_e[0]);
}
/*---------------------------------------------------------------------*/
int ot_scope (int* Nord, double* cpuse, double* work, int* ldin, int* nr)
{   ot_scope_(Nord, cpuse, work, ldin, nr); }
int ot_scope_(int* Nord, double* cpuse, double* work, int* ldin, int* nr)
{ int j, l, lout;
#include "A_SCoPE.h"
  if (A_ShmNum < 0)	return(0);
  ASCOPE = (struct A_scope *)A_ShmAdr[*Nord+1];
  *cpuse = (*ASCOPE).My.CPUse;
  *work  = (*ASCOPE).SCoPE_in.time;
  *nr    = (*ASCOPE).SCoPE_in.Nrho;
  l = (sizeof((*ASCOPE).SCoPE_in)-sizeof(int))/sizeof(double);
  *ldin = l;
  lout = sizeof((*ASCOPE).SCoPE_out)/sizeof(double);
  for (j=0; j < l-1;  j++) *(work+j+1) = *(&((*ASCOPE).SCoPE_in.rho[0])+j);
  for (j=0; j < lout; j++) *(work+j+l) = *(&((*ASCOPE).SCoPE_out.R_0)+j);
  /* Presently: 
     *ldin = l = 8*NC1+1 = 8*NB1+1
     work(0) = time
     work(1:NC1) = rho(1:NC1)
     work(NC1+1:2*NC1) = T_e(1:NC1)
     ...
     work(7*NC1+1:8*NC1) = pr(1:NC1)
     work(l+1:l+lout+1) = SCoPE_output
   */
/*  printf("Size of SCoPE output %d bytes,  %d\n",lout*sizeof(double),lout);
    printf("Size of SCoPE output %d*sizeof(double)\n",lout);
    printf("Exit ot_scope: 1st element %g,  last defined element %g\n",
       (*ASCOPE).SCoPE_out.R_0, *(&((*ASCOPE).SCoPE_out.R_0)+4+25*NC1+*nr));
    printf("Exit ot_scope: 1st element %g,  last defined element %g\n",
        *work,*(work-1+lout-NC1+*nr));
*/
  return(0);
}
/*---------------------------------------------------------------------*/
