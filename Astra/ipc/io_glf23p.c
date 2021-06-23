#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<sys/ipc.h>
#include	<sys/shm.h>
extern char  A_ChNa[][16];	/* Child process ID */
extern int   A_NB1, A_ShmNum, A_ShmL[], A_ShmID[];
extern void *A_ShmAdr[];
#define NC1 A_NB1
#include "A_vars.h"
/*---------------------------------------------------------------------*/
/* This function is called from Astra -> toglf23p -> to_glf23p
   "to_glf23p" fills GLF23P specific data in the dedicated memory segment   
   Similarly, Astra -> otglf23p -> ot_glf23p returns the data calculated
   by the process "glf23p" and stored in shared memory.
*/
int to_glf23p__(int* IS, int* IE,  int* NA1, 
	       double* ER, double* NN, double* NIBM, double* SHEAR,
	       int* N)
{   to_glf23p_(IS, IE, NA1, ER, NN, NIBM, SHEAR, N); }
int to_glf23p (int* IS, int* IE,  int* NA1, 
	       double* ER, double* NN, double* NIBM, double* SHEAR,
	       int* N)
{   to_glf23p_(IS, IE, NA1, ER, NN, NIBM, SHEAR, N); }
int to_glf23p_(int* IS, int* IE,  int* NA1, 
	       double* ER, double* NN, double* NIBM, double* SHEAR,
	       int* N)
{ int i, j, k;
  struct shmid_ds Myshmid_ds;
#include "A_glf_IO.h"
  if (A_ShmNum < 0)	return(0);
  if (shmctl(A_ShmID[*N+1], IPC_STAT, &Myshmid_ds) < 0)
  {   printf(">>> Process # %d: shmctl error >>>\n",*N+1);
      exit(1);
  }
  else if (Myshmid_ds.shm_segsz != A_ShmL[*N+1])
  {   printf(">>> Process No.%d, \"%s\" >>>",*N,&A_ChNa[*N+1][0]);
      // printf(">>> Process No.%d, \"%s\" >>>",*N,&(*IONEUT).My.Path);
      printf(" Size of shared memory mismatch\n");
      printf("    Allocated %d != %d(found)\n",
	     A_ShmL[*N+1], Myshmid_ds.shm_segsz);
      exit(1);
  }
/* printf("\n\nSizeof A_proc_info   %d bytes,   ",sizeof(My));
   printf("Sizeof ShMem[%d]   %d bytes, \n", *N+1, Myshmid_ds.shm_segsz); */
  AVARS = (struct A_vars *)A_ShmAdr[0];
  k = (*AVARS).nrd;
  IOGLF = (struct A_glf_IO *)A_ShmAdr[*N+1];
  (*IOGLF).is = *IS;
  (*IOGLF).ie = *IE;
  (*IOGLF).na1n = *(NA1+9);
  (*IOGLF).na1e = *(NA1+10);
  (*IOGLF).na1i = *(NA1+11);
  for (j=0; j < NC1; j++)
  {					i = k*7;
    IOGLF->vtor[j]   = *(NN+j+i);
    IOGLF->nibm[j]   = *(NIBM+j);
    IOGLF->er[j]     = *(ER+j);		i = k*16;
    IOGLF->ipol[j]   = *(ER+j+i);	i = k*12;
    IOGLF->g11[j]    = *(ER+j+i);
    IOGLF->vrs[j]    = *(ER+j+k);
    IOGLF->shear[j]  = *(SHEAR+j);	i = k*6;
    IOGLF->gradro[j] = *(SHEAR+j-i);
  }
}
/*---------------------------------------------------------------------*/
int ot_glf23p__(int* IS, int* IE, int* N, double* cpuse, double* YY)
{   ot_glf23p_(IS, IE, N, cpuse, YY); }
int ot_glf23p (int* IS, int* IE, int* N, double* cpuse, double* YY)
{   ot_glf23p_(IS, IE, N, cpuse, YY); }
int ot_glf23p_(int* IS, int* IE, int* N, double* cpuse, double* YY)
{ int j, i;
#include "A_glf_IO.h"
  if (A_ShmNum < 0)	return(0);
  AVARS = (struct A_vars *)A_ShmAdr[0];
  IOGLF = (struct A_glf_IO *)A_ShmAdr[*N+1];
  *cpuse = (*IOGLF).My.CPUse;
  for (j=*IS-1; j <= *IE-1; j++)  
    {  i = 1;
/*     YY[j+i] = (*IOGLF).chi[j];		// alternative form    */
       YY[j+i] = IOGLF->chi[j];	   i += AVARS->nrd;	// work(j+1,1)
       YY[j+i] = IOGLF->che[j];	   i += AVARS->nrd;	// work(j+1,2)
       YY[j+i] = IOGLF->dif[j];	   i += AVARS->nrd;	// work(j+1,3)
       YY[j+i] = IOGLF->dph[j];	   i += AVARS->nrd;	// work(j+1,5)
       YY[j+i] = IOGLF->dpl[j];	   i += AVARS->nrd;	// work(j+1,6)
       YY[j+i] = IOGLF->dpr[j];	   i += AVARS->nrd;	// work(j+1,7)
       YY[j+i] = IOGLF->xtb[j];	   i += AVARS->nrd;	// work(j+1,8)
       YY[j+i] = IOGLF->egm[j];	   i += AVARS->nrd;	// work(j+1,9)
       YY[j+i] = IOGLF->gam[j];	   i += AVARS->nrd;	// work(j+1,10)
       YY[j+i] = IOGLF->gm1[j];	   i += AVARS->nrd;	// work(j+1,11)
       YY[j+i] = IOGLF->gm2[j];	   i += AVARS->nrd;	// work(j+1,12)
       YY[j+i] = IOGLF->om1[j];	   i += AVARS->nrd;	// work(j+1,13)
       YY[j+i] = IOGLF->om2[j];	   i += AVARS->nrd;	// work(j+1,14)
       YY[j+i] = IOGLF->vin[j];	   i += AVARS->nrd;	// work(j+1,4)
       YY[j+i] = IOGLF->fr1[j];		 		// not used
    }
  i = AVARS->nrd;
  if ( *IS == 1 ) {
     for (j=0; j <= 15*i; j += i)  YY[j] = 0.;
     /* another option:
        for (j=0; j <= 15*i; j += i)  YY[j] = YY[j+1]; */
     }
  if ( *IE >= AVARS->na1-2 ) {
    //      for (j=0; j <= 15*(AVARS->nrd); j += AVARS->nrd)  YY[j] = YY[j+1];
     }
/*  for (j=0; j <= 15*(AVARS->nrd); j += AVARS->nrd)  YY[j] = 0.;      */
/*  printf("ot_glf23p IO = %d,  %d\n", IOGLF, &(AVARS->ab));           */
}
/*---------------------------------------------------------------------*/
