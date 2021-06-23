/*---------------------------------------------------------------------*/
#include "ipc.head"
/*---------------------------------------------------------------------*/
#include "A_SCoPE.h"
  lS = sizeof(*ASCOPE);
/*---------------------------------------------------------------------*/
#include "ipc.body"
/*---------------------------------------------------------------------*/
/* Set input, run function, get output */
  ASCOPE = (struct A_scope *)ShmAdr;			/* Special  */
  ASCOPE->My = My;					/* is used by ASTRA */
  callscope_(
  /* input is defined in to_scope (ipc/io_scope.c) */
    &(NC1),				/* NB1 */
    &((*ASCOPE).SCoPE_in.Nrho),		/* NA1 */
    &((*ASCOPE).SCoPE_in.time),		/* TIME */
    &((*ASCOPE).SCoPE_in.rho),		/* RHO() */
    &((*ASCOPE).SCoPE_in.T_e),		/* TE() */
    &((*ASCOPE).SCoPE_in.T_i),		/* TI() */
    &((*ASCOPE).SCoPE_in.n_e),		/* NE() */
    &((*ASCOPE).SCoPE_in.n_i),		/* NI() */
    &((*ASCOPE).SCoPE_in.cc),		/* CC() */
    &((*ASCOPE).SCoPE_in.jni),		/* CD() */
    &((*ASCOPE).SCoPE_in.pr),		/* pressure() */
         /* output */
    &((*ASCOPE).SCoPE_out.R_0)
    );
/*  i = (*ASCOPE).SCoPE_in.Nrho; */
/*    printf("Exit scope.c   %g %g %g\n",(*ASCOPE).SCoPE_out.R_0,
      (*ASCOPE).SCoPE_out.B_0,(*ASCOPE).SCoPE_out.coef6[i-2]);    */
/*--------------------------------------------------------------------*/
#include "ipc.tail"
/*--------------------------------------------------------------------*/
