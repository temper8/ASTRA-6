/*---------------------------------------------------------------------*/
#include "ipc.head"
/*---------------------------------------------------------------------*/
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
  lS = sizeof(*IONEUT);
/*---------------------------------------------------------------------*/
#include "ipc.body"
/*---------------------------------------------------------------------*/
  IONEUT = (struct A_neutIO *)ShmAdr;
  IONEUT->My = My;
/* printf("neutex.c:   %d   %d   %d  %g   %g\n",
       (AVARS->na1),
       (AVARS->nab),
       (AVARS->nb1),
       (*AARRS).ne[0],
       (*AARRS).ne[AVARS->na1]);
*/
  /* Call Fortran function */
  a_neutex_(
       /* input */
	  &(AVARS->abc),
          &(AVARS->amj),
          &(AVARS->na1),
          &(AVARS->nab),
          &NC1,
          &((*IONEUT).nncx),
          &(AVARS->nncl),
          &(AVARS->nnwm),
          &(AVARS->encl),
          &(AVARS->enwm),
          &(AARRS->te),
          &(AARRS->ti),
          &(AARRS->ne),
          &(AARRS->ni),
          &(AARRS->amain),
          &(AARRS->ametr),
          &((*IONEUT).snnbm),
       /* output */
          &((*IONEUT).nn),
          &((*IONEUT).tn),
          &((*IONEUT).albpl)
      );
  /* printf("Process # %d finished,   time = %g,   T_n(0) = %g\n",
	 		      My.OrdNr, AVARS->time, AARRS->tn[0]);   */
/*------------------------- Specific part ended ----------------------*/
#include "ipc.tail"
/*---------------------------------------------------------------------*/
