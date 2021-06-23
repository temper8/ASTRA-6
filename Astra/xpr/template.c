/*---------------------------------------------------------------------*/
#include "ipc.head"
/*---------------------------------------------------------------------*/
  struct A_example					/* Special */
  { struct A_proc_info My;
    int	 Size;						/* Special */
    int	 nncx;						/* Special */
    double albpl;					/* Special */
    double snnbm[NC1];					/* Special */
  } *IO;
  lS = sizeof(*IO);
/*---------------------------------------------------------------------*/
#include "ipc.body"
/*---------------------------------------------------------------------*/
/* Set input, run function(sbrbody), get output */
  /* Set input to example(...) */
  IO = (struct A_example *)ShmAdr;			/* Special  */
  IO->My = My;					/* is used by ASTRA */
  /* Define more input data */				/* Special  */
//  getcwd(test,(size_t)32); printf("%d %s\n",My.OrdNr,test);
//  system("pwd");
  /*  example(...); */					/* Special  */
  /* Get more ouput from example(...) */
  /* Define more output data (optional) */		/* Special  */
/*--------------------------------------------------------------------*/
#include "ipc.tail"
/*--------------------------------------------------------------------*/
