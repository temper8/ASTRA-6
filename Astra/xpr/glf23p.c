/*---------------------------------------------------------------------*/
#include "ipc.head"
/*---------------------------------------------------------------------*/
#include "A_glf_IO.h"
lS = sizeof(struct A_glf_IO);  /* Equivalent use: lS = sizeof(*IOGLF); */
/*---------------------------------------------------------------------*/
#include "ipc.body"
/*---------------------------------------------------------------------*/
  IOGLF = (struct A_glf_IO *)ShmAdr;
  IOGLF->My = My;
  /* Call Fortran function */
  a_glf23p_(
       /* input */
          &(IOGLF->is),
          &(IOGLF->ie),
          &(AVARS->na1),
          &(IOGLF->na1n),
          &(IOGLF->na1e),
          &(IOGLF->na1i),
          &(AVARS->hro),
          &(AVARS->hroa),
          &(AVARS->roc),
          &(AVARS->btor),
          &(AVARS->rtor),
          &(AVARS->amj),
          &(AVARS->aim1),
          &(AARRS->ne),
          &(AARRS->te),
          &(AARRS->ni),
          &(AARRS->ti),
          &(AARRS->zef),
          &(AARRS->zim1),
          &(AARRS->amain),
          &(AARRS->mu),
          &(AARRS->rho),
          &(AARRS->ametr),
          &(AARRS->shif),
          &(AARRS->elon),
          &(IOGLF->vtor),	     //		VTOR,	 
          &(IOGLF->er),		     //		ER,	 
          &(IOGLF->nibm),	     //		NIBM,	 
          &(IOGLF->ipol),	     //		IPOL,	 
          &(IOGLF->g11),	     //		G11,	 
          &(IOGLF->vrs),	     //		VRS,	 
          &(IOGLF->gradro),	     //		GRADRO,	 
          &(IOGLF->shear),	     //		SHEAR,	 
       /* output */
          &((*IOGLF).chi),
          &((*IOGLF).che),
          &((*IOGLF).dif),
          &((*IOGLF).vin),
          &((*IOGLF).dph),
          &((*IOGLF).dpl),
          &((*IOGLF).dpr),
          &((*IOGLF).xtb),
          &((*IOGLF).egm),
          &((*IOGLF).gam),
          &((*IOGLF).gm1),
          &((*IOGLF).gm2),
          &((*IOGLF).om1),
          &((*IOGLF).om2),
          &((*IOGLF).fr1)
      );
      /*  printf("Process # %d finished,   No.48 = %g,   No.49 = %g\n",
	  My.OrdNr, (*IOGLF).chi[48], (*IOGLF).chi[49]); */
/*------------------------- Specific part ended ----------------------*/
#include "ipc.tail"
/*--------------------------------------------------------------------*/
