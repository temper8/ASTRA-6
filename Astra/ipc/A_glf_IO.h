  struct A_glf_IO
  {
    struct A_proc_info My;	/* General IO information */
    int	 Size;			/* Control: Size of the Shmem */
    int	 is;			/* Input:   */
    int	 ie;			/* Input:   */
    int	 na1n;			/* Input:   */
    int	 na1e;			/* Input:   */
    int	 na1i;			/* Input:   */
    double vtor[NC1];		/* Input:   */
    double er[NC1];		/* Input:   */
    double nibm[NC1];		/* Input:   */
    double ipol[NC1];		/* Input:   */
    double g11[NC1];		/* Input:   */
    double vrs[NC1];		/* Input:   */
    double gradro[NC1];		/* Input:   */
    double shear[NC1];		/* Input:   */
    double chi[NC1];		/* Output: chi_i  */
    double che[NC1];		/* Output: chi_e  */
    double dif[NC1];		/* Output: diff_i  */
    double vin[NC1];		/* Output: diff_i  */
    double dph[NC1];		/* Output: diff_tor  */
    double dpl[NC1];		/* Output: diff_pll  */
    double dpr[NC1];		/* Output: diff_perp  */
    double xtb[NC1];		/* Output: turb heat exchange  */
    double egm[NC1];		/* Output: egamma  */
    double gam[NC1];		/* Output: gamma_p  */
    double gm1[NC1];		/* Output: lead mode rate  */
    double gm2[NC1];		/* Output: 2nd mode rate  */
    double om1[NC1];		/* Output: lead mode freq  */
    double om2[NC1];		/* Output: 2nd mode freq  */
    double fr1[NC1];		/* Output:   */
  } *IOGLF;
