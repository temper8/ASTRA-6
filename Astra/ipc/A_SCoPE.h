struct A_scope
{ 
  struct A_proc_info My;
  struct toscope
  { double time;
    int    Nrho;
    double rho[NC1];
    double T_e[NC1];
    double T_i[NC1];
    double n_e[NC1];
    double n_i[NC1];
    double cc[NC1];
    double jni[NC1];
    double pr[NC1];
  } SCoPE_in;
  struct otscope
  { double R_0;
    double B_0;
    double Ipl;
    double Rma;
    double Phi_edge;
    double Psi[NC1];
    double cuden[NC1];
    double safac[NC1];
    double rmin[NC1];
    double shif[NC1];
    double elon[NC1];
    double tria[NC1];
    double updn[NC1];
    double Vrho[NC1];
    double g11[NC1];
    double g22[NC1];
    double g33[NC1];
    double Ipol[NC1];
    double Slat[NC1];
    double bb0[NC1];
    double bb02[NC1];
    double b0b2[NC1];
    double bmax[NC1];
    double bmin[NC1];
    double fr_t[NC1];
    double coef1[NC1];
    double coef2[NC1];
    double coef3[NC1];
    double coef4[NC1];
    double coef5[NC1];
    double coef6[NC1];
/*
    int    N_bnd;
    double r_x;
    double z_x;
    double r_bnd[N_bnd];
    double r_bnd[N_bnd];
*/
  } SCoPE_out;
} *ASCOPE;
