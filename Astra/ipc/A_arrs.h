struct A_arrs
{         /* Common block A_VECTORS */
  int NC1;
  int Size;
  double te[NC1];
  double ti[NC1];
  double ne[NC1];
  double fp[NC1];
  double mu[NC1];
  double cu[NC1];
  double zef[NC1];
  double upl[NC1];
  double ni[NC1];
  double zmain[NC1];
  double amain[NC1];
  double vr[NC1];
  double shif[NC1];
  double elon[NC1];
  double tria[NC1];
         /* Common block A_EQUIL */
  double ametr[NC1];
  double rho[NC1];
         /* Common block A_IONS */
  double nn[NC1];
  double tn[NC1];
  double niz1[NC1];
  double niz2[NC1];
  double niz3[NC1];
  double nalf[NC1];
  double nhydr[NC1];
  double ndeut[NC1];
  double ntrit[NC1];
  double nhe3[NC1];
  double zim1[NC1];
  double zim2[NC1];
  double zim3[NC1];
  int CheckWord;
} *AARRS;
