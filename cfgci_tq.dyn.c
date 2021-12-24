/* cfgcisd(tq).c */
/* last updated: 03/06/08 */
#include "cfgci_tq.h"

/* variables private to this file */
static FILE *mcfgfile, *ciftfile;
static int nst, norb, irrep, mxkmax;
static int spin, ngrps, norbx, norbi0, norbi, nnorbx, *mcfg, nmcr[22], nmcrsum[22];
static int ncsfq1tmp, ncsfq2tmp, ncsfxstate, elmt22tot, bunlen, iwrite_elmt;
static int orbsym[MAXNBF+1], orbsymi[MAXNBF+1], nxsym[NSTMAX], sxsym[NSTMAX], exsym[NSTMAX];
static int nvrtx[MAXCFGS], xspwt[MXSPVRTX], spnhead[MXKMAX+1];
static int spvrtxmin[MXKMAX+1], spvrtxmax[MXKMAX+1];
static int jdn[9*MAXVRTX], arcsp[MXSPVRTX*2], jspdn[MXSPVRTX*2], jspup[MXSPVRTX*2];
static int mcfgln[MAXCFGS], ncsfcfg[MXKMAX+1], snglptr[MXKMAX+1];
static int nxprsym[NSTMAX], nxtrsym[NSTMAX], nxqrsym[NSTMAX];
static int ftpt[MXKMAX+1], *gbcsf, *lccsf_11; 
static int bospn[NINTMX+1], dospn[NINTMX+1], delspn[NINTMX+1], nspndif[MXKMAX+3];
static int ovrlap, switch_bra_ket, jcfg, jcfg11;
static int intgrlmo[MAXNBF+1], xorbq2[4];
static int *elmt11cnt, *elmt22cnt, len11, len22;
static double *oneint, *twoint;
static double pbm1[MXKMAX+1], pb01[MXKMAX+1], pb31[MXKMAX+1], sqt, sqt2, isqt, sqrt3, sqrt32;
static double pbm0[MXKMAX+1], pb02[MXKMAX+1], pb12[MXKMAX+1], pb32[MXKMAX+1];
static double sqm0[MXKMAX+1], sq00[MXKMAX+1], sq10[MXKMAX+1], sq20[MXKMAX+1];
static double sq01[MXKMAX+1], sq21[MXKMAX+1];
static double sqm2[MXKMAX+1], sq02[MXKMAX+1], sq12[MXKMAX+1], sq22[MXKMAX+1];
static double headcum[MXKMAX+1][MXKMAX+1], acum[MXKMAX+1][MXKMAX+1], bcum[MXKMAX+1][MXKMAX+1];
static double ccum[MXKMAX+1][MXKMAX+1];
static double *hqq, *hqp, *hqwp;
static double  hdiag, *hdiagtmp, *hdiagtmp2, **hqq2, *hqq2inv, *buff;
static double intk[NINTMX+3], oneintk[MXKMAX+1], exchint[(MXKMAX*(MXKMAX-1))/2], *exchcc, efzc;
static double mtrxelt[((NINTMX+1)>MXKMAX)?(NINTMX+2):(MXKMAX+1)], mtrxel1[NINTMX+2];
static const double zthresh=10.*DBL_EPSILON;
static struct cc_elmt1{
       int icsf;
       int jcsf;
       double elmtone;
       double acum[MXKMAX+1];
} *cc_one, *ptr_cc_one;
static struct cc_elmt2{
       int icsf;
       int jcsf;
       double drct_cc;
       double exch_cc;
} *cc_two, *ptr_cc_two;

void initci(long totalmem) 
{
 extern int my_rank;
 extern int nst, prntflag, invflag, norbi0, norbi, norb, norbx, nnorbx; 
 extern int nxsym[], sxsym[], exsym[], mxkmax, bunlen, iwrite_elmt;
 extern int mcfgln[], irrep, nstates, snglptr[];
 extern int ftpt[], ncsfcfg[], *gbcsf, *lccsf_11, nxprsym[], nxtrsym[], nxqrsym[];
 extern int jspdn[], jspup[], arcsp[], xspwt[], spvrtxmin[], spvrtxmax[], nmcr[], nmcrsum[];
 extern int nmcrt, ngrps, *mcfg, nvrtx[], jdn[], orbsymi[], spin, spnhead[], orbsym[];
 extern int intgrlmo[];
 extern int *elmt11cnt, *elmt22cnt, len11, len22;
 extern long ncsfq1;
 extern double *oneint, *twoint, *exchcc, efzc;
 extern double *hqq, *hqp, *hqwp, **hqq2, *hqq2inv, *buff; 
 extern double pbm1[], pb01[], pb31[], pbm0[], pb02[], pb12[], pb32[], sqt, sqt2, isqt, sqrt3, sqrt32;
 extern double sqm0[], sq00[], sq10[], sq20[], sqm2[], sq02[], sq12[], sq22[], sq21[], sq01[];
 extern FILE *outfile, *mcfgfile, *infofile, *ciftfile;

 int i, ig, ij, iloc, ix, ist, ichk, nnbf, nbf, ncor, len, ae, nvsub, icfg, itmp, ncsfsum;
 int inord[MAXNBF+1], iout[MAXNBF+1], nso[NSTMAX];
 int j, jg, maxnints, nocc, kl, k, l, ii, ik, jk, jl, jj, indx, ijkl;
 int lmax, nnorb, nfzc, ncfgi, mxvrtx, nelx, usesym, nelec, mxkmax2;
 int mcrloc, ncsftmp, ncsftmp2, mult, nspvrtx, nvrtxtmp, spntmp, arccum, nlanes;
 int *jup, spncum[MXSPVRTX];
 int ncsfint, ncc, kst, lst, kx;
 double mem, val2e, val, *buck2e;
 long lenbuf, meminit=0;
 FILE *mointsfile;

 void getbuckij(int i, int j, double *buck2e, FILE *moints);
 int mkhdiagmcr(int nelx, int jup[], int gbcsf[]);
 int diagcc(int nopen, double vec[]);
 
 nbf = getint(infofile, 1);
 nst = getint(infofile, 2);
 nfzc = getint(infofile, 3);
 ncor = getint(infofile, 4);
 norbi0 = ncor + getint(infofile, 5); /* number of internal orbitals */ 
 norbx = getint(infofile, 6); 
 mult = getint(infofile, 8);
 nelec = getint(infofile, 9);
 irrep = getint(infofile, 10);  
 nvsub = getint(infofile, 16);
 ncfgi = getint(infofile, 34);
 usesym = getint(infofile, 41);
 spin = mult-1;                             /* 2 * true spin */ 
 nelec += 2*ncor;
 if (ncor) ngrps = nvsub + 1;
 else ngrps = nvsub;
 norbi = norbi0 + 4;   /* includes four dummy orbitals */ 
 
 if (norbi > NINTMX) {
   (void)fprintf(outfile, "\nERROR: # of internal orbitals exceeds allowed maximum (%d)\n", NINTMX);
   (void)fprintf(outfile, "*** RECOMPILE SOURCE CODE ***\n\n"); 
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
 if (norbx > NEXTMX) {
   (void)fprintf(outfile, "\nERROR: # of external orbitals exceeds allowed maximum (%d)\n", NEXTMX);
   (void)fprintf(outfile, "*** RECOMPILE SOURCE CODE ***\n\n");
   MPI_Abort(MPI_COMM_WORLD, 1); 
 }
 
 if (!usesym) {
   (void)fprintf(outfile, "\n* warning: use of symmetry information suppressed by flag from mcrcfgs *\n"); 
   nst = 1;
   nso[0] = nbf;
 } else {
   ichk = getarray(infofile, 2, nso, nst * sizeof (int)); 
 }
  
 /* >>>>>>>>>>>>>>>>>>>>>>>> construct orbital indexing arrays <<<<<<<<<<<<<<<<<<<<<< */
 ichk = getarray(infofile, 3, inord, nbf * sizeof (int));
 if (ichk != nbf * sizeof (int)) {
   (void)fprintf(outfile, "\n*** ichk = %d for inord read ***\n", ichk);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }

#ifdef DEBUG
 (void)fprintf(outfile, "\n");
 for (i=1; i<=nbf; i++) {
    (void)fprintf(outfile, "inord[%d]=%d\n", i-1, inord[i-1]);
 } 
#endif
  
 for (i=0; i<nbf; i++) iout[inord[i]] = i;
 for (ist=0,ig=0; ist<nst; ist++) {
   for (i=0; i<nso[ist]; i++) orbsym[iout[ig++]+1] = ist;         /* orbsym is (1..nbf) */
 }

#ifdef DEBUG
 for (i=1; i<=nbf; i++) {
    (void)fprintf(outfile, "orbsym[%d]=%d\n", i, orbsym[i]);
 } 
#endif

 /* the symmetry for dummy external orbitals is set to be "0". */
 norb = norbi0 + norbx; 
 for (i=1; i<=norbi0; i++) orbsymi[i] = orbsym[norb+1-i];
 for (; i<=norbi; i++) orbsymi[i] = 0;
 orbsymi[0] = 0;

#ifdef DEBUG
 for (i=1; i<=norbi; i++) {
    (void)fprintf(outfile, "orbsymi[%d]=%d\n", i, orbsymi[i]);
 } 
#endif
 
 for (i=0; i<nst; i++) nxsym[i] = 0;
 for (i=1; i<=norbx; i++) nxsym[orbsym[i]]++;
 sxsym[nst-1] = 1;
 for (i=nst-1; i>0; i--) {
   exsym[i] = sxsym[i] + nxsym[i] - 1;
   sxsym[i-1] = exsym[i] + 1;
 }
 exsym[0] = sxsym[0] + nxsym[0] - 1;
 /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

 /* mapping orbitals from distinct row table to those in integrals */
 for (i=0; i<5; i++) intgrlmo[i] = 0; /* initialize dummy orbitals */
 for (; i<=norbi; i++) {
    intgrlmo[i] = norbx + i - 4;
 }
 
 /* >>>>>>>>>>>>>>>>>>>>>>>>>> process molecular integrals <<<<<<<<<<<<<<<<<<<<<<<<<< */ 
 nnbf = (nbf*(nbf+1))/2;
 lenbuf = (nnbf+1) * sizeof (double);
 oneint = (double *)malloc(lenbuf);
 ichk = getarray(infofile, 6, &oneint[1], nnbf * sizeof (double));
 if (ichk != nnbf * sizeof (double)) {
   (void)fprintf(outfile, "\n*** ichk = %d for oneint read ***\n", ichk);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
 meminit += lenbuf;   
  
 nnorb = (norb*(norb+1))/2;
 lenbuf = (1+(nnorb*(nnorb+1))/2) * sizeof (double);
 twoint = (double *)malloc(lenbuf);
 if (!twoint) {
    (void)fprintf(outfile, "\nAllocation Failure for ERIs\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
 }
 meminit += lenbuf;
 
 for (i=1,efzc=0.; i<=nfzc; i++) {
    ix = i + norb;
    efzc += oneint[((ix)*(ix-1))/2+ix];
 }

 mointsfile = sopen(MOINTSNAME, "rb", infofile);
 if (mointsfile == NULL) {
   (void)fprintf(outfile, "\nError opening file: %s\n", MOINTSNAME);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
  
 nocc = norb + nfzc;
 maxnints = (nocc*(nocc+1))/2; 
  
 buck2e = (double *)malloc(maxnints * sizeof (double));
 if (!buck2e) {
   (void)fprintf(outfile, "\nAllocation Error (buck2e)\n");
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
 
 /* read & process two-electron integrals */
 ijkl = 1;
 for (i=1; i<=norb; i++) {
   for (j=1; j<=i; j++) {
     getbuckij(i, j, buck2e, mointsfile);
     for (k=1,indx=0; k<=i; k++) {
       if (k==i) lmax = j;
       else lmax = k;  
       for (l=1; l<=lmax; l++) {
         val2e = buck2e[indx++];
         twoint[ijkl++] = val2e;  
       }
     }
   }
 }
 
 for (i=norb+1; i<=nocc; i++) {
   for (j=1,jl=0; j<=norb; j++) {
     getbuckij(i, j, buck2e, mointsfile);
     /* only (ij|il) can contribute */
     indx = (i*(i-1))/2;
     for (l=1; l<=j; l++) oneint[++jl] -= buck2e[indx++];
   }
   
   indx = (i*(i-1))/2 + norb;
   for (j=norb+1; j<i; j++) {
     getbuckij(i, j, buck2e, mointsfile);
     /* only (ij|ij) can contribute */
     efzc -= buck2e[indx];
     indx++;
   }

   /* j=i */
   getbuckij(i, i, buck2e, mointsfile);
   for (k=1,indx=0,kl=0; k<=norb; k++) {
     for (l=1; l<=k; l++) oneint[++kl] += 2.*buck2e[indx++];
   }
 
   indx = ((norb+1)*(norb+2))/2 - 1;
   for (k=norb+1; k<i; k++) {
     efzc += 2.*buck2e[indx];
     indx += k + 1;
   }

   /* (ii|ii) */
   indx = (i*(i+1))/2 - 1;
   efzc += buck2e[indx]/2.;

 }

 efzc += efzc;
 free(buck2e);

 /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
  
 /* >>>>>>>>>> construct hamiltonian diagonal elements & indexing arrays <<<<<<<<<<<<<<<< */
 ncsfint = getint(infofile, 12);
 nmcrt = getint(infofile, 32);
 mxkmax = getint(infofile, 33);
 if (nmcrt >= MAXCFGS) {
   (void)fprintf(outfile, "\nERROR: # of macroconfigurations exceeds allowed maximum (%d)\n", MAXCFGS);
   (void)fprintf(outfile, "*** RECOMPILE SOURCE CODE ***\n\n");
   MPI_Abort(MPI_COMM_WORLD, 1); 
 }
 lenbuf = nmcrt * sizeof (int);
 ichk = getarray(infofile, 28, nvrtx, lenbuf);
 lenbuf *= ngrps;
 mcfg = (int *)malloc(lenbuf);
 meminit += lenbuf;
 ichk = getarray(infofile, 29, mcfg, lenbuf);
 /* 30 different macroconfigurations  */
 ichk = getarray(infofile, 30, nmcr, 22 * sizeof (int)); 
 nmcrsum[0] = nmcr[0];
 for (i=1; i<22; i++) 
   nmcrsum[i] += nmcrsum[i-1] + nmcr[i];
 
 if (mxkmax > MXKMAX) {
   (void)fprintf(outfile, "\nmxkmax (%d) exceeds MXKMAX (%d).  CHANGE CODE AND RECOMPILE!\n",
    mxkmax, MXKMAX);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
  
 nspvrtx = 0;
 for (k=0; k<=mxkmax; k++) nspvrtx += (k/2) + 1;
 if (nspvrtx > MXSPVRTX) {
   (void)fprintf(outfile, "\nnspvrtx (%d) exceeds MXSPVRTX (%d).  CHANGE CODE AND RECOMPILE!\n",
    nspvrtx, MXSPVRTX);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
 if (prntflag>1 && my_rank==0) (void)fprintf(outfile, "number of spin vertices = %d\n", nspvrtx);
  
 /* initialize spin DRT */
 for (i=0,ii=0; i<nspvrtx; i++) {
   xspwt[i]=-1;
   jspdn[ii]=-1;
   jspup[ii]=-1;
   arcsp[ii]=-1;
   ii++;
   jspdn[ii]=-1;
   jspup[ii]=-1;
   arcsp[ii]=-1;
   ii++;
 }
 mxvrtx=0;
 spncum[0]=0;
 spvrtxmin[0]=0;
 spvrtxmax[0]=0;
 xspwt[0]=1;
 for (i=1; i<=mxkmax; i++) {
   nvrtxtmp=0;
   spvrtxmin[i] = spvrtxmax[i-1] + 1;
   for (j=spvrtxmin[i-1]; j<=spvrtxmax[i-1]; j++) {
  
     /* spin-increasing arc */
     spntmp=spncum[j]+1;
     for (k=0; k<nvrtxtmp; k++) {
       if (spntmp==spncum[spvrtxmin[i]+k]) {
         /* attach to existing vertex */
         jspdn[2*(spvrtxmin[i]+k)+1] = j;
         jspup[2*j+1] = spvrtxmin[i]+k;
         goto down_arc;
       }
     }
     /* add vertex */
     mxvrtx++;
     nvrtxtmp++;
     spncum[mxvrtx] = spntmp;
     jspdn[2*mxvrtx+1] = j;
     jspup[2*j+1] = mxvrtx;

     down_arc:;
     /* spin-decreasing arc */
     spntmp=spncum[j]-1;
     if (spntmp < 0) continue;
     for (k=0; k<nvrtxtmp; k++) {
       if (spntmp==spncum[spvrtxmin[i]+k]) {
         /* attach to existing vertex */
         jspdn[2*(spvrtxmin[i]+k)] = j;
         jspup[2*j] = spvrtxmin[i]+k;
         goto no_more_arcs;
       }
     }
     /* add vertex */
     mxvrtx++;
     nvrtxtmp++;
     spncum[mxvrtx] = spntmp;
     jspdn[2*mxvrtx] = j;
     jspup[2*j] = mxvrtx;
     
     no_more_arcs:;
   }
   spvrtxmax[i] = spvrtxmax[i-1] + nvrtxtmp;
   
   /* calculate arc weights */
   for (j=spvrtxmin[i]; j<=spvrtxmax[i]; j++) {
     jj = j*2;
     arccum = 0;
     for (k=0; k<=1; k++) {
       if (jspdn[jj+k]!=-1) {
         arcsp[jj+k] = arccum;
         arccum += xspwt[jspdn[jj+k]];
       }
     }
     xspwt[j] = arccum;
   }
 } 
  
 for (i=0; i<=mxkmax; i++) {  
   spntmp = i;
   spnhead[i] = -1;      /* flag */
   for (j=spvrtxmin[i]; j<=spvrtxmax[i]; j++) {
     if (spntmp==spin) {
       spnhead[i] = j;
       break;
     }
     spntmp -= 2;
   }
 }
 
 /* >>>>> precompute segment values <<<<< */
 spntmp = (spin+mxkmax)/2;       /* maximum intermediate spin */
 sqt = sqrt(2.);
 sqt2 = sqt * 0.5;
 isqt = 1./sqt;
 sqrt3 = sqrt(3.);
 sqrt32 = sqrt(3./2.);
 pbm0[0] = pbm0[1] = 0.;
 pbm1[0] = pbm1[1] = 0.;  
 pb01[0] = 0.;
 pb31[0] = sqrt3;
 pb02[0] = 0.;
 pb12[0] = pb01[1] = isqt;
 pb32[0] = sqrt32;
 pb31[1] = sqt;
 pb02[1] = 1./sqrt3;
 pb12[1] = 1./sqrt32;
 pb32[1] = sqrt(4./3.);
 for (i=2; i<=spntmp; i++) {
   pbm0[i] = sqrt((double)(i-1)/i);
   pbm1[i] = sqrt((double)(i-1)/(i+1));
   pb01[i] = sqrt((double)i/(i+1));
   pb31[i] = sqrt((double)(i+3)/(i+1));
   pb02[i] = sqrt((double)i/(i+2));
   pb12[i] = sqrt((double)(i+1)/(i+2));
   pb32[i] = sqrt((double)(i+3)/(i+2));
 }
 sqm0[0] = sqm0[1] = 0.;
 sq00[0] = 0.;
 sq10[0] = 1.;
 sq20[0] = sq10[1] = 1./sqrt3;
 sq01[0] = sq01[1] = 0.;
 sq21[0] = sqrt3/2.;
 sqm2[0] = sqm2[1] = 0.;
 sq02[0] = 0.;
 sq12[0] = 0.;
 sq22[0] = sq12[1] = 1./sqrt32;
 sq21[1] = sqrt(8.)/3.;
 sq00[1] = 1.;
 sq20[1] = sqrt(1./6.);
 sq02[1] = 0.;
 sq22[1] = sqrt(5./6.);
 for (i=2; i<=spntmp; i++) {
   val = sqrt((double)(i-1)*i);
   sqm0[i] = sqt/val;
   sqm2[i] = sqrt((double)(i-2)*(i+1))/val; 
   val = sqrt((double)i*(i+1));
   sq00[i] = sqt/val;
   sq02[i] = sqrt((double)(i-1)*(i+2))/val;
   val = sqrt((double)(i+1)*(i+2));
   sq10[i] = sqt/val;
   sq12[i] = sqrt((double)i*(i+3))/val;
   val = sqrt((double)(i+2)*(i+3));
   sq20[i] = sqt/val;
   sq22[i] = sqrt((double)(i+1)*(i+4))/val;
   sq01[i] = sqrt(i*i-1.)/i;
   val = i + 2;
   sq21[i] = sqrt(val*val-1.)/val;
 }
 
 /* >>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<< */
 
 if ((mcfgfile = sopen(MCFGNAME, "rb", infofile)) == NULL) {
   (void)fprintf(outfile, "\nError opening %s\n", MCFGNAME);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }

 if (nelec > (norbi0+2)) {
   mxkmax2 = 2*norbi0 + 4 - nelec;
 } else mxkmax2 = mxkmax;
 ncsftmp2 = 0;
 for (i=0; i<=mxkmax2; i++) {
   ncsfcfg[i] = 0;
   if (spnhead[i] != -1) {
     ncsfcfg[i] = xspwt[spnhead[i]];  
     ncsftmp2 = (ncsfcfg[i]>ncsftmp2) ? ncsfcfg[i] : ncsftmp2;
   }
 }
 if (prntflag>1 && my_rank==0) {
    (void)fprintf(outfile, 
      "maximum number of csfs per configuration (Q1) = %d\n", ncsftmp2);
    (void)fflush(outfile);
 }

 ncsftmp = ncsftmp2;
 for (; i<=mxkmax; i++) {
   ncsfcfg[i] = 0;
   if (spnhead[i] != -1) {
     ncsfcfg[i] = xspwt[spnhead[i]];        /* csfs per configuration of i open shells */
     ncsftmp = (ncsfcfg[i]>ncsftmp) ? ncsfcfg[i] : ncsftmp;
   }
 }
 if (prntflag>1 && my_rank==0) {
    (void)fprintf(outfile, 
      "maximum number of csfs per configuration (Q2) = %d\n", ncsftmp);
    (void)fflush(outfile);
 }

 
 /* array for evaluation of hwpp */
 hqq = (double *)malloc(nstates*sizeof(double));
 lenbuf = nstates * ncsftmp * sizeof(double);
 hqp = (double *)malloc(lenbuf);
 hqwp = (double *)malloc(lenbuf);
 meminit += 2 * lenbuf;
 lenbuf =  ncsftmp * sizeof(double);
 hdiagtmp  = (double *)malloc(lenbuf);
 hdiagtmp2  = (double *)malloc(lenbuf);
 meminit += 2 * lenbuf;
 if (invflag) {
    lenbuf =  ncsftmp * sizeof(double *);
    hqq2    = (double **)malloc(lenbuf);
    meminit += lenbuf;
    lenbuf =  (((ncsftmp+1) * ncsftmp)/2) * sizeof(double);
    hqq2[0] = (double *)malloc(lenbuf);
    hqq2inv = (double *)malloc(lenbuf);
    buff    = (double *)malloc(lenbuf);
    meminit += 3 * lenbuf;
    for (i=1; i<ncsftmp; i++) {
      hqq2[i] = hqq2[i-1] + i;
    }
 }
 
 /* >>>>>>>>>> external completion tables <<<<<<<<<< */
 /* unpaired two external electrons   (1 1 0 0 )  */
 nnorbx = (norbx*(norbx+1))/2;
 lenbuf = (((mxkmax-spin)/2 + 1) * nnorbx + 1) * sizeof(int);
 lccsf_11 = (int *)malloc(lenbuf);  
 if (!lccsf_11) {
   (void)fprintf(outfile, "\nAllocation failure (lccsf_11)\n");
   MPI_Abort(MPI_COMM_WORLD, 1);
 }
 meminit += lenbuf;

 lccsf_11[0] = 0;
 itmp = 0;
 if (spin%2) iloc = 1; 
 else iloc = 0;
 for (i=(spin/2); i<=(mxkmax/2); i++) {
   ncsftmp = ncsfcfg[2*i+iloc];
   for (ist=0; ist<nst; ist++) {
     ncsfsum = 0;
     for (kst=0; kst<nst; kst++) {
       lst = GRPMUL(ist,kst);
       if (lst == kst) {
         for (k=exsym[kst]; k>=sxsym[kst]; k--) {
           kl = (k*(k+1))/2;
           for (l=k-1; l>=sxsym[lst]; l--) {
             /* unpaired external electrons */
             kl--;
             lccsf_11[kl+itmp] = ncsfsum;
             ncsfsum += ncsftmp;
           }
         }
 
       } else if (lst > kst) {
         /* unpaired external electrons */
         for (k=exsym[kst]; k>=sxsym[kst]; k--) {
           kx = (k*(k-1))/2;
           for (l=exsym[lst]; l>=sxsym[lst]; l--) {
             kl = kx + l + itmp;
             lccsf_11[kl] = ncsfsum;
             ncsfsum += ncsftmp;
           }
         }
       }
     }
   }
   itmp += nnorbx; 
 }

 /*if (prntflag>1) {
   for (i=(spin/2); i<=(mxkmax/2); i++) {
     (void)fprintf(outfile, "two electron with no pair completion table for %d open shells:\n", 2*i+spin%2);
     for (k=1,kl=1; k<=norbx; k++) {
       for (l=1; l<=k; l++) {
         (void)fprintf(outfile, "k=%d, l=%d, lccsf_11[%d] = %d\n", k, l, kl, lccsf_11[kl+i*nnorbx]);
         kl++;
       } 
     }
   } 
   (void)fflush(outfile);
 }*/ 

 /* for pairs of singly occupied external orbitals */
 for (i=0; i<nst; i++) nxprsym[i] = 0;
 for (i=0; i<nst; i++) {
   for (j=0; j<i; j++) {
     /* i > j */
     nxprsym[GRPMUL(i,j)] += nxsym[i] * nxsym[j]; 
   }
   /* i = j */
   nxprsym[0] += (nxsym[i]*(nxsym[i]-1))/2;  
 }

 /* for triples of singly occupied external orbitals */
 for (i=0; i<nst; i++) nxtrsym[i] = 0;
 for (i=0; i<nst; i++) {
   for (j=0; j<i; j++) {
     ij = GRPMUL(i,j);
     for (k=0; k<j; k++) { 
       /* i > j > k */
       nxtrsym[GRPMUL(ij,k)] += nxsym[i] * nxsym[j] * nxsym[k];
     }
     /* i > j = k */
     nxtrsym[i] += nxsym[i] * (nxsym[j]*(nxsym[j]-1))/2;  
     /* i = j > k */
     nxtrsym[j] += nxsym[j] * (nxsym[i]*(nxsym[i]-1))/2;  
   }
   /* i = j = k */
   nxtrsym[i] += (nxsym[i]*(nxsym[i]-1)*(nxsym[i]-2))/6;  
 }

 /* for quadruples of singly occupied external orbitals */
 for (i=0; i<nst; i++) nxqrsym[i] = 0;
 for (i=0; i<nst; i++) {
   for (j=0; j<i; j++) {
     ij = GRPMUL(i,j);
     for (k=0; k<j; k++) {
       ik = GRPMUL(i,k);
       jk = GRPMUL(j,k);
       for (l=0; l<k; l++) { 
         /* i > j > k > l */
         nxqrsym[GRPMUL(ij,GRPMUL(k,l))] += nxsym[i] * nxsym[j] * nxsym[k] * nxsym[l];
       }
       /* i > j > k = l */
       nxqrsym[ij] += nxsym[i] * nxsym[j] * (nxsym[k]*(nxsym[k]-1))/2; 
       /* i > j = k > l */
       nxqrsym[ik] += nxsym[i] * nxsym[k] * (nxsym[j]*(nxsym[j]-1))/2; 
       /* i = j > k > l */
       nxqrsym[jk] += nxsym[j] * nxsym[k] * (nxsym[i]*(nxsym[i]-1))/2; 
     }
     /* i > j = k = l */
     nxqrsym[ij] += nxsym[i] * (nxsym[j]*(nxsym[j]-1)*(nxsym[j]-2))/6;   
     /* i = j > k = l */
     nxqrsym[0]  += ((nxsym[i]*(nxsym[i]-1))/2) * ((nxsym[j]*(nxsym[j]-1))/2);
     /* i = j = k > l */
     nxqrsym[ij] += nxsym[j] * (nxsym[i]*(nxsym[i]-1)*(nxsym[i]-2))/6;  
   }
   /* i = j = k = l */
   nxqrsym[0] += (nxsym[i]*(nxsym[i]-1)*(nxsym[i]-2)*(nxsym[i]-3))/24;   
 }

 /* compute number of coupling coefficients for hamiltonian diagonal */
 ncc = 1;              /* include dummy location */
 itmp = 1;             /* ditto */
 for (i=0; i<=mxkmax; i++) {
   ftpt[i] = 0;
   snglptr[i] = 0;
   if (ncsfcfg[i]) {
     ftpt[i] = ncc;
     if (invflag) {
        ncc += ncsfcfg[i] * (i*(i-1))/2;
     } else {
        ncc += (i*(i-1))/2;
     }
     snglptr[i] = itmp;
     itmp += ncsfcfg[i];
   }
 }
 lenbuf = ncc * sizeof(double);
 exchcc = (double *)malloc(lenbuf);
 meminit += lenbuf;
 if (!invflag) {
    for (i=0; i<ncc; i++) exchcc[i] = 0.0;
 }
 
 for (i=0; i<=mxkmax; i++) {
   if (spnhead[i] != -1) {
     (void)diagcc(i, &exchcc[ftpt[i]]);
   }
 }

 if (!invflag) {
   for (i=0, jg=1; i<=mxkmax; i++) {  /* include dummy location */
     if (ncsfcfg[i]) {
        ncc = (i*(i-1))/2;
        val = 1./ncsfcfg[i];
        for (j=0; j<ncc; j++, jg++) exchcc[jg] *= val;
     }
   }
 }

 /* global position of CSFs */
 lenbuf = ncfgi * sizeof(int);
 gbcsf = (int *)malloc(lenbuf);
 if (!gbcsf) {
    (void)fprintf(outfile, "\nAllocation Failure for gbcsf\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
 }
 meminit += lenbuf;

 /* array for coupling coefficients  */
 /* one-particle-one-hole            */
 /* The order for each block of coupling coefficients is
    elmtone[ncsftmp*ncsftmp2], acum[numa*ncsftmp*ncsftmp2],
    bcum[numb*ncsftmp*ncsftmp2], ccum[numc*ncsftmp*ncsftmp2].  */
 len11 = 3 * norbi0;
 lenbuf = len11 * ncsftmp * ncsftmp2 * sizeof(struct cc_elmt1);
 meminit += lenbuf;
 cc_one = (struct cc_elmt1 *)malloc(lenbuf);
 if (!cc_one) {
    (void)fprintf(outfile, "\nAllocation Failure for coupling coefficients (1) \n");
    MPI_Abort(MPI_COMM_WORLD, 1);
 }
 lenbuf =  len11 * sizeof(int);
 meminit += lenbuf;
 elmt11cnt = (int *)malloc(lenbuf);
 if (!elmt11cnt) {
    (void)fprintf(outfile, "\nAllocation Failure for elmt11cnt[] \n");
    MPI_Abort(MPI_COMM_WORLD, 1);
 }

 /* check memory */
 len22 = 3 * norbi0 * norbi0 * (nelec - 3); 
 lenbuf =  len22 * sizeof(int);
 meminit += lenbuf;
 meminit += ncsfq1 * nstates * sizeof (double); /* memory allocation for CI vectors */ 
 mem = meminit/1.0e9;
 if (totalmem < meminit) {
    (void)fprintf(outfile, "\n***** requires %.2lf gb allocatable memory   *****\n", mem);
    (void)fprintf(outfile, "\n***** increase allocatable memory (MEMORY=#.#GB) *****\n");
    (void)fflush(outfile);
    MPI_Abort(MPI_COMM_WORLD, 1);
 }
 
 /* two-particle-two-hole            */
 lenbuf = len22 * ncsftmp * ncsftmp2 * sizeof(struct cc_elmt2); 
 meminit += lenbuf;
 if (totalmem < meminit) {
    iwrite_elmt = 1;
    mem = totalmem/1.0e9;
    bunlen = (lenbuf + totalmem - meminit) / sizeof (struct cc_elmt2);
    if (bunlen < 20000000) {
       mem = (lenbuf + totalmem - meminit) / 1.0e9;
       mem = 0.5 - mem;
       (void)fprintf(outfile, "\n***** Warning: possible insufficient allocatable memory  ****\n"
             "\n***** Increase allocatable memory (MEMORY=#.#GB) by %.2lfGB *****\n", mem);
       (void)fflush(outfile);
    }
    ciftfile = mpisopen(CIFTNAME, "w+b", infofile, my_rank);
    if (ciftfile == NULL) {
      (void)fprintf(outfile, "\nError opening file: %s\n", CIFTNAME);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  
 } else {
    bunlen = len22 * ncsftmp * ncsftmp2;
    iwrite_elmt = 0;
    mem = meminit/1.0e9;
    if (prntflag && my_rank==0) {
       (void)fprintf(outfile, "\n***** approximate memory usage:  %.4lf gb    *****\n", mem);
       (void)fflush(outfile);
    }
 }
 lenbuf = bunlen * sizeof(struct cc_elmt2); 
 cc_two = (struct cc_elmt2 *)malloc(lenbuf);
 if (!cc_two) {
    (void)fprintf(outfile, "\nAllocation Failure for coupling coefficients (2) \n");
    MPI_Abort(MPI_COMM_WORLD, 1);
 }
 lenbuf =  len22 * sizeof(int);
 elmt22cnt = (int *)malloc(lenbuf);
 if (!elmt22cnt) {
    (void)fprintf(outfile, "\nAllocation Failure for elmt22cnt[] \n");
    MPI_Abort(MPI_COMM_WORLD, 1);
 }
 
 nlanes = 0;
 nelx = 0;
 mcrloc = 0;
 for (ae=0; ae<nmcrsum[0]; ae++) {
  mcrloc += 9*nvrtx[ae]*sizeof(int);
 }
 (void)fseek(mcfgfile, mcrloc, SEEK_SET);

 for (; ae<nmcrsum[2]; ae++) {
   mcfgln[ae] = nlanes; 
   /* read drt corresponding to macroconfiguration aeg */
   len = 9*nvrtx[ae]*sizeof(int);
   if (fread(jdn, len, 1, mcfgfile)!=1) {
     (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   /* derived pointers */
   jup = jdn + 3*nvrtx[ae];

   icfg = mkhdiagmcr(nelx, jup, &gbcsf[nlanes]);
   nlanes += icfg;               /* cumulative number of internal configs (lanes) */
 }

 nelx=1;
 for (; ae<nmcrsum[4]; ae++) {
   mcfgln[ae] = nlanes; 
   /* read drt corresponding to macroconfiguration aeg */
   len = 9*nvrtx[ae]*sizeof(int);
   if (fread(jdn, len, 1, mcfgfile)!=1) {
     (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   /* derived pointers */
   jup = jdn + 3*nvrtx[ae];

   icfg = mkhdiagmcr(nelx, jup, &gbcsf[nlanes]);
   nlanes += icfg;               /* cumulative number of internal configs (lanes) */
 }

 nelx=2;
 for (; ae<nmcrsum[5]; ae++) {
   mcfgln[ae] = nlanes; 
   /* read drt corresponding to macroconfiguration aeg */
   len = 9*nvrtx[ae]*sizeof(int);
   if (fread(jdn, len, 1, mcfgfile)!=1) {
     (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   /* derived pointers */
   jup = jdn + 3*nvrtx[ae];

   icfg = mkhdiagmcr(nelx, jup, &gbcsf[nlanes]);
   nlanes += icfg;               /* cumulative number of internal configs (lanes) */
 }

 nelx=3;
 for (; ae<nmcrsum[6]; ae++) {
   mcfgln[ae] = nlanes; 
   /* read drt corresponding to macroconfiguration aeg */
   len = 9*nvrtx[ae]*sizeof(int);
   if (fread(jdn, len, 1, mcfgfile)!=1) {
     (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
   /* derived pointers */
   jup = jdn + 3*nvrtx[ae];

   icfg = mkhdiagmcr(nelx, jup, &gbcsf[nlanes]);
   nlanes += icfg;               /* cumulative number of internal configs (lanes) */
 }

 if (prntflag>1 && my_rank==0) 
   (void)fprintf(outfile, "\ntotal number of internal configurations (lanes) = %d\n", nlanes);

}

int mkhdiagmcr(int nelx, int jup[], int gbcsf[])
{
 /* calculate diagonal hamiltonian matrix elements between two macroconfigurations */
 extern int norb, norbi, norbx, orbsymi[], irrep, nst, ncsfcfg[];
 extern FILE *outfile;

 int i, icfg, ivrtx, k, kk, kimax;
 int path[NINTMX+1], symcum[NINTMX+1], nstep[NINTMX+1];
 int successful_ascent, nopenb[NINTMX+3];
 int ncsf, ast, bst, extsym;
 static int ncsfsum=0;
#ifdef DEBUG
 extern int nstates, exsym[], sxsym[];
 extern double *civec;
 int ik, istate, icsf, icsfloc, a, b;
#endif

  
 /* generate bra walks in ascending order */
 ivrtx = 0;
 path[0] = ivrtx;
 symcum[0] = 0;
 nopenb[0] = 0;
 icfg = 0;
 i = 0;
 kimax = 2;
 
 for (;;) {
   /* try to ascend (i => i+1) */
   successful_ascent = 0;
   kk = ivrtx*3;
 
   for (k=kimax; k>=0; k--) {
      if (jup[kk+k]>-1) {
         nstep[i+1] = k;
         ivrtx = jup[kk+k];
         path[i+1] = ivrtx;
         if (k==1) {
            symcum[i+1] = GRPMUL(symcum[i],orbsymi[i+1]);
            nopenb[i+1] = nopenb[i] + 1;
         } else {
            symcum[i+1] = symcum[i];
            nopenb[i+1] = nopenb[i];
         }
         kimax = 2;
         i++;
         successful_ascent=1;
         break;
      }
   }
      
   if (successful_ascent) {
      if (i==norbi) {
#ifdef DEBUG
   (void)fprintf(outfile, "\niq1cfg=%d, jcsfbase=%d, ncsf=%d", 
          icfg, ncsfsum, ncsfcfg[nopenb[norbi]]);
    (void)fprintf(outfile, "\n(%d", nstep[1]);
    for (ik=2; ik<=norbi; ik++) {
      (void)fprintf(outfile, " %d", nstep[ik]);
    }
    (void)fprintf(outfile, ")\n");

#endif
         gbcsf[icfg++] = ncsfsum;
         extsym = GRPMUL(symcum[norbi],irrep);
         ncsf = ncsfcfg[nopenb[norbi]]; 
         
         /* presume the order of CSFs from MRCISD has been modified */
         switch (nelx) {
            case 3:
               for (ast=0; ast<nst; ast++) {
                  bst = GRPMUL(extsym,ast);
                  if (bst == ast) {
#ifdef DEBUG
   icsfloc = ncsfsum;
   for (a=exsym[ast]; a>sxsym[ast]; a--) {
      for (b=a-1; b>=sxsym[ast]; b--) {
         for (icsf=0; icsf<ncsf; icsf++) {
            for (istate=0; istate<nstates; istate++) {
               (void)fprintf(outfile, "\n {%d, %d} civec[%d]=%16.12lf", 
                     a, b, icsfloc+icsf, civec[(icsfloc+icsf)*nstates+istate]);
            }
         }
         icsfloc += ncsf; 
      }
   }
#endif
                     ncsfsum += (ncsf*nxsym[ast]*(nxsym[ast]-1))/2;
                  } else if (bst > ast) {
#ifdef DEBUG
   icsfloc = ncsfsum;
   for (a=exsym[ast]; a>=sxsym[ast]; a--) {
      for (b=exsym[bst]; b>=sxsym[bst]; b--) {
         for (icsf=0; icsf<ncsf; icsf++) {
            for (istate=0; istate<nstates; istate++) {
               (void)fprintf(outfile, "\n {%d, %d} civec[%d]=%16.12lf", 
                     a, b, icsfloc+icsf, civec[(icsfloc+icsf)*nstates+istate]);
            }
         }
         icsfloc += ncsf; 
      }
   }
#endif
                     ncsfsum += ncsf*nxsym[ast]*nxsym[bst];
                  }
               }
               break;

            case 2:
#ifdef DEBUG
   icsfloc = ncsfsum;
   if (!extsym) {
      for (a=norbx; a>0; a--) {
         for (icsf=0; icsf<ncsf; icsf++) {
            for (istate=0; istate<nstates; istate++) {
               (void)fprintf(outfile, "\n {%d} civec[%d]=%16.12lf", 
                     a, icsfloc+icsf, civec[(icsfloc+icsf)*nstates+istate]);
            }
         }
         icsfloc += ncsf; 
      }
   }
#endif
               if (!extsym) ncsfsum += ncsf*norbx;
               break;

            case 1:
#ifdef DEBUG
   icsfloc = ncsfsum;
   for (a=exsym[extsym]; a>=sxsym[extsym]; a--) {
      for (icsf=0; icsf<ncsf; icsf++) {
         for (istate=0; istate<nstates; istate++) {
            (void)fprintf(outfile, "\n {%d} civec[%d]=%16.12lf", 
                a, icsfloc+icsf, civec[(icsfloc+icsf)*nstates+istate]);
         }
      }
      icsfloc += ncsf; 
   }
#endif
               ncsfsum += ncsf*nxsym[extsym];
               break;
 
            default:
#ifdef DEBUG
   icsfloc = ncsfsum;
   if (!extsym) {
      for (icsf=0; icsf<ncsf; icsf++) {
         for (istate=0; istate<nstates; istate++) {
            (void)fprintf(outfile, "\n civec[%d]=%16.12lf", 
                  icsfloc+icsf, civec[(icsfloc+icsf)*nstates+istate]);
         }
      }
      icsfloc += ncsf; 
   }
#endif
               if (!extsym) ncsfsum += ncsf;
               break;
         }   /* end switch */ 

         /* re-initialize */
         i = norbi - 2;            /* by-pass a wasted search step */
         if (i>-1) {
            ivrtx = path[i];
            kimax = nstep[i+1] - 1;
         } else {
            break;                 /* normal exit for 1-orbital internal space (rare!) */
         }
      }

   } else {           /* not "succsessful ascent" */
      /* try to descend */
      if (!i) break;               /* normal exit */
      ivrtx = path[--i];
      kimax = nstep[i+1] - 1;
   }
 }

 return icfg;  
}



int diagcc(int nopen, double exchcc[])
/* calculate coupling coefficients for diagonal hamiltonian matrix elements */
/* also maps spin trajectories to those that are singlet coupled for first two electrons */ 
{
 int braspn[MXKMAX+1], bstep[MXKMAX+1];
 int bralev, tmp, icsf, i, j, incflag, nelmt, nopenltr;
 int ij, ir, jr, isngl;
 double val, elmt;
 extern int spin, arcsp[], jspdn[], jspup[], spnhead[], invflag;
 extern FILE *outfile;
 extern double headcum[][MXKMAX+1];
 extern double pbm1[], pb02[], pb31[];
 extern double sq02[], sq12[];

 nelmt = 0; 
 icsf = 0; 
 bralev = nopen;
 if (spin>bralev) return nelmt;
 nopenltr = ((nopen-1)*nopen)/2;
 
 
 braspn[bralev] = spin;
 incflag = 1;
 isngl = 0;

 for (;;) {

   if (!bralev) {
      /* internal trajectory completed */
      
      /* complete the internal loops */
      for (i=1; i<nopen; i++) {
         ir = nopen - i + 1;
         if (bstep[i]==1) val = 1./pb02[braspn[i]];
         else val = -pb02[braspn[i]];
         for (j=i+1; j<=nopen; j++) {
            jr = nopen - j + 1;
            elmt = -headcum[i+1][j]*val/2.;
            if (invflag) {
               ij  = icsf*nopenltr + ((ir-1)*(ir-2))/2 + jr - 1;
               exchcc[ij] = -elmt;         /* J=1 phase factor */
            } else {
               ij  = ((ir-1)*(ir-2))/2 + jr - 1;
               exchcc[ij] -= elmt;         /* J=1 phase factor */
            }
            if (fabs(elmt)>0.) nelmt++;
         }
      }

      icsf++;
            
      /* re-initialize */
      bralev++;
      if (bralev>nopen) break;   /* normal exit for no internal open shells */
      incflag = 0; 
         
   } else if (incflag) { 
      /* try to descend (bralev => bralev-1) with spin increase, if possible */
      tmp = braspn[bralev] + 1;
      if (tmp < bralev) {
         bstep[bralev] = -1;
         if (bralev > 1) {
            headcum[bralev][bralev] = pb31[braspn[bralev]];
            val = sq12[braspn[bralev]];
            for (i=bralev+1; i<=nopen; i++) {
               headcum[bralev][i] = headcum[bralev+1][i] * val;
            }
         }
         braspn[--bralev] = tmp;
         
      } else {
         bstep[bralev] = 1;
         if (bralev > 1) {
            headcum[bralev][bralev] = -pbm1[braspn[bralev]];
            val = sq02[braspn[bralev]];
            for (i=bralev+1; i<=nopen; i++) {
               headcum[bralev][i] = headcum[bralev+1][i] * val;
            }
         }
         braspn[--bralev] = tmp - 2;
      }

   } else if (bstep[bralev]==-1) {
      /* try to descend with spin decrease */
      tmp = braspn[bralev] - 1;
      if (tmp >= 0) {
         bstep[bralev] = 1;
         if (bralev > 1) {
            headcum[bralev][bralev] = -pbm1[braspn[bralev]];
            val = sq02[braspn[bralev]];
            for (i=bralev+1; i<=nopen; i++) {
               headcum[bralev][i] = headcum[bralev+1][i] * val;
            }
         }
         braspn[--bralev] = tmp;
         incflag = 1;
         
      } else {
         /* try to ascend */
         if (++bralev>nopen) break;   /* normal exit */
      }
      
   } else {
      /* try to ascend */
      if (++bralev>nopen) break;   /* normal exit */
   }   
 }

 return nelmt; 
}
 
void mkhw(void)
{

 extern int comm_sz, my_rank;
 extern int nstates, nmcrt, nmcr[], nmcrsum[], nvrtx[], jdn[];
 extern int norbi, ngrps, *mcfg, mcfgln[], prntflag, iter, invflag;
 extern double *energs, *energs2, efzc;
 extern FILE *mcfgfile, *outfile, *scrfile;

 int assignment, readyrank, readytag=0, asmttag=1, Idx_asmtarray;
 int imcr, jmcr, jcnt, jcnttmp, i;
 int imcrloc, imcrtyp, len, nelxb;
 int *jup, **kdn, **ptr_kup, **ptr_arc, **ptr_ketbase, *nelxk;
 static int *asmtarray;//stands for assignments array
 static int mcrloc[15];
 
 int mkq1mcr(int imcr, int jmcrloc, int jstart, int jend, int nsp, 
             int **kdn, int **ptr_kup, int **ptr_arc, 
             int **ptr_ketbase, int jcnt);
 void mcrdrv(int jup[], int **ptr_kup, int **ptr_arc, int jcnt, 
             int **ptr_ketbase, int nelxb, int *nelxk);


 /***************************************************************** 
  *                                                               *
  * The macroconfigurations are further divided into 30 blocks.   *
  * 0----(N,  0)E0     no excitation                              *
  * 1----(N-1, 1)E0    one internal excitation                    *
  * 2----(N-2, 2)E0    two internal excitations                   *
  * 3----(N-1, 0)E1    one external excitation                    *
  * 4----(N-2, 1)E1    one internal, one external excitations     *
  * 5----(N-2, 0)E2    two external excitations                   *
  * 6----(N-2, 0)E11   two external excitations                   *
  * 7----(N-3, 3)E0    three internal excitations                 *
  * 8----(N-4, 4)E0    four internal excitations                  *
  * 9----(N-3, 2)E1    two internal, one external excitations     *
  *10----(N-4, 3)E1    three internal, one external excitations   *
  *11----(N-3, 1)E2    one internal, two external excitations     *
  *12----(N-3, 1)E11   one internal, two external excitations     *
  *13----(N-4, 2)E2    two internal, two external excitations     *
  *14----(N-4, 2)E11   two internal, two external excitations     *
  *15----(N-3, 0)E21   three external excitation                  *
  *16----(N-3, 0)E111  three external excitation                  *
  *17----(N-4, 1)E21   one internal, three external excitations   *
  *18----(N-4, 1)E111  one internal, three external excitations   *
  *19----(N-4, 0)E22   four external excitations                  * 
  *20----(N-4, 0)E211  four external excitations                  * 
  *21----(N-4, 0)E1111 four external excitations                  * 
  *                                                             *
  ***************************************************************/ 
 
 /* contributions solely from frozen orbitals will be cancelled */
 for (i=0; i<nstates; i++) energs[i] -= efzc;
 
 if (invflag) {
   for (i=0; i<nstates; i++) energs2[i] -= efzc; 
 }

 if (iter) {
    if (my_rank != 0) { //only slaves work in later iterations, master assignes jobs in the 1st iter and them rest

       Idx_asmtarray = 0;
       rewind(scrfile);
       (void)fseek(mcfgfile, mcrloc[6], SEEK_SET);
       for (imcrtyp=7; imcrtyp<22; imcrtyp++) {
          switch (imcrtyp) {
             case 21:
                nelxb = 8;
                break;
             case 20:
                nelxb = 7;
                break;
             case 19:
                nelxb = 6;
                break;
             case 18: case 16:
                nelxb = 5;
                break;
             case 17: case 15:
                nelxb = 4;
                break;
             case 14: case 12:
                nelxb = 3;
                break;
             case 13: case 11:
                nelxb = 2;
                break;
             case 10: case 9:
                nelxb = 1;
                break;
             case 8: case 7:
                nelxb = 0;
                break;
          }

          for (imcr=nmcrsum[imcrtyp-1]; imcr<nmcrsum[imcrtyp]; imcr++) {
             /* read drt for imcr from storage */
             len = 9*nvrtx[imcr];
             if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
                (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
                MPI_Abort(MPI_COMM_WORLD, 1);
             }
             jup = jdn + 3*nvrtx[imcr];
    
             if (asmtarray[Idx_asmtarray] == imcr) {
                Idx_asmtarray++;
                mcrdrv(jup, ptr_kup, ptr_arc, 0, ptr_ketbase, nelxb, nelxk);
             } 
          }
       } 
    }
    return;
 } 

 if (my_rank == 0) { //master
    for (imcr=nmcrsum[6]; imcr<nmcrt; imcr++) {
       assignment = imcr;
       MPI_Recv(&readyrank, 1, MPI_INT, MPI_ANY_SOURCE, readytag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       MPI_Send(&assignment, 1, MPI_INT, readyrank, asmttag, MPI_COMM_WORLD);
       if (prntflag>=2) {
          fprintf(outfile, "proc %d is working on macroconfiguration %d in Q space\n", readyrank, assignment);
          fflush(outfile);
       }
    }
    for (i=1; i< comm_sz; i++) {
       assignment = -1;
       MPI_Recv(&readyrank, 1, MPI_INT, MPI_ANY_SOURCE, readytag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       MPI_Send(&assignment, 1, MPI_INT, readyrank, asmttag, MPI_COMM_WORLD);
    }

 }

 else { //slaves
 mcrloc[0] = 0;
 for (imcr=0; imcr<nmcrsum[0]; imcr++) {
  mcrloc[0] += 9*nvrtx[imcr]*sizeof(int);
 }
 for (i=1; i<7; i++) {
  mcrloc[i] = mcrloc[i-1];
  for (; imcr<nmcrsum[i]; imcr++) {
    mcrloc[i] += 9*nvrtx[imcr]*sizeof(int);
  }
 }

 asmtarray   = (int *)malloc(nmcrt * sizeof (int));
 ptr_kup     = (int **)malloc(nmcrsum[6] * sizeof (int *));
 ptr_arc     = (int **)malloc(nmcrsum[6] * sizeof (int *));
 ptr_ketbase = (int **)malloc(nmcrsum[6] * sizeof (int *));
 nelxk       = (int *)malloc(nmcrsum[6] * sizeof (int));
 kdn         = (int **)malloc(nmcrsum[6] * sizeof (int *));
 kdn[0]      = (int *) malloc(nmcrsum[6] * 9 * MAXVRTX * sizeof (int));
 for (imcr=1; imcr<nmcrsum[6]; imcr++)   
   kdn[imcr] = kdn[imcr-1] + 9 * MAXVRTX;

 /* three internal excitations----(N-3, 3)E0 */
 imcrloc = mcrloc[6];
 readyrank = my_rank;
 MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
 MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 Idx_asmtarray = 0;

 for (imcr=nmcrsum[6]; imcr<nmcrsum[7]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[0], nmcrsum[0], nmcrsum[2], 3, kdn, ptr_kup, 
                 ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[3], nmcrsum[3], nmcrsum[4], 3, kdn, ptr_kup, 
                 ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 0, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
 } /* end of three internal excitations----(N-3, 3)E0 */

 /* four internal excitations----(N-4, 4)E0 */
 for (; imcr<nmcrsum[8]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[1], nmcrsum[1], nmcrsum[2], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 0, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of four internal excitations----(N-4, 4)E0 */

 /* two internal, one external excitations----(N-3, 2)E1 */
 for (; imcr<nmcrsum[9]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[0], nmcrsum[0], nmcrsum[2], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[2], nmcrsum[2], nmcrsum[4], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 1, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of two internal, one external excitations----(N-3, 2)E1 */

 /* three internal, one external excitations----(N-4, 3)E1 */
 for (; imcr<nmcrsum[10]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[1], nmcrsum[1], nmcrsum[2], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[3], nmcrsum[3], nmcrsum[4], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 1, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of three internal, one external excitations----(N-4, 3)E1 */

 /* one internal, two external excitations----(N-3, 1)E20 */
 for (; imcr<nmcrsum[11]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[0], nmcrsum[0], nmcrsum[2], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[2], nmcrsum[2], nmcrsum[4], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 2, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of one internal, two external excitations----(N-3, 1)E20 */

 /* one internal, two external excitations----(N-3, 1)E11 */
 for (; imcr<nmcrsum[12]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[0], nmcrsum[0], nmcrsum[2], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[2], nmcrsum[2], nmcrsum[4], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 3, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of one internal, two external excitations----(N-3, 1)E11 */

 /* two internal, two external excitations----(N-4, 2)E20 */
 for (; imcr<nmcrsum[13]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[1], nmcrsum[1], nmcrsum[2], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[3], nmcrsum[3], nmcrsum[4], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 2, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
 } /* end of two internal, two external excitations----(N-4, 2)E2 */

 /* two internal, two external excitations----(N-4, 2)E11 */
 for (; imcr<nmcrsum[14]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[1], nmcrsum[1], nmcrsum[2], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 0;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[3], nmcrsum[3], nmcrsum[4], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 3, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 3, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
 } /* end of two internal, two external excitations----(N-4, 2)E11 */

 /* three external excitations----(N-3, 0)E21 */
 for (; imcr<nmcrsum[15]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[2], nmcrsum[2], nmcrsum[4], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 4, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of three external excitations----(N-3, 0)E21 */

 /* three external excitations----(N-3, 0)E111 */
 for (; imcr<nmcrsum[16]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[2], nmcrsum[2], nmcrsum[4], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 5, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of three external excitations----(N-3, 0)E111*/

 /* one internal, three external excitations----(N-4, 1)E21 */
 for (; imcr<nmcrsum[17]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[3], nmcrsum[3], nmcrsum[4], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 4, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of one internal, three external excitations----(N-4, 1)E21 */

 /* one internal, three external excitations----(N-4, 1)E111 */
 for (; imcr<nmcrsum[18]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[3], nmcrsum[3], nmcrsum[4], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 1;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 2, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 5, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of one internal, three external excitations----(N-4, 1)E111*/

 /* four external excitations----(N-4, 0)E22 */
 for (; imcr<nmcrsum[19]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 6, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of four external excitations----(N-4, 0)E22 */

 /* four external excitations----(N-4, 0)E211 */
 for (; imcr<nmcrsum[20]; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[4], nmcrsum[4], nmcrsum[5], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 2;

  jcnttmp = jcnt;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=jcnttmp; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 7, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
 } /* end of four external excitations----(N-4, 0)E211 */

 /* four external excitations----(N-4, 0)E1111 */
 for (; imcr<nmcrt; imcr++) {
  
  /* read drt for imcr from storage */
  (void)fseek(mcfgfile, imcrloc, SEEK_SET);
  len = 9*nvrtx[imcr];
  if (fread(jdn, len*sizeof(int), 1, mcfgfile)!=1) {
    (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  jup = jdn + 3*nvrtx[imcr];
  imcrloc += len * sizeof(int);
  
  /* The interacting ket macroconfigurations will be put together. */
  jcnt = 0;
  jcnt = mkq1mcr(imcr, mcrloc[5], nmcrsum[5], nmcrsum[6], 1, kdn, ptr_kup, 
          ptr_arc, ptr_ketbase, jcnt);
  for (jmcr=0; jmcr<jcnt; jmcr++) nelxk[jmcr] = 3;

  if (assignment == imcr) {
     asmtarray[Idx_asmtarray] = assignment;
     Idx_asmtarray++;
     mcrdrv(jup, ptr_kup, ptr_arc, jcnt, ptr_ketbase, 8, nelxk);
     MPI_Send(&readyrank, 1, MPI_INT, 0, readytag, MPI_COMM_WORLD);
     MPI_Recv(&assignment, 1, MPI_INT, 0, asmttag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

 } /* end of four external excitations----(N-4, 0)E1111 */

 } //end of else(slave operations)

fflush(outfile);
}
 

void mcrdrv(int jup[], int **ptr_kup, int **ptr_arc, int jcnt, int **ptr_ketbase,
            int nelxb, int *nelxk)
{
 /* calculate contributions to sigma vectors between Q2-macroconfigurations
    and all interacting Q1-macroconfigurations.                             */ 
  
 extern FILE *outfile, *ciftfile;
 extern int norb, norbi0, norbi, norbx, orbsymi[], irrep, nnorbx;
 extern int sxsym[], exsym[], ncsfcfg[], orbsym[];
 extern int elmt22tot, bunlen, iwrite_elmt;
 extern int intgrlmo[], xorbq2[], *elmt11cnt, *elmt22cnt, len11, len22;
 extern int jcfg, jcfg11, ncsfq1tmp, ncsfq2tmp, ncsfxstate, nstates;
 extern double *twoint, *oneint, intk[], exchint[];
 extern double hdiag, *hdiagtmp, *hdiagtmp2;
 extern struct cc_elmt2 *cc_two, *ptr_cc_two;
 extern struct cc_elmt1 *cc_one, *ptr_cc_one;

 int i, j, ij, ivrtx, k, kk, jvrtx, exsp, exsh, kmax, kimax;
 int path[NINTMX+1], symcum[NINTMX+1], nstep[NINTMX+1];
 int kpath[NINTMX+1], ksymcum[NINTMX+1], kstep[NINTMX+1], cfgcum[NINTMX+1];
 int xsp[NINTMX+1], xsh[NINTMX+1];
 int *kupmcr, *arcmcr, *ketbase, ncfgq2;
 int kmcr, ncfgtmp, extsym, extsymk;
 int successful_ascent, nopenb[NINTMX+1], dbra[NINTMX+1], dbmap[MXKMAX+1];
 int nopenk[NINTMX+1], dket[NINTMX+1], dkmap[MXKMAX+1], iorb[4], jorb[4];
 int imo, jmo, mii, mij, mjj, mxj, miiii, miijj, mijji;
 int a, b, c, d, ast, bst, cst, dst, xst, permut; 
 double xint;

 void calchdiag_bw(int nstep[], int dbra[], int dbmap[], int nelxb, int nopenb, 
             double intk, int nxorb);
 void calc_hwpp(int nstep[], int dbra[], int dbmap[], int nopenb, int **ptr_kup, 
             int **ptr_arc, int **ptr_ketbase, int jcnt, int nelxb, int *nelxk, 
             int extsym); 
 int chksym(int nelxb, int nelxk, int extsym, int extsymk, int *permut); 
 void preptwoptwoh(int iorb[], int jorb[], int nopenk, int dket[], int dkmap[],
                   int nopenb, int dbra[], int dbmap[]); 
 void oneponeh(int iint, int jint, int nopenb, int dbra[], int dbmap[],
  int nopenk, int dket[], int dkmap[]); 

 xsp[0] = 0;
 xsh[0] = 0;

 /* generate bra walks in ascending order */
 ivrtx = 0;
 path[0] = ivrtx;
 symcum[0] = 0;
 nopenb[0] = 0;
 i = 0;
 kimax = 2;
 dbra[0] = 0;
 intk[0] = 0.0;
 
 for (;;) {
   /* try to ascend (i => i+1) */
   successful_ascent = 0;
   kk = ivrtx*3;
      
   for (k=kimax; k>=0; k--) {
      if (jup[kk+k]>-1) {
         nstep[i+1] = k;
         ivrtx = jup[kk+k];
         path[i+1] = ivrtx;
         intk[i+1] = intk[i];
         imo = norb-i;
         if (k==1) {
            dbmap[nopenb[i]] = norbi - i;
            symcum[i+1] = GRPMUL(symcum[i],orbsymi[i+1]);
            nopenb[i+1] = nopenb[i] + 1;
            dbra[norbi-i] = 4;
            if (i<norbi0) {
               mii = (imo*(imo+1))/2;
               xint = oneint[mii];
               ij = (nopenb[i]*(nopenb[i]-1))/2;
               for (j=0; j<i; j++) {
                  jmo = norb - j;
                  mxj = (jmo*(jmo-1))/2;
                  mij = mxj + imo;
                  mijji = (mij*(mij+1))/2;
                  mjj = mxj + jmo;
                  miijj = (mjj*(mjj-1))/2 + mii;
                  xint += nstep[j+1]*(twoint[miijj]-0.5*twoint[mijji]);
                  if (nstep[j+1]==1) {                   /* open-shell exchange integrals */
                     exchint[ij++] = twoint[mijji];
                  }
               }
               intk[i+1] += xint;
            }
         } else {
            symcum[i+1] = symcum[i];
            nopenb[i+1] = nopenb[i];
            if (nstep[i+1]) {
               dbra[norbi-i] = 3;
               if (i<norbi0) {
                  mii = (imo*(imo+1))/2;
                  xint = 2 * oneint[mii];
                  miiii = (mii*(mii+1))/2;
                  xint += twoint[miiii];
                  for (j=0; j<i; j++) {
                     jmo = norb - j;
                     mxj = (jmo*(jmo-1))/2;
                     mij = mxj + imo;
                     mijji = (mij*(mij+1))/2;
                     mjj = mxj + jmo;
                     miijj = (mjj*(mjj-1))/2 + mii;
                     xint += nstep[j+1]*(2.*twoint[miijj]-twoint[mijji]);
                  }
                  intk[i+1] += xint;
               }
            } else {
               dbra[norbi-i] = 0;
            }
         }
         kimax = 2;
         i++;
         successful_ascent=1;
         break;
      }
   }
      
   if (successful_ascent) {
      if (i==norbi) { 
         /* check symmetry !!!! */
         extsym = GRPMUL(symcum[norbi], irrep);
         ncsfq2tmp = ncsfcfg[nopenb[norbi]];
         if (!ncsfq2tmp) goto endbra;
         ncsfxstate = ncsfq2tmp * nstates;
         ncfgq2 = chksym(nelxb, -1, extsym, 0, &permut); 
         if (!ncfgq2) goto endbra; 

         /* evaluation of coupling coefficients for two-electron integrals */
         /* initialize the coupling coefficient arrays                     */
         ptr_cc_one = cc_one;
         ptr_cc_two = cc_two;
         for (j=0; j<len22; j++) {
            elmt22cnt[j] = 0;
         }
         for (j=0; j<len11; j++) {
            elmt11cnt[j] = 0;
         }
         jcfg = 0;
         jcfg11 = 0;
         if (iwrite_elmt) {
            elmt22tot = 0;
            rewind(ciftfile);
         }

         /* loop over all viable model macroconfigurations. */
         /* kup => ptr_kup[kmcr] */
         /* arc => ptr_arc[kmcr] */
         /* ketbase => ptr_ketbase[kmcr] */
         for (kmcr=0; kmcr<jcnt; kmcr++) {
            kupmcr = ptr_kup[kmcr];
            arcmcr = ptr_arc[kmcr];
            ketbase = ptr_ketbase[kmcr];

            switch (nelxb) {
               case 6: case 7: case 8:  
                  exsp = 0;
                  exsh = 2;
                  break;
               case 5: 
                  switch (nelxk[kmcr]) {
                     case 1: 
                        exsp = 0;
                        exsh = 2;
                        break;
                     case 2:  
                        exsp = 0;
                        exsh = 1;
                        break;
                     case 3: 
                        exsp = 1;
                        exsh = 2;
                        break;
                  }
                  break;
               case 4: 
                  switch (nelxk[kmcr]) {
                     case 1: 
                        exsp = 0;
                        exsh = 2;
                        break;
                     case 2: case 3: 
                        exsp = 1;
                        exsh = 2;
                        break;
                  }
                  break;
               case 3: 
                  switch (nelxk[kmcr]) {
                     case 0: 
                        exsp = 0;
                        exsh = 2;
                        break;
                     case 1: 
                        exsp = 1;
                        exsh = 2;
                        break;
                     case 2:  
                        exsp = 1;
                        exsh = 1;
                        break;
                     case 3: 
                        exsp = 2;
                        exsh = 2;
                        break;
                  }
                  break;
               case 2: 
                  switch (nelxk[kmcr]) {
                     case 0: 
                        exsp = 0;
                        exsh = 2;
                        break;
                     case 1: 
                        exsp = 1;
                        exsh = 2;
                        break;
                     case 2:  
                        exsp = 2;
                        exsh = 2;
                        break;
                     case 3: 
                        exsp = 1;
                        exsh = 1;
                        break;
                  }
                  break;
               case 1: 
                  switch (nelxk[kmcr]) {
                     case 0: 
                        exsp = 1;
                        exsh = 2;
                        break;
                     case 1: 
                        exsp = 2;
                        exsh = 2;
                        break;
                     case 2: case 3: 
                        exsp = 2;
                        exsh = 1;
                        break;
                  }
                  break;
               case 0:
                  switch (nelxk[kmcr]) {
                     case 0: 
                        exsp = exsh = 2;
                        break;
                     case 1: 
                        exsp = 2;
                        exsh = 1;
                        break;
                  }
                  break;
            }
            /* >>>>>> generate *interacting* ket walks in ascending order <<<<<< */
            jvrtx = 0;
            kpath[0] = jvrtx;
            ksymcum[0] = 0;
            nopenk[0] = 0;
            dket[0] = 0;
            cfgcum[0] = 0;
            kmax = 2;
            j = 0;
           
            for (;;) {
               /* try to ascend (j => j+1) */
               successful_ascent = 0;
               kk = jvrtx*3;

               for (k=kmax; k>=0; k--) {
                  if (kupmcr[kk+k]>-1) {
                     kstep[j+1] = k;
                     /* bra and ket switched here */
                     switch (kstep[j+1] - nstep[j+1]) {
                        case 2:
                           iorb[xsp[j]] = norbi-j;
                           iorb[xsp[j]+1] = norbi-j;
                           xsp[j+1] = xsp[j] + 2;
                           xsh[j+1] = xsh[j];
                           break;
                        case 1:
                           iorb[xsp[j]] = norbi-j;
                           xsp[j+1] = xsp[j] + 1;
                           xsh[j+1] = xsh[j];
                           break;
                        case -1:
                           jorb[xsh[j]] = norbi-j;
                           xsp[j+1] = xsp[j];
                           xsh[j+1] = xsh[j] + 1;
                           break;
                        case -2:
                           jorb[xsh[j]] = norbi-j;
                           jorb[xsh[j]+1] = norbi-j;   
                           xsp[j+1] = xsp[j];
                           xsh[j+1] = xsh[j] + 2;
                           break;
                        default:
                           xsp[j+1] = xsp[j];
                           xsh[j+1] = xsh[j];
                     }
                     /* bra and ket switched here for xsp[] and xsh[] */
                     if (j<norbi0) {
                        if ((xsp[j+1]>exsh)||(xsh[j+1]>exsp)) continue;
                     }
                     cfgcum[j+1] = cfgcum[j] + arcmcr[kk+k];
                     jvrtx = kupmcr[kk+k];
                     kpath[j+1] = jvrtx;
                     if (k==1) {
                        dkmap[nopenk[j]] = norbi - j;
                        ksymcum[j+1] = GRPMUL(ksymcum[j],orbsymi[j+1]);
                        nopenk[j+1] = nopenk[j] + 1;
                        dket[norbi-j] = 4;           /* flag */
                     } else {
                        ksymcum[j+1] = ksymcum[j];
                        nopenk[j+1] = nopenk[j];
                        if (kstep[j+1]) dket[norbi-j] = 3;
                        else dket[norbi-j] = 0;
                     }
                     kmax = 2;
                     j++;
                     successful_ascent=1;
                     break;
                  }
               }

               if (successful_ascent) {
                  if (j==norbi) {
                     extsymk = GRPMUL(ksymcum[norbi], irrep);
                     ncsfq1tmp = ncsfcfg[nopenk[norbi]];
                     if (!ncsfq1tmp) goto endket;
                     ncfgtmp = chksym(nelxb, nelxk[kmcr], extsym, extsymk, &permut); 
                     if (!ncfgtmp) goto endket; 

                     switch (permut) {
                        case 0:
                           if (xsp[norbi]==2) {
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           } else {
                              if (iorb[0] > jorb[0]) {
                                 switch_bra_ket = 0;
                                 oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                    nopenb[norbi], dbra, dbmap);
                              } else {
                                 switch_bra_ket = 1;
                                 oneponeh(jorb[0], iorb[0], nopenb[norbi], dbra, dbmap,
                                    nopenk[norbi], dket, dkmap);
                              }     
                              jcfg11++;
                           }
                           break;

                        case 1:
                           if (xsp[norbi]==2) {
                              if (extsym==extsymk) {
                                 preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                    nopenb[norbi], dbra, dbmap); 
                              }
                           } else {
                              if (extsym==extsymk) {
                                 if (iorb[0] > jorb[0]) {
                                    switch_bra_ket = 0;
                                    oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                       nopenb[norbi], dbra, dbmap);
                                 } else {
                                    switch_bra_ket = 1;
                                    oneponeh(jorb[0], iorb[0], nopenb[norbi], dbra, dbmap,
                                       nopenk[norbi], dket, dkmap);
                                 }     
                                 jcfg11++;
                              }
                              iorb[1] = 3;
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           }
                           break;

                        case 2:
                           if (xsp[norbi]==2) {
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           } else {
                              switch_bra_ket = 0;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              iorb[1] = 3;
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           }
                           break;

                        case 3:
                           if (xsp[norbi]==2) {
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           } else {
                              switch_bra_ket = 0;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              jorb[0] = 4;   
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              iorb[1] = 2;
                              jorb[0] = 4;
                              jorb[1] = 3;
                              dket[4] = 0;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           }
                           break;

                        case 4:
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           if (extsym==extsymk) {
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           }
                           break;

                        case 5:
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[0] = 4;
                           dket[4] = 0;
                           dket[3] = 4;
                           dkmap[nopenk[norbi]-1] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[1] = 3;
                           dket[3] = 0;
                           dket[2] = 4;
                           dkmap[nopenk[norbi]-1] = 2;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           break;

                        case 6:
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           iorb[1] = 3;
                           jorb[1] = 4;
                           dket[4] = 0;
                           dket[3] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           break;

                        case 7:
                           if (xsp[norbi]==2) {
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           } else {
                              switch_bra_ket = 0;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              iorb[1] = 3;
                              jorb[0] = 4;
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           }
                           break;

                        case 8: 
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           iorb[1] = 3;
                           jorb[0] = 4;
                           dket[4] = 0;
                           dket[3] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           iorb[1] = 2;
                           jorb[1] = 3;
                           dket[3] = 0;
                           dket[2] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           break;

                        case 9:   /* change and revert */
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           iorb[1] = 4;
                           dbra[4] = 0;
                           dbra[3] = 4;
                           dbmap[nopenb[norbi]-1] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           dbra[4] = 4;
                           dbra[3] = 0;
                           dbmap[nopenb[norbi]-1] = 4;
                           break;

                        case 10:   /* change and revert */
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           iorb[1] = 4;
                           jorb[1] = 3;
                           dbra[4] = 0;
                           dbra[3] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           dbra[4] = 3;
                           dbra[3] = 0;
                           break;

                        case 11:  /* change and revert */
                           if (xsp[norbi]==2) {
                              if (extsym==extsymk) {
                                 preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                    nopenb[norbi], dbra, dbmap); 
                              }
                           } else {
                              if (extsym==extsymk) {
                                 if (iorb[0] > jorb[0]) {
                                    switch_bra_ket = 0;
                                    oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                       nopenb[norbi], dbra, dbmap);
                                 } else {
                                    switch_bra_ket = 1;
                                    oneponeh(jorb[0], iorb[0], nopenb[norbi], dbra, dbmap,
                                       nopenk[norbi], dket, dkmap);
                                 }     
                                 jcfg11++;
                              }     
                              iorb[1] = 2;
                              jorb[1] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-2] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              iorb[1] = 4;
                              jorb[1] = 2;
                              dket[4] = 4;
                              dket[2] = 0;
                              dkmap[nopenk[norbi]-2] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              dbra[4] = 0;
                              dbra[2] = 4;
                              dbmap[nopenb[norbi]-2] = 3;
                              dbmap[nopenb[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              dbra[4] = 4;
                              dbra[3] = 4;
                              dbra[2] = 0;
                              dbmap[nopenb[norbi]-2] = 4;
                              dbmap[nopenb[norbi]-1] = 3;
                           }
                           break;

                        case 12:   /* change and revert */
                           if (xsp[norbi]==2) {
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 3;
                              dbra[4] = 4;
                              dbra[3] = 3;
                              dbmap[nopenb[norbi]-1] = 4;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              dbra[4] = 3;
                              dbra[3] = 4;
                              dbmap[nopenb[norbi]-1] = 3;
                           } else {
                              switch_bra_ket = 0;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              jorb[0] = 3;
                              dbra[4] = 4;
                              dbra[3] = 3;
                              dbmap[nopenb[norbi]-1] = 4; 
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              jorb[0] = 4;
                              dbra[4] = 3;
                              dbra[3] = 4;
                              dbmap[nopenb[norbi]-1] = 3;
                              iorb[1] = 2;
                              jorb[1] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              iorb[1] = 4;
                              dket[3] = 4;
                              dket[2] = 0;
                              dkmap[nopenk[norbi]-1] = 3;
                              jorb[0] = 3;
                              jorb[1] = 2;
                              dbra[4] = 0;
                              dbra[3] = 3;
                              dbra[2] = 4;
                              dbmap[nopenb[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              iorb[1] = 2;
                              jorb[0] = 4;
                              jorb[1] = 4;
                              dbra[4] = 3;
                              dbra[3] = 4;
                              dbra[2] = 0;
                              dbmap[nopenb[norbi]-1] = 3;
                              dket[4] = 0;
                              dket[3] = 4;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-2] = 3;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              iorb[1] = 3;
                              dbra[3] = 0;
                              dbra[2] = 4;
                              dbmap[nopenb[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              dbra[3] = 4;
                              dbra[2] = 0;
                              dbmap[nopenb[norbi]-1] = 3;
                           }
                           break;

                        case 13:
                           if (xsp[norbi]==2) {
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-2] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              dket[4] = 4;
                              dket[2] = 0;
                              dkmap[nopenk[norbi]-2] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                           } else {
                              switch_bra_ket = 0;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              jorb[0] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              jorb[0] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-2] = 3;
                              dkmap[nopenk[norbi]-1] = 2;
                              oneponeh(iorb[0], jorb[0], nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap);
                              jcfg11++;
                              dket[4] = 4;
                              dket[2] = 0;
                              dkmap[nopenk[norbi]-2] = 4;
                              iorb[1] = 1;
                              jorb[0] = 3;
                              jorb[1] = 2;
                              dket[3] = 0;
                              dket[1] = 4;
                              dkmap[nopenk[norbi]-1] = 1;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[0] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dkmap[nopenk[norbi]-2] = 3;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-2] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              iorb[1] = 4;
                              jorb[0] = 2;
                              jorb[1] = 1;
                              dket[4] = 4;
                              dket[3] = 4;
                              dket[2] = 0;
                              dket[1] = 0;
                              dkmap[nopenk[norbi]-2] = 4;
                              dkmap[nopenk[norbi]-1] = 3;
                              dbra[4] = 0;
                              dbra[1] = 4;
                              dbmap[nopenb[norbi]-3] = 3;
                              dbmap[nopenb[norbi]-2] = 2;
                              dbmap[nopenb[norbi]-1] = 1;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[0] = 3;
                              dket[3] = 0;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              jorb[1] = 2;
                              dket[2] = 0;
                              dket[1] = 4;
                              dkmap[nopenk[norbi]-1] = 1;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                              dbra[4] = 4;
                              dbra[3] = 4;
                              dbra[2] = 4;
                              dbra[1] = 0;
                              dbmap[nopenb[norbi]-3] = 4;
                              dbmap[nopenb[norbi]-2] = 3;
                              dbmap[nopenb[norbi]-1] = 2;
                           }
                           break;

                        case 14:  /* change and revert */
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[1] = 3;
                           dket[3] = 0;
                           dket[2] = 4;
                           dkmap[nopenk[norbi]-1] = 2;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[0] = 3;
                           jorb[1] = 2;
                           dket[3] = 4;
                           dket[2] = 0;
                           dkmap[nopenk[norbi]-1] = 3;
                           dbra[4] = 4;
                           dbra[3] = 3;
                           dbmap[nopenb[norbi]-2] = 4;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[0] = 4;
                           jorb[1] = 2;
                           dket[4] = 0;
                           dket[2] = 4;
                           dkmap[nopenk[norbi]-2] = 3;
                           dkmap[nopenk[norbi]-1] = 2;
                           dbra[3] = 4;
                           dbra[2] = 3;
                           dbmap[nopenb[norbi]-1] = 3;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           dbra[4] = 3;
                           dbra[3] = 4;
                           dbra[2] = 4;
                           dbmap[nopenb[norbi]-2] = 3;
                           dbmap[nopenb[norbi]-1] = 2;
                           if (extsym==extsymk) {
                              jorb[0] = 4;
                              jorb[1] = 4;
                              dket[4] = 0;
                              dket[3] = 4;
                              dket[2] = 4;
                              dkmap[nopenk[norbi]-2] = 3;
                              dkmap[nopenk[norbi]-1] = 2;
                              preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                                 nopenb[norbi], dbra, dbmap); 
                           }
                           break;

                        case 15: 
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[0] = 3;
                           dket[3] = 0;
                           dket[2] = 4;
                           dkmap[nopenk[norbi]-1] = 2;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[1] = 2;
                           dket[2] = 0;
                           dket[1] = 4;
                           dkmap[nopenk[norbi]-1] = 1;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[0] = 4;
                           jorb[1] = 1;
                           dket[4] = 0;
                           dket[3] = 4;
                           dket[2] = 4;
                           dket[1] = 0;
                           dkmap[nopenk[norbi]-2] = 3;
                           dkmap[nopenk[norbi]-1] = 2;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[1] = 2;
                           dket[2] = 0;
                           dket[1] = 4;
                           dkmap[nopenk[norbi]-1] = 1;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           jorb[1] = 3;
                           dket[3] = 0;
                           dket[2] = 4;
                           dkmap[nopenk[norbi]-2] = 2;
                           preptwoptwoh(iorb, jorb, nopenk[norbi], dket, dkmap,
                              nopenb[norbi], dbra, dbmap); 
                           break;

                     } /* switch (permut) */
                     
                     endket: /* no viable CSF in this Q1-configuration */
                     /* re-initialize */
                     j = norbi - 2;            /* by-pass a wasted search step */
                     if (j>-1) {               
                        jvrtx = kpath[j];
                        kmax = kstep[j+1] - 1;   
                     } else {
                        break;                 /* normal exit for 1-orbital internal space (rare!) */
                     }
                  } /* if (j==norbi) */
            
               } else {   
                  /* try to descend */   
                  if (!j) break;               /* normal exit */
                  jvrtx = kpath[--j];
                  kmax = kstep[j+1] - 1;   
               } /* if (successful_ascent) else ... *for ket* */
            }/* for (;;)  *for ket* */   
         } /* loop over kmcr (evaluation of coupling coefficients) */

         if (iwrite_elmt) {
            if (elmt22tot>bunlen) {
               if (elmt22tot%bunlen) {
                  if (fwrite(cc_two, (elmt22tot%bunlen)*sizeof(struct cc_elmt2), 1, ciftfile)!=1){ 
                     (void)fprintf(outfile, "\n\nfwrite error (%s)\n\n", CIFTNAME);
                     MPI_Abort(MPI_COMM_WORLD, 1);
                  }
               } 
            } 
         } 

         /* evaluation of sigma vector */
         switch (nelxb) {

            case 0:
               /* E0                  */
               xorbq2[0] = 0;
               calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 0);
               calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                       jcnt, nelxb, nelxk, extsym); 
               break;

            case 1:
               /* E1                  */
               xorbq2[1] = 0;
               for (a=exsym[extsym]; a>=sxsym[extsym]; a--) {
                  xorbq2[0] = a;
                  calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                  calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                          jcnt, nelxb, nelxk, extsym); 
               }
               break;

            case 2:
               /* E2                  */
               xorbq2[1] = 0;
               for (a=norbx; a>0; a--) {
                  xorbq2[0] = a;
                  calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                  calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                          jcnt, nelxb, nelxk, extsym); 
               }
               break;

            case 3:
               /* E11                 */
               xorbq2[2] = 0;
               for (ast=0; ast<nst; ast++) {
                  bst = GRPMUL(extsym,ast);
                  if (bst == ast) {
                     for (a=exsym[ast]; a>sxsym[ast]; a--) {
                        xorbq2[0] = a;  /* a > b  */
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=a-1; b>=sxsym[ast]; b--) {
                           xorbq2[1] = b;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                   jcnt, nelxb, nelxk, extsym); 
                        }
                     }
                  } else if (bst > ast) {
                     for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                        xorbq2[0] = a;  /* a > b  */
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                           xorbq2[1] = b;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                   jcnt, nelxb, nelxk, extsym); 
                        }
                     }
                  } /* if (bst==ast) else if (bst>ast) */
               } /* loop over ast */
               break;

            case 4:
               /* E21                 */
               xorbq2[2] = 0;
               for (a=exsym[extsym]; a>=sxsym[extsym]; a--) {
                  xorbq2[1] = a;
                  calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                  intk[norbi0+1] = hdiag;
                  for (k=0; k<ncsfq2tmp; k++) 
                     hdiagtmp2[k] = hdiagtmp[k];
                  for (b=norbx; b>a; b--) {
                     xorbq2[0] = b;
                     calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                     calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                             jcnt, nelxb, nelxk, extsym); 
                  }
                  b--;
                  for (; b>0; b--) {
                     xorbq2[0] = b;
                     calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                     calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                             jcnt, nelxb, nelxk, extsym); 
                  }
               }
               break;

            case 5:
               /* E111                */
               xorbq2[3] = 0;
               /* three external open shells */
               /* a > b > c  unless specified otherwise   */
               /* xorbq2[0] > xorbq2[1] > xorbq2[2]       */
               for (ast=0; ast<nst; ast++) {
                  if (ast == extsym) {
                     /* ast = bst = cst = extsym */
                     for (a=exsym[ast]; a>(sxsym[ast]+1); a--) {
                        xorbq2[0] = a;
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=a-1; b>sxsym[ast]; b--) {
                           xorbq2[1] = b;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           intk[norbi0+2] = hdiag;
                           for (c=b-1; c>=sxsym[ast]; c--) {
                              xorbq2[2] = c;
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                              calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                      jcnt, nelxb, nelxk, extsym); 
                           }
                        }
                     }
                     /* ast < bst = cst, ast = extsym */
                     for (bst=ast+1; bst<nst; bst++) {
                        for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                           xorbq2[0] = a;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (b=exsym[bst]; b>sxsym[bst]; b--) {
                              xorbq2[1] = b;
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (c=b-1; c>=sxsym[bst]; c--) {
                                 xorbq2[2] = c;
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                         jcnt, nelxb, nelxk, extsym); 
                              }
                           }
                        }
                     } /* loop over bst */
                  } /* if (ast==extsym) */
                  for (bst=ast+1; bst<nst; bst++) {
                     if (bst == extsym) {
                        /* ast = cst < bst = extsym */
                        /* a > c > b                */
                        for (a=exsym[ast]; a>sxsym[ast]; a--) {
                           xorbq2[0] = a;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (c=a-1; c>=sxsym[ast]; c--) {
                              xorbq2[1] = c;
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                                 xorbq2[2] = b;
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                         jcnt, nelxb, nelxk, extsym); 
                              }
                           }
                        }
                     } /* if (bst==extsym) */
                     cst = GRPMUL(extsym,ast);
                     cst = GRPMUL(bst,cst);
                     if (cst>bst) {
                        for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                           xorbq2[0] = a;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                              xorbq2[1] = b;
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (c=exsym[cst]; c>=sxsym[cst]; c--) {
                                 xorbq2[2] = c;
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                         jcnt, nelxb, nelxk, extsym); 
                              }
                           }
                        }
                     } /* if (cst>bst) */
                  } /* loop over bst */
               } /* loop over ast */
               /* end of three external open shells */
               break;

            case 6:
               /* E220                */
               xorbq2[2] = 0;
               /* xorbq2[0] > xorbq2[1]       */
               /* a > b */
               for (a=norbx; a>1; a--) {
                  xorbq2[0] = a;
                  calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                  intk[norbi0+1] = hdiag;
                  for (k=0; k<ncsfq2tmp; k++) 
                     hdiagtmp2[k] = hdiagtmp[k];
                  for (b=a-1; b>0; b--) {
                     xorbq2[1] = b;
                     calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                     calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                             jcnt, nelxb, nelxk, extsym); 
                  }
               }
               break;

            case 7:
               /* E211                */
               xorbq2[3] = 0;
               /* xorbq2[1] > xorbq2[2]       */
               for (ast=0; ast<nst; ast++) {
                  bst = GRPMUL(extsym,ast);
                  if (bst == ast) {
                     for (a=exsym[ast]; a>sxsym[ast]; a--) {
                        xorbq2[1] = a;  /* a > b  */
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=a-1; b>=sxsym[ast]; b--) {
                           xorbq2[2] = b;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           intk[norbi0+2] = hdiag;
                           for (k=0; k<ncsfq2tmp; k++) 
                              hdiagtmp2[k] = hdiagtmp[k];
                           for (c=norbx; c>0; c--) {
                              if (c==a) continue;
                              if (c==b) continue;
                              xorbq2[0] = c; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                              calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                      jcnt, nelxb, nelxk, extsym); 
                           }
                        }
                     }
                  } else if (bst > ast) {
                     for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                        xorbq2[1] = a;  /* a > b  */
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                           xorbq2[2] = b;
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           intk[norbi0+2] = hdiag;
                           for (k=0; k<ncsfq2tmp; k++) 
                              hdiagtmp2[k] = hdiagtmp[k];
                           for (c=norbx; c>0; c--) {
                              if (c==a) continue;
                              if (c==b) continue;
                              xorbq2[0] = c; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                              calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                      jcnt, nelxb, nelxk, extsym); 
                           }
                        }
                     }
                  } /* if (bst==ast) else if (bst>ast) */
               } /* loop over ast */
               break;

            case 8:
               /* E1111               */
               /* four external open shells */
               /* a > b > c > d  unless specified otherwise   */
               /* xorbq2[0] > xorbq2[1] > xorbq2[2] > xorbq2[3]       */
               if (!extsym) {
                  for (ast=0; ast<nst; ast++) {
                     /* ast = bst = cst = dst */
                     for (a=exsym[ast]; a>(sxsym[ast]+2); a--) {
                        xorbq2[0] = a; 
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=a-1; b>(sxsym[ast]+1); b--) {
                           xorbq2[1] = b; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           intk[norbi0+2] = hdiag;
                           for (c=b-1; c>sxsym[ast]; c--) {
                              xorbq2[2] = c; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                              intk[norbi0+3] = hdiag;
                              for (d=c-1; d>=sxsym[ast]; d--) {
                                 xorbq2[3] = d; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                 calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                         jcnt, nelxb, nelxk, extsym); 
                              }
                           }
                        }
                     } /* end of ast = bst = cst = dst */

                     /* ast = bst < cst = dst */
                     for (cst=ast+1; cst<nst; cst++) {
                        for (a=exsym[ast]; a>sxsym[ast]; a--) {
                           xorbq2[0] = a; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (b=a-1; b>=sxsym[ast]; b--) {
                              xorbq2[1] = b; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (c=exsym[cst]; c>sxsym[cst]; c--) {
                                 xorbq2[2] = c; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 intk[norbi0+3] = hdiag;
                                 for (d=c-1; d>=sxsym[cst]; d--) {
                                    xorbq2[3] = d; 
                                    calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                    calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                            jcnt, nelxb, nelxk, extsym); 
                                 }
                              }
                           }
                        }
                     } /*  end of ast = cst < bst = dst */
                  } /* end of loop over ast */
               } /* end of if(!extsym) */

               for (ast=0; ast<nst; ast++) {
                  bst = GRPMUL(extsym,ast);
                  if (bst > ast) {

                     for (cst=0; cst<ast; cst++) {
                        /* cst = dst < ast < bst */
                        for (c=exsym[cst]; c>sxsym[cst]; c--) {
                           xorbq2[0] = c; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (d=c-1; d>=sxsym[cst]; d--) {
                              xorbq2[1] = d; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                                 xorbq2[2] = a; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 intk[norbi0+3] = hdiag;
                                 for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                                    xorbq2[3] = b; 
                                    calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                    calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                            jcnt, nelxb, nelxk, extsym); 
                                 }
                              }
                           }
                        }
                     } /* cst = dst < ast < bst */

                     /* ast = cst = dst < bst */
                     for (a=exsym[ast]; a>(sxsym[ast]+1); a--) {
                        xorbq2[0] = a; 
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (c=a-1; c>sxsym[ast]; c--) {
                           xorbq2[1] = c; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           intk[norbi0+2] = hdiag;
                           for (d=c-1; d>=sxsym[ast]; d--) {
                              xorbq2[2] = d; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                              intk[norbi0+3] = hdiag;
                              for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                                 xorbq2[3] = b; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                 calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                         jcnt, nelxb, nelxk, extsym); 
                              }
                           }
                        }
                     } /* ast = cst = dst < bst */

                     /* ast < cst = dst < bst */
                     for (cst=ast+1; cst<bst; cst++) {
                        for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                           xorbq2[0] = a; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (c=exsym[cst]; c>sxsym[cst]; c--) {
                              xorbq2[1] = c; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (d=c-1; d>=sxsym[cst]; d--) {
                                 xorbq2[2] = d; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 intk[norbi0+3] = hdiag;
                                 for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                                    xorbq2[3] = b; 
                                    calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                    calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                            jcnt, nelxb, nelxk, extsym); 
                                 }
                              }
                           }
                        }
                     } /* ast < cst = dst < bst */

                     /* ast < bst = cst = dst */
                     for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                        xorbq2[0] = a; 
                        calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                        intk[norbi0+1] = hdiag;
                        for (b=exsym[bst]; b>(sxsym[bst]+1); b--) {
                           xorbq2[1] = b; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                           intk[norbi0+2] = hdiag;
                           for (c=b-1; c>sxsym[bst]; c--) {
                              xorbq2[2] = c; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                              intk[norbi0+3] = hdiag;
                              for (d=c-1; d>=sxsym[bst]; d--) {
                                 xorbq2[3] = d; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                 calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                         jcnt, nelxb, nelxk, extsym); 
                              }
                           }
                        }
                     }  /* ast < bst = cst = dst */

                     /* ast < bst < cst = dst */
                     for (cst=bst+1; cst<nst; cst++) {
                        for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                           xorbq2[0] = a; 
                           calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                           intk[norbi0+1] = hdiag;
                           for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                              xorbq2[1] = b; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                              intk[norbi0+2] = hdiag;
                              for (c=exsym[cst]; c>sxsym[cst]; c--) {
                                 xorbq2[2] = c; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                 intk[norbi0+3] = hdiag;
                                 for (d=c-1; d>=sxsym[cst]; d--) {
                                    xorbq2[3] = d; 
                                    calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                    calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                            jcnt, nelxb, nelxk, extsym); 
                                 }
                              }
                           }
                        }
                     } /* ast < bst < cst = dst */
                  } /* end of if (bst > ast) */

                  xst = GRPMUL(extsym,ast);
                  for (bst=ast+1; bst<nst; bst++) {
                     for (cst=bst+1; cst<nst; cst++) {
                        dst = GRPMUL(cst,GRPMUL(xst,bst));
                        if (dst>cst) {
                           for (a=exsym[ast]; a>=sxsym[ast]; a--) {
                              xorbq2[0] = a; 
                              calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0], 1);
                              intk[norbi0+1] = hdiag;
                              for (b=exsym[bst]; b>=sxsym[bst]; b--) {
                                 xorbq2[1] = b; 
                                 calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+1], 2);
                                 intk[norbi0+2] = hdiag;
                                 for (c=exsym[cst]; c>=sxsym[cst]; c--) {
                                    xorbq2[2] = c; 
                                    calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+2], 3);
                                    intk[norbi0+3] = hdiag;
                                    for (d=exsym[dst]; d>=sxsym[dst]; d--) {
                                       xorbq2[3] = d; 
                                       calchdiag_bw(nstep, dbra, dbmap, nelxb, nopenb[norbi], intk[norbi0+3], 4);
                                       calc_hwpp(nstep, dbra, dbmap, nopenb[norbi], ptr_kup, ptr_arc, ptr_ketbase, 
                                               jcnt, nelxb, nelxk, extsym); 
                                    }
                                 }
                              }
                           }
                        } /* if(dst>cst) */
                     } /* loop over cst */
                  } /* loop over bst */
               } /* loop over ast */
               /* end of four external open shells */

               break;
         } /* switch (nelxb)  */

         /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
            
         endbra: /* no viable CSF in this Q2-configuration */
         /* re-initialize */
         i = norbi - 2;            /* by-pass a wasted search step */
         if (i>-1) {               
            ivrtx = path[i];
            kimax = nstep[i+1] - 1;   
         } else {
            break;                 /* normal exit for 1-orbital internal space (rare!) */
         }

      } /* if (i==norbi) */
      
   } else {   
      /* try to descend */   
      if (!i) break;               /* normal exit */
      ivrtx = path[--i];
      kimax = nstep[i+1] - 1;   
   } /* if (successful_ascent) else ... *for bra* */
 } /* for (;;) *for bra* */
  
}

void twoptwoh(int iint, int jint, int nopenb, int dbra[], int dbmap[],
 int nopenk, int dket[], int dkmap[], 
 double (*segment[])(int, int, int, int, double *))
/* calculate hamiltonian matrix elements between two configurations:
      bra lane has 2 excess particles, and
      ket lane has 2 excess particles; 
   first occurrence of difference at orbital iint,
   last occurrence at orbital jint */ 
{
 extern FILE *outfile;
 extern int norbi, spin, arcsp[], jspdn[], jspup[], spnhead[], delspn[], bospn[], dospn[];
 extern double mtrxelt[], mtrxel1[];
 extern const double zthresh;
 int braspn[MXKMAX+1], bstep[MXKMAX+1], kstep[MXKMAX+1], ketspn[MXKMAX+1];
 int ketwt[MXKMAX+1];
 int bralev, ketlev, ketlevmx, tmp, icsf, i, kvrtx, incflag, ketlevmn;
 int pint, qint, itmp, qintmx, qintmn, ketwttl, spnmn;
 double val, val1;
  
 void x22_null(double elmt0, double elmt1, int iwalk, int jwalk);
 
 bralev = nopenb;
 if (spin>bralev) return;
 ketlev = nopenk;
 if (spin>ketlev) return;
 
 braspn[bralev] = spin;
 itmp = (nopenb) ? dbmap[0] : 0;
 for (i=norbi; i>=itmp; i--) bospn[i] = spin;
 pint = itmp;
 icsf = 0;            /* internal and local */
 incflag = 1;
 dospn[norbi] = spin;
 qintmn = jint - 1;
  
 for (;;) {
      
   if (bralev==0) {
      /* internal bra csf generated */
      
      /* >>>>>> generate *interacting* ket csfs in ascending order <<<<<< */
      incflag = 1;
         
      /* >>>>> extract fixed parts of ket spin vector <<<<< */
      ketlev = 0;
      ketwttl = 0;
      dospn[0] = bospn[0];
      ketspn[0] = braspn[0];
      kvrtx = 0;
      for (i=1; i<jint; i++) {
         if (dbra[i]==1) {
            dket[i] = 1;
            dospn[i] = dospn[i-1] + 1;
            kstep[ketlev] = 1;
            ketlev++;
            ketspn[ketlev] = ketspn[ketlev-1] + 1;
            kvrtx = jspup[2*kvrtx+1];
            ketwttl += arcsp[2*kvrtx+1];
         } else if (dbra[i]==2) {
            dket[i] = 2;
            dospn[i] = dospn[i-1] - 1;
              kstep[ketlev] = -1;
            ketlev++;
            ketspn[ketlev] = ketspn[ketlev-1] - 1;
            kvrtx = jspup[2*kvrtx];
            ketwttl += arcsp[2*kvrtx];
         } else {
            dospn[i] = dospn[i-1];
         }
      }
      ketlevmn = ketlev;             /* minimum accessible level */
      spnmn = ketspn[ketlevmn];      /* spin at minimum accessible level */
      
      ketlev = nopenk;
      ketspn[ketlev] = spin;
      ketwt[ketlev] = ketwttl;
      kvrtx = spnhead[ketlev];
      for (i=norbi; i>iint; i--) {
         if (dbra[i]==1) {
            dket[i] = 1;
            dospn[i-1] = dospn[i] - 1;
            kstep[ketlev] = 1;
            ketlev--;
            ketspn[ketlev] = ketspn[ketlev+1] - 1;
            ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
            kvrtx = jspdn[2*kvrtx+1];
         } else if (dbra[i]==2) {
            dket[i] = 2;
            dospn[i-1] = dospn[i] + 1;
              kstep[ketlev] = -1;
            ketlev--;
            ketspn[ketlev] = ketspn[ketlev+1] + 1;
            ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx];
            kvrtx = jspdn[2*kvrtx];
         } else {
            dospn[i-1] = dospn[i];
         }
      }
      ketlevmx = ketlev;        /* maximum changeable level */
      qintmx = (ketlev>0) ? dkmap[nopenk-ketlev] : 0;
      if (qintmx < jint) qintmx = jint - 1;
      /* >>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<< */
         
      mtrxelt[iint+1] = 1.0;
      mtrxel1[iint+1] = 1.0;   
      for (i=iint; i>qintmx; i--) {
         delspn[i] = dospn[i] - bospn[i];
         val = (*segment[i])(dbra[i], dket[i], delspn[i], dospn[i], &val1);
         mtrxelt[i] = mtrxelt[i+1] * val;
         mtrxel1[i] = mtrxel1[i+1] * val1;
         dospn[i-1] = dospn[i];
      }
         
      qint = qintmx;   
      for (;;) {   
         if (ketlev==ketlevmn) {
            if (ketspn[ketlev]==spnmn) {
               /* internal ket csf generated */
               x22_null(mtrxelt[qintmn+1], mtrxel1[qintmn+1], icsf, ketwt[ketlev]);
            }
             
            /* re-initialize */
            ketlev++;
            if (ketlev>ketlevmx) break;   /* normal exit for no internal open shells */
            kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
            incflag = 0;
            qint = dkmap[nopenk-ketlev];
            
         } else if (incflag) { 
            /* try to descend (ketlev => ketlev-1) with spin increase, if possible */
            tmp = ketspn[ketlev] + 1;
            if (tmp < ketlev) {
               kstep[ketlev] = -1;
               dket[qint] = 2;
               dospn[qint-1] = dospn[qint] + 1;
               ketlev--;
               ketspn[ketlev] = tmp;
               ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx];
               kvrtx = jspdn[2*kvrtx];
               
            } else {
               kstep[ketlev] = 1;
               dket[qint] = 1;
               dospn[qint-1] = dospn[qint] - 1;
               ketlev--;
               ketspn[ketlev] = tmp - 2;
               ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
               kvrtx = jspdn[2*kvrtx+1];
            }
               
            itmp = (ketlev>0) ? dkmap[nopenk-ketlev] : 0;
            if (itmp<qintmn) itmp = qintmn;
            for (i=qint; i>itmp; i--) {
               delspn[i] = dospn[i] - bospn[i];
               val = (*segment[i])(dbra[i], dket[i], delspn[i], dospn[i], &val1);
               mtrxelt[i] = mtrxelt[i+1] * val;
               mtrxel1[i] = mtrxel1[i+1] * val1;
               if (i<qint) dospn[i-1] = dospn[i];
            }
            
            if (fabs(mtrxel1[itmp+1])<zthresh) { 
               if (fabs(mtrxelt[itmp+1])<zthresh) {
                  /* reject descent */
                  incflag = 0;
                  ketlev++;
                  kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
               } else      
                  qint = itmp;
            } else 
               qint = itmp;
            
         } else if (kstep[ketlev]==-1) {
            /* try to descend with spin decrease */
            tmp = ketspn[ketlev] - 1;
            if (tmp >= 0) {
               kstep[ketlev] = 1;
               dket[qint] = 1;
               dospn[qint-1] = dospn[qint] - 1;
               ketlev--;
               ketspn[ketlev] = tmp;
               ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
               kvrtx = jspdn[2*kvrtx+1];
               incflag = 1;
               
               itmp = (ketlev>0) ? dkmap[nopenk-ketlev] : 0;
               if (itmp<qintmn) itmp = qintmn;
               for (i=qint; i>itmp; i--) {
                  delspn[i] = dospn[i] - bospn[i];
                  val = (*segment[i])(dbra[i], dket[i], delspn[i], dospn[i], &val1);
                  mtrxelt[i] = mtrxelt[i+1] * val;
                  mtrxel1[i] = mtrxel1[i+1] * val1;
                  if (i<qint) dospn[i-1] = dospn[i];
               }
               
               if (fabs(mtrxel1[itmp+1])<zthresh) { 
                  if (fabs(mtrxelt[itmp+1])<zthresh) {
                     /* reject descent */
                     incflag = 0;
                     ketlev++;
                     kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
                  } else      
                     qint = itmp;
               } else 
                  qint = itmp;
                  
            } else {
               /* try to ascend */
               ketlev++;
               if (ketlev>ketlevmx) break;   /* normal exit */
               kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
               qint = dkmap[nopenk-ketlev];
            }
               
         } else {
            /* try to ascend */
            ketlev++;
            if (ketlev>ketlevmx) break;   /* normal exit */
            kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
            qint = dkmap[nopenk-ketlev];
         }
      }   
      /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
            
      icsf++;
            
      /* re-initialize */
      bralev++;
      if (bralev>nopenb) break;   /* normal exit for no internal open shells */
      incflag = 0; 
      pint = dbmap[nopenb-bralev];
         
   } else if (incflag) { 
      /* try to descend (bralev => bralev-1) with spin increase, if possible */
      tmp = braspn[bralev] + 1;
      if (tmp < bralev) {
         bstep[bralev] = -1;
         dbra[pint] = 2;
         bralev--;
         braspn[bralev] = tmp;
         
      } else {
         bstep[bralev] = 1;
         dbra[pint] = 1;
         bralev--;
         braspn[bralev] = tmp - 2;
      }
         
      itmp = (bralev>0) ? dbmap[nopenb-bralev] : 0;
      for (i=pint-1; i>=itmp; i--) {
         bospn[i] = braspn[bralev];
      }
      pint = itmp; 
         
   } else if (bstep[bralev]==-1) {
      /* try to descend with spin decrease */
      tmp = braspn[bralev] - 1;
      if (tmp >= 0) {
         bstep[bralev] = 1;
         dbra[pint] = 1;
         bralev--;
         braspn[bralev] = tmp;
         incflag = 1;
         
         itmp = (bralev>0) ? dbmap[nopenb-bralev] : 0;
         for (i=pint-1; i>=itmp; i--) {
            bospn[i] = braspn[bralev];
         }
         pint = itmp; 
         
      } else {
         /* try to ascend */
         bralev++;
         if (bralev>nopenb) break;   /* normal exit */
         pint = dbmap[nopenb-bralev];
      }
      
   } else {
      /* try to ascend */
      bralev++;
      if (bralev>nopenb) break;   /* normal exit */
      pint = dbmap[nopenb-bralev];
   }   
 }
 
}

void x22_null(double elmt0, double elmt1, int iwalk, int jwalk)
{
 /* complete the matrix element for two particles and two holes */
 extern int ovrlap, jcfg, nstates, *elmt22cnt, elmt22tot, bunlen, iwrite_elmt;
 extern struct cc_elmt2 *cc_two, *ptr_cc_two;
 extern FILE *ciftfile;
 extern FILE *outfile;

 /* iwalk:  Q1-space 
    jwalk:  Q2-space  */

 ptr_cc_two->icsf = jwalk*nstates; 
 ptr_cc_two->jcsf = iwalk*nstates; 
 if (ovrlap==1) { 
   ptr_cc_two->exch_cc = elmt0 + elmt1; 
   ptr_cc_two->drct_cc = elmt0 - elmt1; 
 } else if (!ovrlap) {
   ptr_cc_two->exch_cc = elmt0 + elmt1; 
   ptr_cc_two->drct_cc = -2. * elmt0; 
 } else {
   ptr_cc_two->drct_cc = elmt0; 
 }
 elmt22cnt[jcfg]++;
 ptr_cc_two++;
 if (iwrite_elmt) {
   elmt22tot++;
   if (!(elmt22tot%bunlen)) {
      if (fwrite(cc_two, bunlen*sizeof(struct cc_elmt2), 1, ciftfile)) {
          (void)fprintf(outfile, "\n\nfwrite error (%s)\n\n", CIFTNAME);
          MPI_Abort(MPI_COMM_WORLD, 1);
      } 
      ptr_cc_two = cc_two;
   } 
 }

}

void oneponeh(int iint, int jint, int nopenb, int dbra[], int dbmap[],
 int nopenk, int dket[], int dkmap[])
/* calculate hamiltonian matrix elements between two configurations:
		bra lane has 1 excess particle, and
		ket lane has 1 excess particle; 
   excess bra orbital iint > jint */ 
{
 extern FILE *outfile;
 extern int norbi, spin, arcsp[], jspdn[], jspup[], spnhead[], delspn[], bospn[], dospn[];
 extern int nspndif[], snglptr[];
 extern double headcum[][MXKMAX+1], acum[][MXKMAX+1], bcum[][MXKMAX+1], ccum[][MXKMAX+1];
 extern double mtrxelt[];          /* added 08.05.05 to accommodate J=0 for external */
 int braspn[MXKMAX+1], bstep[MXKMAX+1], kstep[MXKMAX+1], ketspn[MXKMAX+1], ketwt[MXKMAX+1];
 int bralev, ketlev, ketlevmx, tmp, icsf, i, j, kvrtx, incflag;
 int pint, qint, itmp, qintmx, opnhd, opntl, opnmd, opnj, accpt;
 int ispnptr, jspnptr, nelxb, nelxk, numa;
 double val, val1;
  
 double wail(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wasl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbir(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbpsr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wcp(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wcpp(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdilsl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdirl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdsrl(int dbra, int dket, int delspn, int ketspn, double *val1);
 void x11_null(double elmt, double acum[], double bcum[], double ccum[], 
  int numa, int numb, int numc, int iwalk, int jwalk);
 
 nelxb = 0;
 nelxk = 0;
 bralev = nopenb+nelxb;
 if (spin>bralev) return;
 ketlevmx = nopenk+nelxk;
 if (spin>ketlevmx) return;
 ispnptr = snglptr[bralev];
 jspnptr = snglptr[ketlevmx];
 
 braspn[bralev] = spin;
 itmp = (nopenb) ? dbmap[0] : 0;
 for (i=norbi; i>=itmp; i--) bospn[i] = spin;
 pint = itmp;
 icsf = 0;
 incflag = 1;
 dospn[norbi] = spin;
 for (i=ketlevmx,opnhd=0,opnmd=nelxk,opntl=nelxk,opnj=0; i>nelxk; i--) {
	itmp = dkmap[nopenk-(i-nelxk)];
	if (itmp > iint) opnhd++;
	else if (itmp < iint) {
		opnmd++;
		if (itmp < jint) opntl++;
		else if (itmp == jint) opnj++;
	}
 }
 nspndif[ketlevmx+1] = 0;
 numa = opnmd-(opntl+opnj);
	
 for (;;) {
		
	if (bralev==nelxb) {
		/* internal bra csf generated */
			 
		/* >>>>>> generate *interacting* ket csfs in ascending order <<<<<< */
		incflag = 1;
			
		ketlev = ketlevmx;
		ketspn[ketlev] = spin;
		ketwt[ketlev] = 0;
		kvrtx = spnhead[ketlev];
		qintmx = qint = (ketlev>nelxk) ? dkmap[nopenk-(ketlev-nelxk)] : 0;
		for (i=norbi; i>qint; i--) {
			dospn[i-1] = dospn[i];
			if (i==iint) {
				delspn[iint] = dospn[iint] - bospn[iint];
				(void)wail(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
				acum[opnmd][opnmd] = val;
			} else if (i==jint) {
				delspn[jint] = dospn[jint] - bospn[jint];
				val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
				ccum[opntl][opntl] = acum[opnmd][opnmd] * val1;
				mtrxelt[opntl] = acum[opnmd][opnmd] * val;          /* added 08.05.05 for J=0 */
				(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
				acum[opnmd][opnmd] *= val;
			}
		}
			
		for (;;) {	
			if (ketlev==nelxk) {
				/* internal ket csf generated */
#ifdef DEBUG
				if (qint) {
					(void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
#endif
					/*x11_null(acum[nelxk][opnmd], &acum[nelxk][opntl+opnj], &bcum[nelxk][ketlevmx-opnhd+1], &ccum[nelxk][nelxk], 
	   	 		 opnmd-(opntl+opnj), opnhd, opntl-nelxk, icsf, ketwt[nelxk]);*/
					x11_null(acum[nelxk][opnmd], &acum[nelxk][opntl+opnj], &bcum[nelxk][ketlevmx-opnhd+1], &ccum[nelxk][nelxk], 
	   	 		 numa, opnhd, opntl, icsf, ketwt[0]);
			 	
				/* re-initialize */
				ketlev++;
				if (ketlev>ketlevmx) break;	/* normal exit for no internal open shells */
				kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
				incflag = 0;
				qint = dkmap[nopenk-(ketlev-nelxk)]; 
#ifdef DEBUG
				if (qint>qintmx) {
					(void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
#endif
				
			} else if (incflag) { 
				/* try to descend (ketlev => ketlev-1) with spin increase, if possible */
				tmp = ketspn[ketlev] + 1;
				if (tmp < ketlev) {
					kstep[ketlev] = -1;
					dket[qint] = 2;
					dospn[qint-1] = dospn[qint] + 1;
					ketlev--;
					ketspn[ketlev] = tmp;
					ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx];
					kvrtx = jspdn[2*kvrtx];
					
				} else {
					kstep[ketlev] = 1;
					dket[qint] = 1;
					dospn[qint-1] = dospn[qint] - 1;
					ketlev--;
					ketspn[ketlev] = tmp - 2;
					ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
					kvrtx = jspdn[2*kvrtx+1];
				}
					
				itmp = (ketlev>nelxk) ? dkmap[nopenk-(ketlev-nelxk)] : 0;
				for (i=qint-1; i>itmp; i--) dospn[i-1] = dospn[i];
				
				if (qint > iint) {
					accpt = 0;
					delspn[qint] = dospn[qint] - bospn[qint];
					(void)wcpp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
					for (j=ketlev+2; j<=ketlevmx; j++) {
						headcum[ketlev+1][j] = headcum[ketlev+2][j] * val;
						if (headcum[ketlev+1][j]!=0.0) accpt++;
					}
					if (nspndif[ketlev+2]) val = 0.;
					else {
						accpt++;
						(void)wdirl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
					}
					headcum[ketlev+1][ketlev+1] = val;
					if (dket[qint]==dbra[qint]) nspndif[ketlev+1] = nspndif[ketlev+2];
					else nspndif[ketlev+1] = nspndif[ketlev+2] + 1;
					
					if (itmp < iint) {
						accpt = 0;
						delspn[iint] = dospn[iint] - bospn[iint];
						if (nspndif[ketlev+1]) val = 0.;
						else {	
							(void)wail(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
							if (val!=0.0) accpt++; 
						}
						acum[opnmd][opnmd] = val;
						(void)wbpsr(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[opnmd][j] = headcum[ketlevmx-opnhd+1][j] * val;
							if (bcum[opnmd][j]!=0.0) accpt++;
						}
					}
					
					if (itmp < jint) {
						accpt = 0;
						delspn[jint] = dospn[jint] - bospn[jint];
						val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
						ccum[opntl][opntl] = acum[opnmd][opnmd] * val1;
						mtrxelt[opntl] = acum[opnmd][opnmd] * val;           /* added 08.05.05 */
						if (ccum[opntl][opntl]!=0.0) accpt++;
						if (mtrxelt[opntl]!=0.0) accpt++;
						(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
						acum[opnmd][opnmd] *= val;
						if (acum[opnmd][opnmd]!=0.0) accpt++;
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[opnmd][j] *= val;
							if (bcum[opnmd][j]!=0.0) accpt++;
						}
					}
					
				} else if (qint == iint) {
					accpt = 0;
					delspn[iint] = dospn[iint] - bospn[iint];
					if (nspndif[ketlev+2]) val = 0.;
					else {
						(void)wail(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
						if (val!=0.0) accpt++;
					}
					acum[opnmd][opnmd] = val;
					(void)wbpsr(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
					for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
						bcum[opnmd][j] = headcum[ketlevmx-opnhd+1][j] * val;
						if (bcum[opnmd][j]!=0.0) accpt++;
					}
					
					if (itmp < jint) {
						accpt = 0;
						delspn[jint] = dospn[jint] - bospn[jint];
						val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
						ccum[opntl][opntl] = acum[opnmd][opnmd] * val1;
						if (ccum[opntl][opntl]!=0.0) accpt++;
						mtrxelt[opntl] = acum[opnmd][opnmd] * val;
						if (mtrxelt[opntl]!=0.0) accpt++;
						(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
						acum[opnmd][opnmd] *= val;
						if (acum[opnmd][opnmd]!=0.0) accpt++;
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[opnmd][j] *= val;
							if (bcum[opnmd][j]!=0.0) accpt++;
						}
					}
				
				} else if (qint > jint) {
					accpt = 0;
					delspn[qint] = dospn[qint] - bospn[qint];
					(void)wcp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
					for (j=ketlev+1; j<=opnmd; j++) {
						acum[ketlev][j] = acum[ketlev+1][j] * val;
						if (acum[ketlev][j]!=0.0) accpt++;
					}
					for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
						bcum[ketlev][j] = bcum[ketlev+1][j] * val;
						if (bcum[ketlev][j]!=0.0) accpt++;
					}
					(void)wdilsl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
					acum[ketlev][ketlev] = acum[ketlev+1][opnmd] * val;
					if (acum[ketlev][ketlev]!=0.0) accpt++;
					
					if (itmp < jint) {
						accpt = 0;
						delspn[jint] = dospn[jint] - bospn[jint];
						val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
						ccum[opntl][opntl] = acum[ketlev][opnmd] * val1;
						if (ccum[opntl][opntl]!=0.0) accpt++;
						mtrxelt[opntl] = acum[ketlev][opnmd] * val;
						if (mtrxelt[opntl]!=0.0) accpt++;
						(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
						for (j=ketlev+1; j<=opnmd; j++) {
							acum[ketlev][j] *= val;
							if (acum[ketlev][j]!=0.0) accpt++;
						}
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[ketlev][j] *= val;
							if (bcum[ketlev][j]!=0.0) accpt++;
						}
					}
						
				} else if (qint == jint) {
					accpt = 0;
					delspn[jint] = dospn[jint] - bospn[jint];
					val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
					ccum[opntl][opntl] = acum[ketlev+1][opnmd] * val1;
					if (ccum[opntl][opntl]!=0.0) accpt++;
					mtrxelt[opntl] = acum[ketlev+1][opnmd] * val;
					if (mtrxelt[opntl]!=0.0) accpt++;
					(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
					for (j=ketlev+1; j<=opnmd; j++) {
						acum[ketlev][j] = acum[ketlev+1][j] * val;
						if (acum[ketlev][j]!=0.0) accpt++;
					}
					for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
						bcum[ketlev][j] = bcum[ketlev+1][j] * val;
						if (bcum[ketlev][j]!=0.0) accpt++;
					}
					
				} else {
					accpt = 0;
					delspn[qint] = dospn[qint] - bospn[qint];
					(void)wdsrl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
					ccum[ketlev][ketlev] = ccum[ketlev+1][opntl] * val;
					if (ccum[ketlev][ketlev]!=0.0) accpt++;
					val = wcpp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val1);
					ccum[ketlev][opntl] = ccum[ketlev+1][opntl] * val1;
					if (ccum[ketlev][opntl]!=0.0) accpt++;
					mtrxelt[ketlev] = mtrxelt[ketlev+1] * val;		
					if (mtrxelt[ketlev]!=0.0) accpt++;
					/* below closed loops */
					if (dbra[qint]==dket[qint]) {
						for (j=opntl+opnj; j<=opnmd; j++) {
							acum[ketlev][j] = acum[ketlev+1][j];
							if (acum[ketlev][j]!=0.0) accpt++;
						}
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[ketlev][j] = bcum[ketlev+1][j];
							if (bcum[ketlev][j]!=0.0) accpt++;
						}
						for (j=ketlev+1; j<opntl; j++) {
							ccum[ketlev][j] = ccum[ketlev+1][j];
							if (ccum[ketlev][j]!=0.0) accpt++;
						}
					} else {
						for (j=opntl+opnj; j<=opnmd; j++) acum[ketlev][j] = 0.;
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) bcum[ketlev][j] = 0.;
						for (j=ketlev+1; j<opntl; j++) ccum[ketlev][j] = 0.;
					}
				}
					
				/* reject descent */ 
				if (!accpt) { 
					incflag = 0;
					ketlev++;
					kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
				} else		
					qint = itmp;
				
			} else if (kstep[ketlev]==-1) {
				/* try to descend with spin decrease */
				tmp = ketspn[ketlev] - 1;
				if (tmp >= 0) {
					kstep[ketlev] = 1;
					dket[qint] = 1;
					dospn[qint-1] = dospn[qint] - 1;
					ketlev--;
					ketspn[ketlev] = tmp;
					ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
					kvrtx = jspdn[2*kvrtx+1];
					incflag = 1;
					
					itmp = (ketlev>nelxk) ? dkmap[nopenk-(ketlev-nelxk)] : 0;
					for (i=qint-1; i>itmp; i--) dospn[i-1] = dospn[i];
					
					/* >>>>> beginning of level evaluation <<<<< */	
					if (qint > iint) {
						accpt = 0;
						delspn[qint] = dospn[qint] - bospn[qint];
						(void)wcpp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
						for (j=ketlev+2; j<=ketlevmx; j++) {
							headcum[ketlev+1][j] = headcum[ketlev+2][j] * val;
							if (headcum[ketlev+1][j]!=0.0) accpt++;
						}
						if (nspndif[ketlev+2]) val = 0.;
						else {
							accpt++;
							(void)wdirl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
						}
						headcum[ketlev+1][ketlev+1] = val;
						if (dket[qint]==dbra[qint]) nspndif[ketlev+1] = nspndif[ketlev+2];
						else nspndif[ketlev+1] = nspndif[ketlev+2] + 1;
						
						if (itmp < iint) {
							accpt = 0;
							delspn[iint] = dospn[iint] - bospn[iint];
							if (nspndif[ketlev+1]) val = 0.;
							else {	
								(void)wail(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
								if (val!=0.0) accpt++; 
							}
							acum[opnmd][opnmd] = val;
							(void)wbpsr(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
							for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
								bcum[opnmd][j] = headcum[ketlevmx-opnhd+1][j] * val;
								if (bcum[opnmd][j]!=0.0) accpt++;
							}
						}
						
						if (itmp < jint) {
							accpt = 0;
							delspn[jint] = dospn[jint] - bospn[jint];
							val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
							ccum[opntl][opntl] = acum[opnmd][opnmd] * val1;
							if (ccum[opntl][opntl]!=0.0) accpt++;
							mtrxelt[opntl] = acum[opnmd][opnmd] * val;
							if (mtrxelt[opntl]!=0.0) accpt++;
							(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
							acum[opnmd][opnmd] *= val;
							if (acum[opnmd][opnmd]!=0.0) accpt++;
							for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
								bcum[opnmd][j] *= val;
								if (bcum[opnmd][j]!=0.0) accpt++;
							}
						}
						
					} else if (qint == iint) {
						accpt = 0;
						delspn[iint] = dospn[iint] - bospn[iint];
						if (nspndif[ketlev+2]) val = 0.;
						else {
							(void)wail(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
							if (val!=0.0) accpt++;
						}
						acum[opnmd][opnmd] = val;
						(void)wbpsr(dbra[iint], dket[iint], delspn[iint], dospn[iint], &val);
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[opnmd][j] = headcum[ketlevmx-opnhd+1][j] * val;
							if (bcum[opnmd][j]!=0.0) accpt++;
						}
							
						if (itmp < jint) {
							accpt = 0;
							delspn[jint] = dospn[jint] - bospn[jint];
							val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
							ccum[opntl][opntl] = acum[opnmd][opnmd] * val1;
							if (ccum[opntl][opntl]!=0.0) accpt++;
							mtrxelt[opntl] = acum[opnmd][opnmd] * val;
							if (mtrxelt[opntl]!=0.0) accpt++;
							(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
							acum[opnmd][opnmd] *= val;
							if (acum[opnmd][opnmd]!=0.0) accpt++;
							for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
								bcum[opnmd][j] *= val;
								if (bcum[opnmd][j]!=0.0) accpt++;
							}
						}
						
					} else if (qint > jint) {
						accpt = 0;
						delspn[qint] = dospn[qint] - bospn[qint];
						(void)wcp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
						for (j=ketlev+1; j<=opnmd; j++) {
							acum[ketlev][j] = acum[ketlev+1][j] * val;
							if (acum[ketlev][j]!=0.0) accpt++;
						}
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[ketlev][j] = bcum[ketlev+1][j] * val;
							if (bcum[ketlev][j]!=0.0) accpt++;
						}
						(void)wdilsl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
						acum[ketlev][ketlev] = acum[ketlev+1][opnmd] * val;
						if (acum[ketlev][ketlev]!=0.0) accpt++;
						
						if (itmp < jint) {
							accpt = 0;
							delspn[jint] = dospn[jint] - bospn[jint];
							val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
							ccum[opntl][opntl] = acum[ketlev][opnmd] * val1;
							if (ccum[opntl][opntl]!=0.0) accpt++;
							mtrxelt[opntl] = acum[ketlev][opnmd] * val;
							if (mtrxelt[opntl]!=0.0) accpt++;
							(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
							for (j=ketlev+1; j<=opnmd; j++) {
								acum[ketlev][j] *= val;
								if (acum[ketlev][j]!=0.0) accpt++;
							}
							for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
								bcum[ketlev][j] *= val;
								if (bcum[ketlev][j]!=0.0) accpt++;
							}
						}
							
					} else if (qint == jint) {
						accpt = 0;
						delspn[jint] = dospn[jint] - bospn[jint];
						val = wbir(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val1);
						ccum[opntl][opntl] = acum[ketlev+1][opnmd] * val1;
						if (ccum[opntl][opntl]!=0.0) accpt++;
						mtrxelt[opntl] = acum[ketlev+1][opnmd] * val;
						if (mtrxelt[opntl]!=0.0) accpt++;
						(void)wasl(dbra[jint], dket[jint], delspn[jint], dospn[jint], &val);
						for (j=ketlev+1; j<=opnmd; j++) {
							acum[ketlev][j] = acum[ketlev+1][j] * val;
							if (acum[ketlev][j]!=0.0) accpt++;
						}
						for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
							bcum[ketlev][j] = bcum[ketlev+1][j] * val;
							if (bcum[ketlev][j]!=0.0) accpt++;
						}
						
					} else {
						accpt = 0;
						delspn[qint] = dospn[qint] - bospn[qint];
						(void)wdsrl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
						ccum[ketlev][ketlev] = ccum[ketlev+1][opntl] * val;
						if (ccum[ketlev][ketlev]!=0.0) accpt++;
						val = wcpp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val1);
						ccum[ketlev][opntl] = ccum[ketlev+1][opntl] * val1;
						if (ccum[ketlev][opntl]!=0.0) accpt++;
						mtrxelt[ketlev] = mtrxelt[ketlev+1] * val;
						if (mtrxelt[ketlev]!=0.0) accpt++;
						/* below closed loops */
						if (dbra[qint]==dket[qint]) {
							for (j=opntl+opnj; j<=opnmd; j++) {
								acum[ketlev][j] = acum[ketlev+1][j];
								if (acum[ketlev][j]!=0.0) accpt++;
							}
							for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) {
								bcum[ketlev][j] = bcum[ketlev+1][j];
								if (bcum[ketlev][j]!=0.0) accpt++;
							}
							for (j=ketlev+1; j<opntl; j++) {
								ccum[ketlev][j] = ccum[ketlev+1][j];
								if (ccum[ketlev][j]!=0.0) accpt++;
							}
						} else {
							for (j=opntl+opnj; j<=opnmd; j++) acum[ketlev][j] = 0.;
							for (j=ketlevmx-opnhd+1; j<=ketlevmx; j++) bcum[ketlev][j] = 0.;
							for (j=ketlev+1; j<opntl; j++) ccum[ketlev][j] = 0.;
						}
					}
					/* >>>>>>>> end of level evaluation <<<<<<<< */
					
					/* reject descent */ 
					if (!accpt) { 
						incflag = 0;
						ketlev++;
						kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
					} else		
						qint = itmp;
					
				} else {
					/* try to ascend */
					ketlev++;
					if (ketlev>ketlevmx) break;	/* normal exit */
					kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
					qint = dkmap[nopenk-(ketlev-nelxk)];
#ifdef DEBUG
					if (qint>qintmx) {
						(void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
						MPI_Abort(MPI_COMM_WORLD, 1);
					}
#endif
				}
					
			} else {
				/* try to ascend */
				ketlev++;
				if (ketlev>ketlevmx) break;	/* normal exit */
				kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
				qint = dkmap[nopenk-(ketlev-nelxk)];
#ifdef DEBUG
				if (qint>qintmx) {
					(void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
#endif
			}
		}	
		/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
				
		icsf++;
				
		/* re-initialize */
		bralev++;
		if (bralev>(nopenb+nelxb)) break;	/* normal exit for no internal open shells */
		incflag = 0; 
		pint = dbmap[nopenb-(bralev-nelxb)];
			
	} else if (incflag) { 
		/* try to descend (bralev => bralev-1) with spin increase, if possible */
		tmp = braspn[bralev] + 1;
		if (tmp < bralev) {
			bstep[bralev] = -1;
			dbra[pint] = 2;
			bralev--;
			braspn[bralev] = tmp;
			
		} else {
			bstep[bralev] = 1;
			dbra[pint] = 1;
			bralev--;
			braspn[bralev] = tmp - 2;
		}
			
		itmp = (bralev>nelxb) ? dbmap[nopenb-(bralev-nelxb)] : 0;
		for (i=pint-1; i>=itmp; i--) {
			bospn[i] = braspn[bralev];
		}
		pint = itmp; 
			
	} else if (bstep[bralev]==-1) {
		/* try to descend with spin decrease */
		tmp = braspn[bralev] - 1;
		if (tmp >= 0) {
			bstep[bralev] = 1;
			dbra[pint] = 1;
			bralev--;
			braspn[bralev] = tmp;
			incflag = 1;
			
			itmp = (bralev>nelxb) ? dbmap[nopenb-(bralev-nelxb)] : 0;
			for (i=pint-1; i>=itmp; i--) {
				bospn[i] = braspn[bralev];
			}
			pint = itmp; 
			
		} else {
			/* try to ascend */
			bralev++;
			if (bralev>(nopenb+nelxb)) break;	/* normal exit */
			pint = dbmap[nopenb-(bralev-nelxb)];
		}
		
	} else {
		/* try to ascend */
		bralev++;
		if (bralev>(nopenb+nelxb)) break;	/* normal exit */
		pint = dbmap[nopenb-(bralev-nelxb)];
	}	
 }
 
}

void x11_null(double elmt, double acum[], double bcum[], double ccum[], 
 int numa, int numb, int numc, int iwalk, int jwalk)
{

 /* complete the matrix element for one-particle-one-hole */
 extern int nstates, switch_bra_ket, jcfg11, *elmt11cnt;
 extern struct cc_elmt1 *ptr_cc_one;

 extern FILE *outfile;

 int k, ik;
 if (switch_bra_ket) {
    ptr_cc_one->icsf = iwalk*nstates;
    ptr_cc_one->jcsf = jwalk*nstates;
 } else {
    ptr_cc_one->icsf = jwalk*nstates;
    ptr_cc_one->jcsf = iwalk*nstates;
 }
 ptr_cc_one->elmtone = elmt;
 for (k=0, ik=0; k<numc; k++, ik++) ptr_cc_one->acum[ik] = ccum[k]; 
 for (k=0; k<numa; k++,ik++) ptr_cc_one->acum[ik] = acum[k];
 for (k=0; k<numb; k++,ik++) ptr_cc_one->acum[ik] = bcum[k];
 elmt11cnt[jcfg11]++;
 ptr_cc_one++;

}

double wail(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of A/L */
{
 extern double pb01[], pb12[];
 double val;
 
 switch (dbra) {
 case 1:
   val = pb01[ketspn];
   break;
 case 2:
   val =-1. / pb12[ketspn];
   break; 
 default:
   val = -1.;
   break;
 }
 *val1 = val;
 return val;
}

double wasl(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of A^L */
{
 double val;
 extern double pb01[], pb12[];
 
 switch (delspn) {
 case 1:
   switch (dbra) {
   case 0:
      if (dket==1) val = 1./pb01[ketspn];
      else val = 0.;
      break;
   case 1:
      val = 0.;
      break;
   default:
      val = 1.;
      break;
   }
   break;
 case -1:
   switch (dbra) {
   case 0:
      if (dket==1) val = 0.;
      else val = -pb12[ketspn];
      break;
   case 1:
      val = 1.;
      break;
   default:
      val = 0.;
      break;
   }
   break;
 default:
   val = 0.;
   break;
 }
 *val1 = val;
 return val;
}
 
double wasr(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of A^R */
{
 extern double pb01[], pb12[];
 double val;
 
 switch (delspn) {
 case 1:
   switch (dbra) {
   case 1:
      val = 0.;
      break;
   case 2:
      val = 1.;
      break;
   default:
      if (dket==1) val = -1./pb01[ketspn];
      else val = 0.;
      break;
   }
   break;
 case -1:
   switch (dbra) {
   case 1:
      val = 1.;
      break;
   case 2:
      val = 0.;
      break;
   default:
      if (dket==1) val = 0.;
      else val = pb12[ketspn];
      break;
   }
   break;
 default:
   val = 0.;
   break;
 }
 *val1 = val;
 return val;
}
 
double wair(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of A/R */
{
 extern double pb01[], pb12[];
 double val;
  
 switch (dbra) {
 case 0:
   val = 1.;
   break;
 case 1:
   val = pb01[ketspn];
   break;
 default:
   val = -1. / pb12[ketspn];
   break;
 }
 *val1 = val;
 return val;
}

double wbil(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value pf B\L */
 extern FILE *outfile;
 extern double isqt, pbm0[], pb01[], pb02[], pb32[], pb12[];
 double val;
 
 switch (delspn) {
 case 1:
   switch (dbra) {
   case 1:
      val = 0.;
      *val1 = -pbm0[ketspn];
      break;
      
   case 2:
      val = -isqt;
      *val1 = isqt / pb02[ketspn];
      break;
      
   default:
      if (dket==1) {
         val = isqt / pb01[ketspn];
         *val1 = -isqt * pbm0[ketspn];
      } else {
         val = 0.;
         *val1 = -1.;
      }
      break;
   }
   break;
 
 case -1:
   switch (dbra) {
   case 1:
      val = -isqt;
      *val1 = -isqt * pb02[ketspn];
      break;
      
   case 2:
      val = 0.;
      *val1 = pb32[ketspn];
      break;
   
   default:
      if (dket==1) {
         val = 0.;
         *val1 = -1.;
      } else {
         val = -isqt * pb12[ketspn];
         *val1 = -isqt * pb32[ketspn];
      }
      break;   
   }
   break;
   
 default:
   val = 0.;
   *val1 = 0.;
   break;   
 }
 return val;
}

double wbsr(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of B^R */
 extern double isqt, pb01[], pb12[], pb02[];
 double val;
 
 switch (delspn) {
 case 0:
   switch (dbra) {
   case 1:
      val = isqt * pb01[ketspn];
      *val1 = -isqt / pb12[ketspn];
      break;
   case 2:
      val = -isqt / pb12[ketspn];
      *val1 = -isqt * pb01[ketspn];
      break;
   default:
      val = -isqt;
      if (dket==1) *val1 = -isqt / pb02[ketspn];
      else *val1 = isqt * pb02[ketspn];
      break;
   }
   break;
   
 case 2:
   val = 0.;
   if (dbra==2) *val1 = -1.;
   else if (dbra==3) {
      if (dket==1) *val1 = -1. / pb01[ketspn];
      else *val1 = 0.;
   } else *val1 = 0.;
   break;
   
 case -2:
   val = 0.;
   if (dbra==1) *val1 = -1.;
   else if (dbra==3) {
      if (dket==2) *val1 = pb12[ketspn];
      else *val1 = 0.;
   } else *val1 = 0.;
   break;
   
 default:
   val = 0.;
   *val1 = 0.;
   break;
 }
 return val;
}

double wbpsr(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of B'^R */
 extern double isqt, pb01[], pb12[], pb02[];
 double val;
 
 switch (delspn) {
 case 0:
   switch (dbra) {
   case 1:
      val = isqt * pb01[ketspn];
      *val1 = isqt / pb12[ketspn];
      break;
   case 2:
      val = -isqt / pb12[ketspn];
      *val1 = isqt * pb01[ketspn];
      break;
   default:
      val = -isqt;
      if (dket==1) *val1 = isqt / pb02[ketspn];
      else *val1 = -isqt * pb02[ketspn];
      break;
   }
   break;
   
 case 2:
   val = 0.;
   if (dbra==2) *val1 = 1.;
   else if (dbra==3) {
      if (dket==1) *val1 = 1. / pb01[ketspn];
      else *val1 = 0.;
   } else *val1 = 0.;
   break;
   
 case -2:
   val = 0.;
   if (dbra==1) *val1 = 1.;
   else if (dbra==3) {
      if (dket==2) *val1 = -pb12[ketspn];
      else *val1 = 0.;
   } else *val1 = 0.;
   break;
   
 default:
   val = 0.;
   *val1 = 0.;
   break;
 }
 return val;
}

double wbsl(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of B^L for J=0 & 1 */
{
 extern FILE *outfile;
 extern double isqt, pb01[], pb12[], pb02[];
 double val;
  
 switch (delspn) {
 case 0:
   switch (dbra) {
   case 0:
      val = isqt;
      if (dket==1) *val1 = -isqt / pb02[ketspn];
      else *val1 = isqt * pb02[ketspn];
      break;
   case 1:
      val = isqt * pb01[ketspn];
      *val1 = isqt / pb12[ketspn]; 
      break;
   case 2:
      val =-isqt / pb12[ketspn];
      *val1 = isqt * pb01[ketspn];
      break;
   }
   break;
   
 case 2:
   val = 0.;
   switch (dbra) {
   case 0:
      if (dket==1) *val1 = -1. / pb01[ketspn];
      else *val1 = 0.;
      break;
   case 2:
      *val1 = 1.;
      break;
   default:
      *val1 = 0.;
      break;
   }
   break;

 case -2:   
   val = 0.;
   switch (dbra) {
   case 0:
      if (dket==2) *val1 = pb12[ketspn];
      else *val1 = 0.;
      break;
   case 1:
      *val1 = 1.;
      break;
   default:
      *val1 = 0.;
      break;   
   }
   break;

 default:
   val = 0.;
   *val1 = 0.;
   break;   
 }
 return val;
}

double wbpsl(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of (B')^L for J=0 & 1 */
{
 extern FILE *outfile;
 extern double isqt, pb01[], pb12[], pb02[];
 double val;
  
 switch (delspn) {
 case 0:
   switch (dbra) {
   case 0:
      val = isqt;
      if (dket==1) *val1 = isqt / pb02[ketspn];
      else *val1 = -isqt * pb02[ketspn];
      break;
   case 1:
      val = isqt * pb01[ketspn];
      *val1 = -isqt / pb12[ketspn]; 
      break;
   case 2:
      val =-isqt / pb12[ketspn];
      *val1 = -isqt * pb01[ketspn];
      break;
   }
   break;
   
 case 2:
   val = 0.;
   switch (dbra) {
   case 0:
      if (dket==1) *val1 = 1. / pb01[ketspn];
      else *val1 = 0.;
      break;
   case 2:
      *val1 = -1.;
      break;
   default:
      *val1 = 0.;
      break;
   }
   break;

 case -2:   
   val = 0.;
   switch (dbra) {
   case 0:
      if (dket==2) *val1 = -pb12[ketspn];
      else *val1 = 0.;
      break;
   case 1:
      *val1 = -1.;
      break;
   default:
      *val1 = 0.;
      break;   
   }
   break;

 default:
   val = 0.;
   *val1 = 0.;
   break;   
 }
 return val;
}

double wbir(int dbra, int dket, int delspn, int ketspn, double *val1)
 /* segment value of B/R for J=0 & 1 */
{
 extern double isqt, pb01[], pb12[], pb02[], pb32[], pbm0[];
 double val;
  
 switch (delspn) {
 case -1:
   switch (dbra) {
   case 0:
      if (dket==1) {
         val = 0.;
         *val1 = -1.;
      } else {
         val = isqt * pb12[ketspn];
         *val1 = -isqt * pb32[ketspn];
      }
      break;
   case 1:
      val = -isqt;
      *val1 = isqt * pb02[ketspn]; 
      break;
   default:
      val = 0.;
      *val1 = -pb32[ketspn];
      break;
   }
   break;
   
 case 1:
   switch (dbra) {
      case 0:
         if (dket==1) {
            val = -isqt / pb01[ketspn];
            *val1 = -isqt * pbm0[ketspn];
         } else {
            val = 0.;
            *val1 = -1.;
         }
         break;
      case 1:
         val = 0.;
         *val1 = pbm0[ketspn];
         break;
      default:
         val = -isqt;
         *val1 = -isqt / pb02[ketspn];
         break;
      }
      break;
      
 default:
   val = 0.;
   *val1 = 0.;
   break;   
 } 
 return val;
}

double wcp(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of C' */
 extern FILE *outfile;
 extern double sq21[], sq01[];
 double val;
 
 if (delspn==-1) {
   switch (dbra) {
   case 1:
      if (dket==1) val = -1.;
      else val = -1./(ketspn+2);
      break;
   case 2:
      if (dket==1) val = 0.;
      else val = -sq21[ketspn];
      break;
   default:
      val = 1.;
      break;
   }
 
 } else if (delspn==1) {
   switch (dbra) {
   case 1:
      if (dket==1) val = -sq01[ketspn];
      else val = 0.;
      break;
   case 2:
      if (dket==1) val = 1./ketspn;
      else val = -1.;
      break;
   default:
      val = 1.;
      break;
   }
   
 } else val = 0.;

 *val1 = val; 
 return val;
} 

double wcpp(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of C" */
 extern FILE *outfile;
 extern double sq02[], sqm2[], sq10[], sq20[], sq12[], sq22[], sq00[], sqm0[];
 double val;
 
 switch (dbra) {
 case 1:
   if (dket==1) {
      switch (delspn) {
      case 0:
         val = 1.;
         *val1 = sq02[ketspn];
         break;
      case 2:
         val = 0.;
         *val1 = sqm2[ketspn];
         break;
      case -2:
         val = 0.;
         *val1 = 1.;
         break;
      default:
         val = 0.;
         *val1 = 0.;
         break;
      }
   } else {
      val = 0.;
      switch (delspn) {
      case 0:
         *val1 = sq10[ketspn];
         break;
      case -2:
         *val1 = sq20[ketspn];
         break;
      default:
         *val1 = 0.;
         break;   
      }
   }
   break;
   
 case 2:
   if (dket==2) {
      switch (delspn) {
      case 0:
         val = 1.;
         *val1 = sq12[ketspn];
         break;
      case 2:
         val = 0.;
         *val1 = 1.;
         break;
      case -2:
         val = 0.;
         *val1 = sq22[ketspn];
         break;
      default:
         val = 0.;
         *val1 = 0.;
         break;
      }
   } else {
      val = 0.;
      switch (delspn) {
      case 0:
         *val1 = -sq00[ketspn];
         break;
      case 2:
         *val1 = -sqm0[ketspn];
         break;
      default:
         *val1 = 0.;
         break;
      }
   }
    break;
   
 default:
   switch (delspn) {
   case 0:
      val = 1.;
      if (ketspn) *val1 = 1.;
      else *val1 = 0.;
      break;
   case 2:
      val = 0.;
      *val1 = 1.;
      break;
   case -2:
      val = 0.;
      *val1 = 1.;
      break;
   default:
      val = 0.;
      *val1 = 0.;
      break;
   }
   break;
 } 
 return val; 
}

double wdill(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D/LL */
 extern double sqt;
 double val;
 
 *val1 = 0.;   
 if (!delspn) val = -sqt;
 else val = 0.;
 return val; 
}

double wdirr(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D/RR */
 extern double sqt;
 double val;
 
 *val1 = 0.;   
 if (!delspn) val = sqt;
 else val = 0.;
 return val; 
}

double wdilsl(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D/L^L */
 extern double pb01[], pb12[], sq01[], sq21[];
 double val;
    
 switch (delspn) {
 case -1:
   switch (dbra) {
   case 1:
      if (dket==1) val = 0.;
      else {
         val = pb12[ketspn];
         val *= -val;
      }
      break;
   case 2:
      if (dket==1) val = 0.;
      else val = sq21[ketspn];   
      break;
   case 3:
      val = -1.;
      break;
   default:
      val = 0.;
      break;   
   }
   break;
      
 case 1:
   switch (dbra) {
   case 1:
      if (dket==1) val = sq01[ketspn];
      else val = 0.;
      break;
   case 2:
      if (dket==1) {
         val = 1. / pb01[ketspn];
         val *= -val;
      } else val = 0.;
      break;
   case 3:
      val = -1.;
      break;
   default:
      val = 0.;
      break;
   }
   break;
   
 default:
   val = 0.;
   break;
 }
 *val1 = val; 
 return val;
}
 
double wdilsr(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D/L^R */
 double val;
  
 *val1 = 0.;
 switch (delspn) {
 case 1:
   val = -1.;
   break;
 case -1:
   val = -1.;
   break;
 default:
   val = 0.;
   break;
 }
 return val;
}

double wdirsl(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D/R^L */
 double val;
  
 *val1 = 0.;
 switch (delspn) {
 case 1:
   val = 1.;
   break;
 case -1:
   val = 1.;
   break;
 default:
   val = 0.;
   break;
 }
 return val;
}

double wdirl(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D/RL */
 extern double sqt, isqt, pbm1[], pb01[], pb12[], pb31[];
 double val;
 
 switch (dbra) {
 case 1:
   if (dket==1) {
      val = -isqt;
      *val1 = -isqt * pbm1[ketspn];
   } else {
      val = 0.;
      *val1 = -pb01[ketspn];
   }
   break;

 case 2:
   if (dket==1) {
      val = 0.;
      *val1 = 1. / pb12[ketspn];
   } else {
      val = -isqt;
      *val1 = isqt * pb31[ketspn];
   }
   break;

 case 3:
   val = -sqt;
   *val1 = 0.;
   break;
   
 default:
   val = 0.;
   *val1 = 0.;
   break;
 }
 return val;
}
 
double wdsll(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D^LL */
 extern double sqt;
 double val;
 
 *val1 = 0.;   
 if (!delspn) val = -sqt;
 else val = 0.;
 return val; 
}

double wdsll2(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D^LL */
 extern double sqt2;
 double val;
 
 *val1 = 0.;   
 if (!delspn) val = -sqt2;
 else val = 0.;
 return val; 
}

double wdsrr(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D^RR */
 extern double sqt;
 double val;
 
 *val1 = 0.;   
 if (!delspn) val = sqt;
 else val = 0.;
 return val; 
}

double wdsrr2(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D^RR */
 extern double sqt2;
 double val;
 
 *val1 = 0.;   
 if (!delspn) val = sqt2;
 else val = 0.;
 return val; 
}

double wdsrl(int dbra, int dket, int delspn, int ketspn, double *val1)
{
 /* segment value of D^RL */
 extern double sqt, isqt, pb01[], pb02[], pb12[];
 double val;
 
 switch (delspn) {
 case 2:
   val = 0.;
   if (dbra==2) {
      if (dket==1) *val1 = 1. / pb01[ketspn];
      else *val1 = 0.;
   } else *val1 = 0.;
   break;
 case 0:
   switch (dbra) {
   case 1:
      if (dket==1) {
         val = isqt;
         *val1 = isqt / pb02[ketspn];
      } else {
         val = 0.;
         *val1 = 0.;
      }
      break;
   case 2:
      if (dket==2) {
         val = isqt;
         *val1 = -isqt * pb02[ketspn];
      } else {
         val = 0.;
         *val1 = 0.;
      }
      break;
   case 3:
      val = sqt;
      *val1 = 0.;
      break;
   default:
      val = 0.;
      *val1 = 0.;
      break;
   }
   break;   
 case -2:
   val = 0.;
   if (dbra==1) {
      if (dket==2) *val1 = -pb12[ketspn];
      else *val1 = 0.;
   } else *val1 = 0.;
   break;
 default:
   val = 0.;
   *val1 = 0.;
   break;
 }
 return val;
}

int chksym(int nelxb, int nelxk, int extsym, int extsymk, int *permut) 
{
  /* check if the Q1-configurations interact with Q2-configuration */

 extern int nxsym[], nxprsym[], nxtrsym[], nxqrsym[];
 extern FILE *outfile;

 *permut = 0;
 switch (nelxk) {
    case 0:
       if (extsymk) return 0;
       break;
    case 1:
       if (!nxsym[extsymk]) return 0;
       switch (nelxb) {
          case 1:
             *permut = 1;
             break;
          case 2:
             *permut = 2;
             break;
          case 3:   /* two open shells */
             *permut = 3;
             break;
          case 4:
             *permut = 4;
             break;
          case 5:
             if (!nxprsym[GRPMUL(extsym, extsymk)]) return 0;
             *permut = 5;  
             break;
       }
       return nxsym[extsymk];
       break;
    case 2:   /* closed shell    */
       if (extsymk) return 0;
       switch (nelxb) {
          case 3:   /* two open shells */
             *permut = 6;
             break;
          case 4:
             *permut = 7;
             break;
          case 5:
             *permut = 8;  
             break;
       }
       return norbx;
       break;
    case 3:   /* two open shells */
       if (!nxprsym[extsymk]) return 0;
       switch (nelxb) {
          case 1:
             if (!nxsym[GRPMUL(extsym, extsymk)]) return 0;
             *permut = 9;
             break;
          case 2:
             *permut = 10;
             break;
          case 3:   /* two open shells */
             if (extsym!=extsymk) {
                if (!nxprsym[GRPMUL(extsym, extsymk)]) return 0;
             }
             *permut = 11;
             break;
          case 4:
             *permut = 12;
             break;
          case 5:
             *permut = 13;
             break;
          case 7:
             if (extsym!=extsymk) {
                if (!nxprsym[GRPMUL(extsym, extsymk)]) return 0;
             }
             *permut = 14;
             break;
          case 8:
             if (!nxprsym[GRPMUL(extsym, extsymk)]) return 0;
             *permut = 15;
             break;
       }
       return nxprsym[extsymk];
       break;

    default:  /* Q2-configurations */
       switch (nelxb) {
          case 0:
             if (extsym) return 0;
             break;
          case 1:
             return nxsym[extsym]; 
             break;
          case 2:   /* closed shell    */
             if (extsym) return 0;
             return norbx; 
             break;
          case 3:   /* two open shells */
             return nxprsym[extsym]; 
             break;
          case 4:
             return nxsym[extsym]*(norbx-1); 
             break;
          case 5:
             return nxtrsym[extsym]; 
             break;
          case 6:
             if (extsym) return 0;
             return (norbx*(norbx-1))/2; 
             break;
          case 7:
             return (norbx-2)*nxprsym[extsym]; 
             break;
          case 8:
             return nxqrsym[extsym]; 
             break;
       }
       break;
 } 

 return 1;

}

/* add the contributtion to Hpp    */

void calc_hwpp(int nstep[], int dbra[], int dbmap[], int nopenb, int **ptr_kup, 
             int **ptr_arc, int **ptr_ketbase, int jcnt, int nelxb, int *nelxk, 
             int extsym) 
{

 extern FILE *outfile, *scrfile, *ciftfile;
 extern int norbi0, norbi, norbx, orbsymi[], spin, switch_bra_ket, irrep, nnorbx;
 extern int sxsym[], exsym[], ncsfcfg[], orbsym[]; 
 extern int mxkmax, iter, mxiter, invflag, prntflag;
 extern int intgrlmo[], xorbq2[], *lccsf_11;
 extern int elmt22tot, bunlen, iwrite_elmt;
 extern int jcfg, jcfg11, ncsfq1tmp, ncsfq2tmp, nstates, ncsfxstate;
 extern double *twoint, *oneint;
 extern double *hqp, *hqwp, *energs, **hwpp, **wwpp, **upp; 
 extern double *energs2, **upp2, **hwpp2, **wwpp2, **hqq2, *hqq2inv;
 extern double *hqq, hdiag, *buff;
 extern struct cc_elmt2 *cc_two, *ptr_cc_two;
 extern struct cc_elmt1 *cc_one, *ptr_cc_one;
#ifdef DEBUG
 extern int norb;
 static int ncsfq2sum=0;
#endif

 int i, j, k, kk, jvrtx, kmax, elmt22sum, bunloc;
 int kpath[NINTMX+1], ksymcum[NINTMX+1], kstep[NINTMX+1], cfgcum[NINTMX+1];
 int xsp[NINTMX+1], xsh[NINTMX+1];
 int *kupmcr, *arcmcr, *ketbase, jcsfbase;
 int kmcr, ncfgtmp, extsymk, kst, abmo, bcmo, akmo, bkmo, jcfgtmp; 
 int successful_ascent, permut, exsp, exsh;
 int nopenk[NINTMX+1], dket[NINTMX+1], dkmap[MXKMAX+1], iorb[4], jorb[4];
 int imo, jmo, kmo, lmo, itmp, iread_elmt;
 int icsf, jcsf, ijcsf, csfloc, jcsf_loc, istate, jstate, icsf_k, jcsf_k; 
 double tmp, hscr[MXSTATE];
 struct cc_elmt2 *ptr22tmp, *ptr22tmp2;

 int chksym(int nelxb, int nelxk, int extsym, int extsymk, int *permut); 
 void getintgrl_one(int iorb, int jorb, int jcsf_loc, int dbra[]); 
 void getintgrl_two(int imo, int jmo, int kmo, int lmo, int jcfgq1, int jcsf_loc); 
 void getintgrl_twox(int imo, int jmo, int kmo, int lmo, int jcfgq1, int jcsf_loc); 
 int  choleski(double *smat, double *sinv, double *buff, int n);

 /* loop over all viable Q1-macroconfigurations. */
 /* kup => ptr_kup[kmcr] */
 /* arc => ptr_arc[kmcr] */
 /* ketbase => ptr_ketbase[kmcr] */
 xsp[0] = 0;
 xsh[0] = 0;
 jcfg = 0;
 jcfg11 = 0;
 if (iwrite_elmt) {
    iread_elmt = 0;
    if (elmt22tot > bunlen) {
       elmt22sum = 0;
       bunloc = 0;
       iread_elmt = 1;
       rewind(ciftfile);
       if (fread(cc_two, bunlen*sizeof(struct cc_elmt2), 1, ciftfile)) {
          (void)fprintf(outfile, "\n\nfread error (%s)\n\n", CIFTNAME); 
          MPI_Abort(MPI_COMM_WORLD, 1);
       } 
   }
 }
 ptr_cc_one = cc_one;
 ptr_cc_two = cc_two;
 for (i=0; i<ncsfxstate; i++) hqp[i] = 0.0;

#ifdef DEBUG
    (void)fprintf(outfile, "\n");
    if (nstep[1]==4) {
      (void)fprintf(outfile, "bra: (1");
    } else {
      (void)fprintf(outfile, "bra: (%d", nstep[1]);
    }
    for (j=2; j<=norbi0; j++) {
      if (nstep[j]==4) {
        (void)fprintf(outfile, " 1");
      } else {
        (void)fprintf(outfile, " %d", nstep[j]);
      }
    }
    switch (nelxb) {
      case 0:
        (void)fprintf(outfile, ")\n");
        break;
      case 1:
        (void)fprintf(outfile, " 1)   ext{%d}\n", xorbq2[0]);
        break;
      case 2:
        (void)fprintf(outfile, " 3)   ext{%d}\n", xorbq2[0]);
        break;
      case 3:
        (void)fprintf(outfile, " 1 1) ext{%d, %d}\n", 
        xorbq2[0], xorbq2[1]);
        break;
      case 4:
        (void)fprintf(outfile, " 2 1) ext{%d, %d}\n", 
        xorbq2[0], xorbq2[1]);
        break;
      case 5:
        (void)fprintf(outfile, " 1 1 1) ext{%d, %d, %d}\n", 
        xorbq2[0], xorbq2[1], xorbq2[2]);
        break;
      case 6:
        (void)fprintf(outfile, " 2 2) ext{%d, %d}\n", 
        xorbq2[0], xorbq2[1]);
        break;
      case 7:
        (void)fprintf(outfile, " 2 1 1) ext{%d, %d, %d}\n", 
        xorbq2[0], xorbq2[1], xorbq2[2]);
        break;
      case 8:
        (void)fprintf(outfile, " 1 1 1 1) ext{%d, %d, %d, %d}\n", 
        xorbq2[0], xorbq2[1], xorbq2[2], xorbq2[3]);
        break;
    }
    fflush(outfile);
#endif

 for (kmcr=0; kmcr<jcnt; kmcr++) {

    kupmcr = ptr_kup[kmcr];
    arcmcr = ptr_arc[kmcr];
    ketbase = ptr_ketbase[kmcr]; /* pointer to the starting CSF in this macroconfiguration */

    switch (nelxb) {
       case 6: case 7: case 8:  
          exsp = 0;
          exsh = 2;
          break;
       case 5: 
          switch (nelxk[kmcr]) {
             case 1: 
                exsp = 0;
                exsh = 2;
                break;
             case 2:  
                exsp = 0;
                exsh = 1;
                break;
             case 3: 
                exsp = 1;
                exsh = 2;
                break;
          }
          break;
       case 4: 
          switch (nelxk[kmcr]) {
             case 1: 
                exsp = 0;
                exsh = 2;
                break;
             case 2: case 3: 
                exsp = 1;
                exsh = 2;
                break;
          }
          break;
       case 3: 
          switch (nelxk[kmcr]) {
             case 0: 
                exsp = 0;
                exsh = 2;
                break;
             case 1: 
                exsp = 1;
                exsh = 2;
                break;
             case 2:  
                exsp = 1;
                exsh = 1;
                break;
             case 3: 
                exsp = 2;
                exsh = 2;
                break;
          }
          break;
       case 2: 
          switch (nelxk[kmcr]) {
             case 0: 
                exsp = 0;
                exsh = 2;
                break;
             case 1: 
                exsp = 1;
                exsh = 2;
                break;
             case 2:  
                exsp = 2;
                exsh = 2;
                break;
             case 3: 
                exsp = 1;
                exsh = 1;
                break;
          }
          break;
       case 1: 
          switch (nelxk[kmcr]) {
             case 0: 
                exsp = 1;
                exsh = 2;
                break;
             case 1: 
                exsp = 2;
                exsh = 2;
                break;
             case 2: case 3: 
                exsp = 2;
                exsh = 1;
                break;
          }
          break;
       case 0:
          switch (nelxk[kmcr]) {
             case 0: 
                exsp = exsh = 2;
                break;
             case 1: 
                exsp = 2;
                exsh = 1;
                break;
          }
          break;
    }
    /* >>>>>> generate *interacting* ket walks in ascending order <<<<<< */
    jvrtx = 0;
    kpath[0] = jvrtx;
    ksymcum[0] = 0;
    cfgcum[0] = 0;
    kmax = 2;
    j = 0;
    nopenk[0] = 0;
    dket[0] = 0;
           
    for (;;) {
       /* try to ascend (j => j+1) */
       successful_ascent = 0;
       kk = jvrtx*3;

       for (k=kmax; k>=0; k--) {
          if (kupmcr[kk+k]>-1) {
             kstep[j+1] = k;
             /* bra and ket switched here for iorb[] and jorb[] */
             switch (kstep[j+1] - nstep[j+1]) {
                case 2:
                   iorb[xsp[j]] = norbi-j;
                   iorb[xsp[j]+1] = norbi-j;
                   xsp[j+1] = xsp[j] + 2;
                   xsh[j+1] = xsh[j];
                   break;
                case 1:
                   iorb[xsp[j]] = norbi-j;
                   xsp[j+1] = xsp[j] + 1;
                   xsh[j+1] = xsh[j];
                   break;
                case -1:
                   jorb[xsh[j]] = norbi-j;
                   xsp[j+1] = xsp[j];
                   xsh[j+1] = xsh[j] + 1;
                   break;
                case -2:
                   jorb[xsh[j]] = norbi-j;
                   jorb[xsh[j]+1] = norbi-j;   
                   xsp[j+1] = xsp[j];
                   xsh[j+1] = xsh[j] + 2;
                   break;
                default:
                   xsp[j+1] = xsp[j];
                   xsh[j+1] = xsh[j];
             }
             /* bra and ket switched here for xsp[] and xsh[] */
             if (j<norbi0) {
                if ((xsp[j+1]>exsh)||(xsh[j+1]>exsp)) continue;
             }
             cfgcum[j+1] = cfgcum[j] + arcmcr[kk+k];
             jvrtx = kupmcr[kk+k];
             kpath[j+1] = jvrtx;
             if (k==1) {
                dkmap[nopenk[j]] = norbi - j;
                ksymcum[j+1] = GRPMUL(ksymcum[j],orbsymi[j+1]);
                nopenk[j+1] = nopenk[j] + 1;
                dket[norbi-j] = 4;           /* flag */
             } else {
                ksymcum[j+1] = ksymcum[j];
                nopenk[j+1] = nopenk[j];
                if (kstep[j+1]) dket[norbi-j] = 3;
                else dket[norbi-j] = 0;
             }
             kmax = 2;
             j++;
             successful_ascent=1;
             break;
          }
       }

       if (successful_ascent) {
          if (j==norbi) {
             jcsfbase = ketbase[cfgcum[norbi]];
             extsymk = GRPMUL(ksymcum[norbi], irrep);
             ncsfq1tmp = ncsfcfg[nopenk[norbi]];
             if (!ncsfq1tmp) goto endket;
             ncfgtmp = chksym(nelxb, nelxk[kmcr], extsym, extsymk, &permut); 
             if (!ncfgtmp) goto endket; 

             if (iread_elmt) {
                elmt22sum = ptr_cc_two - cc_two;
                if ((elmt22sum + 6 * ncsfq1tmp * ncsfq2tmp) > bunlen) {
                   bunloc += elmt22sum;
                   (void)fseek(ciftfile, bunloc * sizeof(struct cc_elmt2), SEEK_SET);
                   itmp = elmt22tot - bunloc;
                   if (itmp>bunlen) {
                      itmp = bunlen*sizeof(struct cc_elmt2); 
                   } else {
                      iread_elmt = 0;
                      itmp *= sizeof(struct cc_elmt2); 
                   }
                   if (fread(cc_two, itmp, 1, ciftfile)!=1){
                       (void)fprintf(outfile, "\n\nfread error (%s)\n\n", CIFTNAME);
                       MPI_Abort(MPI_COMM_WORLD, 1);
                   } 
                   ptr_cc_two = cc_two;
                }
             }
             imo = intgrlmo[iorb[0]];
             jmo = intgrlmo[jorb[0]];
             if (jorb[0]<5) jmo = xorbq2[4-jorb[0]];  
             switch (xsp[norbi]) {
                   case 2:
 
                      /* direct integrals    (00|11) */
                      /* exchange integrals  (01|01) */
                      /* imo >= kmo */
                      kmo = intgrlmo[iorb[1]];
                      lmo = intgrlmo[jorb[1]];
                      if (jorb[1]<5) lmo = xorbq2[4-jorb[1]];  
                      switch (nelxb) {
                         case 0: 
                            switch (nelxk[kmcr]) {
                               case 0: 
                                  /* bra: (N-1, 1)E0  (N-2, 2)E0 
                                     ket: (N-3, 3)E0  (N-4, 4)E0  
                                     imo >= kmo   
                                     jmo >= lmo                     */
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsfbase); 
                                  break;

                               case 1:   
                                  /* bra: (N-2, 1)E1 
                                     ket: (N-3, 3)E0 
                                     kmo => E1        ? sum of integrals? */
                                  jcsf_loc = jcsfbase;
                                  for (kmo=exsym[extsymk]; kmo>=sxsym[extsymk]; kmo--) {
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     jcsf_loc += ncsfq1tmp;
                                     ptr_cc_two -= elmt22cnt[jcfg];
                                  }
                                  ptr_cc_two += elmt22cnt[jcfg];
                                  break;
                            }
                            break;

                         case 1: 
                            switch (nelxk[kmcr]) {
                               case 0: 
                                  /* bra: (N-1, 1)E0  (N-2, 2)E0 
                                     ket: (N-3, 2)E1  (N-4, 3)E1  
                                     imo >= kmo   
                                     jmo >= lmo                     */
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsfbase); 
                                  break;

                               case 1:   
                                  /* bra: (N-1, 0)E1  (N-2, 1)E1   
                                     ket: (N-3, 2)E1  (N-4, 3)E1 */
                                  if (extsym==extsymk) {
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else {
                                     jcfg--;
                                  }
                                  break;

                               case 2:
                                  /* bra: (N-2, 0)E2
                                     ket: (N-3, 2)E1
                                     lmo => E1                    */        
                                  kmo = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  /* bra: (N-2, 0)E11
                                     ket: (N-3, 2)E10
                                     kmo => E(1)1      ? sum of integrals */   
                                  jcfgtmp = jcfg+1;
                                  ptr22tmp = ptr_cc_two;
                                  ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg];
                                  kst = GRPMUL(extsym, extsymk);
                                  for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                                     if (kmo==xorbq2[0]) continue;
                                     if (xorbq2[0]>kmo) {
                                        ptr_cc_two = ptr22tmp;
                                        akmo = IOFF(xorbq2[0]) + kmo;
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else { 
                                        ptr_cc_two = ptr22tmp2;
                                        akmo = IOFF(kmo) + xorbq2[0]; 
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                                     }
                                  }
                                  jcfg++; 
                                  ptr_cc_two = ptr22tmp2 + elmt22cnt[jcfg];
                                  break;
                            }
                            break;

                         case 2: 
                            switch (nelxk[kmcr]) {
                               case 0: 
                                  /* bra: (N-1, 1)E0  (N-2, 2)E0 
                                     ket: (N-3, 1)E2  (N-4, 2)E2
                                     imo >= kmo   
                                     jmo >= lmo                     */
                                  getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsfbase); 
                                  break;

                               case 1: 
                                  /* bra: (N-1, 0)E10  (N-2, 1)E10 
                                     ket: (N-3, 1)E20  (N-4, 2)E20 
                                     lmo => E20                    */        
                                  if (orbsym[xorbq2[0]]==extsymk) {
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  break;

                               case 2:
                                  /* bra: (N-2, 0)E2                               
                                     ket: (N-3, 1)E2  (N-4, 2)E2         */
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  /* bra: (N-2, 0)E11
                                     ket: (N-3, 1)E20
                                     kmo => E1                    */        
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  jcfgtmp = jcfg+1;
                                  ptr22tmp = ptr_cc_two;
                                  ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg];
                                  kst = GRPMUL(orbsym[xorbq2[0]], extsymk);
                                  for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                                     if (kmo==xorbq2[0]) continue;
                                     if (xorbq2[0]>kmo) {
                                        ptr_cc_two = ptr22tmp;
                                        akmo = IOFF(xorbq2[0]) + kmo;
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else {
                                        ptr_cc_two = ptr22tmp2;
                                        akmo = IOFF(kmo) + xorbq2[0]; 
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                                     }
                                  }
                                  jcfg++; 
                                  ptr_cc_two = ptr22tmp2 + elmt22cnt[jcfg];
                                  break;
                            }
                            break;

                         case 3: 
                            switch (nelxk[kmcr]) {
                               case 0: 
                                  /* bra: (N-1, 1)E0  (N-2, 2)E0 
                                     ket: (N-3, 1)E11 (N-4, 2)E11
                                     imo >= kmo   
                                     jmo >= lmo                     */
                                  getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsfbase); 
                                  break;

                               case 1:
                                  /* bra: (N-1, 0)E10  (N-2, 1)E10 
                                     ket: (N-3, 1)E11  (N-4, 2)E11     
                                     lmo => E(1)1   *default*        */        
                                  if (extsymk==orbsym[xorbq2[0]]) {
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++; 
                                  /* bra: (N-1, 0)E01  (N-2, 1)E01 
                                     ket: (N-3, 1)E11  (N-4, 2)E11 
                                     lmo => E1(1)                    */     
                                  if (extsymk==orbsym[xorbq2[1]]) {
                                     lmo = xorbq2[0];
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[1]);
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  break;

                               case 2:
                                  /* bra: (N-2, 0)E20
                                     ket: (N-3, 1)E11 
                                     kmo => E20     *default*      */        
                                  kmo = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  /* bra: (N-2, 0)E02
                                     ket: (N-3, 1)E11 
                                     kmo => E02                    */   
                                  jcfg++;
                                  kmo = xorbq2[1];
                                  lmo = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[1]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  /* bra: (N-2, 0)E11
                                     ket: (N-3, 1)E11   (N-4, 2)E11
                                     *symmetry examined before*         */
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  if (extsym==extsymk) {
                                     abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                     jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else {
                                     jcfg--;
                                  }
                                  break;
                            }
                            break;

                         case 4: 
                            switch (nelxk[kmcr]) {
                               case 1:
                                  /* bra: (N-1, 0)E10  (N-2, 1)E10 
                                     ket: (N-3, 0)E21  (N-4, 1)E21     
                                     jmo =>    *default*  
                                     lmo =>    *default*        */        
                                  if (extsymk==orbsym[xorbq2[0]]) {
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  if (extsym==extsymk) {
                                     /* bra: (N-1, 0)E01  (N-2, 1)E01 
                                        ket: (N-3, 0)E21  (N-4, 2)E21 
                                        jmo = lmo => E2(1)                    */ 
                                     jcfg++;  
                                     lmo = xorbq2[0];
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[1]);
                                     getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  }
                                  break;

                               case 2:
                                  /* bra: (N-2, 0)E20                           
                                     ket: (N-3, 0)E21  (N-4, 1)E21  */
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  /* bra: (N-2, 0)E110                           
                                     ket: (N-3, 0)E210  (N-4, 1)E210   */
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[1]])==extsymk) { 
                                     if (xorbq2[0]>xorbq2[1]) {
                                        abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                        ptr_cc_two += elmt22cnt[jcfg+1];
                                     } else { 
                                        ptr_cc_two += elmt22cnt[jcfg];
                                        abmo = IOFF(xorbq2[1]) + xorbq2[0];
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg+1, jcsf_loc); 
                                     }
                                  } else  ptr_cc_two += elmt22cnt[jcfg] + elmt22cnt[jcfg+1];
                                  jcfg++; 
                                  break;
                            }
                            break;

                         case 5: 
                            switch (nelxk[kmcr]) {
                               case 1:
                                  if (extsymk==orbsym[xorbq2[0]]) {
                                     /* bra: (N-1, 0)E100  (N-2, 1)E100 
                                        ket: (N-3, 0)E111  (N-4, 1)E111     
                                        jmo => E(1)11   *default*  jmo > lmo          
                                        lmo => E(1)11   *default*        */        
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++;  
                                  if (extsymk==orbsym[xorbq2[1]]) {
                                     /* bra: (N-1, 0)E010  (N-2, 1)E010 
                                        ket: (N-3, 0)E111  (N-4, 1)E111     
                                        jmo => E1(1)1     jmo > lmo          
                                        lmo => E1(1)1                    */  
                                     jmo = xorbq2[0];          
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[1]);
                                     getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++; 
                                  if (extsymk==orbsym[xorbq2[2]]) {
                                     /* bra: (N-1, 0)E001  (N-2, 1)E001 
                                        ket: (N-3, 0)E111  (N-4, 1)E111     
                                        imo => E11(1)     imo > kmo          
                                        kmo => E11(1)                    */ 
                                     jmo = xorbq2[0];          
                                     lmo = xorbq2[1];          
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[2]);
                                     getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  break;

                               case 2:
                                  /* bra: (N-2, 0)E200                          
                                     ket: (N-3, 0)E111  (N-4, 1)E111 */
                                  kmo = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  jcfg++; 
                                  kmo = xorbq2[1];
                                  jmo = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[1]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  jcfg++; 
                                  kmo = xorbq2[2];
                                  lmo = xorbq2[1];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[2]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  /* bra: (N-2, 0)E1100   
                                     ket: (N-3, 0)E1110  (N-4, 1)E1110    ( a > b > c)
                                     lmo => E(11)1              */        
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[2]]) { 
                                     abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                     jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++; 
                                  /* bra: (N-2, 0)E1010  
                                     ket: (N-3, 0)E1110  (N-4, 1)E1110    ( a > b > c)
                                     lmo => E(1)1(1)               */        
                                  if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[1]]) { 
                                     lmo=xorbq2[1];
                                     abmo = IOFF(xorbq2[0]) + xorbq2[2];
                                     jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++; 
                                  /* bra: (N-2, 0)E0110 
                                     ket: (N-3, 0)E1110  (N-4, 1)E1110    ( a > b > c)
                                     lmo => E1(11)                  */        
                                  if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[0]]) { 
                                     lmo=xorbq2[0];
                                     abmo = IOFF(xorbq2[1]) + xorbq2[2];
                                     jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } else ptr_cc_two += elmt22cnt[jcfg];
                                  break;
                            }
                            break;

                         case 6: 
                            switch (nelxk[kmcr]) {
                               case 2:
                                  /* bra: (N-2, 0)E20    *default*              
                                     ket: (N-4, 0)E22                 */
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  /* bra: (N-2, 0)E02                           
                                     ket: (N-4, 0)E22                 */
                                  ptr_cc_two -= elmt22cnt[jcfg];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[1]);
                                  jmo = xorbq2[0];
                                  lmo = xorbq2[0];
                                  getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  /* bra: (N-2, 0)E11    *default*              
                                     ket: (N-4, 0)E22    (a > b)      */
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[1]])!=extsymk) {
                                      ptr_cc_two += elmt22cnt[jcfg];
                                      break;
                                  }
                                  abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                  jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                  getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;
                            }
                            break;

                         case 7: 
                            switch (nelxk[kmcr]) {
                               case 2:
                                  /* bra: (N-2, 0)E200   *default*              
                                     ket: (N-4, 0)E211                */
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  break;

                               case 3:
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  ptr22tmp = ptr_cc_two;
                                  if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[1]])==extsymk) { 
                                     /* bra: (N-2, 0)E110   *default*              
                                        ket: (N-4, 0)E211                */
                                     if (xorbq2[0]>xorbq2[1]) { 
                                        abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                        jcsf_loc   = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else {
                                        ptr_cc_two += elmt22cnt[jcfg] + elmt22cnt[jcfg+1];
                                        abmo = IOFF(xorbq2[1]) + xorbq2[0]; 
                                        jcsf_loc   = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_twox(imo, jmo, kmo, lmo, jcfg+2, jcsf_loc); 
                                     }
                                  } 
                                  ptr_cc_two = ptr22tmp + elmt22cnt[jcfg];
                                  ptr22tmp = ptr_cc_two;
                                  jcfg++; 
                                  if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[2]])==extsymk) { 
                                     /* bra: (N-2, 0)E101                          
                                        ket: (N-4, 0)E211                */
                                     lmo = xorbq2[1];
                                     if (xorbq2[0]>xorbq2[2]) {
                                        abmo = IOFF(xorbq2[0]) + xorbq2[2];
                                        jcsf_loc   = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else {
                                     /* lmo and jmo exchanged             */
                                        ptr_cc_two += elmt22cnt[jcfg] + elmt22cnt[jcfg+1];
                                        abmo = IOFF(xorbq2[2]) + xorbq2[0]; 
                                        jcsf_loc   = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_twox(imo, lmo, kmo, jmo, jcfg+2, jcsf_loc); 
                                     }
                                  } 
                                  ptr_cc_two = ptr22tmp + elmt22cnt[jcfg] + elmt22cnt[jcfg+1] + elmt22cnt[jcfg+2];
                                  jcfg += 2; 
                                  if (extsym==extsymk) { 
                                     /* bra: (N-2, 0)E011                          
                                        ket: (N-4, 0)E211                */
                                     jcfg++; 
                                     jmo = xorbq2[0];
                                     lmo = xorbq2[0];
                                     bcmo = IOFF(xorbq2[1]) + xorbq2[2];
                                     jcsf_loc   = jcsfbase + lccsf_11[csfloc+bcmo];
                                     getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  } 
                                  break;
                            }
                            break;

                         case 8: 
                            csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                            if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[1]])==extsymk) { 
                               /* bra: (N-2, 0)E1100                         
                                  ket: (N-4, 0)E1111  a > b > c > d   */
                               abmo = IOFF(xorbq2[0]) + xorbq2[1];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                            } else ptr_cc_two += elmt22cnt[jcfg];
                            jcfg++; 
                            if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[2]])==extsymk) { 
                               /* bra: (N-2, 0)E1010                         
                                  ket: (N-4, 0)E1111                  */
                               jmo = xorbq2[1];
                               abmo = IOFF(xorbq2[0]) + xorbq2[2];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                            } else  ptr_cc_two += elmt22cnt[jcfg];
                            jcfg++; 
                            if (GRPMUL(orbsym[xorbq2[0]], orbsym[xorbq2[3]])==extsymk) { 
                               /* bra: (N-2, 0)E1001                         
                                  ket: (N-4, 0)E1111                  */
                               jmo = xorbq2[1];
                               lmo = xorbq2[2];
                               abmo = IOFF(xorbq2[0]) + xorbq2[3];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                            } else  ptr_cc_two += elmt22cnt[jcfg];
                            jcfg++; 
                            if (GRPMUL(orbsym[xorbq2[1]], orbsym[xorbq2[2]])==extsymk) { 
                               /* bra: (N-2, 0)E0110                         
                                  ket: (N-4, 0)E1111                  */
                               jmo = xorbq2[0];
                               lmo = xorbq2[3];
                               abmo = IOFF(xorbq2[1]) + xorbq2[2];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                            } else  ptr_cc_two += elmt22cnt[jcfg];
                            jcfg++; 
                            if (GRPMUL(orbsym[xorbq2[1]], orbsym[xorbq2[3]])==extsymk) { 
                               /* bra: (N-2, 0)E0101                         
                                  ket: (N-4, 0)E1111                  */
                               jmo = xorbq2[0];
                               lmo = xorbq2[2];
                               abmo = IOFF(xorbq2[1]) + xorbq2[3];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                            } else  ptr_cc_two += elmt22cnt[jcfg];
                            jcfg++; 
                            if (GRPMUL(orbsym[xorbq2[2]], orbsym[xorbq2[3]])==extsymk) { 
                               /* bra: (N-2, 0)E0011                         
                                  ket: (N-4, 0)E1111                  */
                               jmo = xorbq2[0];
                               lmo = xorbq2[1];
                               abmo = IOFF(xorbq2[2]) + xorbq2[3];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_twox(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                            } else  ptr_cc_two += elmt22cnt[jcfg];
                            break;

                      } /* end of switch(nelxb) */
                      jcfg++;
                      break; /* end two-particle-two-hole */ 
                        
                   case 1:
                      switch (nelxb) {
                         case 0:
                            /* bra:  (N-2, 2)E0    
                               ket:  (N-3, 3)E0  */
                            getintgrl_one(iorb[0], jorb[0], jcsfbase, dket); 
                            jcfg11++;
                            break;

                         case 1:
                            switch (nelxk[kmcr]) {
                               case 0:
                                  /* bra:  (N-2, 2)E0    
                                     ket:  (N-3, 2)E1     */
                                  intgrlmo[4] = xorbq2[0];
                                  getintgrl_one(iorb[0], jorb[0], jcsfbase, dket); 
                                  jcfg11++;
                                  break;

                               case 1:
                                  /* bra:  (N-2, 1)E10    
                                     ket:  (N-3, 2)E10  symmetry restricted */
                                  if (extsym==extsymk) {
                                     intgrlmo[4] = xorbq2[0];
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                     jcfg11++;
                                  }
                                  /* bra: (N-2, 1)E01 
                                     ket: (N-3, 2)E10 
                                     kmo => E01        
                                     lmo => E10        */  
                                  jcsf_loc = jcsfbase;
                                  lmo = xorbq2[0];
                                  for (kmo=exsym[extsymk]; kmo>=sxsym[extsymk]; kmo--) {
                                     if (kmo==lmo) {
                                        jcsf_loc += ncsfq1tmp;
                                        continue;
                                     }
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     jcsf_loc += ncsfq1tmp;
                                     ptr_cc_two -= elmt22cnt[jcfg];
                                  }
                                  ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++;
                                  break;
                            }
                            break;

                         case 2:
                            switch (nelxk[kmcr]) {
                               case 1:
                                  /* bra:  (N-2, 1)E10    
                                     ket:  (N-3, 1)E20  */
                                  if (extsymk==orbsym[xorbq2[0]]) {
                                     intgrlmo[4] = xorbq2[0];
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                  } else ptr_cc_one += elmt11cnt[jcfg11];
                                  jcfg11++;
                                  /* bra: (N-2, 1)E01 
                                     ket: (N-3, 1)E20   
                                     jmo == lmo => E20                    
                                     no  exchange integrals        */  
                                  jcsf_loc = jcsfbase;
                                  jmo = xorbq2[0];
                                  lmo = xorbq2[0];
                                  for (kmo=exsym[extsymk]; kmo>=sxsym[extsymk]; kmo--) {
                                     if (kmo==lmo) {
                                        jcsf_loc += ncsfq1tmp;
                                        continue;
                                     }
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     jcsf_loc += ncsfq1tmp;
                                     ptr_cc_two -= elmt22cnt[jcfg];
                                  }
                                  ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++;
                                  break;

                               case 2:
                                  /* bra:  (N-3, 1)E20
                                     ket:  (N-2, 0)E20  */
                                  intgrlmo[4] = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                  jcfg11++;
                                  break;
                            }
                            break;

                         case 3:
                            switch (nelxk[kmcr]) {
                               case 1: 
                                  if (extsymk==orbsym[xorbq2[0]]) {
                                     /* bra:  (N-2, 1)E10  *default   
                                        ket:  (N-3, 1)E11           */
                                     intgrlmo[4] = xorbq2[0];
                                     intgrlmo[3] = xorbq2[1];
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[0]);
                                     getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                  } else ptr_cc_one += elmt11cnt[jcfg11];
                                  jcfg11++;
                                  if (extsymk==orbsym[xorbq2[1]]) {
                                     /* bra:  (N-2, 1)E01   ????reverse????      
                                        ket:  (N-3, 1)E11                      */
                                     jorb[0] = 4;   
                                     dket[4] = 0;
                                     dket[3] = 4;
                                     intgrlmo[4] = xorbq2[0];
                                     intgrlmo[3] = xorbq2[1];
                                     jcsf_loc = jcsfbase + ncsfq1tmp*(exsym[extsymk]-xorbq2[1]);
                                     getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                     dket[4] = 4;
                                     dket[3] = 0;
                                  } else ptr_cc_one += elmt11cnt[jcfg11];
                                  jcfg11++;
                                  /* bra: (N-2, 1)E001 
                                     ket: (N-3, 1)E110   
                                     imo kmo => E11                    */ 
                                  jcsf_loc = jcsfbase;
                                  jmo = xorbq2[0];
                                  lmo = xorbq2[1];
                                  for (kmo=exsym[extsymk]; kmo>=sxsym[extsymk]; kmo--) {
                                     if ((kmo==jmo)||(kmo==lmo)) {
                                        jcsf_loc += ncsfq1tmp;
                                        continue;
                                     }
                                     getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     jcsf_loc += ncsfq1tmp;
                                     ptr_cc_two -= elmt22cnt[jcfg];
                                  }
                                  ptr_cc_two += elmt22cnt[jcfg];
                                  jcfg++;
                                  break;

                               case 3:
                                  /* bra:  (N-2, 0)E11   symmetry restricted  
                                     ket:  (N-3, 1)E11                    */
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  if (extsym==extsymk) {
                                     intgrlmo[4] = xorbq2[0];
                                     intgrlmo[3] = xorbq2[1];
                                     abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                     jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                     getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                     jcfg11++;
                                  }
                                  /* bra: (N-2, 0)E101
                                     ket: (N-3, 1)E110 
                                     kmo => E(1)01    *default*         */
                                  jcfgtmp = jcfg+2;
                                  ptr22tmp = ptr_cc_two;
                                  ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg] + elmt22cnt[jcfg+1];
                                  kst = GRPMUL(orbsym[xorbq2[0]], extsymk);
                                  lmo = xorbq2[1];
                                  for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                                     if (kmo==xorbq2[0]) continue;
                                     if (kmo==xorbq2[1]) continue;
                                     if (xorbq2[0]>kmo) {
                                        ptr_cc_two = ptr22tmp;
                                        akmo = IOFF(xorbq2[0]) + kmo;
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else { 
                                        ptr_cc_two = ptr22tmp2;
                                        akmo = IOFF(kmo) + xorbq2[0]; 
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                                     }
                                  }
                                  /* bra: (N-2, 0)E011
                                     ket: (N-3, 1)E110 
                                     lmo => E0(1)1       *symmetry examined before* */
                                  ptr_cc_two = ptr22tmp + elmt22cnt[jcfg];
                                  jcfg++;  
                                  jcfgtmp = jcfg+2;
                                  ptr22tmp = ptr_cc_two;
                                  ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg] + elmt22cnt[jcfg+1];
                                  kst = GRPMUL(orbsym[xorbq2[1]], extsymk);
                                  lmo = xorbq2[0];
                                  for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                                     if (kmo==xorbq2[0]) continue;
                                     if (kmo==xorbq2[1]) continue;
                                     if (xorbq2[1]>kmo) { 
                                        ptr_cc_two = ptr22tmp;
                                        bkmo = IOFF(xorbq2[1]) + kmo;
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+bkmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else {
                                        ptr_cc_two = ptr22tmp2;
                                        bkmo = IOFF(kmo) + xorbq2[1]; 
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+bkmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                                     }
                                  }
                                  ptr_cc_two = ptr22tmp2 +  elmt22cnt[jcfg+2];
                                  jcfg += 3; 
                                  break;
                            }
                            break;

                         case 4:
                            switch (nelxk[kmcr]) {
                               case 2:
                                  /* bra:  (N-2, 0)E20    
                                     ket:  (N-3, 0)E21  */
                                  intgrlmo[4] = xorbq2[0];
                                  intgrlmo[3] = xorbq2[1];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[0]);
                                  getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                  jcfg11++;
                                  /* ket: (N-2, 0)E02
                                     bra: (N-3, 0)E21  
                                     kmo => E2(1)       imo==kmo      */        
                                  kmo = xorbq2[1];
                                  jmo = xorbq2[0];
                                  lmo = xorbq2[0];
                                  jcsf_loc = jcsfbase + ncsfq1tmp*(norbx-xorbq2[1]);
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                  jcfg++;
                                  break;

                               case 3:
                                  csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                                  if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[0]]) {
                                     if (xorbq2[0]>xorbq2[1]) { 
                                        /* bra:  (N-2, 0)E11    
                                           ket:  (N-3, 0)E21  */
                                        abmo = IOFF(xorbq2[0]) + xorbq2[1];
                                        intgrlmo[4] = xorbq2[0];
                                        intgrlmo[3] = xorbq2[1];
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                        jcfg11++;
                                        ptr_cc_one += elmt11cnt[jcfg11];
                                     } else {
                                        /* bra:  (N-2, 0)E11    
                                           ket:  (N-3, 0)E12   */
                                        /* Change, then revert */
                                        ptr_cc_one += elmt11cnt[jcfg11];
                                        jcfg11++;
                                        jorb[0] = 3;
                                        abmo = IOFF(xorbq2[1]) + xorbq2[0];
                                        intgrlmo[4] = xorbq2[1];
                                        intgrlmo[3] = xorbq2[0];
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                                        getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                                     }
                                     jcfg11++;
                                  } else {
                                     ptr_cc_one += elmt11cnt[jcfg11] + elmt11cnt[jcfg11+1];
                                     jcfg11+=2;
                                  }
                                  /* bra: (N-2, 0)E101
                                     ket: (N-3, 1)E210 
                                     lmo => E(1)01                             */
                                  jcfgtmp = jcfg+1;
                                  ptr22tmp = ptr_cc_two;
                                  ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg];
                                  kst = GRPMUL(orbsym[xorbq2[0]], extsymk);
                                  jmo = xorbq2[0];
                                  lmo = xorbq2[1];
                                  for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                                     if (kmo==xorbq2[0]) continue;
                                     if (kmo==xorbq2[1]) continue;
                                     if (xorbq2[0]>kmo) {
                                        ptr_cc_two = ptr22tmp;
                                        akmo = IOFF(xorbq2[0]) + kmo;
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else {
                                        ptr_cc_two = ptr22tmp2;
                                        akmo = IOFF(kmo) + xorbq2[0]; 
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                                     }
                                  }
                                  ptr_cc_two = ptr22tmp2 + elmt22cnt[jcfg+1];
                                  jcfg += 2;  
                                  ptr22tmp = ptr_cc_two;
                                  ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg];
                                  jcfgtmp = jcfg+1;
                                  /* ket: (N-2, 0)E011
                                     bra: (N-3, 1)E210   (jmo==lmo) 
                                     kmo => E0(1)1                             */
                                  kst = GRPMUL(extsym, extsymk);
                                  jmo = xorbq2[0];
                                  lmo = xorbq2[0];
                                  for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                                     if (kmo==xorbq2[0]) continue;
                                     if (kmo==xorbq2[1]) continue;
                                     if (xorbq2[1]>kmo) {
                                        ptr_cc_two = ptr22tmp;
                                        bkmo = IOFF(xorbq2[1]) + kmo;
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+bkmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                                     } else {
                                        ptr_cc_two = ptr22tmp2;
                                        bkmo = IOFF(kmo) + xorbq2[1]; 
                                        jcsf_loc = jcsfbase + lccsf_11[csfloc+bkmo];
                                        getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                                     }
                                  }
                                  ptr_cc_two = ptr22tmp2 + elmt22cnt[jcfg+1];
                                  jcfg += 2; 
                                  break;
                            }
                            break;

                         case 5:
                            csfloc = ((nopenk[norbi]-spin)/2)*nnorbx;
                            intgrlmo[4] = xorbq2[0];
                            intgrlmo[3] = xorbq2[1];
                            intgrlmo[2] = xorbq2[2];

                            if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[2]]) {
                               /* bra:  (N-2, 0)E110   
                                  ket:  (N-3, 0)E111   */
                               abmo = IOFF(xorbq2[0]) + xorbq2[1];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket);
                            } else ptr_cc_one += elmt11cnt[jcfg11];
                            jcfg11++;
                            if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[1]]) {
                               /* bra:  (N-2, 0)E101    
                                  ket:  (N-3, 0)E111    */
                               jorb[0] = 3;
                               dket[3] = 0;
                               dket[2] = 4;
                               abmo = IOFF(xorbq2[0]) + xorbq2[2];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                            } else ptr_cc_one += elmt11cnt[jcfg11];
                            jcfg11++;
                            if (GRPMUL(extsym, extsymk)==orbsym[xorbq2[0]]) {
                               /* bra:  (N-2, 0)E011   
                                  ket:  (N-3, 0)E111    */
                               jorb[0] = 4;
                               dket[4] = 0;
                               dket[3] = 4;
                               dket[2] = 4;
                               abmo = IOFF(xorbq2[1]) + xorbq2[2];
                               jcsf_loc = jcsfbase + lccsf_11[csfloc+abmo];
                               getintgrl_one(iorb[0], jorb[0], jcsf_loc, dket); 
                            } else ptr_cc_one += elmt11cnt[jcfg11];
                            jcfg11++;
                            dket[4] = 4;
                            dket[3] = 4;
                            dket[2] = 0;
                            /* bra: (N-2, 0)E1001    
                               ket: (N-3, 0)E1110   ( a > b > c)
                               kmo => E(1)001       *symmetry examined before* */
                            ptr22tmp = ptr_cc_two;
                            ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg] + elmt22cnt[jcfg+1] + elmt22cnt[jcfg+2];
                            jcfgtmp = jcfg+3;
                            kst = GRPMUL(orbsym[xorbq2[0]], extsymk);
                            jmo=xorbq2[1];
                            lmo=xorbq2[2];
                            for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                               if (kmo==xorbq2[0]) continue;
                               if (kmo==xorbq2[1]) continue;
                               if (kmo==xorbq2[2]) continue;
                               if (xorbq2[0]>kmo) {
                                  ptr_cc_two = ptr22tmp;
                                  akmo = IOFF(xorbq2[0]) + kmo;
                                  jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                               } else {
                                  ptr_cc_two = ptr22tmp2;
                                  akmo = IOFF(kmo) + xorbq2[0]; 
                                  jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                               }
                            }
                            /* bra: (N-2, 0)E0101    
                               ket: (N-3, 0)E1110   ( a > b > c)
                               kmo => E0(1)01       *symmetry examined before* */
                            ptr_cc_two = ptr22tmp + elmt22cnt[jcfg];
                            jcfg++;
                            jcfgtmp = jcfg+3;
                            ptr22tmp = ptr_cc_two;
                            ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg] + elmt22cnt[jcfg+1] + elmt22cnt[jcfg+2];
                            kst = GRPMUL(orbsym[xorbq2[1]], extsymk);
                            jmo=xorbq2[0];
                            for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                               if (kmo==xorbq2[0]) continue;
                               if (kmo==xorbq2[1]) continue;
                               if (kmo==xorbq2[2]) continue;
                               if (xorbq2[1]>kmo) {
                                  ptr_cc_two = ptr22tmp;
                                  akmo = IOFF(xorbq2[1]) + kmo;
                                  jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                               } else {
                                  ptr_cc_two = ptr22tmp2;
                                  akmo = IOFF(kmo) + xorbq2[1]; 
                                  jcsf_loc = jcsfbase + lccsf_11[csfloc+akmo];
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                               }
                            }
                            /* bra: (N-2, 0)E0011    
                               ket: (N-3, 0)E1110   ( a > b > c)
                               kmo => E00(1)1       *symmetry examined before* */
                            ptr_cc_two = ptr22tmp + elmt22cnt[jcfg];
                            jcfg++;
                            jcfgtmp = jcfg+3;
                            ptr22tmp = ptr_cc_two;
                            ptr22tmp2 = ptr_cc_two + elmt22cnt[jcfg] + elmt22cnt[jcfg+1] + elmt22cnt[jcfg+2];
                            lmo=xorbq2[1];
                            kst = GRPMUL(orbsym[xorbq2[2]], extsymk);
                            for (kmo=exsym[kst]; kmo>=sxsym[kst]; kmo--) {
                               if (kmo==xorbq2[0]) continue;
                               if (kmo==xorbq2[1]) continue;
                               if (kmo==xorbq2[2]) continue;
                               if (xorbq2[2]>kmo) {
                                  ptr_cc_two = ptr22tmp;
                                  akmo = IOFF(xorbq2[2]) + kmo;
                                  jcsf_loc   = jcsfbase + lccsf_11[csfloc+akmo];
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfg, jcsf_loc); 
                               } else {
                                  ptr_cc_two = ptr22tmp2;
                                  akmo = IOFF(kmo) + xorbq2[2]; 
                                  jcsf_loc   = jcsfbase + lccsf_11[csfloc+akmo];
                                  getintgrl_two(imo, jmo, kmo, lmo, jcfgtmp, jcsf_loc); 
                               }
                            }
                            ptr_cc_two = ptr22tmp2 + elmt22cnt[jcfg+3];
                            jcfg += 4;
                            break;
                      } /* end of switch(nelxb) */
                      break; /* end one-particle-one-hole */
                        
                }  /* end switch(xsp) */
                     
             endket: /* no viable CSF in this Q1-configuration */
             /* re-initialize */
             j = norbi - 2;            /* by-pass a wasted search step */
             if (j>-1) {               
                jvrtx = kpath[j];
                kmax = kstep[j+1] - 1;   
             } else {
                break;                 /* normal exit for 1-orbital internal space (rare!) */
             }
          } /* if (j==norbi) */
            
       } else {   
          /* try to descend */   
          if (!j) break;               /* normal exit */
          jvrtx = kpath[--j];
          kmax = kstep[j+1] - 1;   
       } /* if (successful_ascent) else ... *for ket* */
    }/* for (;;)  *for ket* */   
 } /* loop over kmcr *evaluation of sigma vector* */

 if (iter) {
    if (fread(hqp, ncsfxstate*sizeof(double), 1, scrfile)!=1) {
         (void)fprintf(outfile, "\n\nfread error (%s)\n\n", SCRNAME);
         MPI_Abort(MPI_COMM_WORLD, 1);
    }
 } else if (mxiter>1) {
    if (fwrite(hqp, ncsfxstate*sizeof(double), 1, scrfile)!=1) {
         (void)fprintf(outfile, "\n\nfwrite error (%s)\n\n", SCRNAME);
         MPI_Abort(MPI_COMM_WORLD, 1);
    }
     
 }

 /* evaluation of (HW)pp */
 for (istate=0; istate<nstates; istate++) {
    hqq[istate] = energs[istate] - hdiag;
    hqq[istate] = 1.0 / hqq[istate];
 }

 for (icsf=0, icsf_k=0; icsf<ncsfq2tmp; icsf++) {
    for (istate=0; istate<nstates; istate++) 
       hscr[istate] = 0.0;
    for (istate=0; istate<nstates; istate++, icsf_k++) {
       for (jstate=0; jstate<nstates; jstate++) {
          hscr[jstate] += hqp[icsf_k] * upp[istate][jstate];
       }
    }
    for (istate=0; istate<nstates; istate++) {
       for (jstate=0; jstate<istate; jstate++) {
          tmp = hscr[istate] * hscr[jstate];
          hwpp[istate][jstate] += tmp * hqq[jstate];
          hwpp[jstate][istate] += tmp * hqq[istate];
          wwpp[istate][jstate] += tmp * hqq[istate] * hqq[jstate];
       }
       tmp = hscr[istate] * hscr[istate] * hqq[istate];
       hwpp[istate][istate] += tmp; 
       wwpp[istate][istate] += tmp * hqq[istate];
    }
 }




 /* evaluation of (HW)pp with inverse matrix (He2e2-Ep)(-1) */
 if (invflag) {

    /* update hqp for iterative evaluation */
    if (iter) {
       for (icsf=0, icsf_k=0; icsf<ncsfq2tmp; icsf++) {
          for (istate=0; istate<nstates; istate++) 
             hscr[istate] = 0.0;
          for (istate=0; istate<nstates; istate++, icsf_k++) {
             for (jstate=0; jstate<nstates; jstate++) {
                hscr[jstate] += hqp[icsf_k] * upp2[istate][jstate];
             }
          }
          icsf_k -= nstates;
          for (istate=0; istate<nstates; istate++, icsf_k++) {
             hqp[icsf_k] = hscr[istate];
          }
       }
    }

    for (jstate=0; jstate<nstates; jstate++) {

       if (ncsfq2tmp>1) {
          for (icsf=0; icsf<ncsfq2tmp; icsf++) { 
             hqq2[icsf][icsf] = hdiagtmp[icsf] - energs[jstate];
          }

          itmp = choleski(hqq2[0], hqq2inv, buff, ncsfq2tmp);
          if (!itmp) {
             /* skip if inverse matrix failed */
             if (prntflag>1)
                (void)fprintf(outfile, "\ninverse matrix failed!");
             break;
          }

       } else {
          hqq2inv[0] = 1.0/(hdiagtmp[0] - energs[jstate]);
       }

       for (icsf=0; icsf<ncsfq2tmp; icsf++) {
          ijcsf = (icsf * (icsf+1))/2;
          for (jcsf=0, jcsf_k=jstate, tmp=0.0; jcsf<=icsf; jcsf++, ijcsf++) { 
             tmp += hqq2inv[ijcsf] * hqp[jcsf_k]; 
             jcsf_k += nstates;
          }
          ijcsf += icsf;
          for (; jcsf<ncsfq2tmp; jcsf++) { 
             tmp += hqq2inv[ijcsf] * hqp[jcsf_k]; 
             jcsf_k += nstates;
             ijcsf += jcsf + 1;
          }
          icsf_k = icsf * nstates;
          hqwp[icsf_k+jstate] = tmp;
          for (istate=0; istate<jstate; istate++, icsf_k++) {
             hwpp2[istate][jstate] -= hqp[icsf_k] * tmp;
             wwpp2[jstate][istate] += hqwp[icsf_k] * tmp;
          }
          wwpp2[jstate][jstate] += tmp * tmp;
          for (; istate<nstates; istate++, icsf_k++) {
             hwpp2[istate][jstate] -= hqp[icsf_k] * tmp;
          }
       }

    }
 }

#ifdef DEBUG
 for (icsf=0; icsf<ncsfq2tmp; icsf++, ncsfq2sum++) {
    fprintf(outfile, "Hq2q2[%d] = %16.12lf", ncsfq2sum, hdiagtmp[icsf]);
 }
 (void)fprintf(outfile, "\nncsfq2sum =%d,",ncsfq2sum);
 if (ncsfq2tmp==1) 
    fprintf(outfile, "\nhdiag=%16.12lf", hdiag);
#endif

}

void getintgrl_one(int iorb, int jorb, int jcsf_loc, int dket[]) 
{

 extern int nstates, jcfg11, intgrlmo[]; 
 extern double *civec, *hqp, oneintk[]; 
 extern FILE *outfile;
 extern struct cc_elmt1 *ptr_cc_one;

 int imo, jmo, kmo, mij, mik, mkj, mkk, mijkk, mikkj;
 int i, k, icsf_k, jcsf_k, istate, mix, mjx, ijxx, icnt, ntmp;
 double intone, tmp, val;

 imo = intgrlmo[iorb];
 jmo = intgrlmo[jorb];
 mix = IOFF(imo);
 mjx = IOFF(jmo);

 if (iorb > jorb) {  
    /* particle always internal orbital 
       hole might be external orbital     */
    mij = mix + jmo;
    ijxx = IOFF(mij);
    intone = oneint[mij];
    for (k=2,icnt=0; k<jorb; k++) {
       if (!dket[k]) continue;
       kmo = intgrlmo[k];
       mkj = mjx + kmo;
       mik = mix + kmo; 
       mikkj = IOFF(mik) + mkj;
       tmp = twoint[mikkj];
       mkk = IOFF(kmo) + kmo;
       mijkk = ijxx + mkk;
       if (dket[k]==3) intone += 2 * twoint[mijkk] - tmp;
       else if (dket[k]%3) {
          intone += twoint[mijkk] - tmp*0.5;
          oneintk[icnt] = tmp;
          icnt++;
       }   
    }
    kmo = jmo;
    mkk = IOFF(kmo) + kmo;
    mijkk = ijxx + mkk;
    if (dket[jorb]) intone += twoint[mijkk];   
    for (k=jorb+1; k<iorb; k++) {
       if (!dket[k]) continue;
       kmo = intgrlmo[k];
       if (kmo>jmo) mkj = IOFF(kmo) + jmo;
       else mkj = mjx + kmo;
       mik = mix + kmo; 
       mikkj = IOFF(mik) + mkj;
       tmp = twoint[mikkj];   
       mkk = IOFF(kmo) + kmo;
       mijkk = ijxx + mkk;
       if (dket[k]==3) {
          intone += 2 * twoint[mijkk] - tmp;
       } else if (dket[k]%3) {
          intone += twoint[mijkk];
          oneintk[icnt] = tmp;
          icnt++;
       }   
    }
    kmo = imo;
    mkk = IOFF(kmo) + kmo;
    mijkk = IOFF(mkk) + mij;
    if (dket[iorb]==3) intone += twoint[mijkk];
    for (k=iorb+1; k<=norbi; k++) {
       if (!dket[k]) continue;
       kmo = intgrlmo[k];
       mkj = IOFF(kmo) + jmo;
       mik = IOFF(kmo) + imo; 
       mikkj = IOFF(mik) + mkj;
       tmp = twoint[mikkj];   
       mkk = IOFF(kmo) + kmo;
       mijkk = IOFF(mkk) + mij;
       if (dket[k]==3) intone += 2 * twoint[mijkk] - tmp;
       else if (dket[k]%3) {
          intone += twoint[mijkk] - tmp*0.5;
          oneintk[icnt] = tmp;
          icnt++;
       }   
    }

 } else {
    /* particle and hole always internal orbitals */
    mij = mjx + imo;
    ijxx = IOFF(mij);
    intone = oneint[mij];
    for (k=2,icnt=0; k<iorb; k++) {
       if (!dket[k]) continue;
       kmo = intgrlmo[k];
       mik = mix + kmo; 
       mkj = mjx + kmo;
       mikkj = IOFF(mkj) + mik;
       tmp = twoint[mikkj];
       mkk = IOFF(kmo) + kmo;
       mijkk = ijxx + mkk;
       if (dket[k]==3) intone += 2 * twoint[mijkk] - tmp;
       else if (dket[k]%3) {
          intone += twoint[mijkk] - tmp*0.5;
          oneintk[icnt] = tmp;
          icnt++;
       }   
    }
    kmo = imo;
    mkk = IOFF(kmo) + kmo;
    mijkk = ijxx + mkk;
    if (dket[iorb]==3) intone += twoint[mijkk];
    for (k=iorb+1; k<jorb; k++) {
       if (!dket[k]) continue;
       kmo = intgrlmo[k];
       mkj = mjx + kmo;
       mik = IOFF(kmo) + imo;
       mikkj = IOFF(mkj) + mik;
       tmp = twoint[mikkj];   
       mkk = IOFF(kmo) + kmo;
       mijkk = ijxx + mkk;
       if (dket[k]==3) {
          intone += 2 * twoint[mijkk] - tmp;
       } else if (dket[k]%3) {
          intone += twoint[mijkk];
          oneintk[icnt] = tmp;
          icnt++;
       }   
    }
    kmo = jmo;
    mkk = IOFF(kmo) + kmo;
    mijkk = IOFF(mkk) + mij;
    if (dket[jorb]) intone += twoint[mijkk];  
    for (k=jorb+1; k<=norbi; k++) {
       if (!dket[k]) continue;
       kmo = intgrlmo[k];
       mkj = IOFF(kmo) + jmo;
       mik = IOFF(kmo) + imo; 
       mikkj = IOFF(mkj) + mik;
       tmp = twoint[mikkj];   
       mkk = IOFF(kmo) + kmo;
       mijkk = IOFF(mkk) + mij;
       if (dket[k]==3) intone += 2 * twoint[mijkk] - tmp;
       else if (dket[k]%3) {
          intone += twoint[mijkk] - tmp*0.5;
          oneintk[icnt] = tmp;
          icnt++;
       }   
    }
 }   

 jcsf_loc *= nstates;
 ntmp = elmt11cnt[jcfg11];
 for (i=0; i<ntmp; i++) { 
    val = (ptr_cc_one->elmtone) * intone;
    for (k=0; k<icnt; k++) {
       val += (ptr_cc_one->acum[k]) * oneintk[k]; 
    }
    icsf_k = ptr_cc_one->icsf; 
    jcsf_k = (ptr_cc_one->jcsf) + jcsf_loc; 
    for (istate=0; istate<nstates; istate++, icsf_k++, jcsf_k++) {
       hqp[icsf_k] += val * civec[jcsf_k];
    }
    ptr_cc_one++;
 }

}

void getintgrl_two(int imo, int jmo, int kmo, int lmo, int jcfgq1, int jcsf_loc) 
{

 /* jcsfq1:  position of coupling coefficients 
    jcsf_loc:  position of Q1-CSF                    */

 /* The sparseness of coupling coefficients could be explored
    by storing only non-zero elements with indeces.  
    {iwalk[], jwalk[], jcsfq1_loc[], val_drct[], val_exch[] }     */

 extern int nstates, *elmt22cnt;
 extern double *civec, *hqp;
 extern struct cc_elmt2 *ptr_cc_two;
#ifdef DEBUG
 extern FILE *outfile;
#endif

 int mdrct, mexch, i, ntmp, icsf_k, jcsf_k, istate;
 double val, intexch, intdrct;

 mexch = 0;
 if (imo > jmo) { /* iorb[0] (=>imo) is an internal orbital */
    if (imo > kmo) {
       if (jmo > kmo) {
          if (jmo == lmo) {
             mdrct = IOFF(IOFF(imo)+jmo) + IOFF(jmo) + kmo;
          } else {
             mexch = IOFF(IOFF(imo)+lmo) + IOFF(jmo) + kmo;
             if (kmo > lmo) {
                mdrct = IOFF(IOFF(imo)+jmo) + IOFF(kmo) + lmo;
             } else {
                mdrct = IOFF(IOFF(imo)+jmo) + IOFF(lmo) + kmo;
             }
          }
       } else {
          if (jmo == lmo) {
             mdrct = IOFF(IOFF(imo)+jmo) + IOFF(kmo) + jmo;
          } else {
             mexch = IOFF(IOFF(imo)+lmo) + IOFF(kmo) + jmo;
             if (kmo>lmo) {
                mdrct = IOFF(IOFF(imo)+jmo) + IOFF(kmo) + lmo;
             } else {
                mdrct = IOFF(IOFF(imo)+jmo) + IOFF(lmo) + kmo;
             }
          }
       } 
    } else {
       if (jmo > lmo) {
          mdrct = IOFF(IOFF(imo)+jmo) + IOFF(imo) + lmo;
       } else {
          mdrct = IOFF(IOFF(imo)+lmo) + IOFF(imo) + jmo;
       }
    }
 } else {  /* imo < jmo, jmo >= lmo */
    if (jmo > lmo) {
       if (imo > lmo) {
          if (imo > kmo) {
             mexch = IOFF(IOFF(jmo)+kmo) + IOFF(imo) + lmo;
             if (kmo > lmo) {
                mdrct = IOFF(IOFF(jmo)+imo) + IOFF(kmo) + lmo;
             } else {
                mdrct = IOFF(IOFF(jmo)+imo) + IOFF(lmo) + kmo;
             }
          } else {
             mdrct = IOFF(IOFF(jmo)+imo) + IOFF(imo) + lmo;
          }
       } else {
          if (imo > kmo) {
             mexch = IOFF(IOFF(jmo)+kmo) + IOFF(lmo) + imo;
             mdrct = IOFF(IOFF(jmo)+imo) + IOFF(lmo) + kmo;
          } else {
             mdrct = IOFF(IOFF(jmo)+imo) + IOFF(lmo) + imo;
          }
       } 
    } else {
       if (imo > kmo) {
          mdrct = IOFF(IOFF(jmo)+imo) + IOFF(jmo) + kmo;
       } else {
          mdrct = IOFF(IOFF(jmo)+imo) + IOFF(jmo) + imo;
       }
    } 
 }

 jcsf_loc *= nstates;
 ntmp = elmt22cnt[jcfgq1];
 if (mexch) { 
    intdrct = twoint[mdrct];
    intexch = twoint[mexch];
    for (i=0; i<ntmp; i++) { 
       icsf_k = ptr_cc_two->icsf; 
       jcsf_k = (ptr_cc_two->jcsf) + jcsf_loc; 
       val = (ptr_cc_two->drct_cc)*intdrct + (ptr_cc_two->exch_cc)*intexch;
       for (istate=0; istate<nstates; istate++, icsf_k++, jcsf_k++) {
          hqp[icsf_k] += val * civec[jcsf_k];
       }
       ptr_cc_two++;
    }
 } else {
    intdrct = twoint[mdrct];
    for (i=0; i<ntmp; i++) { 
       icsf_k = ptr_cc_two->icsf; 
       jcsf_k = (ptr_cc_two->jcsf) + jcsf_loc; 
       val = (ptr_cc_two->drct_cc)*intdrct;
       for (istate=0; istate<nstates; istate++, icsf_k++, jcsf_k++) {
          hqp[icsf_k] += val * civec[jcsf_k];
       }
       ptr_cc_two++;
    }
 }

}

void getintgrl_twox(int imo, int jmo, int kmo, int lmo, int jcfgq1, int jcsf_loc) 
{

 /* imo, kmo > jmo, lmo
    jcfgq1:  position of coupling coefficients 
    jcsf_loc:  position of Q1-CSF                      
    The sparseness of coupling coefficients could be explored
    by storing only non-zero elements with indeces.  
    {iwalk[], jwalk[], jcsfq1_loc[], val_drct[], val_exch[] }     */

 extern int nstates, *elmt22cnt;
 extern double *civec, *hqp;
 extern struct cc_elmt2 *ptr_cc_two;

 int mdrct, mexch, i, ntmp, icsf_k, jcsf_k, istate;
 double val, intexch, intdrct;

 mexch = 0;
 if (imo > kmo) {
    if (jmo == lmo) {
       mdrct = IOFF(IOFF(imo)+jmo) + IOFF(kmo) + jmo;
    } else {
       mexch = IOFF(IOFF(imo)+lmo) + IOFF(kmo) + jmo;
       mdrct = IOFF(IOFF(imo)+jmo) + IOFF(kmo) + lmo;
    }
 } else {
    if (jmo > lmo) {
       mdrct = IOFF(IOFF(imo)+jmo) + IOFF(imo) + lmo;
    } else {
       mdrct = IOFF(IOFF(imo)+lmo) + IOFF(imo) + jmo;
    }
 }

 jcsf_loc *= nstates;
 ntmp = elmt22cnt[jcfgq1];
 if (mexch) { 
    intdrct = twoint[mdrct];
    intexch = twoint[mexch];
    for (i=0; i<ntmp; i++) { 
       icsf_k = ptr_cc_two->icsf; 
       jcsf_k = (ptr_cc_two->jcsf) + jcsf_loc; 
       val = (ptr_cc_two->drct_cc)*intdrct + (ptr_cc_two->exch_cc)*intexch;
       for (istate=0; istate<nstates; istate++, icsf_k++, jcsf_k++) {
          hqp[icsf_k] += val * civec[jcsf_k];
       }
       ptr_cc_two++;
    }
 } else {
    intdrct = twoint[mdrct];
    for (i=0; i<ntmp; i++) { 
       icsf_k = ptr_cc_two->icsf; 
       jcsf_k = (ptr_cc_two->jcsf) + jcsf_loc; 
       val = (ptr_cc_two->drct_cc)*intdrct;
       for (istate=0; istate<nstates; istate++, icsf_k++, jcsf_k++) {
          hqp[icsf_k] += val * civec[jcsf_k];
       }
       ptr_cc_two++;
    }
 }

}
/* calculate the diagonal elements */
/* BW-type                         */
/* energy contribution from frozen core orbitals only
   considered partially as folded in oneint[]                  */

void calchdiag_bw(int nstep[], int dbra[], int dbmap[], int nelxb, int nopenb, 
       double intk, int nxorb)
{
 /* calculate diagonal hamiltonian matrix elements for Q2-configurations */
 extern int invflag, norb, norbi, norbi0, ncsfq2tmp, ftpt[], xorbq2[];
 extern double *twoint, *oneint, exchint[], *exchcc; 
 extern double hdiag, *hdiagtmp, *hdiagtmp2, **hqq2;
 extern FILE *outfile;
  
 int i, j, ij, k, jg, nopenltr;
 int jmo, maj, mab, mjj, mxj, maajj, majja, mabba, maabb, maaaa;
 int dket[NINTMX+1];
 double xint, tmp, tmp2;
 static int a, b, c, d, maa, mbb, mcc, mdd; 
 
 void nopnoh(int nopen, int dbra[], int dmap[], int dket[]);

 switch (nelxb) {
    case 0:
       /*  (N-3, 3)E0  
           (N-4, 4)E0           */ 
       xint = 0.0;
       break;

    case 1:
       /*  (N-3, 2)E1  
           (N-4, 3)E1           */ 
       a = xorbq2[0];
       maa = (a*(a+1))/2;
       xint = oneint[maa];
       ij = ((nopenb-1)*(nopenb-2))/2;
       for (j=0; j<norbi0; j++) {
          jmo = norb - j;
          maj = IOFF(jmo) + a;
          majja = (maj*(maj+1))/2;
          mjj = (jmo*(jmo+1))/2;
          maajj = IOFF(mjj) + maa;
          xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
          if (nstep[j+1]==1) { /* open-shell exchange integrals */
             exchint[ij++] = twoint[majja];
          }
       }
       break;
                  
    case 2:
       /*  (N-3, 1)E2  
           (N-4, 2)E2           */ 
       a = xorbq2[0];
       maa = (a*(a+1))/2;
       xint = 2 * oneint[maa];
       maaaa = (maa*(maa+1))/2;
       xint += twoint[maaaa];
       for (j=0; j<norbi0; j++) {
          jmo = norb - j;
          maj = IOFF(jmo) + a;
          majja = (maj*(maj+1))/2;
          mjj = (jmo*(jmo+1))/2;
          maajj = IOFF(mjj) + maa;
          xint += nstep[j+1]*(2*twoint[maajj]-twoint[majja]);
       }
       break;

    case 3:
       /*  (N-3, 1)E11 
           (N-4, 2)E11          */ 
       if (nxorb == 1) {
          a = xorbq2[0];
          maa = (a*(a+1))/2;
          xint = oneint[maa];
          ij = ((nopenb-2)*(nopenb-3))/2;
          for (j=0; j<norbi0; j++) {
             jmo = norb - j;
             maj = IOFF(jmo) + a;
             majja = (maj*(maj+1))/2;
             mjj = (jmo*(jmo+1))/2;
             maajj = IOFF(mjj) + maa;
             xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
             if (nstep[j+1]==1) { /* open-shell exchange integrals */
                exchint[ij++] = twoint[majja];
             }
          }
          hdiag = intk + xint;
          return;
       } else {
          a = xorbq2[0];
          b = xorbq2[1];
          maa = (a*(a+1))/2;
          mbb = (b*(b+1))/2;
          xint = oneint[mbb];
          ij = ((nopenb-1)*(nopenb-2))/2;
          for (j=0; j<norbi0; j++) {
             jmo = norb - j;
             maj = IOFF(jmo) + b;
             majja = (maj*(maj+1))/2;
             mjj = (jmo*(jmo+1))/2;
             maajj = IOFF(mjj) + mbb;
             xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
             if (nstep[j+1]==1) { /* open-shell exchange integrals */
                exchint[ij++] = twoint[majja];
             }
          }
          mab = IOFF(a) + b;
          mabba = (mab*(mab+1))/2;
          maabb = IOFF(maa) + mbb;
          xint += twoint[maabb]-0.5*twoint[mabba];
          exchint[ij] = twoint[mabba];
       }
       break;

    case 4:
       /*  (N-3, 0)E21  
           (N-4, 1)E21           */ 
       if (nxorb == 1) {
          b = xorbq2[1];
          mbb = (b*(b+1))/2;
          xint = oneint[mbb];
          ij = ((nopenb-1)*(nopenb-2))/2;
          for (j=0; j<norbi0; j++) {
             jmo = norb - j;
             maj = IOFF(jmo) + b;
             majja = (maj*(maj+1))/2;
             mjj = (jmo*(jmo+1))/2;
             maajj = IOFF(mjj) + mbb;
             xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
             if (nstep[j+1]==1) { /* open-shell exchange integrals */
                exchint[ij++] = twoint[majja];
             }
          }
       } else {
          a = xorbq2[0];
          b = xorbq2[1];
          maa = (a*(a+1))/2;
          xint = 2 * oneint[maa];
          maaaa = (maa*(maa+1))/2;
          xint += twoint[maaaa];
          for (j=0; j<norbi0; j++) {
             jmo = norb - j;
             mxj = (jmo*(jmo-1))/2;
             maj = mxj + a;
             majja = (maj*(maj+1))/2;
             mjj = mxj + jmo;
             maajj = (mjj*(mjj-1))/2 + maa;
             xint += nstep[j+1]*(2*twoint[maajj]-twoint[majja]);
          }
          if (a>b) {
             mab = IOFF(a) + b;
             maabb = IOFF(maa) + mbb;
          } else {
             mab = IOFF(b) + a;
             maabb = IOFF(mbb) + maa;
          }
          mabba = (mab*(mab+1))/2;
          xint += 2*twoint[maabb]-twoint[mabba];

          for (k=0; k<ncsfq2tmp; k++) 
             hdiagtmp[k] = hdiagtmp2[k] + xint;

          hdiag = intk + xint;
          return;
       }
       break;

    case 5:
       /*  (N-3, 0)E111
           (N-4, 1)E111         */ 
       switch (nxorb) {
          case 3:
             a = xorbq2[0];
             b = xorbq2[1];
             c = xorbq2[2];
             mcc = (c*(c+1))/2;
             xint = oneint[mcc];
             ij = ((nopenb-1)*(nopenb-2))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + c;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mcc;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             mab = IOFF(a) + c;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(maa) + mcc;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij++] = twoint[mabba];
             mab = IOFF(b) + c;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(mbb) + mcc;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij] = twoint[mabba];
             break;

          case 2:
             a = xorbq2[0];
             b = xorbq2[1];
             mbb = (b*(b+1))/2;
             xint = oneint[mbb];
             ij = ((nopenb-2)*(nopenb-3))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + b;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mbb;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             mab = IOFF(a) + b;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(maa) + mbb;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij] = twoint[mabba];
             hdiag = intk + xint;
             return;
             break;

          case 1:
             a = xorbq2[0];
             maa = (a*(a+1))/2;
             xint = oneint[maa];
             ij = ((nopenb-3)*(nopenb-4))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + a;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + maa;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             hdiag = intk + xint;
             return;
             break;

       }
       break;

    case 6:
       /*  (N-4, 0)E22            */ 
       if (nxorb == 1) {
          a = xorbq2[0];
          maa = (a*(a+1))/2;
          xint = 2*oneint[maa];
          maaaa = (maa*(maa+1))/2;
          xint += twoint[maaaa];
          for (j=0; j<norbi0; j++) {
             jmo = norb - j;
             maj = IOFF(jmo) + a;
             majja = (maj*(maj+1))/2;
             mjj = (jmo*(jmo+1))/2;
             maajj = IOFF(mjj) + maa;
             xint += nstep[j+1]*(2*twoint[maajj]-twoint[majja]);
          }
       } else {
          a = xorbq2[0];
          b = xorbq2[1];
          mbb = (b*(b+1))/2;
          xint = 2*oneint[mbb];
          maaaa = (mbb*(mbb+1))/2;
          xint += twoint[maaaa];
          for (j=0; j<norbi0; j++) {
             jmo = norb - j;
             maj = IOFF(jmo) + b;
             majja = (maj*(maj+1))/2;
             mjj = (jmo*(jmo+1))/2;
             maajj = IOFF(mjj) + mbb;
             xint += nstep[j+1]*(2*twoint[maajj]-twoint[majja]);
          }
          mab = IOFF(a) + b;
          mabba = (mab*(mab+1))/2;
          maabb = IOFF(maa) + mbb;
          xint += 4*twoint[maabb]-2*twoint[mabba];

          for (k=0; k<ncsfq2tmp; k++) 
             hdiagtmp[k] = hdiagtmp2[k] + xint;

          hdiag = intk + xint;
          return;
       }
       break;
            
    case 7:
       /*  (N-4, 0)E211          */ 
       switch (nxorb) {
          case 3:
             a = xorbq2[0];
             b = xorbq2[1];
             c = xorbq2[2];
             maa = (a*(a+1))/2;
             xint = 2 * oneint[maa];
             maaaa = (maa*(maa+1))/2;
             xint += twoint[maaaa];
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                mxj = (jmo*(jmo-1))/2;
                maj = mxj + a;
                majja = (maj*(maj+1))/2;
                mjj = mxj + jmo;
                maajj = (mjj*(mjj-1))/2 + maa;
                xint += nstep[j+1]*(2*twoint[maajj]-twoint[majja]);
             }
             if (a>b) {
                mab = IOFF(a) + b;
                maabb = IOFF(maa) + mbb;
             } else {
                mab = IOFF(b) + a;
                maabb = IOFF(mbb) + maa;
             }
             mabba = (mab*(mab+1))/2;
             xint += 2*twoint[maabb]-twoint[mabba];
             if (a>c) {
                mab = IOFF(a) + c;
                maabb = IOFF(maa) + mcc;
             } else {
                mab = IOFF(c) + a;
                maabb = IOFF(mcc) + maa;
             }
             mabba = (mab*(mab+1))/2;
             xint += 2*twoint[maabb]-twoint[mabba];

             for (k=0; k<ncsfq2tmp; k++) 
                hdiagtmp[k] = hdiagtmp2[k] + xint;

             hdiag = intk + xint;
             return;
             break;

          case 1:
             b = xorbq2[1];
             mbb = (b*(b+1))/2;
             xint = oneint[mbb];
             ij = ((nopenb-2)*(nopenb-3))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + b;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mbb;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             hdiag = intk + xint;
             return;
             break;

          case 2:
             b = xorbq2[1];
             c = xorbq2[2];
             mcc = (c*(c+1))/2;
             xint = oneint[mcc];
             ij = ((nopenb-1)*(nopenb-2))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + c;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mcc;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             mab = IOFF(b) + c;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(mbb) + mcc;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij] = twoint[mabba];
             break;
       }
       break;

    case 8:
       /*  (N-4, 1)E1111        */ 
       switch (nxorb) {
          case 4:
             a = xorbq2[0];
             b = xorbq2[1];
             c = xorbq2[2];
             d = xorbq2[3];
             mdd = (d*(d+1))/2;
             xint = oneint[mdd];
             ij = ((nopenb-1)*(nopenb-2))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + d;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mdd;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             mab = IOFF(a) + d;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(maa) + mdd;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij++] = twoint[mabba];
             mab = IOFF(b) + d;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(mbb) + mdd;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij++] = twoint[mabba];
             mab = IOFF(c) + d;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(mcc) + mdd;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij] = twoint[mabba];
             break;

          case 3:
             a = xorbq2[0];
             b = xorbq2[1];
             c = xorbq2[2];
             mcc = (c*(c+1))/2;
             xint = oneint[mcc];
             ij = ((nopenb-2)*(nopenb-3))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + c;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mcc;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             mab = IOFF(a) + c;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(maa) + mcc;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij++] = twoint[mabba];
             mab = IOFF(b) + c;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(mbb) + mcc;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij] = twoint[mabba];
             hdiag = intk + xint;
             return;
             break;

          case 2:
             a = xorbq2[0];
             b = xorbq2[1];
             mbb = (b*(b+1))/2;
             xint = oneint[mbb];
             ij = ((nopenb-3)*(nopenb-4))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + b;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + mbb;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             mab = IOFF(a) + b;
             mabba = (mab*(mab+1))/2;
             maabb = IOFF(maa) + mbb;
             xint += twoint[maabb]-0.5*twoint[mabba];
             exchint[ij] = twoint[mabba];
             hdiag = intk + xint;
             return;
             break;

          case 1:
             a = xorbq2[0];
             maa = (a*(a+1))/2;
             xint = oneint[maa];
             ij = ((nopenb-4)*(nopenb-5))/2;
             for (j=0; j<norbi0; j++) {
                jmo = norb - j;
                maj = IOFF(jmo) + a;
                majja = (maj*(maj+1))/2;
                mjj = (jmo*(jmo+1))/2;
                maajj = IOFF(mjj) + maa;
                xint += nstep[j+1]*(twoint[maajj]-0.5*twoint[majja]);
                if (nstep[j+1]==1) { /* open-shell exchange integrals */
                   exchint[ij++] = twoint[majja];
                }
             }
             hdiag = intk + xint;
             return;
             break;
       }

       break;

 }   /* end switch */ 

 hdiag = intk + xint;
 nopenltr = (nopenb*(nopenb-1))/2;
 if (nopenltr) {
    if (invflag) {
       /* initialization of arrays */
       for (i=0; i<ncsfq2tmp; i++) {
          for (j=0; j<i; j++) {
             hqq2[i][j] = 0.0;
          }
          hdiagtmp[i] = hdiag;
       } 

       /* offdiagonal Hamiltonian elements */
       for (i=0; i<=norbi; i++) dket[i] = dbra[i];
       nopnoh(nopenb, dbra, dbmap, dket);

       jg = ftpt[nopenb];
       for (k=0, tmp2=0.0; k<ncsfq2tmp; k++) {
          for (j=0,tmp=0.0; j<nopenltr; j++) {
             tmp += exchint[j] * exchcc[jg++];
          }
          hdiagtmp[k] += tmp;
          tmp2 += tmp;
       }   
       hdiag += tmp2/ncsfq2tmp;
    } else {
       jg = ftpt[nopenb];
       for (j=0; j<nopenltr; j++) {
          /* use averaged exchcc[] */
          hdiag += exchint[j] * exchcc[jg++]; 
       }
    }
 } else {
    hdiagtmp[0] = hdiag;
 }

}

int mkq1mcr(int imcr, int jmcrloc, int jstart, int jend, int nsp, int **kdn, 
             int **ptr_kup, int **ptr_arc, int **ptr_ketbase, int jcnt)
{

 extern int my_rank;
 extern int prntflag, ngrps;
 extern int nvrtx[], *mcfg, mcfgln[];
 extern FILE *mcfgfile, *outfile;

 int i, jmcr, len, xsp, ipt, jpt;

 (void)fseek(mcfgfile, jmcrloc, SEEK_SET);
 for (jmcr=jstart; jmcr<jend; jmcr++) {  
    len = 9*nvrtx[jmcr]*sizeof(int);
    if (fread(kdn[jcnt], len, 1, mcfgfile)!=1) {
       (void)fprintf(outfile, "\n\nfread error (%s)\n\n", MCFGNAME);
       MPI_Abort(MPI_COMM_WORLD, 1);
    }
    xsp = 0;
    for (i=0; i<ngrps; i++) {
       ipt = i + imcr*ngrps;
       jpt = i + jmcr*ngrps;
       if (mcfg[ipt]>mcfg[jpt]) xsp += mcfg[ipt] - mcfg[jpt];
    }
    if (xsp<nsp) {
      if (prntflag>=2 && my_rank==0) {
        (void)fprintf(outfile, 
          "macroconfigurations #%d and #%d interact\n", imcr, jmcr);
        (void)fflush(outfile);
      }
      ptr_kup[jcnt] = kdn[jcnt] + 3*nvrtx[jmcr];
      ptr_arc[jcnt] = kdn[jcnt] + 6*nvrtx[jmcr];
      ptr_ketbase[jcnt] = gbcsf + mcfgln[jmcr];
      jcnt++;
    }
  }
  
  return jcnt;

}

void preptwoptwoh(int iorb[], int jorb[], int nopenb, int dbra[], int dbmap[],
  int nopenk, int dket[], int dkmap[])
{

 extern int jcfg, ovrlap;

 int iint, jint, k;

 void twoptwoh(int iint, int jint, int nopenb, int dbra[], int dbmap[],
  int nopenk, int dket[], int dkmap[],
  double (*segment[])(int, int, int, int, double *));
 double (*segment[NINTMX+2])(int, int, int, int, double *);
 
 double wail(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wair(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wasl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wasr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbil(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbir(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbpsl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbpsr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbsl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wbsr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wcp(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wcpp(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdill(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdirr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdilsr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdirsl(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdsll(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdsll2(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdsrr(int dbra, int dket, int delspn, int ketspn, double *val1);
 double wdsrr2(int dbra, int dket, int delspn, int ketspn, double *val1);

 if (iorb[0] > jorb[0]) {
    iint = iorb[0];
    if (iorb[0] > iorb[1]) {
       segment[iorb[0]] = wail; 
       if (jorb[0] > iorb[1]) {
          for (k=iorb[0]-1; k>jorb[0]; k--) segment[k] = wcp;
          if (jorb[0] > jorb[1]) {
             ovrlap = 0;
             segment[jorb[0]] = wbir; 
             if (jorb[1] > iorb[1]) {
                /* Paldus & Boyle c1 & c3 */ 
                for (k=jorb[0]-1; k>jorb[1]; k--) segment[k] = wcpp; 
                segment[jorb[1]] = wbpsl; 
                for (k=jorb[1]-1; k>iorb[1]; k--) segment[k] = wcp; 
                segment[iorb[1]] = wasr; 
                jint = iorb[1];
             } else {
                /* Paldus & Boyle d1 & c5 */ 
                for (k=jorb[0]-1; k>iorb[1]; k--) segment[k] = wcpp;
                segment[iorb[1]] = wbpsr;
                for (k=iorb[1]-1; k>jorb[1]; k--) segment[k] = wcp;
                segment[jorb[1]] = wasl;
                jint = jorb[1];
             }
          } else {
             /* Paldus & Boyle c2 */ 
             ovrlap = -1;               /* adjacent */
             segment[jorb[0]] = wdirsl;
             for (k=jorb[0]-1; k>iorb[1]; k--) segment[k] = wcp;
             segment[iorb[1]] = wasr;
             jint = iorb[1];
          }

       } else {
          for (k=iorb[0]-1; k>iorb[1]; k--) segment[k] = wcp;
          segment[iorb[1]] = wbil;
          for (k=iorb[1]-1; k>jorb[0]; k--) segment[k] = wcpp;
          if (jorb[0] > jorb[1]) {
             /* Paldus & Boyle d3 & d5 */
             ovrlap = 1;
             segment[jorb[0]] = wbpsl; 
             for (k=jorb[0]-1; k>jorb[1]; k--) segment[k] = wcp;
             segment[jorb[1]] = wasl;

          } else {
             /* Paldus & Boyle d4 */
             ovrlap = -1;               /* only one integral */
             segment[jorb[0]] = wdsll;
          }
          jint = jorb[1];
       } 

    } else {
       ovrlap = -1;                  /* only one integral */
       segment[iorb[0]] = wdill;
       for (k=iorb[0]-1; k>jorb[0]; k--) segment[k] = wcpp;
       if (jorb[0] > jorb[1]) {
          /* Paldus & Boyle d6 */
          segment[jorb[0]] = wbsl; 
          for (k=jorb[0]-1; k>jorb[1]; k--) segment[k] = wcp;
          segment[jorb[1]] = wasl;

       } else {
          /* Paldus & Boyle d7 */
          segment[jorb[0]] = wdsll2; 
          /* note: factor of 2 accounts for absence of permutational
          symmetry present in other coupling coefficients */
       }
       jint = jorb[1];
    }

 } else {

    iint = jorb[0];
    if (jorb[0] > jorb[1]) {
       segment[jorb[0]] = wair;
       if (iorb[0] > jorb[1]) {
          for (k=jorb[0]-1; k>iorb[0]; k--) segment[k] = wcp;
          if (iorb[0] > iorb[1]) {
             ovrlap = 0;
             segment[iorb[0]] = wbil;
             if (iorb[1] > jorb[1]) {
                /* Paldus & Boyle b1 & b3 */
                for (k=iorb[0]-1; k>iorb[1]; k--) segment[k] = wcpp;
                segment[iorb[1]] = wbpsr;
                for (k=iorb[1]-1; k>jorb[1]; k--) segment[k] = wcp;
                segment[jorb[1]] = wasl;
                jint = jorb[1];
             } else {
                /* Paldus & Boyle a1 & b5 */ 
                for (k=iorb[0]-1; k>jorb[1]; k--) segment[k] = wcpp;
                segment[jorb[1]] = wbpsl;
                for (k=jorb[1]-1; k>iorb[1]; k--) segment[k] = wcp;
                segment[iorb[1]] = wasr;
                jint = iorb[1];
             }
          } else {
             /* Paldus & Boyle b2 */ 
             ovrlap = -1;               /* adjacent */
             segment[iorb[0]] = wdilsr;
             for (k=iorb[0]-1; k>jorb[1]; k--) segment[k] = wcp;
             segment[jorb[1]] = wasl;
             jint = jorb[1];
          }

       } else {
          for (k=jorb[0]-1; k>jorb[1]; k--) segment[k] = wcp;
          segment[jorb[1]] = wbir;
          for (k=jorb[1]-1; k>iorb[0]; k--) segment[k] = wcpp;
          if (iorb[0] > iorb[1]) {
             /* Paldus & Boyle a3 & a5 */
             ovrlap = 1;
             segment[iorb[0]] = wbpsr; 
             for (k=iorb[0]-1; k>iorb[1]; k--) segment[k] = wcp;
             segment[iorb[1]] = wasr;
          } else {
             /* Paldus & Boyle a4 */
             ovrlap = -1;                 /* only one integral */
             segment[iorb[0]] = wdsrr; 
          }
          jint = iorb[1];
       } 

    } else {
       ovrlap = -1;                      /* only one integral */
       segment[jorb[0]] = wdirr;
       for (k=jorb[0]-1; k>iorb[0]; k--) segment[k] = wcpp;
       if (iorb[0] > iorb[1]) {
          /* Paldus & Boyle a6 */
          segment[iorb[0]] = wbsr; 
          for (k=iorb[0]-1; k>iorb[1]; k--) segment[k] = wcp;
          segment[iorb[1]] = wasr;
       } else {
          /* Paldus & Boyle a7 */
          segment[iorb[0]] = wdsrr2;
          /* note: factor of 2 accounts for absence of permutational
          symmetry present in other coupling coefficients */
       }
       jint = iorb[1];
    }
 }
 twoptwoh(iint, jint, nopenb, dbra, dbmap, nopenk, dket, dkmap, segment); 
 jcfg++;

}

void nopnoh(int nopen, int dbra[], int dmap[], int dket[])
/* calculate offdiagonal hamiltonian matrix elements */ 
{
 int braspn[MXKMAX+1], bstep[MXKMAX+1], kstep[MXKMAX+1], ketspn[MXKMAX+1];
 int ketwt[MXKMAX+1];
 int bralev, ketlev, ketlevmx, tmp, icsf, i, j, kvrtx, incflag, nelx;
 int pint, qint, itmp, qintmx, accpt, headmin;
 int jtmp, ipt, jpt, ij, ispnptr, jspnptr;
 double val, val1, sum;
 extern int norbi, spin, arcsp[], jspdn[], jspup[], spnhead[], delspn[], bospn[], dospn[];
 extern int nspndif[], norbx, snglptr[];
 extern FILE *outfile;
 extern double headcum[][MXKMAX+1], exchint[], **hqq2;
 /* extern double *twoint */
  
 double wcpp(int dbra, int dket, int delspn, int ketspn, double *val);
 double wdirl(int dbra, int dket, int delspn, int ketspn, double *val);
 double wdsrl(int dbra, int dket, int delspn, int ketspn, double *val);

 nelx = 0;
 
 bralev = nopen+nelx;
 if (spin>bralev) return;
 ketlevmx = nopen+nelx;
 ispnptr = snglptr[bralev];
 jspnptr = snglptr[ketlevmx];
 
 braspn[bralev] = spin;
 itmp = (nopen) ? dmap[0] : 0;
 for (i=norbi; i>=itmp; i--) bospn[i] = spin;
 pint = itmp;
 icsf = 0;
 incflag = 1;
 dospn[norbi] = spin;
 nspndif[ketlevmx+1] = 0;
   
 for (;;) {
      
   if (bralev==nelx) {
      /* internal bra csf generated */
          
      /* >>>>>> generate *interacting* ket csfs in ascending order <<<<<< */
      incflag = 1;
         
      ketlev = ketlevmx;
      headmin = ketlevmx + 1;
      ketspn[ketlev] = spin;
      ketwt[ketlev] = 0;
      kvrtx = spnhead[ketlev];
      qintmx = qint = (ketlev>nelx) ? dmap[nopen-(ketlev-nelx)] : 0;
      for (i=norbi; i>qint; i--) dospn[i-1] = dospn[i];
         
      for (;;) {   
         if (ketlev==nelx) {
            
            if (icsf==ketwt[0]) break;
               
            /* internal ket csf generated */
#ifdef DEBUG   
            if (qint) {
               (void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
               MPI_Abort(MPI_COMM_WORLD, 1);
            }
#endif
            
            /* >>>>> complete internal loops <<<<< */
            for (i=nelx+1,sum=0.; i<ketlevmx; i++) {
               if (ketspn[i-1]==braspn[i-1]) {
                  ipt = nopen-i+nelx;
                  itmp = dmap[nopen-(i-nelx)];
                  delspn[itmp] = dospn[itmp] - bospn[itmp];
                  (void)wdsrl(dbra[itmp], dket[itmp], delspn[itmp], dospn[itmp], &val1);
                  itmp += norbx;
                  for (j=(i+1>headmin)?i+1:headmin; j<=ketlevmx; j++) {
                     val = headcum[i+1][j] * val1;
                     jpt = nopen-j+nelx;
                     jtmp = dmap[nopen-(j-nelx)] + norbx;
                     ij = (ipt*(ipt-1))/2 + jpt;
                     sum += val * exchint[ij];
                  }
               } else break;
            }
            /*??? lexical order ???*/
            if (icsf > ketwt[0]) 
               hqq2[icsf][ketwt[0]] += sum;
            else 
               hqq2[ketwt[0]][icsf] += sum;
               
            /* >>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<< */ 
                
            /* re-initialize */
            ketlev++;
            if (ketlev>ketlevmx) break;   /* normal exit for no internal open shells */
            kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
            incflag = 0;
            qint = dmap[nopen-(ketlev-nelx)];
#ifdef DEBUG 
            if (qint>qintmx) {
               (void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
               MPI_Abort(MPI_COMM_WORLD, 1);
            }
#endif
            
         } else if (incflag) { 
            /* try to descend (ketlev => ketlev-1) with spin increase, if possible */
            tmp = ketspn[ketlev] + 1;
            if (tmp < ketlev) {
               kstep[ketlev] = -1;
               dket[qint] = 2;
               dospn[qint-1] = dospn[qint] + 1;
               ketlev--;
               ketspn[ketlev] = tmp;
               ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx];
               kvrtx = jspdn[2*kvrtx];
               
            } else {
               kstep[ketlev] = 1;
               dket[qint] = 1;
               dospn[qint-1] = dospn[qint] - 1;
               ketlev--;
               ketspn[ketlev] = tmp - 2;
               ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
               kvrtx = jspdn[2*kvrtx+1];
            }
               
            itmp = (ketlev>nelx) ? dmap[nopen-(ketlev-nelx)] : 0;
            for (i=qint-1; i>itmp; i--) dospn[i-1] = dospn[i];
            
            accpt = 0;
            delspn[qint] = dospn[qint] - bospn[qint];
            if (nspndif[ketlev+2]) {
               if (headmin<ketlev+2) headmin = ketlev + 2;
               tmp = headmin;
            } else {
               accpt++;
               headmin = ketlev + 1;
               tmp = headmin + 1;
               (void)wdirl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
               headcum[ketlev+1][ketlev+1] = val;
            }
            (void)wcpp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
            for (j=tmp; j<=ketlevmx; j++) {
               headcum[ketlev+1][j] = headcum[ketlev+2][j] * val;
               if (headcum[ketlev+1][j]!=0.0) accpt++;
            }
            
            if (dket[qint]==dbra[qint]) nspndif[ketlev+1] = nspndif[ketlev+2];
            else nspndif[ketlev+1] = nspndif[ketlev+2] + 1;
               
            /* reject descent */
            if (ketwt[ketlev]>icsf) { 
               incflag = 0;
               ketlev++;
               kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
            } else if (!ketlev) {
               if (ketwt[ketlev]==icsf) {
                  incflag = 0;
                  ketlev++;
                  kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
               } else qint = itmp; 
            } else      
               qint = itmp; 
            
         } else if (kstep[ketlev]==-1) {
            /* try to descend with spin decrease */
            tmp = ketspn[ketlev] - 1;
            if (tmp >= 0) {
               kstep[ketlev] = 1;
               dket[qint] = 1;
               dospn[qint-1] = dospn[qint] - 1;
               ketlev--;
               ketspn[ketlev] = tmp;
               ketwt[ketlev] = ketwt[ketlev+1] + arcsp[2*kvrtx+1];
               kvrtx = jspdn[2*kvrtx+1];
               incflag = 1;
               
               itmp = (ketlev>nelx) ? dmap[nopen-(ketlev-nelx)] : 0;
               for (i=qint-1; i>itmp; i--) dospn[i-1] = dospn[i];
               
               /* >>>>> beginning of level evaluation <<<<< */   
               accpt = 0;
               delspn[qint] = dospn[qint] - bospn[qint];
               if (nspndif[ketlev+2]) {
                  if (headmin<ketlev+2) headmin = ketlev + 2;
                  tmp = headmin;
               } else {
                  accpt++;
                  headmin = ketlev + 1;
                  tmp = headmin + 1;
                  (void)wdirl(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
                  headcum[ketlev+1][ketlev+1] = val;
               }
               (void)wcpp(dbra[qint], dket[qint], delspn[qint], dospn[qint], &val);
               for (j=tmp; j<=ketlevmx; j++) {
                  headcum[ketlev+1][j] = headcum[ketlev+2][j] * val;
                  if (headcum[ketlev+1][j]!=0.0) accpt++;
               }
                  
               if (dket[qint]==dbra[qint]) nspndif[ketlev+1] = nspndif[ketlev+2];
               else nspndif[ketlev+1] = nspndif[ketlev+2] + 1;
               /* >>>>>>>> end of level evaluation <<<<<<<< */
               
               /* reject descent */
               if (ketwt[ketlev]>icsf) { 
                  incflag = 0;
                  ketlev++;
                  kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
               } else if (!ketlev) {
                  if (ketwt[ketlev]==icsf) {
                     incflag = 0;
                     ketlev++;
                     kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
                  } else qint = itmp;     
               } else qint = itmp; 
               
            } else {
               /* try to ascend */
               ketlev++;
               if (ketlev>ketlevmx) break;   /* normal exit */
               kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
               qint = dmap[nopen-(ketlev-nelx)];
#ifdef DEBUG
               if (qint>qintmx) {
                  (void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
                  MPI_Abort(MPI_COMM_WORLD, 1);
               }
#endif
            }
               
         } else {
            /* try to ascend */
            ketlev++;
            if (ketlev>ketlevmx) break;   /* normal exit */
            kvrtx = jspup[2*kvrtx+(kstep[ketlev]+1)/2];
            qint = dmap[nopen-(ketlev-nelx)];
#ifdef DEBUG
            if (qint>qintmx) {
               (void)fprintf(outfile, "algorithm error: qint=%d\n", qint);
               MPI_Abort(MPI_COMM_WORLD, 1);
            }
#endif
         }
      }   
      /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
            
      icsf++;
            
      /* re-initialize */
      bralev++;
      if (bralev>(nopen+nelx)) break;   /* normal exit for no internal open shells */
      incflag = 0; 
      pint = dmap[nopen-(bralev-nelx)];
         
   } else if (incflag) { 
      /* try to descend (bralev => bralev-1) with spin increase, if possible */
      tmp = braspn[bralev] + 1;
      if (tmp < bralev) {
         bstep[bralev] = -1;
         dbra[pint] = 2;
         bralev--;
         braspn[bralev] = tmp;
         
      } else {
         bstep[bralev] = 1;
         dbra[pint] = 1;
         bralev--;
         braspn[bralev] = tmp - 2;
      }
         
      itmp = (bralev>nelx) ? dmap[nopen-(bralev-nelx)] : 0;
      for (i=pint-1; i>=itmp; i--) bospn[i] = braspn[bralev];
      pint = itmp; 
         
   } else if (bstep[bralev]==-1) {
      /* try to descend with spin decrease */
      tmp = braspn[bralev] - 1;
      if (tmp >= 0) {
         bstep[bralev] = 1;
         dbra[pint] = 1;
         bralev--;
         braspn[bralev] = tmp;
         incflag = 1;
         
         itmp = (bralev>nelx) ? dmap[nopen-(bralev-nelx)] : 0;
         for (i=pint-1; i>=itmp; i--) bospn[i] = braspn[bralev];
         pint = itmp; 
         
      } else {
         /* try to ascend */
         bralev++;
         if (bralev>(nopen+nelx)) break;   /* normal exit */
         pint = dmap[nopen-(bralev-nelx)];
      }
      
   } else {
      /* try to ascend */
      bralev++;
      if (bralev>(nopen+nelx)) break;   /* normal exit */
      pint = dmap[nopen-(bralev-nelx)];
   }   
 }
 
}

