// this version uses mpi version of scratch files. it also expects imcr's to
// be seperated to different processors dynaimcally


/***********************************************************************
 *                                                                     *
 * main program for MRCISD(TQ) based on CFGCI                          *
 * March 9, 2008 by WJ                                                 *
 *                                                                     *
 **********************************************************************/

#include "cfgci_tq.h" 

int comm_sz, my_rank; //MPI related variables
FILE *outfile, *infofile, *scrfile;
int prntflag, mxdim, nmcrt, nstates, iter, mxiter=1, invflag=0;
long ncsft, ncsfq1, ncsfm;
double **upp, **hwpp, **wwpp, *energs, *civec;
double **upp2, **hwpp2, **wwpp2, *energs2;

int main(void)
{
 extern int comm_sz, my_rank;
 extern int prntflag, nstates, iter, mxiter;
 extern long ncsft, ncsfq1, ncsfm;
 extern double **upp, **hwpp, **wwpp, **hwpp2, **wwpp2, *energs, *energs2;
 extern FILE *outfile, *infofile, *scrfile;

 int ichk, i, j, k, ij, nroots, nnstates, nstates2;
 int conv=0;
 long ncsfs[3], totalmem;
 double repnuc, cpu_time, runtime;
 double diff, etol=1.e-12, temp, sum, sum2, mem;
 double **upptemp, **upptemp2, **uvecs, **vpp; 
 double **wpp, **wpp2, **wppinv, **wppinv2, **wqp;
 double **hwpp_sum, **wwpp_sum, **wwpp_sum2, **hwpp_sum2;
 double *heff, *eigval, *old_energs, *old_energs2, *energs0;
 char inname[LENNAME], outname[LENNAME], *cptr;
 char *keyword_ptr, line[LENLINE], str[LENLINE]; 
 char *asctime(const struct tm *tmptr);
 struct tm *tmptr;
 time_t lt;
 size_t memtmp;
 FILE *infile, *cifile;
  
 void initci(long totalmem);
 void mkhw(void);
 void multiply(double **a, double **b, double **c, int n);
 void multiplyplus(double **a, double **b, double **c, int n);

 
 MPI_Init(NULL, NULL);
 MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
 MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

 if (comm_sz == 1) {
  (void)printf("\nError: This is the parallel version of tqcorr.exe, please use more than one processor or use the serial code\n" );
  MPI_Abort(MPI_COMM_WORLD, 1);
 }

 cpu_time = (double)clock()/CLOCKS_PER_SEC;
 if ((infofile = fopen(INFONAME, "r+b")) == NULL) {
  (void)printf("\nError opening %s\n", INFONAME);
  MPI_Abort(MPI_COMM_WORLD, 1);
 }
  
 getname(infofile, inname, LENNAME);
 
 (void)strcpy(outname, inname);
 cptr = strstr(outname, ".inp");
 (void)strcpy(cptr, ".log");
  
 if ((outfile = fopen(outname, "a")) == NULL) {
  (void)printf("\nError opening %s\n", outname);
  MPI_Abort(MPI_COMM_WORLD, 1);
 }

 if (my_rank == 0) { 
    
   /* data from infofile */
   prntflag = getint(infofile, 36);    /* set default printing level */
   nstates = getint(infofile, 31);
   repnuc  = getreal(infofile, 1);
   mem     = getreal(infofile, 5);
   getarray(infofile, 60, ncsfs, 3 * sizeof(long int)); /* for all space */
   ncsft =  ncsfs[0];
   ncsfq1 = ncsfs[1];
   ncsfm  = ncsfs[2];
   totalmem = (long) (1.e9 * mem);
    
   if ((infile = fopen(inname, "r")) == NULL) {
    (void)printf("\nError opening input file\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
   }
   (void)strcpy(line, KEY);
   str[0] = '\0';
   ichk = locate(line, infile);
   if (ichk) {
     do {
       (void) fgets(line, sizeof(line), infile);
       
       keyword_ptr = strstr(line, "NSTATES=");
       if (keyword_ptr != NULL) {
         keyword_ptr += 8;
         i = sscanf(keyword_ptr,"%d", &nroots); 
       }
       
       keyword_ptr = strstr(line, "MXITER=");
       if (keyword_ptr != NULL) {
         keyword_ptr += 7;
         i = sscanf(keyword_ptr,"%d", &mxiter);  
       }
       
       keyword_ptr = strstr(line, "INVERSE=");
       if (keyword_ptr != NULL) {
         keyword_ptr += 8;
         i = sscanf(keyword_ptr,"%d", &invflag);
       }
       
       keyword_ptr = strstr(line, "ETOL=");
       if (keyword_ptr != NULL) {
         keyword_ptr += 5;
         i = sscanf(keyword_ptr,"%lf", &etol);
       }
       
       keyword_ptr = strstr(line, "PRINT=");
       if (keyword_ptr != NULL) {
         keyword_ptr += 6;
         i = sscanf(keyword_ptr,"%d", &prntflag);
       }
       
     } while ((keyword_ptr=strstr(line,"#END")) == NULL);
   }
   (void)fclose(infile);

   if (prntflag) {
      (void)fprintf(outfile,
        "\n********** Configuration Interaction program (TQ)");
      (void)fprintf(outfile,"**********\n");
      lt = time(NULL);
      tmptr = localtime(&lt);
      (void)fprintf(outfile, "%s", asctime(tmptr));
      if (!ichk) {
        (void)fprintf(outfile, "\ninput not found => using defaults\n");
        (void)fflush(outfile);
      }
      (void)fprintf(outfile, "\ntotal number of processors (including master proc): %d\n", comm_sz);
      (void)fprintf(outfile, "\ndimension of model space: %ld", ncsfm);
      (void)fprintf(outfile, "\ndimension of M+Q1 space: %ld", ncsfq1+ncsfm);
      (void)fprintf(outfile, "\ndimension of M+Q1+Q2 space: %ld", ncsft);
      (void)fprintf(outfile, "\nprint flag = %d", prntflag);
      (void)fprintf(outfile, "\nnroots = %d\n\n", nstates);
      (void)fflush(outfile);
   }

   if (nstates<1) {
      (void)fprintf(outfile, "\nNo converged root from MRCISD!\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
   }
   mxdim = nstates;
   nnstates = (nstates*(nstates+1))/2;
 }


 MPI_Bcast(&prntflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&nstates, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&nnstates, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&mxiter, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&invflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
 MPI_Bcast(&ncsfq1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
 MPI_Bcast(&totalmem, 1, MPI_LONG, 0, MPI_COMM_WORLD);

 if (mxiter>1) {
   if ((scrfile = mpisopen(SCRNAME, "w+b", infofile, my_rank)) == NULL) {
     (void)printf("\nProcessor %d: Error opening %s\n", my_rank, SCRNAME);
     MPI_Abort(MPI_COMM_WORLD, 1);
   }
 }
 
 energs  = (double *)malloc(nstates * sizeof(double));
 getarray(infofile, 27, energs, nstates * sizeof(double)); /* for all space */
 
 if (my_rank == 0) {
    energs0 = (double *)malloc(nstates * sizeof(double));
    old_energs = (double *)malloc(nstates * sizeof(double));
    eigval = (double *)malloc(nstates * sizeof(double));
    for (i=0; i<nstates; i++)  {
      old_energs[i] = energs[i];
      energs0[i] = energs[i];
    }
   
    if (prntflag) {
      (void)fprintf(outfile, "\n MRCISD energies:\n");
      (void) fprintf(outfile, "\n\nState       Energy      ");
      (void) fprintf(outfile,   "\n (I)        /a.u./      ");
      (void) fprintf(outfile,   "\n.....  ..................");
      for (i=0; i<nstates; i++) 
        (void) fprintf(outfile,"\n  %d     %+.12lf",
          i+1, old_energs[i]+repnuc);
      (void)fprintf(outfile, "\n\n");
      (void) fflush(outfile);
    }
 }   

 civec  = (double *)malloc(ncsfq1 * nstates * sizeof(double));
 if (!civec) {
   (void)fprintf(outfile, "\n\nmemory allocation error for civec in processor %d\n\n", my_rank);
   MPI_Abort(MPI_COMM_WORLD, 1);
 }

 if (my_rank == 0) {
    cifile = sopen(CINAME, "rb", infofile);
    if (cifile == NULL) {
      (void)fprintf(outfile, "\nError opening file: %s\n", CINAME);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
   
    (void)fseek(cifile, ncsfm*nstates*sizeof(double), SEEK_SET);
    if (fread(civec, ncsfq1*nstates*sizeof(double), 1, cifile)!=1) {
      (void)fprintf(outfile, "\n\nfread error (%s)\n\n", CINAME);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
   
 }


 MPI_Bcast(civec, ncsfq1*nstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);

 initci(totalmem);  

 hwpp    = dmatrix0(nstates, nstates); 
 upp     = dmatrix0(nstates, nstates); 
 wwpp    = (double **)malloc(nstates * sizeof(double *));
 memtmp = nnstates * sizeof(double);
 memtmp += 4 * nstates * sizeof(double);
 if (sizeof(double) > sizeof(int)) memtmp += nstates * sizeof(double);
 else memtmp += nstates * sizeof(int);
 wwpp[0] = (double *)malloc(memtmp);
 for (i=1; i<nstates; i++) {
   wwpp[i] = wwpp[i-1] + i;
 }
 for (i=0; i<nstates; i++) {
   upp[i][i] = 1.0;
 }

 /* prepare the arrays for inverse matrix (He2e2-Ep)(-1)  */
 if (invflag) {
   hwpp2    = dmatrix0(nstates, nstates); 
   upp2     = dmatrix0(nstates, nstates); 
   wwpp2    = (double **)malloc(nstates * sizeof(double *));
   wwpp2[0] = (double *)malloc(memtmp);
   for (i=1; i<nstates; i++) {
     wwpp2[i] = wwpp2[i-1] + i;
   }
   energs2  = (double *)malloc(nstates * sizeof(double));
   for (i=0; i<nstates; i++) {
     energs2[i] = energs[i];
     upp2[i][i] = 1.0;
   }
 }

    hwpp_sum= dmatrix0(nstates, nstates);
    wwpp_sum= (double **)malloc(nstates * sizeof(double *));
    wwpp_sum[0] = (double *)malloc(memtmp);
    for (i=1; i<nstates; i++) {
      wwpp_sum[i] = wwpp_sum[i-1] + i;
    }
    if (invflag) {
       hwpp_sum2= dmatrix0(nstates, nstates);
       wwpp_sum2= (double **)malloc(nstates * sizeof(double *));
       wwpp_sum2[0] = (double *)malloc(memtmp);
       for (i=1; i<nstates; i++) {
         wwpp_sum2[i] = wwpp_sum2[i-1] + i;
       }
    }

 if (my_rank == 0) {

    wpp     = dmatrix0(nstates, nstates); 
    wppinv  = dmatrix0(nstates, nstates); 
    wqp     = dmatrix0(nstates, nstates); 
    upptemp = dmatrix0(nstates, nstates); 
    uvecs   = dmatrix0(nstates, nstates); 
    vpp     = dmatrix0(nstates, nstates); 
    heff    = (double *)malloc(memtmp);
    for (i=0; i<nstates; i++) {
      wppinv[i][i] = 1.0;
      upptemp[i][i] = 1.0;
    }
   
    /* prepare the arrays for inverse matrix (He2e2-Ep)(-1)  */
    if (invflag) {
      upptemp2 = dmatrix0(nstates, nstates); 
      wpp2     = dmatrix0(nstates, nstates); 
      wppinv2  = dmatrix0(nstates, nstates); 
      old_energs2 = (double *)malloc(nstates * sizeof(double));
      for (i=0; i<nstates; i++) {
        old_energs2[i] = energs[i];
        wppinv2[i][i] = 1.0;
        upptemp2[i][i] = 1.0;
      }
    }

 }

 nstates2 = nstates*nstates;

 /* start TQ iteration */
 for (iter=0; iter<mxiter && (!conv); iter++) { 
   for (i=0; i<nstates2; i++) {
     hwpp[0][i] = 0.0;
     hwpp_sum[0][i] = 0.0;
   }
   for (i=0; i<nnstates; i++) {
     wwpp[0][i] = 0.0;
     wwpp_sum[0][i] = 0.0;
   }

   if (invflag) {
     for (i=0; i<nstates2; i++) {
       hwpp2[0][i] = 0.0;
       hwpp_sum2[0][i] = 0.0;
     }
     for (i=0; i<nnstates; i++) {
       wwpp2[0][i] = 0.0;
       wwpp_sum2[0][i] = 0.0;
     }
   }

   MPI_Bcast(&iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(energs, nstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   for (i = 0; i < nstates; i++) {
      MPI_Bcast(upp[i], nstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
   if (invflag) {
      MPI_Bcast(energs2, nstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for (i = 0; i < nstates; i++) {
         MPI_Bcast(upp2[i], nstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
   }


   MPI_Barrier(MPI_COMM_WORLD);

   mkhw(); 

   for (i = 0; i < nstates; i++) {
      MPI_Reduce(hwpp[i], hwpp_sum[i], nstates, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(wwpp[i], wwpp_sum[i], i+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   }

   if (invflag) {
      for (i = 0; i < nstates; i++) {
         MPI_Reduce(hwpp2[i], hwpp_sum2[i], nstates, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         MPI_Reduce(wwpp2[i], wwpp_sum2[i], i+1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      }
   }
 
   MPI_Barrier(MPI_COMM_WORLD);

   if (my_rank == 0) { //copy data from sum to regular so the following code stay the saem
   
      for (i = 0; i < nstates; i++) {
         for (j = 0; j <= i; j++) {
            wwpp[i][j] = wwpp_sum[i][j];
         }
         for (j = 0; j < nstates; j++) {
            hwpp[i][j] = hwpp_sum[i][j];
         }
      }
      if (invflag) {
         for (i = 0; i < nstates; i++) {
            for (j = 0; j <= i; j++) {
               wwpp2[i][j] = wwpp_sum2[i][j];
            }
            for (j = 0; j < nstates; j++) {
               hwpp2[i][j] = hwpp_sum2[i][j];
            }
         }
      }
   }


   if (my_rank == 0) { //all operations in the current iteration can be done by the master.
      if (prntflag) {
        (void)fprintf(outfile, "\n******* iteration #%d ********* \n", iter+1);
        (void)fprintf(outfile, "\n (HW)pp matrix for MRCISD(TQ):\n");
        for (i=0; i<nstates; i++) {
          for (j=0; j<nstates; j++) {
            (void)fprintf(outfile, " (HW)pp[%d][%d]=%.12f", i, j, hwpp[i][j]);
          }
          (void)fprintf(outfile, "\n");
        }
        for (i=0; i<nstates; i++) {
          (void)fprintf(outfile, " Lower[%d][%d]=%.12f\n", 
           i, i, old_energs[i]+hwpp[i][i]+repnuc);
        }
   
        (void)fprintf(outfile, "\n (WW)pp matrix for MRCISD(TQ):\n");
        for (i=0; i<nstates; i++) {
          for (j=0; j<=i; j++) {
            (void)fprintf(outfile, " WW[%d][%d]=%.12f", i, j, wwpp[i][j]);
          }
          (void)fprintf(outfile, "\n");
        }
        for (i=0; i<nstates; i++) {
          (void)fprintf(outfile, " Upper[%d][%d]=%.12f\n", i, i, 
                old_energs[i]+hwpp[i][i]/(1.+wwpp[i][i])+repnuc);
        }
      }
      (void)fflush(outfile);
      
      /* =============================================================== *
       * Update  (HW)pp                                                  *
       *         (HW)pp = Wpp(-1)*(HW)pp                                 *
       * scratch use of vpp[][]                                              *
       * =============================================================== */
   
       multiply(wppinv, hwpp, vpp, nstates);
   
      /* =============================================================== *
       * update matrices Wpp and Wpp(-1)                                 *
       * =============================================================== */
   
      for (i=0; i<nstates; i++) eigval[i] = 0.0;
      if (nstates==1) {
        uvecs[0][0] = 1.0;
        eigval[0] = sqrt(1.-wwpp[0][0]);
      } else {
        for (i=0; i<nstates; i++) {
          for (j=0; j<i; j++) wwpp[i][j] *= -1.;
          wwpp[i][i] = 1.-wwpp[i][i];
        }
        jacobi2dm(wwpp[0], eigval, uvecs, nstates, 0);
        for (i=0; i<nstates;  i++) eigval[i] = sqrt(eigval[i]);
      }
   
      for (i=0; i<nstates; i++) {
        for (j=0; j<=i; j++) {
          sum = 0.0;
          sum2 = 0.0;
          for (k=0; k<nstates; k++) {
            temp = uvecs[i][k]*uvecs[j][k];
            sum += temp*eigval[k];
            sum2 += temp/eigval[k];
          }
          wpp[i][j] = sum;
          wppinv[i][j] = sum2;  
          if(i>j) {
            wpp[j][i] = sum;
            wppinv[j][i] = sum2; 
          }
        }
      } 
   
      /* =============================================================== *
       * Construct Heff                                                  *
       * =============================================================== */
   
      /* scratch use of upp[][] to get (hpw)pp */
      for (i=0; i<nstates; i++) {
        for (j=0; j<=i; j++) {
          sum = 0.0;
          for (k=0; k<nstates; k++) {
            temp = upptemp[k][i]*upptemp[k][j];
            sum += temp*energs0[k];
          }
          upp[i][j] = sum;
          if(i>j) {
            upp[j][i] = sum;
          }
        }
      }

      /* add (hpw)pp to (hqw)pp, scratch use of vpp[][] */
      multiplyplus(upp, wpp, vpp, nstates);
   
      multiply(wppinv, vpp, hwpp, nstates);
   
      for (i=0, ij=0; i<nstates; i++) {
        for (j=0; j<i; j++, ij++) {
          heff[ij] = 0.5 * (hwpp[i][j] + hwpp[j][i]);
        }
        heff[ij++] = hwpp[i][i];
      }
   
      if (prntflag) {
        (void)fprintf(outfile, "\n (Heff)pp matrix for MRCISD(TQ):\n");
        for (i=0,ij=0; i<nstates; i++) {
          for (j=0; j<=i; j++,ij++) {
            (void)fprintf(outfile, " (Heff)pp[%d]=%.12f", ij, heff[ij]);
          }
          (void)fprintf(outfile, "\n");
        }
      }
      (void)fflush(outfile);
   
      /* =============================================================== *
       * Diagonalization of Heff and update Upp and Vpp                  *
       * =============================================================== */
   
      for (i=0; i<nstates; i++) eigval[i] = 0.0;
      if (nstates == 1) {
        vpp[0][0] = 1.0;
        energs[0] = heff[0];
      } else {
        jacobi2dm(heff, energs, vpp, nstates, 0);
      }
   
      if (prntflag) {
        (void) fprintf(outfile, "\n\nState       Energy      ");
        (void) fprintf(outfile,   "\n (I)        /a.u./      ");
        (void) fprintf(outfile,   "\n.....  ..................");
        for (i=0; i<nstates; i++) 
          (void) fprintf(outfile,"\n  %d     %+.12lf",
            i+1, energs[i]+repnuc);
        (void) fflush(outfile);
        /* Test only */
        /*(void) fprintf(outfile, "\n\nState       Energy      ");
        (void) fprintf(outfile,   "\n (I)        /a.u./      ");
        (void) fprintf(outfile,   "\n.....  ..................");*/
        /* Test only */
      }
   
      for (i=0, diff=0.0; i<nstates; i++) {
        temp = old_energs[i] - energs[i];
        /* Test only */
        /*if (prntflag) {
          (void) fprintf(outfile,"\n  %d     %+16.14e", i+1, temp);
          (void)fflush(outfile);
        }*/
        /* Test only */
        diff += temp * temp;
        old_energs[i] = energs[i];
      }
      diff = sqrt(diff);
   
      if (prntflag) {
        (void)fprintf(outfile, "\n Convergence test: %16.12lf\n", diff);
        (void)fflush(outfile);
      }
   
      if ((!invflag)&&(diff<etol)) conv = 1;
   }

   MPI_Bcast(&conv, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (prntflag>1) prntflag = 1;

   if (my_rank == 0) { 
      /* update Upp and save it to Upptemp */
      multiply(upptemp, vpp, upp, nstates);
      for (i=0; i<nstates2; i++) upptemp[0][i] = upp[0][i]; 
   
      /* update (UW)pp to Upp */
      multiply(upptemp, wpp, upp, nstates);
   }


   if (my_rank == 0) { 
      if (invflag) {
   
          (void)fprintf(outfile, "\n******* inverse matrix ********* \n");
          (void)fprintf(outfile, "\n (HW)pp matrix for MRCISD(TQ):\n");
          for (i=0; i<nstates; i++) {
            for (j=0; j<nstates; j++) {
              (void)fprintf(outfile, " (HW)pp[%d][%d]=%.12f", i, j, hwpp2[i][j]);
            }
            (void)fprintf(outfile, "\n");
          }
          for (i=0; i<nstates; i++) {
            (void)fprintf(outfile, " Lower[%d][%d]=%.12f\n", 
             i, i, old_energs2[i]+hwpp2[i][i]+repnuc);
          }
   
          (void)fprintf(outfile, "\n (WW)pp matrix for MRCISD(TQ):\n");
          for (i=0; i<nstates; i++) {
            for (j=0; j<=i; j++) {
              (void)fprintf(outfile, " WW[%d][%d]=%.12f", i, j, wwpp2[i][j]);
            }
            (void)fprintf(outfile, "\n");
          }
          for (i=0; i<nstates; i++) {
            (void)fprintf(outfile, " Upper[%d][%d]=%.12f\n", i, i, 
                  old_energs2[i]+hwpp2[i][i]/(1.+wwpp2[i][i])+repnuc);
          }
   
      
         /* =============================================================== *
          * Update  (HW)pp                                                  *
          *         (HW)pp = Wpp(-1)*(HW)pp                                 *
          * scratch use of vpp[][]                                          *
          * =============================================================== */
   
          multiply(wppinv2, hwpp2, vpp, nstates);
   
        /* =============================================================== *
         * Calculation of Wpp and Wpp(-1)                                  *
         * =============================================================== */
   
        for (i=0; i<nstates; i++) eigval[i] = 0.0;
        if (nstates==1) {
          uvecs[0][0] = 1.0;
          eigval[0]   = sqrt(1.-wwpp2[0][0]);
        } else {
          for (i=0; i<nstates; i++) {
            for (j=0; j<i; j++) wwpp2[i][j] *= -1.;
            wwpp2[i][i] = 1.-wwpp2[i][i];
          }
          jacobi2dm(wwpp2[0], eigval, uvecs, nstates, 0);
          for (i=0; i<nstates;  i++) eigval[i] = sqrt(eigval[i]);
        }
   
        for (i=0; i<nstates; i++) {
          for (j=0; j<=i; j++) {
            sum = 0.0;
            sum2 = 0.0;
            for (k=0; k<nstates; k++) {
              temp = uvecs[i][k]*uvecs[j][k];
              sum += temp*eigval[k];
              sum2 += temp/eigval[k];
            }
            wpp2[i][j] = sum;
            wppinv2[i][j] = sum2;  
            if(i>j) {
              wpp2[j][i] = sum;
              wppinv2[j][i] = sum2; 
            }
          }
        } 
   
        /* =============================================================== *
         * Construct Heff                                                  *
         * =============================================================== */
   
        /* scratch use of upp2[][] to get (hpw)pp */
        for (i=0; i<nstates; i++) {
          for (j=0; j<=i; j++) {
            sum = 0.0;
            for (k=0; k<nstates; k++) {
              temp = upptemp2[k][i]*upptemp2[k][j];
              sum += temp*energs0[k];
            }
            upp2[i][j] = sum;
            if (i>j) {
              upp2[j][i] = sum;
            }
          }
        }
   
        /* add (hpw)pp to (hqw)pp */
        multiplyplus(upp2, wpp2, vpp, nstates);
   
        multiply(wppinv2, vpp, hwpp2, nstates);
   
        for (i=0, ij=0; i<nstates; i++) {
          for (j=0; j<i; j++, ij++) {
            heff[ij] = 0.5 * (hwpp2[i][j] + hwpp2[j][i]);
          }
          heff[ij++] = hwpp2[i][i];
        }
   
        if (prntflag) {
          (void)fprintf(outfile, "\n ****** with inverse matrix ******* \n");
          (void)fprintf(outfile, "\n (Heff)pp matrix for MRCISD(TQ):\n");
          for (i=0,ij=0; i<nstates; i++) {
            for (j=0; j<=i; j++,ij++) {
              (void)fprintf(outfile, " (Heff)pp[%d]=%.12f", ij, heff[ij]);
            }
            (void)fprintf(outfile, "\n");
          }
          (void)fflush(outfile);
        }
   
        /* =============================================================== *
         * Diagonalization of Heff and update Upp and Vpp                  *
         * =============================================================== */
   
        for (i=0; i<nstates; i++) eigval[i] = 0.0;
        if (nstates == 1) {
          vpp[0][0] = 1.0;
          energs2[0] = heff[0];
        } else {
          jacobi2dm(heff, energs2, vpp, nstates, 0);
        }
   
        if (prntflag) {
   
          (void) fprintf(outfile, "\n\nState       Energy      ");
          (void) fprintf(outfile,   "\n (I)        /a.u./      ");
          (void) fprintf(outfile,   "\n.....  ..................");
          for (i=0; i<nstates; i++) 
            (void) fprintf(outfile,"\n  %d     %+.12lf",
              i+1, energs2[i]+repnuc);
          (void) fflush(outfile);
        }
   
        for (i=0, diff=0.0; i<nstates; i++) {
          diff += (energs2[i] - old_energs2[i]) * (energs2[i] - old_energs2[i]);
          old_energs2[i] = energs2[i];
        }
        diff = sqrt(diff);
   
        if (prntflag) {
          (void)fprintf(outfile, "\n Convergence test: %16.12lf\n", diff);
          (void)fflush(outfile);
        }
        if (diff<etol) conv = 1;

      }
   }

   MPI_Bcast(&conv, 1, MPI_INT, 0, MPI_COMM_WORLD);

   if (prntflag>1) prntflag = 1;

   if (my_rank == 0) {
      if (invflag) { 
        /* update Upp and save it to Upptemp */
        multiply(upptemp2, vpp, upp2, nstates);
        for (i=0; i<nstates2; i++) upptemp2[0][i] = upp2[0][i]; 
   
        /* update (UW)pp to Upp */
        multiply(upptemp2, wpp2, upp2, nstates);
   
      } /* end of if (invflag) */
      if ((!iter)&&prntflag) {
        runtime = (double)clock()/CLOCKS_PER_SEC - cpu_time;
        (void)fprintf(outfile, 
         " \n\n first iteration cpu time MRCISD(TQ) = %.2lf s\n\n", runtime);
        lt    = time(NULL);
        tmptr = localtime(&lt);
        (void)fprintf(outfile, "%s", asctime(tmptr));
        (void)fflush(outfile);
      }
   }  
 } /* end of iterations */

 if (my_rank == 0) {
    if (iter==mxiter) 
      (void) fprintf(outfile,"\n MRCISD(TQ) does not converge in %d iterations!", iter);
   
    /* ................................... *
       Print energies of the MRCISD(TQ) states
     * ................................... */
   
    /*if (prntflag) {
   
      (void) fprintf(outfile, "\n\nState       Energy      ");
      (void) fprintf(outfile,   "\n (I)        /a.u./      ");
      (void) fprintf(outfile,   "\n.....  ..................");
      for (i=0; i<nstates; i++) 
        (void) fprintf(outfile,"\n  %d     %+.12lf",
          i+1, energs[i]+repnuc);
   
      (void) fflush(outfile);
    }*/
    
    if (prntflag) {
      runtime = (double)clock()/CLOCKS_PER_SEC - cpu_time;
      (void)fprintf(outfile, " \n\ntotal cpu time MRCISD(TQ) = %.2lf s\n\n",
        runtime);
      lt    = time(NULL);
      tmptr = localtime(&lt);
      (void)fprintf(outfile, "%s", asctime(tmptr));
      (void)fprintf(outfile, "\n**********    "
        "End of MRCISD(TQ) (CFGCI variant)    **********\n\n");
    }
 }

 (void)fclose(outfile);
 (void)fclose(infofile);

 MPI_Finalize();
 return 0;
}


/* C = A * B */
void multiply(double **a, double **b, double **c, int n) {

 int i, j, k;
 double temp;

 for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
     temp = 0.0;
     for (k=0; k<n; k++) {
       temp += a[i][k] * b[k][j]; 
     }
     c[i][j] = temp;
   }
 } 

}

/* C += A * B */
void multiplyplus(double **a, double **b, double **c, int n) {
 int i, j, k;
 double temp;

 for (i=0; i<n; i++) {
   for (j=0; j<n; j++) {
     temp = 0.0;
     for (k=0; k<n; k++) {
       temp += a[i][k] * b[k][j]; 
     }
     c[i][j] += temp;
   }
 } 

}

