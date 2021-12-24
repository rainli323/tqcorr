#include "cfgci_tq.h"

int choleski(double smat[], double amat[], double bmat[], int dim)
{
  /* ****************************************************************** *
   * The routine inverses a given symmetric positive definite matrix    *
   * S = A*A(+) of a dimension dim, using the algorithm of Choleski.    *
   * The low triangle of the inverse matrix S**(-1) = B(+)*B, where     *
   * B = A**(-1), is constructed in the array sinv[dim*(dim+1)/2]       *
   * (amat[])                                                           *
   * Input data:                                                        *
   *   dim                 - dimension of the S matrix                  *
   *                       - (dim > 1)                                  *
   *   smat[dim*(dim+1)/2] - the low triangle of the S matrix           *
   *   buf[dim*(dim+1)/2]  - a buffer array                             *
   *   (bmat[dim*(dim+1)/2])                                            *
   *                                                                    *
   * Output data:                                                       *
   *   the returning value = 0 if S is not positive definite and        *
   *                           has not been inverted                    * 
   *                         1 if S is positive definite and the        *
   *                           low triangle of the S**(-1) matrix       *
   *                           has been constructed (=> sinv[ddim])     *
   *                                                   (amat[])         *
   *                                                                    *
   * ****************************************************************** */
  int i, ii_in, ij,
      j, jj_in, jj,
      m, mj_in, mi, mj;
  double sum,   temp;
  
  /* ......................................... *
     Construct low triangle A matrix => amat[]
   * ......................................... */
  jj_in = 0; /* jj_in = j*(j+1)/2 */
  for (j=0; j<dim; j++) {
    jj = jj_in + j;
    sum = smat[jj];
    for (m=0; m<j; m++) { temp = amat[jj_in+m]; sum -= temp*temp; }
    if(sum < 1.0e-10) return(0);
    amat[jj] = temp = sqrt(sum);
    ii_in = ++jj; /* ii_in = i(i+1)/2 */
    for (i=j+1; i<dim; i++,ii_in+=i) {
      ij = ii_in + j;
      sum = smat[ij];
      for (m=0; m<j; m++) sum -= amat[ii_in+m]*amat[jj_in+m];
      amat[ij] = sum/temp;
    } /* loop over i */
    jj_in = jj;
  } /* loop over j */ 
   
  /* ............................................... *
     Construct low triangle A**(-1) matrix => bmat[] 
   * ............................................... */
  jj = 0; /* jj = j*(j+1)/2 */
  for (j=0; j<dim; j++) {
    jj += j;
    bmat[jj] = 1.0/amat[jj];
    mj_in = jj;
    ii_in = ++jj; /* ii_in = i(i+1)/2 */
    for (i=j+1; i<dim; i++,ii_in+=i) {
      mj  = mj_in; 
      sum = 0.0;
      for (m=j; m<i; m++,mj+=m) sum -= amat[ii_in+m]*bmat[mj];
      bmat[ii_in+j] = sum/amat[ii_in+i];
    } /* loop over i */
  } /* loop over j */
   
  /* ...................................................... *
     Construct low triangle of the S**(-1) matrix => sinv[]
   * ...................................................... */
  ii_in = 0; /* ii_in = i(i+1)/2 */
  for (i=0; i<dim; i++,ii_in+=i) {
    for (j=0; j<=i; j++) {
      mi  = ii_in + i;
      mj  = ii_in + j;
      sum = 0.0;
      for (m=i; m<dim; m++,mi+=m,mj+=m) sum += bmat[mi]*bmat[mj];
      amat[ii_in+j] = sum;
    } /* loop over i */
  } /* loop over j */
  
  return(1);
} /* int choleski */
  

  
