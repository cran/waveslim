#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

/*************************************************************************/

void dwt(double *Vin, int *M, int *L, double *h, double *g, 
	 double *Wout, double *Vout)
{

  int n, t, u;

  for(t = 0; t < *M/2; t++) {
    u = 2 * t + 1;
    Wout[t] = h[0] * Vin[u];
    Vout[t] = g[0] * Vin[u];
    for(n = 1; n < *L; n++) {
      u -= 1;
      if(u < 0) u = *M - 1;
      Wout[t] += h[n] * Vin[u];
      Vout[t] += g[n] * Vin[u];
    } 
  }
}

/*************************************************************************/

void idwt(double *Win, double *Vin, int *M, int *L, double *h, double *g, 
	  double *Xout)
{

  int i, j, l, t, u;
  int m = -2, n = -1;

  for(t = 0; t < *M; t++) {
    m += 2;
    n += 2;
    u = t;
    i = 1;
    j = 0;
    Xout[m] = h[i] * Win[u] + g[i] * Vin[u];
    Xout[n] = h[j] * Win[u] + g[j] * Vin[u];
    if(*L > 2) {
      for(l = 1; l < *L/2; l++) {
	u += 1;
	if(u >= *M) u = 0;
	i += 2;
	j += 2;
	Xout[m] += h[i] * Win[u] + g[i] * Vin[u];
	Xout[n] += h[j] * Win[u] + g[j] * Vin[u];
      }
    }
  }
}

/*************************************************************************/

void modwt(double *Vin, int *N, int *j, int *L, double *ht, double *gt, 
	   double *Wout, double *Vout)
{

  int k, n, t;

  for(t = 0; t < *N; t++) {
    k = t;
    Wout[t] = ht[0] * Vin[k];
    Vout[t] = gt[0] * Vin[k];
    for(n = 1; n < *L; n++) {
      k -= (int) pow(2.0, (double) *j - 1.0);
      if(k < 0) k += *N;
      Wout[t] += ht[n] * Vin[k];
      Vout[t] += gt[n] * Vin[k];
    }
  }

}

/*************************************************************************/

void imodwt(double *Win, double *Vin, int *N, int *j, int *L, 
	    double *ht, double *gt, double *Vout)
{

  int k, n, t;

  for(t = 0; t < *N; t++) {
    k = t;
    Vout[t] = (ht[0] * Win[k]) + (gt[0] * Vin[k]);
    for(n = 1; n < *L; n++) {
      k += (int) pow(2.0, (double) *j - 1.0);
      if(k >= *N) k -= *N;
      Vout[t] += (ht[n] * Win[k]) + (gt[n] * Vin[k]);
    }
  }
}

/***************************************************************************
 ***************************************************************************
   This DWT algorithm is shifted to the left by one in order to match with 
   the interval boundary conditions.
 ***************************************************************************
 ***************************************************************************/

void dwt_shift(double *Vin, int *M, int *L, double *h, double *g, 
	       double *Wout, double *Vout)
{

  int n, t, u;

  for(t = 0; t < *M/2; t++) {
    /* u = 2 * t + 1; */
    u = 2 * t + 2;
    Wout[t] = h[0] * Vin[u];
    Vout[t] = g[0] * Vin[u];
    for(n = 1; n < *L; n++) {
      u -= 1;
      if(u < 0) u = *M - 1;
      Wout[t] += h[n] * Vin[u];
      Vout[t] += g[n] * Vin[u];
    } 
  }
}

/***************************************************************************
 ***************************************************************************
   shifted iDWT
 ***************************************************************************
 ***************************************************************************/

void idwt_shift(double *Win, double *Vin, int M, int L, double *h, 
		double *g, double *Xout)
{

  int i, j, l, t, u;
  int m = -2, n = -1;

  for(t = 0; t < M; t++) {
    m += 2;
    n += 2;
    u = t;
    i = 1;
    j = 0;
    Xout[m] = h[i] * Win[u] + g[i] * Vin[u];
    Xout[n] = h[j] * Win[u] + g[j] * Vin[u];
    if(L > 2) {
      for(l = 1; l < L/2; l++) {
	u += 1;
	if(u >= M) u = 0;
	i += 2;
	j += 2;
	Xout[m] += h[i] * Win[u] + g[i] * Vin[u];
	Xout[n] += h[j] * Win[u] + g[j] * Vin[u];
      }
    }
  }
}

/***************************************************************************
 ***************************************************************************
   2D DWT 
 ***************************************************************************
 ***************************************************************************/

void two_D_dwt(double *X, int *M, int *N, int *L, double *h, double *g, 
	       double *LL, double *HL, double *LH, double *HH)
{
  int i, j, k;
  double *data, *Wout, *Vout, *Low, *High;

  /*
   *  Perform one-dimensional DWT on columns (length M).
   */

  Wout = (double *) malloc((*M) * sizeof(double));
  Vout = (double *) malloc((*M) * sizeof(double));

  /*
   *  Create temporary "matrices" to store DWT of columns.
   */
  Low = (double *) malloc((*N*(*M/2)) * sizeof(double));
  High = (double *) malloc((*N*(*M/2)) * sizeof(double));
  
  for(i = 0; i < *N; i++) {
    /*
     *  Must take column from X and place into vector for DWT.
     */
    data = (double *) malloc((*M) * sizeof(double));
    for(j = 0; j < *M; j++) {
      /* printf("X[%d][%d] = %f\n", i, j, X[i*(*M)+j]); */
      data[j] = X[i*(*M)+j];
    }
    /*
     *  Perform DWT and read into temporary matrices.
     */
    dwt(data, M, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *M/2; k++) {
      Low[i*(*M/2)+k] = Vout[k]; 
      /* printf("Low[%d][%d] = %f\n", i, k, Low[i*(*M/2)+k]); */
      High[i*(*M/2)+k] = Wout[k];
      /* printf("High[%d][%d] = %f\n", i, k, High[i*(*M/2)+k]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  /*
   *  Perform one-dimensional DWT on rows (length N).
   */

  Wout = (double *) malloc((*N) * sizeof(double));
  Vout = (double *) malloc((*N) * sizeof(double));

  for(i = 0; i < (int) *M/2; i++) {
    /*
     *  Must take row from "Low" and place into vector for DWT.
     */
    data = (double *) malloc((*N) * sizeof(double));
    for(j = 0; j < *N; j++) {
      /* printf("Low[%d][%d] = %f\n", i, j, Low[i+j*(*M/2)]); */
      data[j] = Low[i+j*(*M/2)];
    }
    /*
     *  Perform DWT and read into final "Low" matrices.
     */
    dwt(data, N, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *N/2; k++) {
      LL[i*(*N/2)+k] = Vout[k]; 
      /* printf("LL[%d][%d] = %f\n", i, k, LL[i*(*N/2)+k]); */
      HL[i*(*N/2)+k] = Wout[k];
      /* printf("LH[%d][%d] = %f\n", i, k, HL[i*(*N/2)+k]); */
    }
    free(data);

    /*
     *  Must take row from "High" and place into vector for DWT.
     */
    data = (double *) malloc((*N) * sizeof(double));
    for(j = 0; j < *N; j++) {
      /* printf("High[%d][%d] = %f\n", j, i, High[i+j*(*M/2)]); */
      data[j] = High[i+j*(*M/2)];
    }
    /*
     *  Perform DWT and read into final "High" matrices.
     */
    dwt(data, N, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *N/2; k++) {
      LH[i*(*N/2)+k] = Vout[k]; 
      /* printf("HL[%d][%d] = %f\n", i, k, LH[i*(*N/2)+k]); */
      HH[i*(*N/2)+k] = Wout[k];
      /* printf("HH[%d][%d] = %f\n", i, k, HH[i*(*N/2)+k]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  free(Low);
  free(High);

}

/***************************************************************************
 ***************************************************************************
   2D DWT using shifted DWT
 ***************************************************************************
 ***************************************************************************/

void two_D_dwt_shift(double *X, int *M, int *N, int *L, double *h, 
		     double *g, double *LL, double *LH, double *HL, 
		     double *HH)
{
  int i, j, k;
  double *data, *Wout, *Vout, *Low, *High;

  /*
   *  Perform one-dimensional DWT on columns (length M).
   */

  Wout = (double *) malloc((*M) * sizeof(double));
  Vout = (double *) malloc((*M) * sizeof(double));

  /*
   *  Create temporary "matrices" to store DWT of columns.
   */
  Low = (double *) malloc((*N*(*M/2)) * sizeof(double));
  High = (double *) malloc((*N*(*M/2)) * sizeof(double));
  
  for(i = 0; i < *N; i++) {
    /*
     *  Must take column from X and place into vector for DWT.
     */
    data = (double *) malloc((*M) * sizeof(double));
    for(j = 0; j < *M; j++) {
      /* printf("X[%d][%d] = %f\n", i, j, X[i*(*M)+j]); */
      data[j] = X[i*(*M)+j];
    }
    /*
     *  Perform DWT and read into temporary matrices.
     */
    dwt_shift(data, M, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *M/2; k++) {
      Low[i*(*M/2)+k] = Vout[k]; 
      /* printf("Low[%d][%d] = %f\n", i, k, Low[i*(*M/2)+k]); */
      High[i*(*M/2)+k] = Wout[k];
      /* printf("High[%d][%d] = %f\n", i, k, High[i*(*M/2)+k]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  /*
   *  Perform one-dimensional DWT on rows (length N).
   */

  Wout = (double *) malloc((*N) * sizeof(double));
  Vout = (double *) malloc((*N) * sizeof(double));

  for(i = 0; i < (int) *M/2; i++) {
    /*
     *  Must take row from "Low" and place into vector for DWT.
     */
    data = (double *) malloc((*N) * sizeof(double));
    for(j = 0; j < *N; j++) {
      /* printf("Low[%d][%d] = %f\n", i, j, Low[i+j*(*M/2)]); */
      data[j] = Low[i+j*(*M/2)];
    }
    /*
     *  Perform DWT and read into final "Low" matrices.
     */
    dwt_shift(data, N, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *N/2; k++) {
      LL[i*(*N/2)+k] = Vout[k]; 
      /* printf("LL[%d][%d] = %f\n", i, k, LL[i*(*N/2)+k]); */
      LH[i*(*N/2)+k] = Wout[k];
      /* printf("LH[%d][%d] = %f\n", i, k, LH[i*(*N/2)+k]); */
    }
    free(data);

    /*
     *  Must take row from "High" and place into vector for DWT.
     */
    data = (double *) malloc((*N) * sizeof(double));
    for(j = 0; j < *N; j++) {
      /* printf("High[%d][%d] = %f\n", j, i, High[i+j*(*M/2)]); */
      data[j] = High[i+j*(*M/2)];
    }
    /*
     *  Perform DWT and read into final "High" matrices.
     */
    dwt_shift(data, N, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *N/2; k++) {
      HL[i*(*N/2)+k] = Vout[k]; 
      /* printf("HL[%d][%d] = %f\n", i, k, HL[i*(*N/2)+k]); */
      HH[i*(*N/2)+k] = Wout[k];
      /* printf("HH[%d][%d] = %f\n", i, k, HH[i*(*N/2)+k]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  free(Low);
  free(High);

}

/***************************************************************************
 ***************************************************************************
   printdvec()
 ***************************************************************************
 ***************************************************************************/

void printdvec(double *v, int n)
{ 
  int i;

 for(i=0;i<=n-1;i++) printf("%f ",v[i]);
 printf("\n");
}

/***************************************************************************
 ***************************************************************************
   2D DWT only executed on the columns of input matrix (using shifted DWT)
 ***************************************************************************
 ***************************************************************************/

void dwt_columns(double *X, int *M, int *N, int *L, double *h, double *g,
		 double *Low, double *High)
{
  int i, j, k;
  double *data, *Wout, *Vout;

  /*
   *  Perform one-dimensional DWT on columns (length M).
   */

  Wout = (double *) malloc((*M/2) * sizeof(double));
  Vout = (double *) malloc((*M/2) * sizeof(double));
  data = (double *) malloc((*M) * sizeof(double));

  for(i = 0; i < *N; i++) {
    /*
     *  Must take column from X and place into vector for DWT.
     */
    for(j = 0; j < *M; j++)
      data[j] = X[i*(*M)+j];
    /* printdvec(data, *M-1); */
    /*
     *  Perform DWT and read into temporary matrices.
     */
    /* dwt(data, M, L, h, g, Wout, Vout); */
    dwt_shift(data, M, L, h, g, Wout, Vout);
    for(k = 0; k < (int) *M/2; k++) {
      Low[i*(*M/2)+k] = Vout[k]; 
      High[i*(*M/2)+k] = Wout[k];
    }
  }

  free(Wout);
  free(Vout);
  free(data);

}

/***************************************************************************
 ***************************************************************************
   2D iDWT
 ***************************************************************************
 ***************************************************************************/

void two_D_idwt(double *LL, double *HL, double *LH, double *HH, int *M, 
		int *N, int *L, double *h, double *g, double *image)
{
  int i, j, k;
  double *Win, *Vin, *Low, *High, *Xout;

  Low = (double *) malloc((*M)*2*(*N) * sizeof(double));
  High = (double *) malloc((*M)*2*(*N) * sizeof(double));

  Win = (double *) malloc((*N) * sizeof(double));
  Vin = (double *) malloc((*N) * sizeof(double));
  Xout = (double *) malloc(2*(*N) * sizeof(double));
  
  for(i = 0; i < *M; i++) {
    /*
     *  Must take row from LL and HL and place into vectors for iDWT.
     */
    for(j = 0; j < *N; j++) {
      Win[j] = HL[i*(*N)+j];
      Vin[j] = LL[i*(*N)+j];
    }
    
    idwt(Win, Vin, N, L, h, g, Xout);
    
    for(k = 0; k < 2*(*N); k++) {
      Low[i+k*(*M)] = Xout[k];
      /* printf("Low[%d][%d] = %f\n", k, i, Low[i+k*(*M)]); */
    }

    /*
     *  Must take row from LH and HH and place into vectors for iDWT.
     */
    for(j = 0; j < *N; j++) {
      Win[j] = HH[i*(*N)+j];
      Vin[j] = LH[i*(*N)+j];
    }

    idwt(Win, Vin, N, L, h, g, Xout);

    for(k = 0; k < 2*(*N); k++) {
      High[i+k*(*M)] = Xout[k];
      /* printf("High[%d][%d] = %f\n", k, i, High[i+k*(*M)]); */
    }

  }

  free(Vin);
  free(Win);
  free(Xout);

  Vin = (double *) malloc((*M) * sizeof(double));
  Win = (double *) malloc((*M) * sizeof(double));
  Xout = (double *) malloc(2*(*M) * sizeof(double));

  for(i = 0; i < 2*(*N); i++) {
    /*
     *  Must take columns from High and Low and place into vectors for iDWT.
     */
    for(k = 0; k < *M; k++) {
      Vin[k] = Low[i*(*M)+k];
      Win[k] = High[i*(*M)+k];
    }

    idwt(Win, Vin, M, L, h, g, Xout);

    for(j = 0; j < 2*(*M); j++) 
      image[i*2*(*M)+j] = Xout[j];

  }

  free(Vin);
  free(Win);
  free(Xout);

  free(Low);
  free(High);

}

/***************************************************************************
 ***************************************************************************
   2D iDWT using shifted DWT
 ***************************************************************************
 ***************************************************************************/

void two_D_idwt_shift(double *LL, double *LH, double *HL, double *HH, 
		      int *M, int *N, int *L, double *h, double *g, 
		      double *image)
{
  int i, j, k;
  double *Win, *Vin, *Low, *High, *Xout;

  Low = (double *) malloc((*M)*2*(*N) * sizeof(double));
  High = (double *) malloc((*M)*2*(*N) * sizeof(double));

  Win = (double *) malloc((*N) * sizeof(double));
  Vin = (double *) malloc((*N) * sizeof(double));
  Xout = (double *) malloc(2*(*N) * sizeof(double));
  
  for(i = 0; i < *M; i++) {
    /*
     *  Must take row from LL and LH and place into vectors for iDWT.
     */
    for(j = 0; j < *N; j++) {
      Win[j] = LH[i*(*M)+j];
      Vin[j] = LL[i*(*M)+j];
    }
    
    idwt_shift(Win, Vin, *N, *L, h, g, Xout);
    
    for(k = 0; k < 2*(*N); k++) 
      Low[i+k*(*M)] = Xout[k];

    /*
     *  Must take row from HL and HH and place into vectors for iDWT.
     */
    for(j = 0; j < *N; j++) {
      Win[j] = HH[i*(*M)+j];
      Vin[j] = HL[i*(*M)+j];
    }

    idwt_shift(Win, Vin, *N, *L, h, g, Xout);

    for(k = 0; k < 2*(*N); k++) 
      High[i+k*(*M)] = Xout[k];

  }

  free(Vin);
  free(Win);
  free(Xout);

  Vin = (double *) malloc((*M) * sizeof(double));
  Win = (double *) malloc((*M) * sizeof(double));
  Xout = (double *) malloc(2*(*M) * sizeof(double));

  for(i = 0; i < 2*(*N); i++) {
    /*
     *  Must take columns from High and Low and place into vectors for iDWT.
     */
    for(k = 0; k < *M; k++) {
      Vin[k] = Low[i*(*M)+k];
      Win[k] = High[i*(*M)+k];
    }

    idwt_shift(Win, Vin, *M, *L, h, g, Xout);

    for(j = 0; j < 2*(*M); j++) 
      image[i+j*2*(*N)] = Xout[j];

  }

  free(Vin);
  free(Win);
  free(Xout);

  free(Low);
  free(High);

}

/***************************************************************************
 ***************************************************************************
   2D iDWT only executed on the columns (using shifted iDWT)
 ***************************************************************************
 ***************************************************************************/

void idwt_columns(double *Low, double *High, int *M, int *N, int *L, 
		  double *h, double *g, double *X)
{
  int i, j, k;
  double *out, *Win, *Vin;

  /*
   *  Perform one-dimensional DWT on columns (length M).
   */

  Win = (double *) malloc((*M) * sizeof(double));
  Vin = (double *) malloc((*M) * sizeof(double));
  out = (double *) malloc(2*(*M) * sizeof(double));

  for(i = 0; i < *N; i++) {
    /*
     *  Must take column from Low and High and place into vectors for iDWT.
     */
    for(j = 0; j < *M; j++) {
      Vin[j] = Low[i*(*M)+j];
      Win[j] = High[i*(*M)+j];
    }
    /* printdvec(Vin, *M-1); printdvec(Win, *M-1); */
    /*
     *  Perform iDWT and read into a temporary matrix.
     */
    /* idwt(Win, Vin, M, L, h, g, out); */
    idwt_shift(Win, Vin, *M, *L, h, g, out);
    for(k = 0; k < (int) 2*(*M); k++)
      X[i*2*(*M)+k] = out[k]; 
  }

  free(Win);
  free(Vin);
  free(out);

}

/***************************************************************************
 ***************************************************************************
   2D MODWT 
 ***************************************************************************
 ***************************************************************************/

void two_D_modwt(double *X, int *M, int *N, int *J, int *L, double *h, 
		 double *g, double *LL, double *LH, double *HL, double *HH)
{
  int i, j, k;
  double *data, *Wout, *Vout, *Low, *High;

  /*
   *  Perform one-dimensional MODWT on columns (length M).
   */

  Wout = (double *) malloc((*M) * sizeof(double));
  Vout = (double *) malloc((*M) * sizeof(double));

  /*
   *  Create temporary "matrices" to store MODWT of columns.
   */
  Low = (double *) malloc((*N*(*M)) * sizeof(double));
  High = (double *) malloc((*N*(*M)) * sizeof(double));
  
  for(i = 0; i < *N; i++) {
    /*
     *  Must take column from X and place into vector for MODWT.
     */
    data = (double *) malloc((*M) * sizeof(double));
    for(j = 0; j < *M; j++) {
      /* printf("X[%d][%d] = %f\n", i, j, X[i*(*M)+j]); */
      data[j] = X[i*(*M)+j];
    }
    /*
     *  Perform MODWT and read into temporary matrices.
     */
    modwt(data, M, J, L, h, g, Wout, Vout);
    for(k = 0; k < *M; k++) {
      Low[i*(*M)+k] = Vout[k]; 
      /* printf("Low[%d][%d] = %f\n", i, k, Low[i*(*M)+k]); */
      High[i*(*M)+k] = Wout[k];
      /* printf("High[%d][%d] = %f\n", i, k, High[i*(*M)+k]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  /*
   *  Perform one-dimensional MODWT on rows (length N).
   */

  Wout = (double *) malloc((*N) * sizeof(double));
  Vout = (double *) malloc((*N) * sizeof(double));

  for(i = 0; i < *M; i++) {
    /*
     *  Must take row from "Low" and place into vector for DWT.
     */
    data = (double *) malloc((*N) * sizeof(double));
    for(j = 0; j < *N; j++) {
      /* printf("Low[%d][%d] = %f\n", i, j, Low[i+j*(*M)]); */
      data[j] = Low[i+j*(*M)];
    }
    /*
     *  Perform MODWT and read into final "Low" matrices.
     */
    modwt(data, N, J, L, h, g, Wout, Vout);
    for(k = 0; k < *N; k++) {
      LL[i*(*N)+k] = Vout[k]; 
      /* printf("LL[%d][%d] = %f\n", i, k, LL[i*(*N)+k]); */
      LH[i*(*N)+k] = Wout[k];
      /* printf("LH[%d][%d] = %f\n", i, k, LH[i*(*N)+k]); */
    }
    free(data);

    /*
     *  Must take row from "High" and place into vector for MODWT.
     */
    data = (double *) malloc((*N) * sizeof(double));
    for(j = 0; j < *N; j++) {
      /* printf("High[%d][%d] = %f\n", j, i, High[i+j*(*M)]); */
      data[j] = High[i+j*(*M)];
    }
    /*
     *  Perform MODWT and read into final "High" matrices.
     */
    modwt(data, N, J, L, h, g, Wout, Vout);
    for(k = 0; k < *N; k++) {
      HL[i*(*N)+k] = Vout[k]; 
      /* printf("HL[%d][%d] = %f\n", i, k, HL[i*(*N)+k]); */
      HH[i*(*N)+k] = Wout[k];
      /* printf("HH[%d][%d] = %f\n", i, k, HH[i*(*N)+k]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  free(Low);
  free(High);

}

/***************************************************************************
 ***************************************************************************
   2D iMODWT
 ***************************************************************************
 ***************************************************************************/

void two_D_imodwt(double *LL, double *LH, double *HL, double *HH, int *M, 
		  int *N, int *J, int *L, double *h, double *g, 
		  double *image)
{
  int i, j, k;
  double *Win, *Vin, *Low, *High, *Xout;

  Low = (double *) malloc((*M)*(*N) * sizeof(double));
  High = (double *) malloc((*M)*(*N) * sizeof(double));

  Win = (double *) malloc((*N) * sizeof(double));
  Vin = (double *) malloc((*N) * sizeof(double));
  Xout = (double *) malloc((*N) * sizeof(double));
  
  for(i = 0; i < *M; i++) {
    /*
     *  Must take row from LL and LH and place into vectors for iMODWT.
     */
    for(j = 0; j < *N; j++) {
      Win[j] = LH[i*(*M)+j];
      Vin[j] = LL[i*(*M)+j];
    }
    
    imodwt(Win, Vin, N, J, L, h, g, Xout);
    
    for(k = 0; k < *N; k++) 
      Low[i+k*(*M)] = Xout[k];

    /*
     *  Must take row from HL and HH and place into vectors for iMODWT.
     */
    for(j = 0; j < *N; j++) {
      Win[j] = HH[i*(*M)+j];
      Vin[j] = HL[i*(*M)+j];
    }

    imodwt(Win, Vin, N, J, L, h, g, Xout);

    for(k = 0; k < *N; k++) 
      High[i+k*(*M)] = Xout[k];

  }

  free(Vin);
  free(Win);
  free(Xout);

  Vin = (double *) malloc((*M) * sizeof(double));
  Win = (double *) malloc((*M) * sizeof(double));
  Xout = (double *) malloc((*M) * sizeof(double));

  for(i = 0; i < *N; i++) {
    /*
     *  Must take columns from High and Low and place into vectors for iMODWT.
     */
    for(k = 0; k < *M; k++) {
      Vin[k] = Low[i*(*M)+k];
      Win[k] = High[i*(*M)+k];
    }

    imodwt(Win, Vin, M, J, L, h, g, Xout);

    for(j = 0; j < *M; j++) 
      image[i+j*(*N)] = Xout[j];

  }

  free(Vin);
  free(Win);
  free(Xout);

  free(Low);
  free(High);

}

/***************************************************************************
 ***************************************************************************
   3D DWT 
 ***************************************************************************
 ***************************************************************************/

void three_D_dwt(double *X, int *NX, int *NY, int *NZ, int *L, 
		 double *h, double *g, double *LLL, double *HLL, 
		 double *LHL, double *LLH, double *HHL, double *HLH, 
		 double *LHH, double *HHH)
{
  int i, j, k, l;
  int printall = 0;
  double *data, *Wout, *Vout, *Xl, *Xh, *Yll, *Ylh, *Yhl, *Yhh;

  /*
    printf("Original Data (N = %d)...\n", *NX * (*NY) * (*NZ));
    printdvec(X, *NX * (*NY) * (*NZ));
  */

  /*
   *  Perform one-dimensional DWT on first dimension (length NX).
   */
  Wout = (double *) malloc((*NX) * sizeof(double));
  Vout = (double *) malloc((*NX) * sizeof(double));

  /*
   *  Create temporary "hyperrectangles" to store DWT of X-dimension.
   */
  Xl = (double *) malloc((*NZ*(*NY)*(*NX/2)) * sizeof(double));
  Xh = (double *) malloc((*NZ*(*NY)*(*NX/2)) * sizeof(double));
  
  for(i = 0; i < *NZ*(*NY); i++) {
    /*
     *  Must take column from X-dimension and place into vector for DWT.
     */
    data = (double *) malloc((*NX) * sizeof(double));
    for(j = 0; j < *NX; j++) {
      /* printf("X[%d][%d] = %f\n", i, j, X[i*(*M)+j]); */
      data[j] = X[i*(*NX)+j];
    }
    /*
     *  Perform DWT and read into temporary matrices.
     */
    dwt(data, NX, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NX/2; j++) {
      Xl[i*(*NX/2)+j] = Vout[j]; 
      /* printf("Low[%d][%d] = %f\n", i, j, Low[i*(*M/2)+j]); */
      Xh[i*(*NX/2)+j] = Wout[j];
      /* printf("High[%d][%d] = %f\n", i, j, High[i*(*M/2)+j]); */
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  /*
    printf("X Low...\n");
    printdvec(Xl, (*NX/2) * (*NY) * (*NZ));
    printf("X High...\n");
    printdvec(Xh, (*NX/2) * (*NY) * (*NZ));
  */

  /*
   *  Perform one-dimensional DWT on second dimension (length NY).
   */
  Wout = (double *) malloc((*NY) * sizeof(double));
  Vout = (double *) malloc((*NY) * sizeof(double));
  data = (double *) malloc((*NY) * sizeof(double));

  /*
   *  Create temporary "hyperrectangles" to store DWT of X-dimension.
   */
  Yll = (double *) malloc((*NZ*(*NY/2)*(*NX/2)) * sizeof(double));
  Ylh = (double *) malloc((*NZ*(*NY/2)*(*NX/2)) * sizeof(double));
  Yhl = (double *) malloc((*NZ*(*NY/2)*(*NX/2)) * sizeof(double));
  Yhh = (double *) malloc((*NZ*(*NY/2)*(*NX/2)) * sizeof(double));

  k = 0;
  l = 0;
  for(i = 0; i < *NZ * (int) *NX/2; i++) {
    /* 
     * Must adjust for 3D array structure.
     *   k: vertical dimension (Z) adjustment when reading in data
     *   l: vertical dimension (Z) adjustment when writing wavelet coeffs.
     */
    /* printf("fmod(%d, %d = %f\n", i, (int) *NX/2, fmod(i, (int) *NX/2)); */
    if(i > 0 && fmod(i, (int) *NX/2) == 0.0) {
      k = k + (*NY - 1) * ((int) *NX/2);
      l = l + ((int) *NY/2 - 1) * ((int) *NX/2);
    }
    /* printf("i = %d\tk = %d\tl = %d\n", i, k, l); */
    /*
     *  Must take row from "Xl" and place into vector for DWT.
     */
    for(j = 0; j < *NY; j++)
      data[j] = Xl[i + j * ((int) *NX/2) + k];
    /*
     *  Perform DWT and read into temporary "Yll" and "Yhl" hyperrectangles.
     */
    dwt(data, NY, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NY/2; j++) {
      Yll[i + j * ((int) *NX/2) + l] = Vout[j]; 
      Yhl[i + j * ((int) *NX/2) + l] = Wout[j];
      if(printall == 1)
	printf("Y.LL[%d][%d] = %f\nY.HL[%d][%d] = %f\n", i, j, 
	       Yll[i*((int) *NX/2)+j+k], i, j, Yhl[i*((int) *NX/2)+j+k]);
    }

    /*
     *  Must take row from "Xh" and place into vector for DWT.
     */
    data = (double *) malloc((*NY) * sizeof(double));
    for(j = 0; j < *NY; j++)
      data[j] = Xh[i + j * ((int) *NX/2) + k];
    /*
     *  Perform DWT and read into temporary "Yhl" and "Yhh" hyperrectangles.
     */
    dwt(data, NY, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NY/2; j++) {
      Ylh[i + j * ((int) *NX/2) + l] = Vout[j]; 
      Yhh[i + j * ((int) *NX/2) + l] = Wout[j];
      if(printall == 1)
	printf("Y.LH[%d][%d] = %f\nY.HH[%d][%d] = %f\n", i, j, 
	       Ylh[i*((int) *NX/2)+j+k], i, j, Ylh[i*((int) *NX/2)+j+k]);
    }
  }

  free(Wout);
  free(Vout);
  free(data);

  free(Xl);
  free(Xh);

  /*
    printf("Y Low-Low...\n");
    printdvec(Yll, (*NX/2) * (*NY/2) * (*NZ));
    printf("Y High-Low...\n");
    printdvec(Yhl, (*NX/2) * (*NY/2) * (*NZ));
    printf("Y Low-High...\n");
    printdvec(Ylh, (*NX/2) * (*NY/2) * (*NZ));
    printf("Y High-High...\n");
    printdvec(Yhh, (*NX/2) * (*NY/2) * (*NZ));
  */

  /*
   *  Perform one-dimensional DWT on third dimension (length NZ).
   */
  Wout = (double *) malloc((*NZ) * sizeof(double));
  Vout = (double *) malloc((*NZ) * sizeof(double));
  data = (double *) malloc((*NZ) * sizeof(double));

  for(i = 0; i < (int) *NY/2 * (int) *NX/2; i++) {
    /*
     *  Must take vertical column from "Yll" and place into vector for DWT.
     */
    for(j = 0; j < *NZ; j++)
      data[j] = Yll[i+j*((int) *NY/2 * (int) *NX/2)];
    /*
     *  Perform DWT and read into final "LLL" and "HLL" hyperrectangles.
     */
    dwt(data, NZ, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NZ/2; j++) {
      LLL[i+j*((int) *NY/2 * (int) *NX/2)] = Vout[j]; 
      HLL[i+j*((int) *NY/2 * (int) *NX/2)] = Wout[j];
      if(printall == 1)
	printf("LLL[%d][%d] = %f\nHLL[%d][%d] = %f\n", i, j, 
	       LLL[i+j*((int) *NY/2 * (int) *NX/2)], i, j, 
	       HLL[i+j*((int) *NY/2 * (int) *NX/2)]);
    }
    /*
     *  Must take row from "Yhl" and place into vector for DWT.
     */
    data = (double *) malloc((*NZ) * sizeof(double));
    for(j = 0; j < *NZ; j++) {
      /* printf("High[%d][%d] = %f\n", j, i, High[i+j*(*M/2)]); */
      data[j] = Yhl[i+j*((int) *NY/2 * (int) *NX/2)];
    }
    /*
     *  Perform DWT and read into final "LHL" and "HHL" hyperrectangles.
     */
    dwt(data, NZ, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NZ/2; j++) {
      LHL[i+j*((int) *NY/2 * (int) *NX/2)] = Vout[j]; 
      /* printf("HL[%d][%d] = %f\n", i, j, LH[i*(*N/2)+j]); */
      HHL[i+j*((int) *NY/2 * (int) *NX/2)] = Wout[j];
      /* printf("HH[%d][%d] = %f\n", i, j, HH[i*(*N/2)+j]); */
    }
    /*
     *  Must take row from "Ylh" and place into vector for DWT.
     */
    data = (double *) malloc((*NZ) * sizeof(double));
    for(j = 0; j < *NZ; j++) {
      /* printf("High[%d][%d] = %f\n", j, i, High[i+j*(*M/2)]); */
      data[j] = Ylh[i+j*((int) *NY/2 * (int) *NX/2)];
    }
    /*
     *  Perform DWT and read into final "LLH" and "HLH" hyperrectangles.
     */
    dwt(data, NZ, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NZ/2; j++) {
      LLH[i+j*((int) *NY/2 * (int) *NX/2)] = Vout[j]; 
      /* printf("HL[%d][%d] = %f\n", i, j, LH[i*(*N/2)+j]); */
      HLH[i+j*((int) *NY/2 * (int) *NX/2)] = Wout[j];
      /* printf("HH[%d][%d] = %f\n", i, j, HH[i*(*N/2)+j]); */
    }
    /*
     *  Must take row from "Yhh" and place into vector for DWT.
     */
    data = (double *) malloc((*NZ) * sizeof(double));
    for(j = 0; j < *NZ; j++) {
      /* printf("High[%d][%d] = %f\n", j, i, High[i+j*(*M/2)]); */
      data[j] = Yhh[i+j*((int) *NY/2 * (int) *NX/2)];
    }
    /*
     *  Perform DWT and read into final "LHH" and "HHH" hyperrectangles.
     */
    dwt(data, NZ, L, h, g, Wout, Vout);
    for(j = 0; j < (int) *NZ/2; j++) {
      LHH[i+j*((int) *NY/2 * (int) *NX/2)] = Vout[j]; 
      /* printf("HL[%d][%d] = %f\n", i, j, LH[i*(*N/2)+j]); */
      HHH[i+j*((int) *NY/2 * (int) *NX/2)] = Wout[j];
      /* printf("HH[%d][%d] = %f\n", i, j, HH[i*(*N/2)+j]); */
    }
  }

  free(Wout);
  free(Vout);
  free(data);

  free(Yll);
  free(Ylh);
  free(Yhl);
  free(Yhh);
}

/***************************************************************************
 ***************************************************************************
   3D iDWT
 ***************************************************************************
 ***************************************************************************/

void three_D_idwt(double *LLL, double *HLL, double *LHL, double *LLH, 
		  double *HHL, double *HLH, double *LHH, double *HHH, 
		  int *NX, int *NY, int *NZ, int *L, double *h, 
		  double *g, double *image)
{
  int i, j, k, l;
  int printall = 0;
  double *Win, *Vin, *Xl, *Xh, *Yll, *Ylh, *Yhl, *Yhh, *Xout;

  /*
   *  Create temporary "hyperrectangles" to store iDWT of Z-dimension.
   */
  Yll = (double *) malloc((2*(*NZ)*(*NY)*(*NX)) * sizeof(double));
  Ylh = (double *) malloc((2*(*NZ)*(*NY)*(*NX)) * sizeof(double));
  Yhl = (double *) malloc((2*(*NZ)*(*NY)*(*NX)) * sizeof(double));
  Yhh = (double *) malloc((2*(*NZ)*(*NY)*(*NX)) * sizeof(double));

  Win = (double *) malloc((*NZ) * sizeof(double));
  Vin = (double *) malloc((*NZ) * sizeof(double));
  Xout = (double *) malloc(2*(*NZ) * sizeof(double));
  
  for(i = 0; i < *NY * (*NX); i++) {
    /*
     *  Must take row from LLL and HLL and place into vectors for iDWT.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HLL[i + j * (*NY) * (*NX)];
      Vin[j] = LLL[i + j * (*NY) * (*NX)];
    }
    idwt(Win, Vin, NZ, L, h, g, Xout);
    for(j = 0; j < 2 * (*NZ); j++)
      Yll[i + j * (*NY) * (*NX)] = Xout[j];

    /*
     *  Must take row from LHL and HHL and place into vectors for iDWT.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HHL[i + j * (*NY) * (*NX)];
      Vin[j] = LHL[i + j * (*NY) * (*NX)];
    }
    idwt(Win, Vin, NZ, L, h, g, Xout);
    for(j = 0; j < 2 * (*NZ); j++)
      Yhl[i + j * (*NY) * (*NX)] = Xout[j];

    /*
     *  Must take row from LLH and HLH and place into vectors for iDWT.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HLH[i + j * (*NY) * (*NX)];
      Vin[j] = LLH[i + j * (*NY) * (*NX)];
    }
    idwt(Win, Vin, NZ, L, h, g, Xout);
    for(j = 0; j < 2 * (*NZ); j++)
      Ylh[i + j * (*NY) * (*NX)] = Xout[j];

    /*
     *  Must take row from LHH and HHH and place into vectors for iDWT.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HHH[i + j * (*NY) * (*NX)];
      Vin[j] = LHH[i + j * (*NY) * (*NX)];
    }
    idwt(Win, Vin, NZ, L, h, g, Xout);
    for(j = 0; j < 2 * (*NZ); j++)
      Yhh[i + j * (*NY) * (*NX)] = Xout[j];
  }

  free(Vin);
  free(Win);
  free(Xout);

  /*
    printf("Y Low-Low...\n");
    printdvec(Yll, (*NX) * (*NY) * 2 * (*NZ));
    printf("Y High-Low...\n");
    printdvec(Yhl, (*NX) * (*NY) * 2 * (*NZ));
    printf("Y Low-High...\n");
    printdvec(Ylh, (*NX) * (*NY) * 2 * (*NZ));
    printf("Y High-High...\n");
    printdvec(Yhh, (*NX) * (*NY) * 2 * (*NZ));
  */

  Xl = (double *) malloc((2*(*NZ)*2*(*NY)*(*NX)) * sizeof(double));
  Xh = (double *) malloc((2*(*NZ)*2*(*NY)*(*NX)) * sizeof(double));

  Vin = (double *) malloc((*NY) * sizeof(double));
  Win = (double *) malloc((*NY) * sizeof(double));
  Xout = (double *) malloc(2*(*NY) * sizeof(double));

  k = 0;
  l = 0;
  for(i = 0; i < 2 * (*NZ) * (*NX); i++) {
    /* 
     * Must adjust for 3D array structure.
     *   k: vertical dimension (Z) adjustment when reading in data
     *   l: vertical dimension (Z) adjustment when writing wavelet coeffs.
     */
    if(i > 0 && fmod(i, *NX) == 0.0) {
      k = k + (*NY - 1) * (*NX);
      l = l + (2 * (*NY) - 1) * (*NX);
    }
    /* printf("k = %d \t l = %d\n", k, l); */

    /*
     *  Must take columns from Yll and Yhl and place into vectors for iDWT.
     */
    for(j = 0; j < *NY; j++) {
      Vin[j] = Yll[i + j * (*NX) + k];
      Win[j] = Yhl[i + j * (*NX) + k];
    }
    idwt(Win, Vin, NY, L, h, g, Xout);
    for(j = 0; j < 2 * (*NY); j++) 
      Xl[i + j * (*NX) + l] = Xout[j];
    /*
     *  Must take columns from Ylh and Yhh and place into vectors for iDWT.
     */
    for(j = 0; j < *NY; j++) {
      Vin[j] = Ylh[i + j * (*NX) + k];
      Win[j] = Yhh[i + j * (*NX) + k];
    }
    idwt(Win, Vin, NY, L, h, g, Xout);
    for(j = 0; j < 2 * (*NY); j++) 
      Xh[i + j * (*NX) + l] = Xout[j];
  }

  /*
    printf("X Low...\n");
    printdvec(Xl, (*NX) * 2 * (*NY) * 2 * (*NZ));
    printf("X High...\n");
    printdvec(Xh, (*NX) * 2 * (*NY) * 2 * (*NZ));
  */

  free(Vin);
  free(Win);
  free(Xout);

  free(Yll);
  free(Ylh);
  free(Yhl);
  free(Yhh);

  Vin = (double *) malloc((*NX) * sizeof(double));
  Win = (double *) malloc((*NX) * sizeof(double));
  Xout = (double *) malloc(2*(*NX) * sizeof(double));

  for(i = 0; i < 2 * (*NZ) * 2 * (*NY); i++) {
    /*
     *  Must take columns from Xl and Xh and place into vectors for iDWT.
     */
    for(j = 0; j < *NX; j++) {
      Vin[j] = Xl[i * (*NX) + j];
      Win[j] = Xh[i * (*NX) + j];
    }
    idwt(Win, Vin, NX, L, h, g, Xout);
    for(j = 0; j < 2 * (*NX); j++)
      image[i * 2 * (*NX) + j] = Xout[j];
  }

  free(Vin);
  free(Win);
  free(Xout);

  free(Xl);
  free(Xh);
}

/***************************************************************************
 ***************************************************************************
   3D MODWT
 ***************************************************************************
 ***************************************************************************/

void three_D_modwt(double *X, int *NX, int *NY, int *NZ, int *J, int *L, 
		   double *h, double *g, double *LLL, double *HLL, 
		   double *LHL, double *LLH, double *HHL, double *HLH, 
		   double *LHH, double *HHH)
{
  int i, j, k;
  double *data, *Wout, *Vout, *Xl, *Xh, *Yll, *Ylh, *Yhl, *Yhh;

  /*
    printf("Original Data (N = %d)...\n", *NX * (*NY) * (*NZ));
    printdvec(X, *NX * (*NY) * (*NZ));
  */

  /*
   *  Perform one-dimensional MODWT on first dimension (length NX).
   */
  Wout = (double *) malloc((*NX) * sizeof(double));
  Vout = (double *) malloc((*NX) * sizeof(double));

  /*
   *  Create temporary "hyperrectangles" to store MODWT of X-dimension.
   */
  Xl = (double *) malloc((*NZ*(*NY)*(*NX)) * sizeof(double));
  Xh = (double *) malloc((*NZ*(*NY)*(*NX)) * sizeof(double));
  
  for(i = 0; i < *NZ*(*NY); i++) {
    /*
     *  Must take column from X-dimension and place into vector for DWT.
     */
    data = (double *) malloc((*NX) * sizeof(double));
    for(j = 0; j < *NX; j++) {
      data[j] = X[i * (*NX) + j];
    }
    /*
     *  Perform MODWT and read into temporary matrices.
     */
    modwt(data, NX, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NX; j++) {
      Xl[i * (*NX) + j] = Vout[j]; 
      Xh[i * (*NX) + j] = Wout[j];
    }
    free(data);
  }

  free(Wout);
  free(Vout);

  /*
    printf("X Low...\n");
    printdvec(Xl, (*NX) * (*NY) * (*NZ));
    printf("X High...\n");
    printdvec(Xh, (*NX) * (*NY) * (*NZ));
  */

  /*
   *  Perform one-dimensional MODWT on second dimension (length NY).
   */
  Wout = (double *) malloc((*NY) * sizeof(double));
  Vout = (double *) malloc((*NY) * sizeof(double));
  data = (double *) malloc((*NY) * sizeof(double));

  /*
   *  Create temporary "hyperrectangles" to store MODWT of X-dimension.
   */
  Yll = (double *) malloc((*NZ*(*NY)*(*NX)) * sizeof(double));
  Ylh = (double *) malloc((*NZ*(*NY)*(*NX)) * sizeof(double));
  Yhl = (double *) malloc((*NZ*(*NY)*(*NX)) * sizeof(double));
  Yhh = (double *) malloc((*NZ*(*NY)*(*NX)) * sizeof(double));

  k = 0;
  for(i = 0; i < *NZ * (*NX); i++) {
    /* 
     * Must adjust for 3D array structure.
     *   k: vertical dimension (Z) adjustment when reading in data
     *   l: vertical dimension (Z) adjustment when writing wavelet coeffs.
     */
    if(i > 0 && fmod(i, *NX) == 0.0)
      k = k + (*NY - 1) * (*NX);
    /*
     *  Must take row from "Xl" and place into vector for DWT.
     */
    for(j = 0; j < *NY; j++)
      data[j] = Xl[i + j * (*NX) + k];
    /*
     *  Perform MODWT and read into temporary "Yll" and "Yhl" hyperrectangles.
     */
    modwt(data, NY, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NY; j++) {
      Yll[i + j * (*NX) + k] = Vout[j]; 
      Yhl[i + j * (*NX) + k] = Wout[j];
    }
    /*
     *  Must take row from "Xh" and place into vector for DWT.
     */
    data = (double *) malloc((*NY) * sizeof(double));
    for(j = 0; j < *NY; j++)
      data[j] = Xh[i + j * (*NX) + k];
    /*
     *  Perform MODWT and read into temporary "Yhl" and "Yhh" hyperrectangles.
     */
    modwt(data, NY, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NY; j++) {
      Ylh[i + j * (*NX) + k] = Vout[j]; 
      Yhh[i + j * (*NX) + k] = Wout[j];
    }
  }

  free(Wout);
  free(Vout);
  free(data);

  free(Xl);
  free(Xh);

  /*
    printf("Y Low-Low...\n");
    printdvec(Yll, (*NX) * (*NY) * (*NZ));
    printf("Y High-Low...\n");
    printdvec(Yhl, (*NX) * (*NY) * (*NZ));
    printf("Y Low-High...\n");
    printdvec(Ylh, (*NX) * (*NY) * (*NZ));
    printf("Y High-High...\n");
    printdvec(Yhh, (*NX) * (*NY) * (*NZ));
  */

  /*
   *  Perform one-dimensional MODWT on third dimension (length NZ).
   */
  Wout = (double *) malloc((*NZ) * sizeof(double));
  Vout = (double *) malloc((*NZ) * sizeof(double));
  data = (double *) malloc((*NZ) * sizeof(double));

  for(i = 0; i < *NY * (*NX); i++) {
    /*
     *  Must take vertical column from "Yll" and place into vector for MODWT.
     */
    for(j = 0; j < *NZ; j++)
      data[j] = Yll[i + j * (*NY) * (*NX)];
    /*
     *  Perform MODWT and read into final "LLL" and "HLL" hyperrectangles.
     */
    modwt(data, NZ, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NZ; j++) {
      LLL[i + j * (*NY) * (*NX)] = Vout[j]; 
      HLL[i + j * (*NY) * (*NX)] = Wout[j];
    }
    /*
     *  Must take row from "Yhl" and place into vector for MODWT.
     */
    data = (double *) malloc((*NZ) * sizeof(double));
    for(j = 0; j < *NZ; j++) {
      data[j] = Yhl[i + j * (*NY) * (*NX)];
    }
    /*
     *  Perform MODWT and read into final "LHL" and "HHL" hyperrectangles.
     */
    modwt(data, NZ, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NZ; j++) {
      LHL[i + j * (*NY) * (*NX)] = Vout[j]; 
      HHL[i + j * (*NY) * (*NX)] = Wout[j];
    }
    /*
     *  Must take row from "Ylh" and place into vector for MODWT.
     */
    data = (double *) malloc((*NZ) * sizeof(double));
    for(j = 0; j < *NZ; j++) {
      data[j] = Ylh[i + j * (*NY) * (*NX)];
    }
    /*
     *  Perform MODWT and read into final "LLH" and "HLH" hyperrectangles.
     */
    modwt(data, NZ, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NZ; j++) {
      LLH[i + j * (*NY) * (*NX)] = Vout[j]; 
      HLH[i + j * (*NY) * (*NX)] = Wout[j];
    }
    /*
     *  Must take row from "Yhh" and place into vector for MODWT.
     */
    data = (double *) malloc((*NZ) * sizeof(double));
    for(j = 0; j < *NZ; j++) {
      data[j] = Yhh[i + j * (*NY) * (*NX)];
    }
    /*
     *  Perform MODWT and read into final "LHH" and "HHH" hyperrectangles.
     */
    modwt(data, NZ, J, L, h, g, Wout, Vout);
    for(j = 0; j < *NZ; j++) {
      LHH[i + j * (*NY) * (*NX)] = Vout[j]; 
      HHH[i + j * (*NY) * (*NX)] = Wout[j];
    }
  }

  free(Wout);
  free(Vout);
  free(data);

  free(Yll);
  free(Ylh);
  free(Yhl);
  free(Yhh);
}

/***************************************************************************
 ***************************************************************************
   3D iMODWT
 ***************************************************************************
 ***************************************************************************/

void three_D_imodwt(double *LLL, double *HLL, double *LHL, double *LLH, 
		    double *HHL, double *HLH, double *LHH, double *HHH, 
		    int *NX, int *NY, int *NZ, int *J, int *L, double *h, 
		    double *g, double *image)
{
  int i, j, k;
  double *Win, *Vin, *Xl, *Xh, *Yll, *Ylh, *Yhl, *Yhh, *Xout;

  /*
   *  Create temporary "hyperrectangles" to store imodwt of Z-dimension.
   */
  Yll = (double *) malloc(((*NZ)*(*NY)*(*NX)) * sizeof(double));
  Ylh = (double *) malloc(((*NZ)*(*NY)*(*NX)) * sizeof(double));
  Yhl = (double *) malloc(((*NZ)*(*NY)*(*NX)) * sizeof(double));
  Yhh = (double *) malloc(((*NZ)*(*NY)*(*NX)) * sizeof(double));

  Win = (double *) malloc((*NZ) * sizeof(double));
  Vin = (double *) malloc((*NZ) * sizeof(double));
  Xout = (double *) malloc((*NZ) * sizeof(double));
  
  for(i = 0; i < *NY * (*NX); i++) {
    /*
     *  Must take row from LLL and HLL and place into vectors for imodwt.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HLL[i + j * (*NY) * (*NX)];
      Vin[j] = LLL[i + j * (*NY) * (*NX)];
    }
    imodwt(Win, Vin, NZ, J, L, h, g, Xout);
    for(j = 0; j < *NZ; j++)
      Yll[i + j * (*NY) * (*NX)] = Xout[j];

    /*
     *  Must take row from LHL and HHL and place into vectors for imodwt.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HHL[i + j * (*NY) * (*NX)];
      Vin[j] = LHL[i + j * (*NY) * (*NX)];
    }
    imodwt(Win, Vin, NZ, J, L, h, g, Xout);
    for(j = 0; j < *NZ; j++)
      Yhl[i + j * (*NY) * (*NX)] = Xout[j];

    /*
     *  Must take row from LLH and HLH and place into vectors for imodwt.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HLH[i + j * (*NY) * (*NX)];
      Vin[j] = LLH[i + j * (*NY) * (*NX)];
    }
    imodwt(Win, Vin, NZ, J, L, h, g, Xout);
    for(j = 0; j < *NZ; j++)
      Ylh[i + j * (*NY) * (*NX)] = Xout[j];

    /*
     *  Must take row from LHH and HHH and place into vectors for imodwt.
     */
    for(j = 0; j < *NZ; j++) {
      Win[j] = HHH[i + j * (*NY) * (*NX)];
      Vin[j] = LHH[i + j * (*NY) * (*NX)];
    }
    imodwt(Win, Vin, NZ, J, L, h, g, Xout);
    for(j = 0; j < *NZ; j++)
      Yhh[i + j * (*NY) * (*NX)] = Xout[j];
  }

  free(Vin);
  free(Win);
  free(Xout);

  /*
    printf("Y Low-Low...\n");
    printdvec(Yll, (*NX) * (*NY) * (*NZ));
    printf("Y High-Low...\n");
    printdvec(Yhl, (*NX) * (*NY) * (*NZ));
    printf("Y Low-High...\n");
    printdvec(Ylh, (*NX) * (*NY) * (*NZ));
    printf("Y High-High...\n");
    printdvec(Yhh, (*NX) * (*NY) * (*NZ));
  */

  Xl = (double *) malloc(((*NZ)*(*NY)*(*NX)) * sizeof(double));
  Xh = (double *) malloc(((*NZ)*(*NY)*(*NX)) * sizeof(double));

  Vin = (double *) malloc((*NY) * sizeof(double));
  Win = (double *) malloc((*NY) * sizeof(double));
  Xout = (double *) malloc((*NY) * sizeof(double));

  k = 0;
  for(i = 0; i < (*NZ) * (*NX); i++) {
    /* 
     * Must adjust for 3D array structure.
     *   k: vertical dimension (Z) adjustment when reading in data
     *   l: vertical dimension (Z) adjustment when writing wavelet coeffs.
     */
    if(i > 0 && fmod(i, *NX) == 0.0)
      k = k + (*NY - 1) * (*NX);
    /*
     *  Must take columns from Yll and Yhl and place into vectors for imodwt.
     */
    for(j = 0; j < *NY; j++) {
      Vin[j] = Yll[i + j * (*NX) + k];
      Win[j] = Yhl[i + j * (*NX) + k];
    }
    imodwt(Win, Vin, NY, J, L, h, g, Xout);
    for(j = 0; j < (*NY); j++) 
      Xl[i + j * (*NX) + k] = Xout[j];
    /*
     *  Must take columns from Ylh and Yhh and place into vectors for imodwt.
     */
    for(j = 0; j < *NY; j++) {
      Vin[j] = Ylh[i + j * (*NX) + k];
      Win[j] = Yhh[i + j * (*NX) + k];
    }
    imodwt(Win, Vin, NY, J, L, h, g, Xout);
    for(j = 0; j < (*NY); j++) 
      Xh[i + j * (*NX) + k] = Xout[j];
  }

  /*
    printf("X Low...\n");
    printdvec(Xl, (*NX) * (*NY) * (*NZ));
    printf("X High...\n");
    printdvec(Xh, (*NX) * (*NY) * (*NZ));
  */

  free(Vin);
  free(Win);
  free(Xout);

  free(Yll);
  free(Ylh);
  free(Yhl);
  free(Yhh);

  Vin = (double *) malloc((*NX) * sizeof(double));
  Win = (double *) malloc((*NX) * sizeof(double));
  Xout = (double *) malloc((*NX) * sizeof(double));

  for(i = 0; i < (*NZ) * (*NY); i++) {
    /*
     *  Must take columns from Xl and Xh and place into vectors for imodwt.
     */
    for(j = 0; j < *NX; j++) {
      Vin[j] = Xl[i * (*NX) + j];
      Win[j] = Xh[i * (*NX) + j];
    }
    imodwt(Win, Vin, NX, J, L, h, g, Xout);
    for(j = 0; j < (*NX); j++)
      image[i * (*NX) + j] = Xout[j];
  }

  free(Vin);
  free(Win);
  free(Xout);

  free(Xl);
  free(Xh);
}

