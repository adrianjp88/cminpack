/*      driver for hybrd1 example. */


#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <minpack.h>
#define real __minpack_real__

void fcn(const int *n, const real *x, real *fvec, int *iflag);

int main()
{
  int j, n, info, lwa;
  real tol, fnorm;
  real x[9], fvec[9], wa[180];
  int one=1;

  n = 9;

/*      the following starting values provide a rough solution. */

  for (j=1; j<=9; j++)
    {
      x[j-1] = -1.;
    }

  lwa = 180;

/*      set tol to the square root of the machine precision. */
/*      unless high solutions are required, */
/*      this is the recommended setting. */

  tol = sqrt(__minpack_func__(dpmpar)(&one));
  __minpack_func__(hybrd1)(&fcn, &n, x, fvec, &tol, &info, wa, &lwa);
  fnorm = __minpack_func__(enorm)(&n, fvec);

  printf("     final L2 norm of the residuals %15.7g\n", (double)fnorm);
  printf("     exit parameter                 %10i\n", info);
  printf("     final approximates solution\n");
  for (j=1; j<=n; j++) printf("%s%15.7g",j%3==1?"\n     ":"", (double)x[j-1]);
  printf("\n");

  return 0;
}

void fcn(const int *n, const real *x, real *fvec, int *iflag)
{
/*      subroutine fcn for hybrd1 example. */

  int k;
  real one=1, temp, temp1, temp2, three=3, two=2, zero=0;
  assert(*n == 9);
  (void)iflag;

  for (k=1; k <= *n; k++)
    {
      temp = (three - two*x[k-1])*x[k-1];
      temp1 = zero;
      if (k != 1) temp1 = x[k-1-1];
      temp2 = zero;
      if (k != *n) temp2 = x[k+1-1];
      fvec[k-1] = temp - temp1 - two*temp2 + one;
    }
  return;
}
