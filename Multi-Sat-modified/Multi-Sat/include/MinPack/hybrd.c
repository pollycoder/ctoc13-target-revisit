/* hybrd.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "cminpack.h"
#include <math.h>
#include "cminpackP.h"


int hybrd(int(* fcnnn)(int n, const double *x, double *fvec, int iflag, const double* para), int n, double *x, double *	fvec, 
		  const double* para, double xtol, int maxfev, int ml, int mu, 
	double epsfcn, double *diag, int mode, double
	factor, int nprint, int *nfev, double *
	fjac, int ldfjac, double *r, int lr, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4)
{
    /* Initialized data */

#define p1 .1
#define p5 .5
#define p001 .001
#define p0001 1e-4

    /* System generated locals */
    int fjac_dim1, fjac_offset, i1;
    double d1, d2;

    /* Local variables */
    int i, j, l, jm1, iwa[1];
    double sum;
    int sing;
    int iter = 0;
    double temp;
    int msum, iflag;
    double delta = 0.;
    int jeval;
    int ncsuc;
    double ratio;
    double fnorm;
    double pnorm, xnorm = 0., fnorm1;
    int nslow1, nslow2;
    int ncfail;
    double actred, epsmch, prered;
    int info;

/*     ********** */

/*     subroutine hybrd */

/*     the purpose of hybrd is to find a zero of a system of */
/*     n nonlinear functions in n variables by a modification */
/*     of the powell hybrid method. the user must provide a */
/*     subroutine which calculates the functions. the jacobian is */
/*     then calculated by a forward-difference approximation. */

/*     the subroutine statement is */

/*       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn, */
/*                        diag,mode,factor,nprint,info,nfev,fjac, */
/*                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4) */

/*     where */

/*       fcn is the name of the user-supplied subroutine which */
/*         calculates the functions. fcn must be declared */
/*         in an external statement in the user calling */
/*         program, and should be written as follows. */

/*         subroutine fcn(n,x,fvec,iflag) */
/*         integer n,iflag */
/*         double precision x(n),fvec(n) */
/*         ---------- */
/*         calculate the functions at x and */
/*         return this vector in fvec. */
/*         --------- */
/*         return */
/*         end */

/*         the value of iflag should not be changed by fcn unless */
/*         the user wants to terminate execution of hybrd. */
/*         in this case set iflag to a negative integer. */

/*       n is a positive integer input variable set to the number */
/*         of functions and variables. */

/*       x is an array of length n. on input x must contain */
/*         an initial estimate of the solution vector. on output x */
/*         contains the final estimate of the solution vector. */

/*       fvec is an output array of length n which contains */
/*         the functions evaluated at the output x. */

/*       xtol is a nonnegative input variable. termination */
/*         occurs when the relative error between two consecutive */
/*         iterates is at most xtol. */

/*       maxfev is a positive integer input variable. termination */
/*         occurs when the number of calls to fcn is at least maxfev */
/*         by the end of an iteration. */

/*       ml is a nonnegative integer input variable which specifies */
/*         the number of subdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         ml to at least n - 1. */

/*       mu is a nonnegative integer input variable which specifies */
/*         the number of superdiagonals within the band of the */
/*         jacobian matrix. if the jacobian is not banded, set */
/*         mu to at least n - 1. */

/*       epsfcn is an input variable used in determining a suitable */
/*         step length for the forward-difference approximation. this */
/*         approximation assumes that the relative errors in the */
/*         functions are of the order of epsfcn. if epsfcn is less */
/*         than the machine precision, it is assumed that the relative */
/*         errors in the functions are of the order of the machine */
/*         precision. */

/*       diag is an array of length n. if mode = 1 (see */
/*         below), diag is internally set. if mode = 2, diag */
/*         must contain positive entries that serve as */
/*         multiplicative scale factors for the variables. */

/*       mode is an integer input variable. if mode = 1, the */
/*         variables will be scaled internally. if mode = 2, */
/*         the scaling is specified by the input diag. other */
/*         values of mode are equivalent to mode = 1. */

/*       factor is a positive input variable used in determining the */
/*         initial step bound. this bound is set to the product of */
/*         factor and the euclidean norm of diag*x if nonzero, or else */
/*         to factor itself. in most cases factor should lie in the */
/*         interval (.1,100.). 100. is a generally recommended value. */

/*       nprint is an integer input variable that enables controlled */
/*         printing of iterates if it is positive. in this case, */
/*         fcn is called with iflag = 0 at the beginning of the first */
/*         iteration and every nprint iterations thereafter and */
/*         immediately prior to return, with x and fvec available */
/*         for printing. if nprint is not positive, no special calls */
/*         of fcn with iflag = 0 are made. */

/*       info is an integer output variable. if the user has */
/*         terminated execution, info is set to the (negative) */
/*         value of iflag. see description of fcn. otherwise, */
/*         info is set as follows. */

/*         info = 0   improper input parameters. */

/*         info = 1   relative error between two consecutive iterates */
/*                    is at most xtol. */

/*         info = 2   number of calls to fcn has reached or exceeded */
/*                    maxfev. */

/*         info = 3   xtol is too small. no further improvement in */
/*                    the approximate solution x is possible. */

/*         info = 4   iteration is not making good progress, as */
/*                    measured by the improvement from the last */
/*                    five jacobian evaluations. */

/*         info = 5   iteration is not making good progress, as */
/*                    measured by the improvement from the last */
/*                    ten iterations. */

/*       nfev is an integer output variable set to the number of */
/*         calls to fcn. */

/*       fjac is an output n by n array which contains the */
/*         orthogonal matrix q produced by the qr factorization */
/*         of the final approximate jacobian. */

/*       ldfjac is a positive integer input variable not less than n */
/*         which specifies the leading dimension of the array fjac. */

/*       r is an output array of length lr which contains the */
/*         upper triangular matrix produced by the qr factorization */
/*         of the final approximate jacobian, stored rowwise. */

/*       lr is a positive integer input variable not less than */
/*         (n*(n+1))/2. */

/*       qtf is an output array of length n which contains */
/*         the vector (q transpose)*fvec. */

/*       wa1, wa2, wa3, and wa4 are work arrays of length n. */

/*     subprograms called */

/*       user-supplied ...... fcn */

/*       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1, */
/*                            qform,qrfac,r1mpyq,r1updt */

/*       fortran-supplied ... dabs,dmax1,dmin1,min0,mod */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
    /* Parameter adjustments */
    --wa4;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --diag;
    --fvec;
    --x;
    fjac_dim1 = ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;
    --r;

    /* Function Body */

/*     epsmch is the machine precision. */

    epsmch = dpmpar(1);

    info = 0;
    iflag = 0;
    *nfev = 0;

/*     check the input parameters for errors. */

    if (n <= 0 || xtol < 0. || maxfev <= 0 || ml < 0 || mu < 0 ||
	    factor <= 0. || ldfjac < n || lr < n * (n + 1) / 2) {
	goto TERMINATE;
    }
    if (mode == 2) {
        for (j = 1; j <= n; ++j) {
            if (diag[j] <= 0.) {
                goto TERMINATE;
            }
        }
    }

/*     evaluate the function at the starting point */
/*     and calculate its norm. */
	if(nprint>0)
		printf("%5s%6s%9s%11s\n", "iter.","nfev","f(x)","step");

    iflag = fcnnn(n, &x[1], &fvec[1], 1, para);
    *nfev = 1;
    if (iflag < 0) {
	goto TERMINATE;
    }
    fnorm = enorm(n, &fvec[1]);

/*     determine the number of calls to fcn needed to compute */
/*     the jacobian matrix. */

/* Computing MIN */
    i1 = ml + mu + 1;
    msum = min(i1,n);

/*     initialize iteration counter and monitors. */

    iter = 1;
    ncsuc = 0;
    ncfail = 0;
    nslow1 = 0;
    nslow2 = 0;

/*     beginning of the outer loop. */

    for (;;) {
        jeval = TRUE_;

/*        calculate the jacobian matrix. */

        iflag = fdjac1(fcnnn, n, &x[1], &fvec[1], para, &fjac[fjac_offset], ldfjac,
                       ml, mu, epsfcn, &wa1[1], &wa2[1]);
        *nfev += msum;
        if (iflag < 0) {
            goto TERMINATE;
        }

/*        compute the qr factorization of the jacobian. */

        qrfac(n, n, &fjac[fjac_offset], ldfjac, FALSE_, iwa, 1,
              &wa1[1], &wa2[1], &wa3[1]);

/*        on the first iteration and if mode is 1, scale according */
/*        to the norms of the columns of the initial jacobian. */

        if (iter == 1) {
            if (mode != 2) {
                for (j = 1; j <= n; ++j) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.) {
                        diag[j] = 1.;
                    }
                }
            }

/*        on the first iteration, calculate the norm of the scaled x */
/*        and initialize the step bound delta. */

            for (j = 1; j <= n; ++j) {
                wa3[j] = diag[j] * x[j];
            }
            xnorm = enorm(n, &wa3[1]);
            delta = factor * xnorm;
            if (delta == 0.) {
                delta = factor;
            }
        }

/*        form (q transpose)*fvec and store in qtf. */

        for (i = 1; i <= n; ++i) {
            qtf[i] = fvec[i];
        }
        for (j = 1; j <= n; ++j) {
            if (fjac[j + j * fjac_dim1] != 0.) {
                sum = 0.;
                for (i = j; i <= n; ++i) {
                    sum += fjac[i + j * fjac_dim1] * qtf[i];
                }
                temp = -sum / fjac[j + j * fjac_dim1];
                for (i = j; i <= n; ++i) {
                    qtf[i] += fjac[i + j * fjac_dim1] * temp;
                }
            }
        }

/*        copy the triangular factor of the qr factorization into r. */

        sing = FALSE_;
        for (j = 1; j <= n; ++j) {
            l = j;
            jm1 = j - 1;
            if (jm1 >= 1) {
                for (i = 1; i <= jm1; ++i) {
                    r[l] = fjac[i + j * fjac_dim1];
                    l = l + n - i;
                }
            }
            r[l] = wa1[j];
            if (wa1[j] == 0.) {
                sing = TRUE_;
            }
        }

/*        accumulate the orthogonal factor in fjac. */

        qform(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);

/*        rescale if necessary. */

        if (mode != 2) {
            for (j = 1; j <= n; ++j) {
                /* Computing MAX */
                d1 = diag[j], d2 = wa2[j];
                diag[j] = max(d1,d2);
            }
        }

/*        beginning of the inner loop. */

        for (;;) {

/*           if requested, call fcn to enable printing of iterates. */

            if (nprint > 0) {
                iflag = 0;
                if ((iter - 1) % nprint == 0) {
                    iflag = fcnnn(n, &x[1], &fvec[1], 0, para);					
					printf("%4d%7d%12.3e%12.3e\n",iter,*nfev,enorm(n, &fvec[1]), delta);
                }
                if (iflag < 0) {
                    goto TERMINATE;
                }
            }

/*           determine the direction p. */

            dogleg(n, &r[1], lr, &diag[1], &qtf[1], delta, &wa1[1], &wa2[1], &wa3[1]);

/*           store the direction p and x + p. calculate the norm of p. */

            for (j = 1; j <= n; ++j) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            pnorm = enorm(n, &wa3[1]);

/*           on the first iteration, adjust the initial step bound. */

            if (iter == 1) {
                delta = min(delta,pnorm);
            }

/*           evaluate the function at x + p and calculate its norm. */

            iflag = fcnnn(n, &wa2[1], &wa4[1], 1, para);
            ++(*nfev);
            if (iflag < 0) {
                goto TERMINATE;
            }
            fnorm1 = enorm(n, &wa4[1]);

/*           compute the scaled actual reduction. */

            actred = -1.;
            if (fnorm1 < fnorm) {
                /* Computing 2nd power */
                d1 = fnorm1 / fnorm;
                actred = 1. - d1 * d1;
            }

/*           compute the scaled predicted reduction. */

            l = 1;
            for (i = 1; i <= n; ++i) {
                sum = 0.;
                for (j = i; j <= n; ++j) {
                    sum += r[l] * wa1[j];
                    ++l;
                }
                wa3[i] = qtf[i] + sum;
            }
            temp = enorm(n, &wa3[1]);
            prered = 0.;
            if (temp < fnorm) {
                /* Computing 2nd power */
                d1 = temp / fnorm;
                prered = 1. - d1 * d1;
            }

/*           compute the ratio of the actual to the predicted */
/*           reduction. */

            ratio = 0.;
            if (prered > 0.) {
                ratio = actred / prered;
            }

/*           update the step bound. */

            if (ratio < p1) {
                ncsuc = 0;
                ++ncfail;
                delta = p5 * delta;
            } else {
                ncfail = 0;
                ++ncsuc;
                if (ratio >= p5 || ncsuc > 1) {
                    /* Computing MAX */
                    d1 = pnorm / p5;
                    delta = max(delta,d1);
                }
                if (fabs(ratio - 1.) <= p1) {
                    delta = pnorm / p5;
                }
            }

/*           test for successful iteration. */

            if (ratio >= p0001) {

/*           successful iteration. update x, fvec, and their norms. */

                for (j = 1; j <= n; ++j) {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                    fvec[j] = wa4[j];
                }
                xnorm =enorm(n, &wa2[1]);
                fnorm = fnorm1;
                ++iter;
            }

/*           determine the progress of the iteration. */

            ++nslow1;
            if (actred >= p001) {
                nslow1 = 0;
            }
            if (jeval) {
                ++nslow2;
            }
            if (actred >= p1) {
                nslow2 = 0;
            }

/*           test for convergence. */

            if (delta <= xtol * xnorm || fnorm == 0.) {
                info = 1;
            }
            if (info != 0) {
                goto TERMINATE;
            }

/*           tests for termination and stringent tolerances. */

            if (*nfev >= maxfev) {
                info = 2;
            }
            /* Computing MAX */
            d1 = p1 * delta;
            if (p1 * max(d1,pnorm) <= epsmch * xnorm) {
                info = 3;
            }
            if (nslow2 == 5) {
                info = 4;
            }
            if (nslow1 == 10) {
                info = 5;
            }
            if (info != 0) {
                goto TERMINATE;
            }

/*           criterion for recalculating jacobian approximation */
/*           by forward differences. */

            if (ncfail == 2) {
                goto TERMINATE_INNER_LOOP;
            }

/*           calculate the rank one modification to the jacobian */
/*           and update qtf if necessary. */

            for (j = 1; j <= n; ++j) {
                sum = 0.;
                for (i = 1; i <= n; ++i) {
                    sum += fjac[i + j * fjac_dim1] * wa4[i];
                }
                wa2[j] = (sum - wa3[j]) / pnorm;
                wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
                if (ratio >= p0001) {
                    qtf[j] = sum;
                }
            }

/*           compute the qr factorization of the updated jacobian. */

            r1updt(n, n, &r[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing);
            r1mpyq(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
            r1mpyq(1, n, &qtf[1], 1, &wa2[1], &wa3[1]);

/*           end of the inner loop. */

            jeval = FALSE_;
        }
TERMINATE_INNER_LOOP:
        ;
/*        end of the outer loop. */

    }
TERMINATE:

/*     termination, either normal or user imposed. */

    if (iflag < 0) {
		info = iflag;
		return info;
    }
    if (nprint > 0 && iflag > 0)
	{
		fcnnn(n, &x[1], &fvec[1], 0, para);
		printf("%4d%7d%12.3e%12.3e\n", iter, *nfev, enorm(n, &fvec[1]), delta);
    }
    return info;

/*     last card of subroutine hybrd. */

} /* hybrd_ */

