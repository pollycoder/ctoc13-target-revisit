#ifndef _J2_LAMBERT_
#define _J2_LAMBERT_

#include "MinpackSolver.h"
#include "Lambert.h"
#include "J2propagation.h"
#include "OrbitMath.h"

// Shooting function adapting to MinPack
// x[6]: v1[3]
// fevc: ||r2[3] - r2_real[3]||
// para[7]: r1[3], r2[3], tof
int shooting_func_lambert(int n, const double* x, double* fevc, int iflag, const double* para);
int solve_lambert(const double* R1, const double* R2, const double& tof, double* v1vec, double* v2vec);


// Lambert solver for problems with J2 perturbation in less than 2 days
// Use two-body Lambert solver to get initial value, then shoot to the precise result
void J2Lambert_short(double* v1vec, double* v2vec, double& a, double& e, const double* R1, const double* R2, const double& TOF,
   const double& mu, int way, int N, int branch, int Maxiter, double tol);

// Test J2 Lambert
void test_lambert();


#endif