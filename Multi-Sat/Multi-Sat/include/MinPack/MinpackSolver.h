#ifndef MINPACKSOLVER_H
#define MINPACKSOLVER_H
#include "stdio.h"
#include "cminpack.h"

namespace MinpackSolver{
//Define a header file for fast integration
typedef int(*DBSolveFun)(int n, const double *x, double *fvec, int iflag, const double* para);
typedef int (*SolveJac)(int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, const double* para);
class MinpackSolver{
private:
	int dim;
	double *work;
	double xtol;
	int nprint;
	int maxfevno;
	double fvecnorm;
	double *fjac;
public:
	MinpackSolver(const int n, const int nprint = 1);
	int Solve(DBSolveFun, double*, double*, const double*);
	int Solve(SolveJac, double *, double *, const double *);
	void setDim(const int n);
	void setXtol(const double tol);
	void setNprint(const int nprint);
	void setMaxfevno(const int maxfevno);
	~MinpackSolver();
};

}
typedef MinpackSolver::MinpackSolver MSolver;
#endif
