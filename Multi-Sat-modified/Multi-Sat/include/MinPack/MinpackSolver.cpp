#include "MinpackSolver.h"
#include <iostream>
#include "stdlib.h"

namespace MinpackSolver{
//Init
	MinpackSolver::MinpackSolver(const int n, const int nprint_) {
	dim = n;
	xtol = 1e-10;
	nprint = nprint_;
	maxfevno = 2000;
	fvecnorm = 0;
	work = new double[n*(3*n + 13)/2];
	if(work == NULL){
		std::cout << "MinpackSolver allocate error\n";
		exit(0);
	}
	fjac = NULL;
}
int MinpackSolver::Solve(DBSolveFun fun, double *x, double *fvec, const double *para) {
	int flag = hybrd1(fun, dim, x, fvec, para, work, xtol, nprint, maxfevno);
	return flag;
}
int MinpackSolver::Solve(SolveJac fun, double *x, double *fvec, const double *para){
	if(fjac == NULL){
		fjac = new double[dim*dim];
		if(fjac == NULL){
			std::cout << "MinpackSolver allocate error\n";
			exit(0);
		}
	}
	int flag = hybrj1(fun, dim, x, fvec, para, fjac, work, xtol, nprint, maxfevno);
	return flag;
}

void MinpackSolver::setDim(const int n) {
	if(dim < n){
		delete[] work;
		work = new double[n*(3*n + 13)/2];
		if(work == NULL){
			std::cout << "MinpackSolver allocate error\n";
			exit(0);
		}
	}
	dim = n;
}

void MinpackSolver::setXtol(const double tol) {
	xtol = tol;
}
void MinpackSolver::setNprint(const int nprint_) {
	nprint = nprint_;
}
void MinpackSolver::setMaxfevno(const int maxfevno_) {
	maxfevno = maxfevno_;
}
MinpackSolver::~MinpackSolver(){
	delete[] work;
	if(fjac)
		delete[] fjac;
}
}
