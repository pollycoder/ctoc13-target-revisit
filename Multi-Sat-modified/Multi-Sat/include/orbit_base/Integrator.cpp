#include "Integrator.h"
#include "ODE45.h"
#include "RungeKutta.h"
#include <iostream>
#include "stdlib.h"

//using namespace ODE45_;
namespace TigerIntegrator{

//Init
ODE45Integrator::ODE45Integrator(int n){
	dim = n;
	RelTol = 1e-10;
	EventTol = 1e-10;
	AbsTol = new double[n];
	work = new double[10*n];
	if(AbsTol == NULL || work == NULL){
		cout << "ODE allocate error\n";
		exit(0);
	}
	for(int i = 0;i < n;i++)	AbsTol[i] = 1e-10;
	MaxIter = 100000;
    MaxStep = -1;
    InitialStep = -1;
	
}
//int ODE45Integrator::Integrate(ODE45Fun fun, double *state, const double *para, double tin, double tout, FILE *fp){
//	int odeflag = ode45(fun, state, para, tin, tout, dim, NumPoint, work, AbsTol, RelTol, MaxIter, 0, MaxStep, InitialStep, fp);
//	return odeflag;
//}
//
//int ODE45Integrator::Integrate(ODE45Fun fun, double *state, const double *para, EventFun eventfun, const double *eventPara, int &eventFlag, const int eventDir, double t0, double &tf, FILE *fp){
//	int odeflag = ode45(fun, eventfun, state, para, eventPara, EventTol, eventFlag, eventDir, t0, tf, dim, NumPoint, work, AbsTol, RelTol, MaxIter, 0, -1, -1, fp);
//	return odeflag;
//}
//
//int ODE45Integrator::Integrate(ODE45Fun fun, double *state, const double *para, EventFun eventfun, const double *eventPara, int &eventFlag, const int eventDir, double t0, double &tf, FILE *fp, double InitialStep, bool output) {
//	int odeflag = ode45(fun, eventfun, state, para, eventPara, EventTol, eventFlag, eventDir, t0, tf, dim, NumPoint, work, AbsTol, RelTol, MaxIter, 0, -1, InitialStep, fp, output);
//	return odeflag;
//}


void ODE45Integrator::setDim(int n){
	if(dim < n){
		delete[] AbsTol;
		delete[] work;
		AbsTol = new double[n];
		work = new double[10*n];
		if(AbsTol == NULL || work == NULL){
			cout << "ODE allocate error\n";
			exit(0);
		}
		for(int i = 0;i < n;i++)	AbsTol[i] = 1e-10;
	}
	dim = n;
}
void ODE45Integrator::setRelTol(double tol){
	RelTol = tol;
}
void ODE45Integrator::setEventTol(double tol){
	EventTol = tol;
}
void ODE45Integrator::setAbsTol(double *tol){
	for(int i = 0;i < dim;i++)	AbsTol[i] = tol[i];
}
void ODE45Integrator::setAbsTol(double tol){
	for(int i = 0;i < dim;i++)	AbsTol[i] = tol;
}
void ODE45Integrator::setMaxIter(int maxiter){
	MaxIter = maxiter;
}
void ODE45Integrator::setMaxStep(double step){
	MaxStep = step;
}
void ODE45Integrator::setInitialStep(double step){
	InitialStep = step;
}
ODE45Integrator::~ODE45Integrator(){
	delete[] AbsTol;
	delete[] work;
}

RKF78Integrator::RKF78Integrator(int n){
	dim = n;
	RelTol = 1e-10;
	EventTol = 1e-10;
	AbsTol = new double[n];
	work = new double[14*n];
	if(AbsTol == NULL || work == NULL){
		cout << "ODE allocate error\n";
		exit(0);
	}
	for(int i = 0;i < n;i++)	AbsTol[i] = 1e-10;
	RKFmaxnfe = 10000;
}
int RKF78Integrator::Integrate(RKF78Fun fun, double *state, const double *para, double tin, double tout, FILE *fp){
	int RKFiflag = 1;
	int steps = RK::RKF78(fun, state, para, tin, tout, dim, RelTol, AbsTol, RKFiflag, work, RKFmaxnfe, fp);
	return RKFiflag;
}
/*
int RKF78Integrator::Integrate(RKF78Fun fun, double *state, const double *para, EventFun eventfun, const double *eventPara, int eventFlag, const int eventDir, double &t0, double tout, FILE *fp){
	int RKFiflag = 1;
	int steps = RK::RKF78(fun, eventfun, state, para, eventPara, EventTol, eventFlag, eventDir, t0, tout, dim, RelTol, AbsTol, RKFiflag, work, RKFmaxnfe, fp);
	return RKFiflag;
	
}
*/
void RKF78Integrator::setDim(int n){
	if(dim < n){
		delete[] AbsTol;
		delete[] work;
		AbsTol = new double[n];
		work = new double[14*n];
		if(AbsTol == NULL || work == NULL){
			cout << "ODE allocate error\n";
			exit(0);
		}
		for(int i = 0;i < n;i++)	AbsTol[i] = 1e-10;
	}
	dim = n;
}
void RKF78Integrator::setRelTol(double tol){
	RelTol = tol;
}
void RKF78Integrator::setAbsTol(double *tol){
	for(int i = 0;i < dim;i++)	AbsTol[i] = tol[i];
}
void RKF78Integrator::setAbsTol(double tol){
	for(int i = 0;i < dim;i++)	AbsTol[i] = tol;
}
void RKF78Integrator::setMaxFe(int maxfe){
	RKFmaxnfe = maxfe;
}
RKF78Integrator::~RKF78Integrator(){
	delete[] AbsTol;
	delete[] work;
}
}
