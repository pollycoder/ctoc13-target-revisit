#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include "stdio.h"

namespace TigerIntegrator{
//Define a header file for fast integration
typedef int (*ODE45Fun)(double, const double *, double *, const double *);  
typedef void (*RKF78Fun)(double, const double *, double *, const double *); 
typedef int (*EventFun)(double t, const double *x, const double *eventPara, double &fval, double &fgrad, int calcflag);
class ODE45Integrator{
private:
	int dim;
	double RelTol;
	double *AbsTol;
	double *work;
	int MaxIter;
	int NumPoint;
	double EventTol;
    double MaxStep;
    double InitialStep;
public:
	ODE45Integrator(int n);
	int Integrate(ODE45Fun, double*, const double*, double, double, FILE *);
	int Integrate(ODE45Fun, double *x, const double *para, EventFun, const double *eventPara, int &eventFlag, const int eventDir, 
double t0, double &tf, FILE *);
	int Integrate(ODE45Fun, double *x, const double *para, EventFun, const double *eventPara, int &eventFlag, const int eventDir,
		double t0, double &tf, FILE *, double, bool output = false);
	void setDim(int n);
	void setRelTol(double tol);
	void setAbsTol(double *tol);
	void setAbsTol(double tol);
	void setEventTol(double tol);
	void setMaxIter(int);
    void setMaxStep(double);
    void setInitialStep(double);
	~ODE45Integrator();
};
class RKF78Integrator{
	private:
	int dim;
	double RelTol;
	double *AbsTol;
	double EventTol;
	int RKFmaxnfe;
	double *work;
public:
	RKF78Integrator(int n);
	int Integrate(RKF78Fun , double *, const double *, double , double , FILE *);
	int Integrate(RKF78Fun, double *x, const double *para, EventFun, const double *eventPara, int eventFlag, const int eventDir, double &t0, double tf, FILE *);
	void setDim(int n);
	void setRelTol(double tol);
	void setAbsTol(double *tol);
	void setAbsTol(double tol);
	void setEventTol(double tol);
	void setMaxFe(int);
	~RKF78Integrator();
	
};
}
typedef TigerIntegrator::ODE45Integrator ODE45;
#endif
