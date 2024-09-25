#pragma once

#include <vector>
#include "ODE45.h"


int Propagate(double mass, const double* rv0, double t1, double* rv1);

void PropagateRkf78(double mass, bool is_drag, const double* rv0, double t0,double t1, double* rv1);

class OrbitPropagator
{
public:
	void SetParam(double mass, bool is_drag);  // 先设参数、启动积分器
	void FlyTo(double tout, double *rv);  // 轨道外推
	int iflag() { return iflag_; }  // 如有需要，查看积分是否成功

private:
	long int iflag_;
	double t_;
	long int neqn_ = 6;
	double relerr_ = 1e-10;
	double abserr_ = 1e-10;
	double work_[226];
	long int iwork_[5];
};

int RhsFnforODE45(double t, const double* rv, double* drv, const double* para);
void PropagateODE45(double mass, bool is_drag, const double* rv0, double t0, double t1, double* rv1);