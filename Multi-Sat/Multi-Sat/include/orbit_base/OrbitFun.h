#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_
#include <stdio.h>
#include"Constant.h"
#include "MyVector.h"

/**************************************************************************************************************************************/
/****************************************************经典轨道根数三种角度关系**********************************************************/
/**************************************************************************************************************************************/
// E2f 根据偏近点角和偏心率求真近点角
double E2f(int& flag, double E, double e);
// E2M 根据偏近点角和偏心率求平近点角
double E2M(int& flag, double E, double e);
// f2E 根据真近点角和偏心率求偏近点角
double f2E(int& flag, double f, double e);
// M2E 根据平近点角和偏心率求偏近点角
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = 1.0e-14);
// f0dt2ft 根据初始真近点角和演化时间求最终真近点角
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter = 100, double epsilon = 1.0e-14);
// f0ft2dt 根据初始真近点角和最终真近点角求演化时间
//double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441800e+14, int MaxIter = 100, double epsilon = 1.0e-14);
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441500e+14);
/**************************************************************************************************************************************/
/*******************************************************春分点轨道根数角度关系*********************************************************/
/**************************************************************************************************************************************/
double L0dt2Lt(int& flag, double L0, double dt, const double* ee, double mu);
double L0Lt2dt(int& flag, double L0, double Lt, const double* ee, double mu);
/**************************************************************************************************************************************/
/*********************************************经典轨道根数、直角坐标、改进春分点轨道根数的转换*****************************************/
/**************************************************************************************************************************************/
// coe2rv 根据经典轨道根数求惯性直角坐标系下的位置和速度分量
void coe2rv(int& flag, double* rv, const double* coe, double mu);
// rv2coe 根据地心惯性直角坐标系下的位置和速度分量求经典轨道根数
void rv2coe(int& flag, double* coe, const double* RV, double mu);
// coe2ee 根据经典轨道根数求改进春分点轨道根数，对轨道倾角180度奇异
void coe2ee(int& flag, double* ee, const double* coe, double mu);
// ee2coe 根据改进春分点根数求经典轨道根数
void ee2coe(int& flag, double* coe, const double* ee, double mu);
// ee2rv 根据改进春分点轨道根数求地心惯性直角坐标系下的位置和速度分量
void ee2rv(int& flag, double* rv, const double* ee, double mu);
// rv2ee 根据地心惯性直角坐标系下的位置和速度分量求改进春分点轨道根数，对轨道倾角180度时奇异
void rv2ee(int& flag, double* ee, const double* RV, double mu);
/**************************************************************************************************************************************/
/********************************计算地固系中位置R0处为出发点的单位向量za与地面相交时交点的经纬度**************************************/
/**************************************************************************************************************************************/
//根据卫星在地固系下的位置(x,y,z)(m),光轴在地固系下的单位分量求出光轴指向地面点的经度(°),纬度(°),高度(m)
MyVec UAxis2Geod(MyVec R0, MyVec za);
/**************************************************************************************************************************************/
/*************************************************地固坐标系与大地坐标系的转换*********************************************************/
/**************************************************************************************************************************************/
// Geo2Fix 大地坐标Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)向地固坐标X, Y, Z,(m)转换.
void Geo2Fix(int& flag, double* Fix, const double* Geo);
// Fix2Geo 地固坐标X, Y, Z,(m)向大地坐标Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)转换.
void Fix2Geo(int& flag, double* Geo, const double* Fix);
//地固坐标X, Y, Z,(m)向大地坐标Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)转换.
//http://images.google.cn/imgres?imgurl=http://upload.wikimedia.org/math/a/2/2/a227314ce288bfdf0ad8099b6847aa97.png&imgrefurl=
//http://en.wikipedia.org/wiki/Geodetic_system&h=422&w=456&sz=12&hl=zh-CN&start=18&um=1&tbnid=NCbrKim-hJNzIM:&tbnh=118&tbnw=
//128&prev=/images%3Fq%3Dgeodetic%2Bgeocentric%2Blongitude%26um%3D1%26complete%3D1%26hl%3Dzh-CN%26newwindow%3D1%26sa%3DN
MyVec Fixed2Geodetic(const MyVec& FixedR);
//大地坐标Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)向地固坐标X, Y, Z,(m)转换.
//http://images.google.cn/imgres?imgurl=http://upload.wikimedia.org/math/a/2/2/a227314ce288bfdf0ad8099b6847aa97.png&imgrefurl=
//http://en.wikipedia.org/wiki/Geodetic_system&h=422&w=456&sz=12&hl=zh-CN&start=18&um=1&tbnid=NCbrKim-hJNzIM:&tbnh=118&tbnw=
//128&prev=/images%3Fq%3Dgeodetic%2Bgeocentric%2Blongitude%26um%3D1%26complete%3D1%26hl%3Dzh-CN%26newwindow%3D1%26sa%3DN
MyVec Geodetic2Fixed(const MyVec& Geodetic);
/**************************************************************************************************************************************/
/************************************************已知初值和时间，求末端状态************************************************************/
/**************************************************************************************************************************************/
//根据初始时刻状态coe0求末端时刻dt的状态coe1，按二体推进。若计算成功,flag返回1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu = 3.98600441500e+14);
//根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu = 3.98600441500e+14);
//根据初始时刻状态ee0求末端时刻dt的状态ee1，按二体推进。若计算成功,flag返回1
void ee02eef(int& flag, double* ee1, const double* ee0, double dt, double mu = 3.98600441500e+14);
//积分求解,根据初始时刻t0的状态rv0求末端时刻t1的状态rv1，按二体推进。若计算成功,flag返回1
void rv02rvf_file(int& flag, double* rv1, const double* rv0, double dt, FILE* fp, double mu = 3.98600441500e+14);

/**************************************************************************************************************************************/
/*********************************************************考虑j2摄动的轨道根数转换*****************************************************/
/**************************************************************************************************************************************/
void coe2soe(const double* coe, double* soe);
void soe2coe(const double* soe, double* coe);
//输入为平均经典轨道根数，输出为瞬时轨道根数
void M2O(const double* soe_m, double* soe_o, double j2 = J2);
void O2M(const double* soe_o, double* soe_m, double j2 = J2);
int myfun(int n, const double* x, double* fvec, int iflag, const double* para);
/**************************************************************************************************************************************/
/************************************************************考虑j2摄动的轨道递推******************************************************/
/**************************************************************************************************************************************/
void j2mcoe02mcoef(const double* me0, const double dt, double* mef);
void j2ocoe02ocoef(double* ocoef, const double* ocoe0, const double dt);
void j2rv02rvf(const double* rv0, const double dt, double* rvf);

int dynODE_j2(double t, const double* y, double* yp, const double* para);
/**************************************************************************************************************************************/
/*********************************************************考虑j2的Lambert问题近似解****************************************************/
/**************************************************************************************************************************************/
//void j2Lambert(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
			  //const double* rv0, const double* rv1, double t, double mu=3.98600441800e+14);

/**************************************************************************************************************************************/
/************************************************************lambert问题***************************************************************/
/**************************************************************************************************************************************/
//Lambert问题求解
void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
	double tf, const double* unith, int& flag, double mu = 3.98600441500e+14, int way = 0, int N = 0, int branch = 0,
	int Maxiter = 60, double tol = 1.0e-12);
//根据初始状态rv0和末端状态rv1及转移时间t，求特征速度最小的双脉冲转移轨道。
//返回：dv0,初始脉冲矢量;dv1，末端脉冲矢量;Mdv0，dv0的幅值；Mdv1,dv1的幅值；N，最优圈次;branch，最优分枝;flag，1表示计算成功
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double mu = 3.98600441500e+14);

void LambEval2(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double GM);
/**************************************************************************************************************************************/
/*********************************************************考虑j2的Lambert问题近似解****************************************************/
/**************************************************************************************************************************************/

int Fun_J2Lambert_Approximation(int n, const double* x, double* fvec, int iflag, const double* para);

void J2Lambert_Approximation1(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double mu);

void J2Lambert_Approximation2(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double mu);
/**************************************************************************************************************************************/
/*********************************************************考虑j2的Lambert问题精确解****************************************************/
/**************************************************************************************************************************************/

void j2mcoe02mcoef_eJ2(const double* me0, const double dt, const double J2, double* mef);


#endif