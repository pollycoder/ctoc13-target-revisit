#ifndef _ORBITFUN_H_
#define _ORBITFUN_H_
#include <stdio.h>
#include"Constant.h"
#include "MyVector.h"

/**************************************************************************************************************************************/
/****************************************************�������������ֽǶȹ�ϵ**********************************************************/
/**************************************************************************************************************************************/
// E2f ����ƫ����Ǻ�ƫ������������
double E2f(int& flag, double E, double e);
// E2M ����ƫ����Ǻ�ƫ������ƽ�����
double E2M(int& flag, double E, double e);
// f2E ���������Ǻ�ƫ������ƫ�����
double f2E(int& flag, double f, double e);
// M2E ����ƽ����Ǻ�ƫ������ƫ�����
double M2E(int& flag, double M, double e, int MaxIter = 100, double epsilon = 1.0e-14);
// f0dt2ft ���ݳ�ʼ�����Ǻ��ݻ�ʱ��������������
double f0dt2ft(int& flag, double f0, double dt, double a, double e, double mu, int MaxIter = 100, double epsilon = 1.0e-14);
// f0ft2dt ���ݳ�ʼ�����Ǻ��������������ݻ�ʱ��
//double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441800e+14, int MaxIter = 100, double epsilon = 1.0e-14);
double f0ft2dt(int& flag, double f0, double ft, double a, double e, double mu = 3.98600441500e+14);
/**************************************************************************************************************************************/
/*******************************************************���ֵ��������Ƕȹ�ϵ*********************************************************/
/**************************************************************************************************************************************/
double L0dt2Lt(int& flag, double L0, double dt, const double* ee, double mu);
double L0Lt2dt(int& flag, double L0, double Lt, const double* ee, double mu);
/**************************************************************************************************************************************/
/*********************************************������������ֱ�����ꡢ�Ľ����ֵ���������ת��*****************************************/
/**************************************************************************************************************************************/
// coe2rv ���ݾ��������������ֱ������ϵ�µ�λ�ú��ٶȷ���
void coe2rv(int& flag, double* rv, const double* coe, double mu);
// rv2coe ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ����󾭵�������
void rv2coe(int& flag, double* coe, const double* RV, double mu);
// coe2ee ���ݾ�����������Ľ����ֵ����������Թ�����180������
void coe2ee(int& flag, double* ee, const double* coe, double mu);
// ee2coe ���ݸĽ����ֵ�����󾭵�������
void ee2coe(int& flag, double* coe, const double* ee, double mu);
// ee2rv ���ݸĽ����ֵ�����������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
void ee2rv(int& flag, double* rv, const double* ee, double mu);
// rv2ee ���ݵ��Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ�����Ľ����ֵ����������Թ�����180��ʱ����
void rv2ee(int& flag, double* ee, const double* RV, double mu);
/**************************************************************************************************************************************/
/********************************����ع�ϵ��λ��R0��Ϊ������ĵ�λ����za������ཻʱ����ľ�γ��**************************************/
/**************************************************************************************************************************************/
//���������ڵع�ϵ�µ�λ��(x,y,z)(m),�����ڵع�ϵ�µĵ�λ�����������ָ������ľ���(��),γ��(��),�߶�(m)
MyVec UAxis2Geod(MyVec R0, MyVec za);
/**************************************************************************************************************************************/
/*************************************************�ع�����ϵ��������ϵ��ת��*********************************************************/
/**************************************************************************************************************************************/
// Geo2Fix �������Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)��ع�����X, Y, Z,(m)ת��.
void Geo2Fix(int& flag, double* Fix, const double* Geo);
// Fix2Geo �ع�����X, Y, Z,(m)��������Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)ת��.
void Fix2Geo(int& flag, double* Geo, const double* Fix);
//�ع�����X, Y, Z,(m)��������Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)ת��.
//http://images.google.cn/imgres?imgurl=http://upload.wikimedia.org/math/a/2/2/a227314ce288bfdf0ad8099b6847aa97.png&imgrefurl=
//http://en.wikipedia.org/wiki/Geodetic_system&h=422&w=456&sz=12&hl=zh-CN&start=18&um=1&tbnid=NCbrKim-hJNzIM:&tbnh=118&tbnw=
//128&prev=/images%3Fq%3Dgeodetic%2Bgeocentric%2Blongitude%26um%3D1%26complete%3D1%26hl%3Dzh-CN%26newwindow%3D1%26sa%3DN
MyVec Fixed2Geodetic(const MyVec& FixedR);
//�������Altitude(m), Geodetic Latitude(rad), Geodetic Longitude(rad)��ع�����X, Y, Z,(m)ת��.
//http://images.google.cn/imgres?imgurl=http://upload.wikimedia.org/math/a/2/2/a227314ce288bfdf0ad8099b6847aa97.png&imgrefurl=
//http://en.wikipedia.org/wiki/Geodetic_system&h=422&w=456&sz=12&hl=zh-CN&start=18&um=1&tbnid=NCbrKim-hJNzIM:&tbnh=118&tbnw=
//128&prev=/images%3Fq%3Dgeodetic%2Bgeocentric%2Blongitude%26um%3D1%26complete%3D1%26hl%3Dzh-CN%26newwindow%3D1%26sa%3DN
MyVec Geodetic2Fixed(const MyVec& Geodetic);
/**************************************************************************************************************************************/
/************************************************��֪��ֵ��ʱ�䣬��ĩ��״̬************************************************************/
/**************************************************************************************************************************************/
//���ݳ�ʼʱ��״̬coe0��ĩ��ʱ��dt��״̬coe1���������ƽ���������ɹ�,flag����1
void coe02coef(int& flag, double* coe1, const double* coe0, double dt, double mu = 3.98600441500e+14);
//���ݳ�ʼʱ��t0��״̬rv0��ĩ��ʱ��t1��״̬rv1���������ƽ���������ɹ�,flag����1
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu = 3.98600441500e+14);
//���ݳ�ʼʱ��״̬ee0��ĩ��ʱ��dt��״̬ee1���������ƽ���������ɹ�,flag����1
void ee02eef(int& flag, double* ee1, const double* ee0, double dt, double mu = 3.98600441500e+14);
//�������,���ݳ�ʼʱ��t0��״̬rv0��ĩ��ʱ��t1��״̬rv1���������ƽ���������ɹ�,flag����1
void rv02rvf_file(int& flag, double* rv1, const double* rv0, double dt, FILE* fp, double mu = 3.98600441500e+14);

/**************************************************************************************************************************************/
/*********************************************************����j2�㶯�Ĺ������ת��*****************************************************/
/**************************************************************************************************************************************/
void coe2soe(const double* coe, double* soe);
void soe2coe(const double* soe, double* coe);
//����Ϊƽ�����������������Ϊ˲ʱ�������
void M2O(const double* soe_m, double* soe_o, double j2 = J2);
void O2M(const double* soe_o, double* soe_m, double j2 = J2);
int myfun(int n, const double* x, double* fvec, int iflag, const double* para);
/**************************************************************************************************************************************/
/************************************************************����j2�㶯�Ĺ������******************************************************/
/**************************************************************************************************************************************/
void j2mcoe02mcoef(const double* me0, const double dt, double* mef);
void j2ocoe02ocoef(double* ocoef, const double* ocoe0, const double dt);
void j2rv02rvf(const double* rv0, const double dt, double* rvf);

int dynODE_j2(double t, const double* y, double* yp, const double* para);
/**************************************************************************************************************************************/
/*********************************************************����j2��Lambert������ƽ�****************************************************/
/**************************************************************************************************************************************/
//void j2Lambert(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
			  //const double* rv0, const double* rv1, double t, double mu=3.98600441800e+14);

/**************************************************************************************************************************************/
/************************************************************lambert����***************************************************************/
/**************************************************************************************************************************************/
//Lambert�������
void lambert(double* v1, double* v2, double& a, double& e, const double* R1, const double* R2,
	double tf, const double* unith, int& flag, double mu = 3.98600441500e+14, int way = 0, int N = 0, int branch = 0,
	int Maxiter = 60, double tol = 1.0e-12);
//���ݳ�ʼ״̬rv0��ĩ��״̬rv1��ת��ʱ��t���������ٶ���С��˫����ת�ƹ����
//���أ�dv0,��ʼ����ʸ��;dv1��ĩ������ʸ��;Mdv0��dv0�ķ�ֵ��Mdv1,dv1�ķ�ֵ��N������Ȧ��;branch�����ŷ�֦;flag��1��ʾ����ɹ�
void LambEval(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double mu = 3.98600441500e+14);

void LambEval2(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double GM);
/**************************************************************************************************************************************/
/*********************************************************����j2��Lambert������ƽ�****************************************************/
/**************************************************************************************************************************************/

int Fun_J2Lambert_Approximation(int n, const double* x, double* fvec, int iflag, const double* para);

void J2Lambert_Approximation1(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double mu);

void J2Lambert_Approximation2(int& flag, double* dv0, double* dv1, double& Mdv0, double& Mdv1, int& N, int& branch,
	const double* rv0, const double* rv1, double t, double mu);
/**************************************************************************************************************************************/
/*********************************************************����j2��Lambert���⾫ȷ��****************************************************/
/**************************************************************************************************************************************/

void j2mcoe02mcoef_eJ2(const double* me0, const double dt, const double J2, double* mef);


#endif