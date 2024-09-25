#include "shootfunc.h"


#include "Mat3D.h"
#include "OrbitFun.h"

/****************************************************************************
* ������   : shoot_func()
* ��  ��   : ��к���
* �� ��    : x: 0:dvx 1:dvy 2:dvz
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-11:�������tf�µĵع�ϵ��λʸ�� 12:tf
* ��    �� : fevc: 0:0 1:0.0 2:0 3:0
****************************************************************************/
int shoot_func(int n, const double* x, double* fevc, int iflag, const double* para)
{
	double rv0[6] = { para[0],para[1],para[2], para[3] + x[0],para[4] + x[1], para[5] + x[2] };
	double t0 = para[8], tf = para[12];
	//double para_prop[2] = { para[6] * exp(-sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) * 1000.0 / 9.80665 / 300.0),para[7] };
	double para_prop[2] = { para[6], para[7] };// ���������嵼�µ������仯

	double node_RF_norm[3] = { para[9], para[10], para[11] };

	double* abstol = new double[6];
	double* newwork = new double[6 * 10];
	int num;
	for (int i = 0; i < 6; i++)  abstol[i] = 1e-10;
	ode45(RhsFnforODE45, rv0, para_prop, t0, tf, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1);
	//rv0��ʱΪ��ǰλ���ٶ�
	double norm2 = V_Norm2(rv0, 3);
	for (int i = 0; i < 3; i++) fevc[i] = node_RF_norm[i] * norm2 - rv0[i];

	delete[]abstol;
	delete[]newwork;
	return 1;
}

/****************************************************************************
* ������   : shoot_func_2()
* ��  ��   : ��к���
* �� ��    : x: 0:dv
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-10:����㾭γ�ȣ�rad�� 11:tf 12-14:dv�ķ���
* ��    �� : fevc: 0:���Ȳ�
****************************************************************************/
int shoot_func_2(int n, const double* x, double* fevc, int iflag, const double* para)
{
	double norm_dv[3] = { para[12],para[13] ,para[14] };
	double rv0[6] = { para[0],para[1],para[2], para[3] + x[0] * norm_dv[0],para[4] + x[0] * norm_dv[1], para[5] + x[0] * norm_dv[2] };
	double t0 = para[8], tf = para[11];
	double para_prop[2] = { para[6] * exp(-fabs(x[0] * 1000.0) / g0 / Isp) ,para[7] };
	double lumda = para[9];
	double phi = para[10];

	double* abstol = new double[6];
	double* newwork = new double[6 * 10];
	int num;

	for (int i = 0; i < 6; i++)  abstol[i] = 1e-10;
	ode45(RhsFnforODE45, rv0, para_prop, t0, tf, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1); //rv0�ǵ�ǰʱ�̵�λ���ٶ�
	int flag_coe;
	double coe[6];
	rv2coe(flag_coe, coe, rv0, mu * 1e-9);  //���������

	double u0 = coe[4] + coe[5]; u0 = NiceAngle(u0);//���//��Χ��0 ��2pi
	double u1 = asin(sin(phi) / sin(coe[2]));
	double u2 = DPI - u1;
	//��Χ�������pi 
	double du1 = NiceAngle(u1 - u0);
	double du2 = NiceAngle(u2 - u0);
	if (du1 > DPI) du1 = du1 - 2 * DPI;
	if (du2 > DPI) du2 = du2 - 2 * DPI;

	//startʱ�̵�λ���ٶ�
	double r[3] = { rv0[0],rv0[1],rv0[2] };
	double v[3] = { rv0[3],rv0[4],rv0[5] };

	//Ѱ�����µ�γ�ȵ�һ����Ŀ���غϵ�����ʱ��
	double _omega_ = V_Norm2(v, 3) / V_Norm2(r, 3);
	double t1 = du1 / _omega_;
	double t2 = du2 / _omega_;
	double t = fabs(t1) < fabs(t2) ? t1 : t2;
	double  t_access = tf + t;

	ode45(RhsFnforODE45, rv0, para_prop, tf, t_access, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1); //rv0�ǵ�ǰʱ�̵�λ���ٶ�

	MyVec RF(rv0[0], rv0[1], rv0[2]);
	double alpha_G = alpha_G0 + OmegaEarth * t_access; //rad
	Mat3D Rz(cos(alpha_G), sin(alpha_G), 0, -sin(alpha_G), cos(alpha_G), 0, 0, 0, 1);
	RF = Rz * RF;  //������ڵع�ϵ������

	double lumda0 = atan2(RF[1], RF[0]);

	//rv0��ʱΪ��ǰλ���ٶ�
	fevc[0] = lumda0 - lumda;
	//para[11] = t_access;
	delete[]abstol;
	delete[]newwork;
	return 1;
}


void acc_validate(const double* x, const double* para,  double& t_access, double * rvf)
{
	double norm_dv[3] = { para[12],para[13] ,para[14] };
	double rv0[6] = { para[0], para[1],para[2], para[3] + x[0] * norm_dv[0],para[4] + x[0] * norm_dv[1], para[5] + x[0] * norm_dv[2] };
	double t0 = para[8], tf = para[11];
	double para_prop[2] = { para[6] * exp(-fabs(x[0] * 1000.0) / g0 / Isp) ,para[7] };
	double lumda = para[9];
	double phi = para[10];

	double * abstol = new double[6];
	double * newwork = new double[6 * 10];
	int num;

	for (int i = 0; i < 6; i++)  abstol[i] = 1e-10;
	ode45(RhsFnforODE45, rv0, para_prop, t0, tf, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1); //rv0�ǵ�ǰʱ�̵�λ���ٶ�
	int flag_coe;
	double coe[6];
	rv2coe(flag_coe, coe, rv0, mu * 1e-9);  //���������

	double u0 = coe[4] + coe[5]; u0 = NiceAngle(u0);//���//��Χ��0 ��2pi
	double u1 = asin(sin(phi) / sin(coe[2]));
	double u2 = DPI - u1;
	//��Χ�������pi 
	double du1 = NiceAngle(u1 - u0);
	double du2 = NiceAngle(u2 - u0);
	if (du1 > DPI) du1 = du1 - 2 * DPI;
	if (du2 > DPI) du2 = du2 - 2 * DPI;

	//startʱ�̵�λ���ٶ�
	double r[3] = { rv0[0],rv0[1],rv0[2] };
	double v[3] = { rv0[3],rv0[4],rv0[5] };

	//Ѱ�����µ�γ�ȵ�һ����Ŀ���غϵ�����ʱ��
	double _omega_ = V_Norm2(v, 3) / V_Norm2(r, 3);
	double t1 = du1 / _omega_;
	double t2 = du2 / _omega_;
	double t = fabs(t1) < fabs(t2) ? t1 : t2;
	t_access = tf + t;

	ode45(RhsFnforODE45, rv0, para_prop, tf, t_access, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1); //rv0�ǵ�ǰʱ�̵�λ���ٶ�

	memcpy(rvf, rv0, 6 * sizeof(double));
	delete[] abstol;
	delete[] newwork;
	
}
/****************************************************************************
* ������   : solve_shoot_dv()
* ��  ��   : ͨ����з�������ڴ�������ʵdv��tf
* ������� :  RVΪ���ǰλ���ٶȣ�mΪû�б��ǰ�����������Ϊ���������������
			  dvΪ��ʼԤ���dv�����Ϊ����dv����tfΪ��ʼԤ���tf�����ٸı�
****************************************************************************/
//int solve_shoot_dv(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double& tf)
//{
//	MyVec RF;
//	double alpha_G = alpha_G0 + OmegaEarth * tf; //rad
//	Mat3D Rz(cos(alpha_G), sin(alpha_G), 0, -sin(alpha_G), cos(alpha_G), 0, 0, 0, 1);
//	Mat3D Rzt = Rz.Trans();
//	RF = Rzt * point.PosRF;
//
//	double rf[3] = { RF[0], RF[1], RF[2] };
//	double normrf = V_Norm2(rf, 3);
//	double rf_norm[3] = { rf[0] / normrf, rf[1] / normrf, rf[2] / normrf };
//
//	double mass = m;
//
//	double para[13] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5],mass, 1.0, t0 , rf_norm[0],rf_norm[1],rf_norm[2],tf };
//	double x[3] = { dv[0],dv[1],dv[2] };
//	double fx[3];
//	MSolver A(3, 0);
//	//A.setMaxfevno(3); //��3��
//	int flag = A.Solve(shoot_func, x, fx, para);
//
//	for (int i = 0; i < 3; i++) dv[i] = x[i];
//	double dv_norm = V_Norm2(dv, 3);
//	//m = mass * exp(-fabs(dv_norm * 1000.0) / g0 / Isp);// ��ʱ���������嵼�µ������仯
//	return 1;
//}

/****************************************************************************
* ������   : solve_trans_dv()
* ��  ��   : ͨ����з�������ڴ�������ʵdv��tf
* ������� :  RVΪ���ǰλ���ٶȣ�mΪû�б��ǰ�����������Ϊ���������������
			  dvΪ��ʼԤ���dv�����Ϊ����dv����tfΪ��ʼԤ���tf���ı�
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-10:����㾭γ�ȣ�rad�� 11:tf 12-14:dv�ķ���
****************************************************************************/
//int solve_trans_dv(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double tf)
//{
//	double norm_dv[3] = { dv[0] / V_Norm2(dv,3), dv[1] / V_Norm2(dv,3), dv[2] / V_Norm2(dv,3) };
//	double mass = m;
//	double para[15] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5],mass, 1.0, t0 , point.Pos.longitude, point.Pos.latitude,tf, norm_dv[0],norm_dv[1],norm_dv[2] };
//	double x[1] = { V_Norm2(dv,3) };
//	double fx[1];
//	MSolver A(1, 0);
//	A.setMaxfevno(3); //��3��
//	int flag = A.Solve(shoot_func_2, x, fx, para);
//	if (fabs(fx[0]) > 1e-3)
//	{
//		return 0;
//	}
//
//
//	for (int i = 0; i < 3; i++) dv[i] = x[0] * norm_dv[i];
//	m = mass * exp(-fabs(x[0] * 1000.0) / g0 / Isp);
//	return 1;
//}
int solve_trans_dv(const double & longitude, const double & latitude, const double* RV, double& m, double* dv, const double t0, double tf)
{
	double norm_dv[3] = { dv[0] / V_Norm2(dv,3), dv[1] / V_Norm2(dv,3), dv[2] / V_Norm2(dv,3) };
	double mass = m;
	double para[15] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5],mass, 1.0, t0 , longitude, latitude,tf, norm_dv[0],norm_dv[1],norm_dv[2] };
	double x[1] = { V_Norm2(dv,3) };
	double fx[1];
	MSolver A(1, 0);
	A.setMaxfevno(20); //��3��
	int flag = A.Solve(shoot_func_2, x, fx, para);
	if (fabs(fx[0]) > 1e-3)
	{
		return 0;
	}


	for (int i = 0; i < 3; i++) dv[i] = x[0] * norm_dv[i];
	m = mass * exp(-fabs(x[0] * 1000.0) / g0 / Isp);
	return 1;
}

int acc_validate_total(const double& longitude, const double& latitude, const double* RV, double m0, double* dv, const double t0, double & tf, double *rvf, double & mf)
{
	double norm_dv[3] = { dv[0] / V_Norm2(dv,3), dv[1] / V_Norm2(dv,3), dv[2] / V_Norm2(dv,3) };
	double mass = m0;
	double para[15] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5], mass, 1.0, t0 , longitude, latitude,tf, norm_dv[0],norm_dv[1],norm_dv[2] };

	double x[1] = { V_Norm2(dv,3) };
	
	acc_validate(x, para, tf, rvf);

	mf = mass * exp(-fabs(x[0] * 1000.0) / g0 / Isp);
	return 1;
}
/****************************************************************************
* ������   : shoot_func_ToRSAT
* ��  ��   : ��к���,���м���
* �� ��    : x: 0:dvx 1:dvy 2:dvz
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-11:rvf 12:tf
* ��    �� : fevc: 0:0 1:0.0 2:0 3:0
****************************************************************************/
int shoot_func_ToRSAT(int n, const double* x, double* fevc, int iflag, const double* para)
{
	double rv0[6] = { para[0],para[1],para[2], para[3] + x[0],para[4] + x[1], para[5] + x[2] };
	double t0 = para[8], tf = para[12];
	double para_prop[2] = { para[6]*exp(-sqrt(x[0]* x[0]+ x[1] * x[1] + x[2] * x[2])*1000.0/9.80665/300.0),para[7] };
	//double para_prop[2] = { para[6],para[7] };
	double rvf[3] = { para[9], para[10], para[11] };

	double* abstol = new double[6];
	double* newwork = new double[6 * 10];
	int num;
	for (int i = 0; i < 6; i++)  abstol[i] = 1e-10;
	ode45(RhsFnforODE45, rv0, para_prop, t0, tf, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1);
	//rv0��ʱΪ��ǰλ���ٶ�
	for (int i = 0; i < 3; i++) fevc[i] = rvf[i] - rv0[i];

	delete[]abstol;
	delete[]newwork;
	return 1;
}

/****************************************************************************
* ������   : solve_shoot_dv_ToRSAT()
* ��  ��   : ͨ����з������ת�����м���λ�õ�׼ȷ���ٶ�����
* ������� :  RVΪ���ǰλ���ٶȣ�RfΪĿ��λ�ã�mΪû�б��ǰ�����������Ϊ���������������
			  dvΪ��ʼԤ���dv�����Ϊ����dv����tfΪ��ʼԤ���tf�����ٸı�
****************************************************************************/
int solve_shoot_dv_ToRSAT(bool& flag, const double* RV, double* Rf, double& m, double* dv, const double t0, double& tf)
{
	double mass = m;

	double para[13] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5],mass, 1.0, t0 , Rf[0],Rf[1],Rf[2],tf };
	double x[3] = { dv[0],dv[1],dv[2] };
	double fx[3];
	MSolver A(3, 0);
	A.setMaxfevno(20); //��10��
	A.Solve(shoot_func_ToRSAT, x, fx, para);

	for (int i = 0; i < 3; i++) dv[i] = x[i];
	double dv_norm = V_Norm2(dv, 3);
	m = mass * exp(-fabs(dv_norm * 1000.0) / g0 / Isp);
	if (fx[0] * fx[0] + fx[1] * fx[1] + fx[2] * fx[2] > 0.0001)
		flag = false;
	else
		flag = true;
	return 1;
}

int solve_shoot_dv_ToRSAT2(bool& flag, const double* RV, double* Rf, double& m, double* dv, const double t0, double& tf)
{
	double mass = m;

	double para[13] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5],mass, 1.0, t0 , Rf[0],Rf[1],Rf[2],tf };
	double x[3] = { dv[0],dv[1],dv[2] };
	double fx[3];

	for (int iter = 0; iter < 100; iter++)
	{
		para[7] = 0.01 * (iter + 1);
		MSolver A(3, 0);
		//A.setMaxfevno(20); //��10��
		A.Solve(shoot_func_ToRSAT, x, fx, para);
	}

	for (int i = 0; i < 3; i++) dv[i] = x[i];
	double dv_norm = V_Norm2(dv, 3);
	m = mass * exp(-fabs(dv_norm * 1000.0) / g0 / Isp);
	if (fx[0] * fx[0] + fx[1] * fx[1] + fx[2] * fx[2] > 0.0001)
		flag = false;
	else
		flag = true;
	return 1;
}

int shoot_func3(int n, const double* x, double* fevc, int iflag, const double* para)
{
	double rv0[6] = { para[0],para[1],para[2], para[3] + x[0],para[4] + x[1], para[5] + x[2] };
	double t0 = para[8], tf = para[12];
	double para_prop[2] = { para[6] * exp(-sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) * 1000.0 / 9.80665 / 300.0),para[7] };
	double node_RF_norm[3] = { para[9], para[10], para[11] };

	double* abstol = new double[6];
	double* newwork = new double[6 * 10];
	int num;
	for (int i = 0; i < 6; i++)  abstol[i] = 1e-10;
	ode45(RhsFnforODE45, rv0, para_prop, t0, tf, 6, num, newwork, abstol, 1e-12, 0, 100000, 0.1);
	//rv0��ʱΪ��ǰλ���ٶ�
	double norm2 = 1.0;
	for (int i = 0; i < 3; i++) fevc[i] = node_RF_norm[i] * norm2 - rv0[i];

	delete[]abstol;
	delete[]newwork;
	return 1;
}

//int solve_shoot_dv3(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double& tf)
//{
//	MyVec RF;
//	double alpha_G = alpha_G0 + OmegaEarth * tf; //rad
//	Mat3D Rz(cos(alpha_G), sin(alpha_G), 0, -sin(alpha_G), cos(alpha_G), 0, 0, 0, 1);
//	Mat3D Rzt = Rz.Trans();
//	RF = Rzt * point.PosRF;
//
//	double rf[3] = { RF[0], RF[1], RF[2] };
//	double normrf = 1.0;
//	double rf_norm[3] = { rf[0] / normrf, rf[1] / normrf, rf[2] / normrf };
//
//	double mass = m;
//
//	double para[13] = { RV[0], RV[1],RV[2],RV[3],RV[4],RV[5],mass, 1.0, t0 , rf_norm[0],rf_norm[1],rf_norm[2],tf };
//	double x[3] = { dv[0],dv[1],dv[2] };
//	double fx[3];
//	MSolver A(3, 0);
//	//A.setMaxfevno(3); //��3��
//	int flag = A.Solve(shoot_func3, x, fx, para);
//
//	for (int i = 0; i < 3; i++) dv[i] = x[i];
//	double dv_norm = V_Norm2(dv, 3);
//	m = mass * exp(-fabs(dv_norm * 1000.0) / g0 / Isp);
//	return 1;
//}