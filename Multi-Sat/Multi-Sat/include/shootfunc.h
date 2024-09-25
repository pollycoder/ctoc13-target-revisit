/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: shootfunc.h
* ���ݼ�������еõ���ʵ�ٶ�����
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2020-10-04    ����      �����ļ�
****************************************************************************/

#ifndef _SHOOTFUNC_H
#define _SHOOTFUNC_H

#include "MinpackSolver.h"
#include "propagate.h"
#include "Constant.h"
#include "OrbitMath.h"
//#include "MyPoint.h"

using namespace MinpackSolver;

/****************************************************************************
* ������   : shoot_func()
* ��  ��   : ��к���
* �� ��    : x: 0:dvx 1:dvy 2:dvz
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-11:�������tf�µĵع�ϵ��λʸ�� 12:tf
* ��    �� : fevc: 0:0 1:0.0 2:0 3:0
****************************************************************************/
int shoot_func(int n, const double* x, double* fevc, int iflag, const double* para);

/****************************************************************************
* ������   : shoot_func_2()
* ��  ��   : ��к���
* �� ��    : x: 0:dv
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-10:����㾭γ�ȣ�rad�� 11:tf 12-14:dv�ķ���
* ��    �� : fevc: 0:���Ȳ�
****************************************************************************/
int shoot_func_2(int n, const double* x, double* fevc, int iflag, const double* para);

void acc_validate(const double* x, const double* para, double& t_access, double* rvf);
int acc_validate_total(const double& longitude, const double& latitude, const double* RV, double m0, double* dv, const double t0, double& tf, double* rvf, double& mf);

/****************************************************************************
* ������   : solve_shoot_dv()
* ��  ��   : ͨ����з�������ڴ�������ʵdv��tf
* ������� :  RVΪ���ǰλ���ٶȣ�mΪû�б��ǰ�����������Ϊ���������������
              dvΪ��ʼԤ���dv�����Ϊ����dv����tfΪ��ʼԤ���tf�����ٸı�
****************************************************************************/
//int solve_shoot_dv(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double& tf);

/****************************************************************************
* ������   : solve_trans_dv()
* ��  ��   : ͨ����з�������ڴ�������ʵdv��tf
* ������� :  RVΪ���ǰλ���ٶȣ�mΪû�б��ǰ�����������Ϊ���������������
			  dvΪ��ʼԤ���dv�����Ϊ����dv����tfΪ��ʼԤ���tf���ı�
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-10:����㾭γ�ȣ�rad�� 11:tf 12-14:dv�ķ���
****************************************************************************/
//int solve_trans_dv(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double tf);
int solve_trans_dv(const double& longitude, const double& latitude, const double* RV, double& m, double* dv, const double t0, double tf);

/****************************************************************************
* ������   : shoot_func_ToRSAT
* ��  ��   : ��к���,���м���
* �� ��    : x: 0:dvx 1:dvy 2:dvz
* ����     : para: 0-5:rv0 6:mass 7:if_drag(����Ϊ1.0 ������Ϊ0.0) 8:t0 9-11:rvf 12:tf
* ��    �� : fevc: 0:0 1:0.0 2:0 3:0
****************************************************************************/
int shoot_func_ToRSAT(int n, const double* x, double* fevc, int iflag, const double* para);

/****************************************************************************
* ������   : solve_shoot_dv_ToRSAT()
* ��  ��   : ͨ����з������ת�����м���λ�õ�׼ȷ���ٶ�����
* ������� :  RVΪ���ǰλ���ٶȣ�RfΪĿ��λ�ã�mΪû�б��ǰ�����������Ϊ���������������
			  dvΪ��ʼԤ���dv�����Ϊ����dv����tfΪ��ʼԤ���tf�����ٸı�
****************************************************************************/
int solve_shoot_dv_ToRSAT(bool& flag, const double* RV, double* Rf, double& m, double* dv, const double t0, double& tf);

int solve_shoot_dv_ToRSAT2(bool& flag, const double* RV, double* Rf, double& m, double* dv, const double t0, double& tf);

int shoot_func3(int n, const double* x, double* fevc, int iflag, const double* para);

//int solve_shoot_dv3(const MyPoint& point, const double* RV, double& m, double* dv, const double t0, double& tf);

#endif