#ifndef SINGLE_IMPLUSE_H
#define SINGLE_IMPLUSE_H

#include<math.h>
#include<vector>
#include"single_impluse.h"
#include"Constant.h"
#include"OrbitFun.h"
#include"J2propagation.h"
#include"OrbitMath.h"
#include"visibility_4_targets.h"
#include"PSO_DE_NLOPT.h"
#include "J2Lambert.h"

/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���շ� jiaoyf20@mails.tsinghua.edu.cn
* �ļ���: single_impluse.h
* ���ݼ��������ڵ���Ŀ��ͺ���Ŀ�꣬��������J2�㶯�ĵ�������ٻ����㷨��
*			�����ݿ����Խ��оֲ��Ż�������ѧģ����ATK�ײ�ģ�����
*			�ο��Ÿ����ģ�Analytical Approximate Solutions to Ground Track 
*			Adjustment for Responsive Space
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2020-09-24    ���շ�     �������ļ�
* 01b		2020-09-26    ��         �������
* 01c		2020-09-28    ��         �����������
* 02a       2020-10-07    ��         �����۲⶯Ŀ��
* 02b		2020-10-08    ��         �޸���һ������
* 03a		2024-09-28	  ��ܲ		 �ع��˳�������CTOC13����
* 04a		2024-10-31	  ��		 ���Ӻ���Ŀ��۲�
****************************************************************************/


//���������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		lambda:		Ŀ���ĵ��澭�ȣ�rad
//		phi:		Ŀ���ĵ���γ�ȣ�rad
//		Day:		�趨�����������㷶Χ(Day-1)~Day
//		branch:		������0�򽵽���1
//		sign:		��Ҫ��dv�ķ��ţ�1����-1��
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, const int sign);


//���������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		lambda:		Ŀ���ĵ��澭�ȣ�rad
//		phi:		Ŀ���ĵ���γ�ȣ�rad
//		Day:		�趨�����������㷶Χ(Day-1)~Day
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch);


//��Ժ���Ŀ��ĵ��������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		lambda:		Ŀ���ĵ��澭�ȣ�rad
//		phi:		Ŀ���ĵ���γ�ȣ�rad
//		Day:		�趨�����������㷶Χ(Day-1)~Day
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp_ship(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch);


//CTOC13�����������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		lambda:		Ŀ���ĵ��澭�ȣ�rad
//		phi:		Ŀ���ĵ���γ�ȣ�rad
//		dt:			ָ���ķ���ʱ��
//�����
//		flag:		����ɹ�����1�����򷵻�0
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const double dt,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, const int sign);


//CTOC13������J2Lambert���⣬���е�����������
//TODO����Ҫ�ĳ��������߶ȵ�������Ա�Lambert�Ż����Ż��߶Ⱥ�ʵ�ʸ��ǵ�Ŀ���
//���Ÿ������������ƣ��������岢���ǽ���������һ������߶ȷ�Χ�ڽ���ѡ������������ܹ۲⵽����Ŀ���
//���룺
//		RV0[6]����ʼλ���ٶ�
//		t0����ʼʱ��
//		tf���ն�ʱ��
//�����
//		dv[3]�����յ�����
//		RVf[6]���ն�λ���ٶ�
//		flag�����ɹ�����1�����ʧ�ܷ���0    
//		h�����ղ��õĸ߶�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
void single_imp_J2Lambert(double* dv, double* RVf, int& flag, const double* RV0, const double& t0, const double& tf, const double& lambda0, const double& phi0, const double& h);


//CTOC13������J2Lambert���Ż�����Ŀ���Ĺ���߶�
//�Ż����������㶯���
//		h������߶�
//		(lambda, phi)��ʵ�ʾ����ĵ��澭γ�ȣ����Ŷ�ʱ�䣬�ٽ�һ���Ŷ���γ��
//		tf��ʵ�ʵ���ĩʱ��
//�Ż�ָ�꣺
//		dv�������С
//������
//		t0����ʼʱ��
//		RV0[6]����ʼλ���ٶ�
//		tf���Ż�ǰ����ĩʱ��
//		h0���Ż�ǰ�Ĺ���߶�
//		target_id��Ҫ�۲��Ŀ����ţ�0-20��
//Լ����
//		��֤Ŀ��ɼ�
//		����߶�200-1000km
//��ʽ��PSO��Ӧ��Obj_func��ȫһ�£�����ֻ��Ҫ��ȡ���յ��Ż�����X
//��ͬ��֮ǰÿ�ζ�Ҫ��дget_value�ٰ�װ������Ϊ��childnode����������˼·������
double obj_func_shooting(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const int& target_id, const int& NR, const int& branch);

#endif // !SINGLE_IMPLUSE_H