#ifndef SINGLE_IMPLUSE_H
#define SINGLE_IMPLUSE_H
//#include"MyPoint.h"

/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���շ� jiaoyf20@mails.tsinghua.edu.cn
* �ļ���: single_impluse.h
* ���ݼ��������ڵ���Ŀ�꣬��������J2�㶯�ĵ�������ٻ����㷨
*			�ο��Ÿ����ģ�Analytical Approximate Solutions to Ground Track 
*			Adjustment for Responsive Space
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01a       2020-09-24    ���շ�     �������ļ�
* 01b		2020-09-26    ��         �������
* 01c		2020-09-28    ��         �����������
* 02a       2020-10-07    ��         �����۲⶯Ŀ��
* 02b		2020-10-08    ��         �޸���һ������
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
	int& flag0, double& mf, double& tf, double* dv, int& NR);


//���ٶȼӳɵķ���
//�ƶ�Ŀ��ĵ��������(ָ����������С�����)
//���룺
//		m0:			��ʼ������kg
//		t0:			��ʼ�����壩ʱ�̣�s
//		rv0[6]:		���ǳ�ʼλ���ٶȣ�km��km/s
//		ID:			��Ŀ���ţ�201~205
//		Day:		�趨��������Ϊ��֤����ͨ��ȡ1���2��
//		branch:		������0�򽵽���1
//�����
//		flag:		����ɹ�����1��ʧ�ܷ���0����������ʱ�䷵��-1
//		mf:			����������������kg
//		tf:			����ʱ����s
//		dv[3]:		�����ٶ�������km/s
//		NR:			ת��Ȧ��
//		rvt[6]		�������������NULL����
void single_imp2(const double m0, const double t0, const double* rv0, const int ID, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, double* rvt);

//��Ŀ����صĳ�����γ�ȡ����㵽30~390�ľ��ȱ߽硢����ʱ��
const double m_lat[5] = { 37.7000, -38.7300, 40.6000, -22.6866, -6.3085 };
const double m_long[5][2] = { 140.9629, 237.3835, 286.6089, 177.9531, 286.7668, 350.7447, 374.3732, 317.9802, 105.2459, 39.0247 };
const double oneway[5] = { 8.0 * 86400.0, 8.0 * 86400.0, 5.0 * 86400.0, 6.0 * 86400.0, 7.0 * 86400.0 };

//����ʵʱ��Ŀ��λ�ã�rad
int GetMov(double t, int ID, double& lat, double& lon);


//�����ӳ���gamma����λrad
int cal_gamma(double& gamma, const double* robs, int tar_ID, double t);

//����̫���߶Ƚǣ��Ƕȵ�λrad��ʱ���Ǵ�alpha_G0����ǰʱ�̵�����
int cal_sunangle(double& sunangle, double t, int tar_ID);

//���ת��ʱ������������С��Ŀ��Ϊ�۲��
int shootfun1(int n, const double* x, double* fvec, int iflag, const double* para);

#endif // !SINGLE_IMPLUSE_H