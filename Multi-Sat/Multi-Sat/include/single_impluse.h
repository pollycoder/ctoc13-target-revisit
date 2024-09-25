#ifndef SINGLE_IMPLUSE_H
#define SINGLE_IMPLUSE_H
//#include"MyPoint.h"

/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 焦艺菲 jiaoyf20@mails.tsinghua.edu.cn
* 文件名: single_impluse.h
* 内容简述：对于地面目标，考虑线性J2摄动的单脉冲快速机动算法
*			参考张刚论文，Analytical Approximate Solutions to Ground Track 
*			Adjustment for Responsive Space
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2020-09-24    焦艺菲     创建该文件
* 01b		2020-09-26    焦         调试完成
* 01c		2020-09-28    焦         增加输出参数
* 02a       2020-10-07    焦         机动观测动目标
* 02b		2020-10-08    焦         修改了一处错误
****************************************************************************/


//单脉冲机动(指定天数的最小脉冲解)
//输入：
//		m0:			初始质量，kg
//		t0:			初始（脉冲）时刻，s
//		rv0[6]:		卫星初始位置速度，km，km/s
//		lambda:		目标点的地面经度，rad
//		phi:		目标点的地面纬度，rad
//		Day:		设定的天数，计算范围(Day-1)~Day
//		branch:		升交段0或降交段1
//		sign:		需要的dv的符号：1正，-1负
//输出：
//		flag:		计算成功返回1，否则返回0
//		mf:			机动后卫星质量，kg
//		tf:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, const int sign);


//单脉冲机动(指定天数的最小脉冲解)
//输入：
//		m0:			初始质量，kg
//		t0:			初始（脉冲）时刻，s
//		rv0[6]:		卫星初始位置速度，km，km/s
//		lambda:		目标点的地面经度，rad
//		phi:		目标点的地面纬度，rad
//		Day:		设定的天数，计算范围(Day-1)~Day
//输出：
//		flag:		计算成功返回1，否则返回0
//		mf:			机动后卫星质量，kg
//		tf:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR);


//角速度加成的方法
//移动目标的单脉冲机动(指定天数的最小脉冲解)
//输入：
//		m0:			初始质量，kg
//		t0:			初始（脉冲）时刻，s
//		rv0[6]:		卫星初始位置速度，km，km/s
//		ID:			动目标编号，201~205
//		Day:		设定的天数，为保证精度通常取1天或2天
//		branch:		升交段0或降交段1
//输出：
//		flag:		计算成功返回1，失败返回0，超出单程时间返回-1
//		mf:			机动后卫星质量，kg
//		tf:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
//		rvt[6]		不输出，传参用NULL即可
void single_imp2(const double m0, const double t0, const double* rv0, const int ID, const int Day,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, double* rvt);

//动目标相关的常数：纬度、换算到30~390的经度边界、单程时间
const double m_lat[5] = { 37.7000, -38.7300, 40.6000, -22.6866, -6.3085 };
const double m_long[5][2] = { 140.9629, 237.3835, 286.6089, 177.9531, 286.7668, 350.7447, 374.3732, 317.9802, 105.2459, 39.0247 };
const double oneway[5] = { 8.0 * 86400.0, 8.0 * 86400.0, 5.0 * 86400.0, 6.0 * 86400.0, 7.0 * 86400.0 };

//计算实时动目标位置，rad
int GetMov(double t, int ID, double& lat, double& lon);


//计算视场角gamma，单位rad
int cal_gamma(double& gamma, const double* robs, int tar_ID, double t);

//计算太阳高度角：角度单位rad，时间是从alpha_G0到当前时刻的秒数
int cal_sunangle(double& sunangle, double t, int tar_ID);

//打靶转移时间和切向脉冲大小，目标为观测角
int shootfun1(int n, const double* x, double* fvec, int iflag, const double* para);

#endif // !SINGLE_IMPLUSE_H