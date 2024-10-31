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
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 焦艺菲 jiaoyf20@mails.tsinghua.edu.cn
* 文件名: single_impluse.h
* 内容简述：对于地面目标和海上目标，考虑线性J2摄动的单脉冲快速机动算法，
*			并根据可视性进行局部优化，动力学模型用ATK底层模型替代
*			参考张刚论文，Analytical Approximate Solutions to Ground Track 
*			Adjustment for Responsive Space
* 文件历史：
* 版本号     日期         作者       说明
* 01a       2020-09-24    焦艺菲     创建该文件
* 01b		2020-09-26    焦         调试完成
* 01c		2020-09-28    焦         增加输出参数
* 02a       2020-10-07    焦         机动观测动目标
* 02b		2020-10-08    焦         修改了一处错误
* 03a		2024-09-28	  周懿		 重构了程序，适配CTOC13丙题
* 04a		2024-10-31	  周		 增加海上目标观测
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
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch);


//针对海上目标的单脉冲机动(指定天数的最小脉冲解)
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
void single_imp_ship(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const int Day,
	int& flag0, double& mf, double& tf, double* dv, const int& NR, const int& branch);


//CTOC13：单脉冲机动(指定天数的最小脉冲解)
//输入：
//		m0:			初始质量，kg
//		t0:			初始（脉冲）时刻，s
//		rv0[6]:		卫星初始位置速度，km，km/s
//		lambda:		目标点的地面经度，rad
//		phi:		目标点的地面纬度，rad
//		dt:			指定的飞行时间
//输出：
//		flag:		计算成功返回1，否则返回0
//		mf:			机动后卫星质量，kg
//		tf:			飞行时长，s
//		dv[3]:		脉冲速度增量，km/s
//		NR:			转移圈数
void single_imp(const double m0, const double t0, const double* rv0, const double lambda0, const double phi0, const double dt,
	int& flag0, double& mf, double& tf, double* dv, int& NR, const int branch, const int sign);


//CTOC13：利用J2Lambert问题，进行单次脉冲修正
//TODO：需要改成输入轨道高度的情况，以便Lambert优化器优化高度和实际覆盖的目标点
//与张刚论文作用类似，但是脉冲并不是近似切向，在一个轨道高度范围内近似选择脉冲最低且能观测到地面目标的
//输入：
//		RV0[6]：初始位置速度
//		t0：初始时刻
//		tf：终端时刻
//输出：
//		dv[3]：最终的脉冲
//		RVf[6]：终端位置速度
//		flag：求解成功返回1，求解失败返回0    
//		h：最终采用的高度                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
void single_imp_J2Lambert(double* dv, double* RVf, int& flag, const double* RV0, const double& t0, const double& tf, const double& lambda0, const double& phi0, const double& h);


//CTOC13：利用J2Lambert，优化覆盖目标点的轨道高度
//优化变量（被摄动项）：
//		h：轨道高度
//		(lambda, phi)：实际经过的地面经纬度，先扰动时间，再进一步扰动经纬度
//		tf：实际的终末时刻
//优化指标：
//		dv：脉冲大小
//参数：
//		t0：初始时刻
//		RV0[6]：初始位置速度
//		tf：优化前的终末时刻
//		h0：优化前的轨道高度
//		target_id：要观测的目标序号（0-20）
//约束：
//		保证目标可见
//		轨道高度200-1000km
//格式与PSO对应的Obj_func完全一致，我们只需要获取最终的优化变量X
//不同于之前每次都要先写get_value再包装，这是为了childnode函数看起来思路更清晰
double obj_func_shooting(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const int& target_id, const int& NR, const int& branch);

#endif // !SINGLE_IMPLUSE_H