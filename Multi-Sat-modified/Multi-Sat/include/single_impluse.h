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




//CTOC13: 动力学模型替代
//优化切向脉冲，使之仍然保持可见且尽量最小（局部优化）
//优化变量：
//		dv_perturb[3]：当次的脉冲扰动（三个分量），扰动量为±0.01km/s
//		tf_perturb：飞行时间的扰动，扰动量为±0.5h
//参数：
//		t0：初始时刻
//		tf：终末时刻 - 以single_impulse计算的结果为准（需要再叠加上t0）
//		rv0：初始位置速度
//		lambda0：目标的经度
//		phi0：目标的纬度
//		time_stamp：指定的时间戳 - 以ATK需要补充观测的时间点为准
//		dv0：切向脉冲的初值
//输出：
//		dv[3]：当次的脉冲
//		tf：优化之后实际的终末时刻
//		rvf：优化之后实际的rvf
//约束：
//		可见性：不可见则惩罚
double obj_func(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
void get_score_data(const std::vector<double>& X, const double* para, double* dv, double* rvf, double& tf);


//CTOC13：从一个目标向下一个目标打靶观测
//不要求精确过点，在视场范围内即可，由single_impulse提供初值，nlopt进行局部优化
//输入：
//		t0：初始时刻
//		rv0：初始位置速度
//		time_stamp：指定的打靶时刻，实际可以上下浮动0.5h
//		lambda0：下一个目标点的经度
//		phi0：下一个目标点的纬度
//输出：
//		tf：打靶之后的终末时间
//		dv[3]：打靶之后的脉冲向量
//		rvf：观测下一个目标时卫星的地心惯性系位置速度，供下一次打靶填补观测缺漏使用
void shooting_target2target(const double t0, const double* rv0, const double time_stamp, const double lambda0, const double phi0, double& tf, double* dv, double* rv, const int& branch);


//CTOC13：获取用来搜索的最佳初始轨道根数
//计算无机动情况下，时刻表上所有点平均经纬度误差最小的初始轨道根数
//输入：
//		gap_list[2][26]：时刻表
//	    X：优化变量
//输出：
//		coe0[6]：最优初始轨道根数
//		average_err：平均经纬度误差（指标）
void optpara_2_real_mission(const std::vector<double>& X, const int sat_num, std::vector<std::vector<double>>& coe_chief);
double obj_func_initial_coe(const std::vector<double>& X, std::vector<double>& grad, void* f_data);
void get_initial_coe(const std::vector<double>& X, double* coe0, double& average_err);


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
void single_imp(double* dv, double* RVf, int& flag, const double* RV0, const double& t0, const double& tf, const double& lambda0, const double& phi0, const double& h);


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
void obs_shooting(int& flag, double* dv, double& tf, double* RVf, const double& t0, const double* RV0, const double& h0, const int& target_id);

#endif // !SINGLE_IMPLUSE_H