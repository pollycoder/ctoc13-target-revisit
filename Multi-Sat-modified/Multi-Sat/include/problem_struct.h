/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: problem_struct.h
* 内容简述：关于解结构的类
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-05-12    张众      创建文件
****************************************************************************/
#ifndef Problem_struct
#define Problem_struct

#include <vector>

#include "Constant.h"


/****************************************************************************
* 结构体名 : OutputResult
* 功  能   : 最终结果的输出信息,一行内容
****************************************************************************/
class OutputResult
{
public:
	int action_;		//标志位，‘0’表示轨道参数初始行，‘1’表示施加脉冲，‘2’表示观测、
	double time_acc_;   //可见时刻
	double rv_acc_[6];  //可见时刻位置速度
	double dv_[3];		//速度增量
	int point_id_;		//若action_=0或1，则0占位，若动作标志位为2，则为观测对应的目标编号1~20（CTOC13对应0~25的gap编号）
	
	OutputResult()
	{
		action_ = 0;
		time_acc_ = 0.0;
		for (int i = 0; i < 3; i++)	dv_[i] = 0.0;
		point_id_ = 0;
	}
	OutputResult(int _id)
	{
		action_ = 0;
		time_acc_ = 0;
		for (int i = 0; i < 3; i++)	dv_[i] = 0.0;
		point_id_ = _id;
	}

	OutputResult(int _id, double time, double* rv)
	{
		action_ = 2;
		time_acc_ = time;
		for (int i = 0; i < 3; i++)	dv_[i] = 0.0;
		for (int i = 0; i < 6; i++) rv_acc_[i] = rv[i];
		point_id_ = _id; 
	}
	OutputResult(double time, double* RV, double* dv)
	{
		action_ = 1;
		time_acc_ = time;
		memcpy(rv_acc_, RV, 6 * sizeof(double));
		memcpy(dv_, dv, 3 * sizeof(double));
		point_id_ = 0;
	}
	void operator =(const OutputResult& OutputResult_result);    //等号重载
};

/****************************************************************************
* 结构体名 : Node_char
* 功  能   : 记录一个节点关于问题的信息
****************************************************************************/
class Node_problem
{
public:
	//TODO: 问题的信息，包含卫星的质量，时间，高度，位置速度，仅只当前节点内容
	std::vector<OutputResult> node_info_;

	Node_problem()
	{
	}
	void operator =(const Node_problem& result);    //等号重载
};



/****************************************************************************
* 结构体名 : Optimization_index
* 功  能   : 优化指标的结构体
****************************************************************************/
struct Optimization_index
{
public:
	//TODO;确定最终指标：时间，燃料消耗，观测质量
	int observed_num_;  //空间碎片清理个数（CTOC13对应已观测的次数）
	double total_impulse_;   //总速度增量
	double time_cost;   //总观测时间
	double height_aver; //总观测时间
	int total_ship_num;
	double time_min;	// 最小树的时间

	Optimization_index()
	{
		observed_num_ = 0;
		total_impulse_ = 0.0;
		height_aver = 0.0;
		time_cost = 0.0;
		total_ship_num = 0;
		time_min = 0.0;
	}
	Optimization_index(Optimization_index & a): observed_num_(a.observed_num_),total_impulse_(a.total_impulse_), time_cost(a.time_cost),
		height_aver(a.height_aver), total_ship_num(a.total_ship_num), time_min(a.time_min){}
};

/****************************************************************************
* 结构体名 : Solution_one
* 功  能   : 单个卫星的结果
****************************************************************************/
class Solution_one
{
public:
	std::vector<OutputResult> node_info_;  //一颗卫星的节点信息
	void operator =(const Solution_one& a); //符号重载
};


/****************************************************************************
* 结构体名 : Solution
* 功  能   : 所有结果
****************************************************************************/
class Solution
{
public:
	Solution_one  solution_[TreeNum];      //所有的卫星的结果
	double total_impulse_;                    //总的速度增量
	int total_observed_num_;               //总观测个数
	double time_cost_;                      //总观测时间
	double height_aver_;                    //总观测时间
	
	Solution():total_impulse_(0.0), total_observed_num_(0), time_cost_(0), height_aver_(0){}
	
	Solution(const Solution& a) : total_impulse_(a.total_impulse_), total_observed_num_(a.total_observed_num_), time_cost_(a.time_cost_),
	height_aver_(a.height_aver_)
	{
		for (int i = 0; i< TreeNum; i++)
			solution_[i] = a.solution_[i];	
	}
	
	void operator =(const Solution& a);   //符号重载
};
#endif

