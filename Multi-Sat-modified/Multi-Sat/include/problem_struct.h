/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: problem_struct.h
* ���ݼ��������ڽ�ṹ����
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01       2021-05-12    ����      �����ļ�
****************************************************************************/
#ifndef Problem_struct
#define Problem_struct

#include <vector>

#include "Constant.h"


/****************************************************************************
* �ṹ���� : OutputResult
* ��  ��   : ���ս���������Ϣ,һ������
****************************************************************************/
class OutputResult
{
public:
	int action_;		//��־λ����0����ʾ���������ʼ�У���1����ʾʩ�����壬��2����ʾ�۲⡢
	double time_acc_;   //�ɼ�ʱ��
	double rv_acc_[6];  //�ɼ�ʱ��λ���ٶ�
	double dv_[3];		//�ٶ�����
	double m_;			//������������ʩ�����ٶ�����������ʩ���ٶ������������
	int point_id_;		//��action_=0��1����0ռλ����������־λΪ2����Ϊ�۲��Ӧ��Ŀ����1~21
	int branch_;        //0Ϊ������1Ϊ����
	OutputResult()
	{
		action_ = 0;
		time_acc_ = 0.0;
		for (int i = 0; i < 3; i++)	dv_[i] = 0.0;
		point_id_ = 0;
		branch_ = 1;
	}
	OutputResult(int _id)
	{
		action_ = 0;
		time_acc_ = 0;
		for (int i = 0; i < 3; i++)	dv_[i] = 0.0;
		point_id_ = _id;
		branch_ = 1;
	}

	OutputResult(int _id, double time, double* rv, double m, int branch)
	{
		action_ = 2;
		time_acc_ = time;
		for (int i = 0; i < 3; i++)	dv_[i] = 0.0;
		for (int i = 0; i < 6; i++) rv_acc_[i] = rv[i];
		point_id_ = _id;
		m_ = m;
		branch_ = branch;
	}
	OutputResult(double time, double* RV, double* dv, double m, int branch)
	{
		action_ = 1;
		time_acc_ = time;
		memcpy(rv_acc_, RV, 6 * sizeof(double));
		memcpy(dv_, dv, 3 * sizeof(double));
		m_ = m;
		point_id_ = 0;
		branch_ = branch;
	}
	void operator =(const OutputResult& OutputResult_result);    //�Ⱥ�����
};

/****************************************************************************
* �ṹ���� : Node_char
* ��  ��   : ��¼һ���ڵ�����������Ϣ
****************************************************************************/
class Node_problem
{
public:
	//TODO: �������Ϣ���������ǵ�������ʱ�䣬�߶ȣ�λ���ٶȣ���ֻ��ǰ�ڵ�����
	std::vector<OutputResult> node_info_;

	Node_problem()
	{
	}
	void operator =(const Node_problem& result);    //�Ⱥ�����
};



/****************************************************************************
* �ṹ���� : Optimization_index
* ��  ��   : �Ż�ָ��Ľṹ��
****************************************************************************/
struct Optimization_index
{
public:
	//TODO;ȷ������ָ�꣺ʱ�䣬ȼ�����ģ��۲�����
	int observed_num_;  //�ռ���Ƭ���������CTOC13��Ӧ�����gap������
	double total_mass_;   //���ٶ�����
	double time_cost;   //�ܹ۲�ʱ��
	double height_aver; //�ܹ۲�ʱ��

	Optimization_index()
	{
		observed_num_ = 0;
		total_mass_ = 0.0;
		height_aver = 0.0;
		time_cost = 0.0;
	}
	Optimization_index(Optimization_index & a): observed_num_(a.observed_num_),total_mass_(a.total_mass_), time_cost(a.time_cost),
		height_aver(a.height_aver){}
};

/****************************************************************************
* �ṹ���� : Solution_one
* ��  ��   : �������ǵĽ��
****************************************************************************/
class Solution_one
{
public:
	std::vector<OutputResult> node_info_;  //һ�����ǵĽڵ���Ϣ
	void operator =(const Solution_one& a); //��������
};


/****************************************************************************
* �ṹ���� : Solution
* ��  ��   : ���н��
****************************************************************************/
class Solution
{
public:
	Solution_one  solution_[TreeNum];      //���е����ǵĽ��
	double total_mass_;                    //�ܵ��ٶ�����
	int total_observed_num_;               //�ܹ۲����
	double time_cost_;                      //�ܹ۲�ʱ��
	double height_aver_;                    //�ܹ۲�ʱ��
	
	Solution():total_mass_(0.0), total_observed_num_(0), time_cost_(0), height_aver_(0){}
	
	Solution(const Solution& a) : total_mass_(a.total_mass_), total_observed_num_(a.total_observed_num_), time_cost_(a.time_cost_),
	height_aver_(a.height_aver_)
	{
		for (int i = 0; i< TreeNum; i++)
			solution_[i] = a.solution_[i];	
	}
	
	void operator =(const Solution& a);   //��������
};
#endif

