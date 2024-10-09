#include "problem_struct.h"

void OutputResult::operator =(const OutputResult& OutputResult_result)
{
	action_ = OutputResult_result.action_;
	time_acc_ = OutputResult_result.time_acc_;
	memcpy(rv_acc_, OutputResult_result.rv_acc_, 6 * sizeof(double));
	memcpy(dv_, OutputResult_result.dv_, 3 * sizeof(double));
	point_id_ = OutputResult_result.point_id_;
}

void Node_problem::operator=(const Node_problem& result)
{
	node_info_.resize(result.node_info_.size());
	for(int i =0 ; i< result.node_info_.size(); i++)
	{
		node_info_[i] = result.node_info_[i];
	}
}

void Solution_one::operator=(const Solution_one& a)
{
	node_info_.clear();
	for (int i= 0; i< a.node_info_.size();i ++)
	{
		node_info_.push_back(a.node_info_[i]);
	}
}

void Solution::operator=(const Solution& a)
{
	total_impulse_ = a.total_impulse_;
	total_observed_num_ = a.total_observed_num_;
	for (int i = 0; i < TreeNum; i++)
	{
		solution_[i] = a.solution_[i];
	}
}
