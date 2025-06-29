/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: multitree_beam.h
* 内容简述：多树结构的类声明
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-05-12    张众      创建文件
* 02	   2024-09-30	 周懿	   修改以适配CTOC13丙题
****************************************************************************/
#ifndef _BEAM_H
#define _BEAM_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <omp.h>
#include <iomanip>

#include "problem_struct.h"
#include "PSO_DE_NLOPT.h"
#include "max_reseetime.h"

/****************************************************************************
* 类  名   : Node
* 功  能   : 节点类，用于存储当前卫星信息
****************************************************************************/
class Node
{
public:
	Node* parent_;                       //父节点
	std::vector<Node*> child_;           //子节点
	Node_problem problem_;               //问题信息的结构体
	unsigned int key_;                   //特征码，用于区分父节点下不同的子节点
	unsigned int inTNCcounter_;          //统计该节点在TNC中的个数

	omp_lock_t lock;
	//构造函数
	Node()                               
	{
		parent_ = NULL; key_ = 21; inTNCcounter_ = 0; 
		omp_init_lock(&lock);
	}
	
	//构造函数
	Node(Node* e, const Node_problem& a) 
	{    
		parent_ = e;                     //指向父节点
		problem_ = a;

		e->child_.push_back(this);       //父指向子节点

		int t_60_int = static_cast<int>(round(a.node_info_.back().time_acc_ / 60.0));
		this->key_ = a.node_info_.back().point_id_*1e4 + t_60_int;     //该节点的特征码赋值
		inTNCcounter_ = 0;
		omp_init_lock(&lock);
	}

	//析构函数，从下而上删除，即如果该点有子节点，不能删除
	~Node()                             
	{
		if (this->child_.size() == 0)
		{
			//查出该节点在 (父亲的孩子节点）的哪个位置，删掉它
			std::vector<Node*>& parent_child = (this->parent_)->child_; //临时变量
			int value = this->key_;

			//omp_set_lock(&this->parent_->lock);
			
			auto iter = std::find_if(parent_child.begin(), parent_child.end(),
				[value](Node* temp) {return temp->key_ == value; });
			parent_child.erase(iter);
			
			//omp_unset_lock(&this->parent_->lock);
		}
		else
		{
			std::cout << "存在子节点却删除" << std::endl;
		}
		omp_destroy_lock(&lock);
	}

	void return_node_sequence(std::vector<Node*>& solution_one_node);            //返回节点序列，不包含根节点
	void getback_problem(std::vector<Node*>& solution_one_node, Solution_one & temp_solution_one);        //回溯函数，返回一个节点直至根节点的信息
	
	//void nlopt2solution(Solution_one& temp, double& dv);                         //通过nlopt优化时间，并返回至解结构体
	void update_node(std::vector<Node*> &solution_one_node, Solution_one &temp); //根据解将所有节点的时间和速度增量更新


	//void node_pso(Solution_one& temp, double& dv);                               //用于对一条解进行pso的时间优化
	//double node_nlopt();                                                         //用于对一个解进行nlopt的含初值优化
	
};

/****************************************************************************
* 类  名   : Tree
* 功  能   : 在树上存在所有节点
****************************************************************************/
class Tree
{
public:
	Node* root_;
	double dv_max_;
	
	Tree()
	{
		dv_max_ = 0.0;
		root_ = new Node();
	}
	~Tree()
	{
		if (root_) delete root_;
	}

	std::vector<Node*> ExpandNode(Node*, const int* visited, const std::vector<Node_problem>& problem_); //扩展函数
};

/****************************************************************************
* 类  名   : TNC
* 功  能   : 包含多个节点，表达当前状态
****************************************************************************/
class TNC
{
public:
	Node *tnc_[TreeNum];               //数组指针，一个数组，里面存放的是节点指针
	Optimization_index op_index_;

	TNC()
	{
		for (int i = 0; i< TreeNum; i++)
		tnc_[i] = nullptr;
	}
	
	TNC(Node* a[TreeNum])
	{
		for (int i = 0; i < TreeNum; i++)
			tnc_[i] = a[i];

		//对应节点的TNCnumbe+1
		for (int i = 0; i< TreeNum; i++)
		{
			omp_set_lock(&a[i]->lock);
			a[i]->inTNCcounter_++;
			omp_unset_lock(&a[i]->lock);
		}
	}
	
	TNC(std::vector<Node*> a)
	{
		for (int i = 0; i < TreeNum; i++)
			tnc_[i] = a[i];

		//对应节点的TNCnumbe+1
		for (int i = 0; i < TreeNum; i++)
		{
			//omp_set_lock(&a[i]->lock);
			a[i]->inTNCcounter_++;
			//omp_unset_lock(&a[i]->lock);
		}
	}
	
	//拷贝构造函数
	TNC(const TNC & old)
	{
		op_index_ = old.op_index_;
		for (int i = 0; i< TreeNum; i++)
		{
			tnc_[i] = old.tnc_[i];
			//old.tnc_[i] = nullptr;
		}
		std::cout << "copy_construct " << std::endl;
	}

	//移动构造函数
	TNC(TNC && old) noexcept 
	: op_index_(old.op_index_)
	{
		for (int i = 0; i < TreeNum; i++)
		{
			tnc_[i] = old.tnc_[i];
			old.tnc_[i] = nullptr;
		}
	}

	TNC& operator=(const TNC& old) {
		
		op_index_ = old.op_index_;
		for (int i = 0; i < TreeNum; i++)
		{
			tnc_[i] = old.tnc_[i];
			//old.tnc_[i] = nullptr;
		}
		return *this;
	}

	TNC& operator=(TNC&& a)  noexcept
	{
		op_index_ = a.op_index_;
		for (int i = 0; i < TreeNum; i++)
		{
			tnc_[i] = a.tnc_[i];
			a.tnc_[i] = nullptr;
		}
		return *this;
	}
	
	~TNC()
	{

		if( tnc_[0])  //不为空
		{
			//对应节点的TNCnumbe-1
			for (int i = 0; i < TreeNum; i++)
			{
				omp_set_lock(&tnc_[i]->lock);
				tnc_[i]->inTNCcounter_--;
				omp_unset_lock(&tnc_[i]->lock);
			}
		}
	}

	

	Solution GetBack();  //获取TNC中所有节点的回溯点信息

	void Calculate_op_index(); //计算最优指标:清理个数和总速度增量
};

/****************************************************************************
* 类  名   : MultiTree
* 功  能   : 多树框架的使用
****************************************************************************/
class MultiTree
{
public:
	int layer_;                                                          //层数
	int W_;                                                              //集束宽度
	int n_;                                                              //初始节点宽度
	int b_;                                                              //一个tnc扩展后筛选
	double dv_max_;
	bool ifFinished_ = false;											// 重访任务是否完成
	//int range;

	
	Tree multi_tree_[TreeNum];										     //存放各个树

	Solution    result_now_;											 //当前代最优解
	Solution    result_all_;											 //全局最优解

	MultiTree(int beamwidth, int n, int b, double dvmax) :
	layer_(0), W_(beamwidth), n_(n),b_(b),dv_max_(dvmax){} //构造函数

	void Expansion_one_TNC(const TNC & tnc, std::vector<TNC> & newTNCs); //扩展一个TNC

	void Expansion(std::vector<TNC> & expandinglist);                    //扩展函数

	static bool SortTNC(const TNC& a, const TNC& b);                     //排序依据
	static bool EqualTNC(const TNC& a, const TNC& b);
	void unique_remove(std::vector<TNC>& expandinglist);

	void beam_discard(std::vector<TNC>& expandinglist);                  //排序并取前W_个

	void Initialize(std::vector<TNC>& expandinglist);                    //初始化首次扩展表
	
	void Run();                                                          //运行主函数

	void RecordBestResult(std::vector<TNC>& expandinglist,std::ofstream &fout0);   //记录最好解

	void Remove(Node* a);                                                //从树中删除该节点及其父节点，前提是不存在其余子节点
	
	void Traverse(Node* node, std::vector<Node*>& last_layer_nodes);          //遍历树中所有节点，返回底部需要删除的节点

	void delete_redundant_nodes();                                       //删除所有多余节点
};
#endif
