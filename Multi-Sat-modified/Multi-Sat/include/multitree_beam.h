/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: multitree_beam.h
* ���ݼ����������ṹ��������
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01       2021-05-12    ����      �����ļ�
* 02	   2024-09-30	 ��ܲ	   �޸�������CTOC13����
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
* ��  ��   : Node
* ��  ��   : �ڵ��࣬���ڴ洢��ǰ������Ϣ
****************************************************************************/
class Node
{
public:
	Node* parent_;                       //���ڵ�
	std::vector<Node*> child_;           //�ӽڵ�
	Node_problem problem_;               //������Ϣ�Ľṹ��
	unsigned int key_;                   //�����룬�������ָ��ڵ��²�ͬ���ӽڵ�
	unsigned int inTNCcounter_;          //ͳ�Ƹýڵ���TNC�еĸ���

	omp_lock_t lock;
	//���캯��
	Node()                               
	{
		parent_ = NULL; key_ = 21; inTNCcounter_ = 0; 
		omp_init_lock(&lock);
	}
	
	//���캯��
	Node(Node* e, const Node_problem& a) 
	{    
		parent_ = e;                     //ָ�򸸽ڵ�
		problem_ = a;

		e->child_.push_back(this);       //��ָ���ӽڵ�

		int t_60_int = static_cast<int>(round(a.node_info_.back().time_acc_ / 60.0));
		this->key_ = a.node_info_.back().point_id_*1e4 + t_60_int;     //�ýڵ�������븳ֵ
		inTNCcounter_ = 0;
		omp_init_lock(&lock);
	}

	//�������������¶���ɾ����������õ����ӽڵ㣬����ɾ��
	~Node()                             
	{
		if (this->child_.size() == 0)
		{
			//����ýڵ��� (���׵ĺ��ӽڵ㣩���ĸ�λ�ã�ɾ����
			std::vector<Node*>& parent_child = (this->parent_)->child_; //��ʱ����
			int value = this->key_;

			//omp_set_lock(&this->parent_->lock);
			
			auto iter = std::find_if(parent_child.begin(), parent_child.end(),
				[value](Node* temp) {return temp->key_ == value; });
			parent_child.erase(iter);
			
			//omp_unset_lock(&this->parent_->lock);
		}
		else
		{
			std::cout << "�����ӽڵ�ȴɾ��" << std::endl;
		}
		omp_destroy_lock(&lock);
	}

	void return_node_sequence(std::vector<Node*>& solution_one_node);            //���ؽڵ����У����������ڵ�
	void getback_problem(std::vector<Node*>& solution_one_node, Solution_one & temp_solution_one);        //���ݺ���������һ���ڵ�ֱ�����ڵ����Ϣ
	
	//void nlopt2solution(Solution_one& temp, double& dv);                         //ͨ��nlopt�Ż�ʱ�䣬����������ṹ��
	void update_node(std::vector<Node*> &solution_one_node, Solution_one &temp); //���ݽ⽫���нڵ��ʱ����ٶ���������


	//void node_pso(Solution_one& temp, double& dv);                               //���ڶ�һ�������pso��ʱ���Ż�
	//double node_nlopt();                                                         //���ڶ�һ�������nlopt�ĺ���ֵ�Ż�
	
};

/****************************************************************************
* ��  ��   : Tree
* ��  ��   : �����ϴ������нڵ�
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

	std::vector<Node*> ExpandNode(Node*, const int* visited, const std::vector<Node_problem>& problem_); //��չ����
};

/****************************************************************************
* ��  ��   : TNC
* ��  ��   : ��������ڵ㣬��ﵱǰ״̬
****************************************************************************/
class TNC
{
public:
	Node *tnc_[TreeNum];               //����ָ�룬һ�����飬�����ŵ��ǽڵ�ָ��
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

		//��Ӧ�ڵ��TNCnumbe+1
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

		//��Ӧ�ڵ��TNCnumbe+1
		for (int i = 0; i < TreeNum; i++)
		{
			//omp_set_lock(&a[i]->lock);
			a[i]->inTNCcounter_++;
			//omp_unset_lock(&a[i]->lock);
		}
	}
	
	//�������캯��
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

	//�ƶ����캯��
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

		if( tnc_[0])  //��Ϊ��
		{
			//��Ӧ�ڵ��TNCnumbe-1
			for (int i = 0; i < TreeNum; i++)
			{
				omp_set_lock(&tnc_[i]->lock);
				tnc_[i]->inTNCcounter_--;
				omp_unset_lock(&tnc_[i]->lock);
			}
		}
	}

	

	Solution GetBack();  //��ȡTNC�����нڵ�Ļ��ݵ���Ϣ

	void Calculate_op_index(); //��������ָ��:������������ٶ�����
};

/****************************************************************************
* ��  ��   : MultiTree
* ��  ��   : ������ܵ�ʹ��
****************************************************************************/
class MultiTree
{
public:
	int layer_;                                                          //����
	int W_;                                                              //�������
	int n_;                                                              //��ʼ�ڵ���
	int b_;                                                              //һ��tnc��չ��ɸѡ
	double dv_max_;
	bool ifFinished_ = false;											// �ط������Ƿ����
	//int range;

	
	Tree multi_tree_[TreeNum];										     //��Ÿ�����

	Solution    result_now_;											 //��ǰ�����Ž�
	Solution    result_all_;											 //ȫ�����Ž�

	MultiTree(int beamwidth, int n, int b, double dvmax) :
	layer_(0), W_(beamwidth), n_(n),b_(b),dv_max_(dvmax){} //���캯��

	void Expansion_one_TNC(const TNC & tnc, std::vector<TNC> & newTNCs); //��չһ��TNC

	void Expansion(std::vector<TNC> & expandinglist);                    //��չ����

	static bool SortTNC(const TNC& a, const TNC& b);                     //��������
	static bool EqualTNC(const TNC& a, const TNC& b);
	void unique_remove(std::vector<TNC>& expandinglist);

	void beam_discard(std::vector<TNC>& expandinglist);                  //����ȡǰW_��

	void Initialize(std::vector<TNC>& expandinglist);                    //��ʼ���״���չ��
	
	void Run();                                                          //����������

	void RecordBestResult(std::vector<TNC>& expandinglist,std::ofstream &fout0);   //��¼��ý�

	void Remove(Node* a);                                                //������ɾ���ýڵ㼰�丸�ڵ㣬ǰ���ǲ����������ӽڵ�
	
	void Traverse(Node* node, std::vector<Node*>& last_layer_nodes);          //�����������нڵ㣬���صײ���Ҫɾ���Ľڵ�

	void delete_redundant_nodes();                                       //ɾ�����ж���ڵ�
};
#endif
