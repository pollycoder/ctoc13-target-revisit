#include "multitree_beam.h"

#include "OrbitMath.h"
#include "single_impluse.h"


bool sort_by_out(const OutputResult& a, const OutputResult& b)
{
	bool result = false;
	if (fabs(a.time_acc_ - b.time_acc_) < 10e-5)
	{
		return a.action_ < b.action_;
	}

	if (a.time_acc_ < b.time_acc_)
	{
		result = true;
	}
	return result;
}


/****************************************************************************
* 函数名   : SortNode(int a,int b)
* 功  能   : 排序函数：根据起始节点编号进行h指标从大到小排序
* 输 入    : 1.节点编号 2.节点编号
****************************************************************************/
inline bool MultiTree::SortTNC(const TNC& a, const TNC& b)
{
	double h_a, h_b;
	h_a = 1.0 / pow(a.op_index_.time_cost, 1);
	h_b = 1.0 / pow(b.op_index_.time_cost, 1);
	return  h_a > h_b;
}

/****************************************************************************
* 函数名   : EqualTNC(const TNC& a, const TNC& b)
* 功  能   : 判断两个tnc序列是否相等，速度增量极其接近
****************************************************************************/
inline bool MultiTree::EqualTNC(const TNC& a, const TNC& b)
{
	//if (fabs(a.op_index_.total_mass_ - b.op_index_.total_mass_) < 1.0e-10)
	//{
	//	//得到两个序列
	//	return true;
	//}
	
	return  false;
}

/****************************************************************************
* 函数名   : void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
* 功  能   : 去重函数(之前已排序过，得到前W个是去重的)
****************************************************************************/
void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
{
	std::vector<TNC> new_expandinglist;
	new_expandinglist.reserve(W_);

	for (int i = 0; i<expandinglist.size()-1; i++) //每次删除后大小都会变
	{
		if ( !EqualTNC(expandinglist[i],expandinglist[i+1])) //相等
		{
			new_expandinglist.push_back(std::move(expandinglist[i]));
		}
		if (new_expandinglist.size() >= W_) break;
	}

	expandinglist.clear();
	expandinglist.reserve(new_expandinglist.size());
	std::move(new_expandinglist.begin(), new_expandinglist.end(), back_inserter(expandinglist));
}
/****************************************************************************
* 函数名   : beam_discard(std::vector<TNC>& expandinglist)
* 功  能   : 排序函数，并取前W_个
****************************************************************************/
void MultiTree::beam_discard(std::vector<TNC>& expandinglist)
{
	//按照h从大到小排序
	sort(expandinglist.begin(), expandinglist.end(), MultiTree::SortTNC);

	unique_remove(expandinglist);

}

/****************************************************************************
* 函数名   : Initialize(std::vector<TNC>& expandinglist)
* 功  能   : 初始化待扩展节点
****************************************************************************/
void MultiTree::Initialize(std::vector<TNC>& expandinglist)
{
	layer_ = 0;
	
	//构建每颗树的初始轨道参数
	for (int i= 0; i< TreeNum; i++)
	{
		OutputResult temp_output;
		temp_output.branch_ = 0;
		temp_output.m_ = 100000.0;
		temp_output.rv_acc_;
		temp_output.action_ = 0;
		double sats_rv0[6]; int flag;
		coe2rv(flag, sats_rv0, sats_coe0[i], mu_km_s);
		memcpy(temp_output.rv_acc_, sats_rv0, 6 * sizeof(double));

		multi_tree_[i].root_->problem_.node_info_.push_back(temp_output);
	}

	//最后将tnc完成
	std::vector<Node*> old_expand_nodes;
	for (int i=0; i< TreeNum; i++)
	{
		old_expand_nodes.push_back(multi_tree_[i].root_);
	}

	TNC temp(old_expand_nodes);
	expandinglist.push_back(std::move(temp));
	

}


/****************************************************************************
* 函数名   : Tree::getback(Node* node)
* 功  能   : 回溯函数，返回一个节点直至根节点的信息
****************************************************************************/
inline void Node::getback_problem(std::vector<Node*>& solution_one_node, Solution_one & temp_solution_one)
{

	//最后一个是最开始的节点
	for (int i = solution_one_node.size() - 1; i > -1; i--)
	{
		for(int j =0; j < solution_one_node[i]->problem_.node_info_.size();j++)
		temp_solution_one.node_info_.push_back(solution_one_node[i]->problem_.node_info_[j]);
	}
}


inline void Node::return_node_sequence(std::vector<Node*>& solution_one_node)
{
	Node* tempnode = this;
	for (; tempnode->parent_; tempnode = tempnode->parent_)
	{
		solution_one_node.push_back(tempnode);
	}
	solution_one_node.push_back(tempnode);
}


inline void Node::update_node(std::vector<Node*> &solution_one_node, Solution_one &temp)
{
	for (int i = solution_one_node.size() - 1; i > -1; i--) //反向循环，从根节点开始
	{
		int j = solution_one_node.size() - 1 - i; //表示对称的位置，即根节点时表示为0
		Node_problem & now_node_problem = solution_one_node[i]->problem_; //当前节点的问题结构体引用
	}
}

/****************************************************************************
* 函数名   : Remove()
* 功  能   : 删除一个节点及其父节点，前提是该节点没有子节点
****************************************************************************/
inline void MultiTree::Remove(Node* a)
{
	Node* tempnode = a;

	while (true)
	{
		Node* tempfaternode = tempnode->parent_;

		if (tempnode->child_.size() == 0 && tempnode->inTNCcounter_ == 0 && tempnode->parent_)
		{
			delete tempnode;
		}
		else
		{
			break;
		}
		tempnode = tempfaternode;
	}
}

void MultiTree::Traverse(Node* node, std::vector<Node*>& last_layer_nodes)
{
	if (node->child_.size() == 0 && node->inTNCcounter_ == 0)
	{
		last_layer_nodes.push_back(node);
		return;
	}

	for (int i = 0; i < node->child_.size(); i++)
	{
		Traverse(node->child_[i], last_layer_nodes);
	}
}

void MultiTree::delete_redundant_nodes()
{
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < TreeNum; i++)
	{
		std::vector<Node*> last_layer_nodes;
		Traverse(multi_tree_[i].root_, last_layer_nodes);
		for (int j = 0; j < last_layer_nodes.size(); j++)
		{
			Remove(last_layer_nodes[j]);
		}
	}
}

/****************************************************************************
* 函数名   : GetBack()
* 功  能   : 返回TNC中所有节点的全部信息
****************************************************************************/
Solution TNC::GetBack()
{
	Solution temp;
	for (int i =0; i< TreeNum; i++)
	{
		std::vector<Node*> solution_one_node;
	    tnc_[i]->return_node_sequence(solution_one_node);             
		tnc_[i]->getback_problem(solution_one_node, temp.solution_[i]);//得到碎片序列
	}

	temp.total_mass_ = op_index_.total_mass_;
	temp.total_observed_num_ = op_index_.observed_num_;

	return temp;
}

/****************************************************************************
* 函数名   : Calculate_op_index()
* 功  能   : 计算指标结构体
****************************************************************************/
void TNC::Calculate_op_index()
{
	//最优指标
	double total_mass = 0.0;
	int  total_number = 0;
	Solution solution_temp = GetBack();
	double height_total = 0.0;
	
	double time_last = -1.0e10;
	for (int i= 0; i< TreeNum; i++)
	{
		if(solution_temp.solution_[i].node_info_.size() == 0) continue;
		
		double time_temp = solution_temp.solution_[i].node_info_.back().time_acc_; //最终时间
		if (time_temp > time_last)  time_last = time_temp;

		for(int j = 0; j < solution_temp.solution_[i].node_info_.size(); j++)      //观测个数
		{
			if (solution_temp.solution_[i].node_info_[j].action_ == 2)
			{
				total_number++;

				height_total += V_Norm2(solution_temp.solution_[i].node_info_[j].rv_acc_, 3) - Re_km;
			}
		}

		total_mass += solution_temp.solution_[i].node_info_.back().m_;             //最终质量
	}

	op_index_.time_cost = time_last;
	op_index_.height_aver = height_total / total_number;
	op_index_.total_mass_ = total_mass;
	op_index_.observed_num_ = total_number;
}

/****************************************************************************
* 函数名   : vector<Node*> Tree::ExpandNode(Node* node, int* visited_debris)
* 功  能   : 扩展一个节点
* 输 入    : 该节点node，以及已经观测的序列visited，观测过的为1，否则为0
* 输    出 : vector<Node*> 一个扩展的列表
****************************************************************************/
std::vector<Node*> Tree::ExpandNode(Node* node, const int* visited, const std::vector<Node_problem> & problem_) //基于问题扩展子节点,本棵树中已经重复的节点就不需要重复了
{

	std::vector<Node*> expandnodes;                 //扩展子节点

	
	for (int i = 0; i < problem_.size(); i++)
	{
		if (visited[problem_[i].node_info_.back().point_id_] > 0 || V_Norm2(problem_[i].node_info_[0].dv_,3) > dv_max_) continue; //已经访问过 大于dv_max

		int id = problem_[i].node_info_.back().point_id_;
		
		auto ifchild = std::find_if(node->child_.begin(), node->child_.end(),
			[id](Node* a) {return a->key_ == id; }); //判断该点是否在其子节点内

		if (ifchild == node->child_.end()) //不存在子节点
		{
			Node* temp_node = new Node(node, problem_[i]);
			expandnodes.push_back(temp_node);           //放入新节点
		}
		else //存在即放入已有子节点
		{
			expandnodes.push_back(*ifchild);
		}
	}

	return expandnodes;
}


/****************************************************************************
* 函数名   : Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
* 功  能   : 扩展一个TNC
* 输 入    : const TNC& tnc 原TNC
* 输    出 : std::vector<TNC>& newTNCs  新扩展的tncs
****************************************************************************/
inline  void children_nodes(Node* node, const int* visited, std::vector<Node_problem>& child_node_problems)
{
	for (int j = 0; j < GapNum; j++)
	{
		// 升交和降交都要算，branch没法在创建节点时指定，只能在这里修改，保持child_node_problems内的节点正确即可
		for (int branch = 0; branch < 2; branch++) 
		{									
			if (visited[j] > 0)
			{
				continue;
			}
			double mass0, t0, rv0[6], lambda0, phi0, tf, dv[3];

			mass0 = node->problem_.node_info_.back().m_;
			memcpy(rv0, node->problem_.node_info_.back().rv_acc_, 6 * sizeof(double));
			t0 = node->problem_.node_info_.back().time_acc_;

			// 经纬度重新指定：适配visit_gap
			double next_target_geodetic[2];
			int id = static_cast<int>(visit_gap[j][0]);
			get_target_geogetic(id, visit_gap[j][1], next_target_geodetic);
			lambda0 = next_target_geodetic[1];
			phi0 = next_target_geodetic[0];

			double norm_r = V_Norm2(rv0, 3);      //轨道高度约束
			if (norm_r - Re_km < 200.0) continue;

			double mass_end = 1000.0;
			double rvf[6];

			// 打靶到下一个目标，适配visit_gap
			shooting_target2target(t0, rv0, visit_gap[j][1], lambda0, phi0, tf, dv, rvf, branch);

			//TODO: 将结果输出
			Node_problem temp;		// temp_D代表降交情况，两种都要push back
			OutputResult temp_out, temp_out2; //机动的信息，机动后的信息

			temp_out.action_ = 1;
			for (int i = 0; i < 3; i++)  rv0[3 + i] += dv[i];
			memcpy(temp_out.rv_acc_, rv0, 6 * sizeof(double));
			memcpy(temp_out.dv_, dv, 3 * sizeof(double));
			temp_out.m_ = mass_end;
			temp_out.time_acc_ = t0;
			temp_out.point_id_ = 0;
			temp_out.branch_ = branch;


			temp_out2.action_ = 2;
			temp_out2.time_acc_ = tf;
			memcpy(temp_out2.rv_acc_, rvf, 6 * sizeof(double));
			for (int i = 0; i < 3; i++) temp_out2.dv_[i] = 0.0;
			temp_out2.m_ = mass_end;
			temp_out2.point_id_ = j;
			temp_out2.branch_ = branch;

			temp.node_info_.push_back(temp_out);
			temp.node_info_.push_back(temp_out2);

			child_node_problems.push_back(temp);
		}
	}
}

/****************************************************************************
* 函数名   : Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
* 功  能   : 扩展一个TNC
* 输 入    : const TNC& tnc 原TNC
* 输    出 : std::vector<TNC>& newTNCs  新扩展的tncs 
****************************************************************************/
void MultiTree::Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
{
	//考虑的是所有节点都可以扩展，但扩展后两个tnc是有可能相同的
	newTNCs.clear();
	std::vector<TNC> delete_tnc;

	int visited[GapNum]{};                     //该tnc的已观测序列，初始化为0，观测后为1


	for (int j = 0; j < TreeNum; j++)             //将访问过的所有节点按目标位置+1
	{
		std::vector<Node*> solution_one_node;
		tnc.tnc_[j]->return_node_sequence(solution_one_node);
		Solution_one temp;
		tnc.tnc_[j]->getback_problem(solution_one_node, temp);
		for (int k = 0; k < temp.node_info_.size(); k++)
		{
			if(temp.node_info_[k].action_ == 2)
				visited[temp.node_info_[k].point_id_] ++;
		}
	}

	double time[TreeNum]; double time_min = 1.e20; int counter = 0; //扩展时间最小的那颗树
	for (int i = 0; i < TreeNum; i++)
	{
		time[i] = tnc.tnc_[i]->problem_.node_info_.back().time_acc_;

		if (time[i] < time_min - 1.0e-3) {
			time_min = time[i]; counter = i;
		}
	}

	for (int i = 0; i < TreeNum; i++) //按树扩展
	{
		if (i != counter) continue;

		std::vector<Node_problem> child_node_problems;
		children_nodes(tnc.tnc_[i], visited, child_node_problems);     //排除掉那些很大速度增量的碎片 TODO: 需要改

		omp_set_lock(&tnc.tnc_[i]->lock);
		std::vector<Node*> leafnode = multi_tree_[i].ExpandNode(tnc.tnc_[i], visited, child_node_problems); //该节点扩展的子节点
		omp_unset_lock(&tnc.tnc_[i]->lock);
		
		for (int j = 0; j < leafnode.size(); j++) //该节点扩展的子节点信息
		{
			Node* temp[TreeNum];
			for (int i = 0; i < TreeNum; i++) temp[i] = tnc.tnc_[i];
			temp[i] = leafnode[j];
			TNC temp_tnc(temp);
			newTNCs.push_back(std::move(temp_tnc));
		}
		
	}

	//扩展后，获得新TNC的指标
	for (int i = 0; i < newTNCs.size(); i++)
	{
		newTNCs[i].Calculate_op_index();
	}

	//先筛选
	sort(newTNCs.begin(), newTNCs.end(), MultiTree::SortTNC);
	//赋值到expand列表
	if (newTNCs.size() > b_)
	{
		newTNCs.erase(newTNCs.begin() + b_, newTNCs.end());
	}
}

/****************************************************************************
* 函数名   : Expansion(std::vector<TNC>& expandinglist) 扩展阶段函数
* 功  能   : 根据扩展信息得到新的扩展列表
****************************************************************************/
void MultiTree::Expansion(std::vector<TNC>& expandinglist)
{
	layer_++;                                             //当前层数
	std::vector< std::vector<TNC>> new_childnode_private;               //局部私有变量
	
#pragma omp parallel
	{
		auto nthreads = omp_get_num_threads();
		auto id = omp_get_thread_num();

#pragma omp single
		{
			new_childnode_private.resize(nthreads);
		}

		//所有需要扩展的节点
#pragma omp for schedule(dynamic) 
		for (int i = 0; i < expandinglist.size(); i++)
		{
			std::vector<TNC> newTNCs;
			Expansion_one_TNC(expandinglist[i], newTNCs);
			move(newTNCs.begin(), newTNCs.end(), back_inserter(new_childnode_private[id]));
		}
#pragma omp single
		{
			expandinglist.clear();
			for (auto& buffer : new_childnode_private) {
				move(buffer.begin(), buffer.end(), back_inserter(expandinglist));
				buffer.clear();
			}
		}
	}

}


/****************************************************************************
* 函数名   : Run()
* 功  能   : 多树搜索的运行总函数
****************************************************************************/
void MultiTree::Run()
{
	for (int i = 0; i < TreeNum; i++)  multi_tree_[i].dv_max_ = dv_max_; 

	std::ofstream fout0("../output_result/result.txt");
	
	std::vector<TNC> expandinglist;
	Initialize( expandinglist);                         //初始化首次扩展表
	
	while (expandinglist.size() > 0 && layer_ < GapNum) //非空
	{
		Expansion(expandinglist);

		beam_discard(expandinglist);                    //排序并取前W_个

		delete_redundant_nodes();                       //删除多余的子节点

		RecordBestResult(expandinglist,fout0);			//记录最好信息

		if (layer_ == GapNum)
		{
			std::vector<TNC> x;
			RecordBestResult(x, fout0);			         //记录最好信息
		}
	}

	
}

void MultiTree::RecordBestResult(std::vector<TNC>& expandinglist, std::ofstream& fout1)
{
	if (expandinglist.size() > 0)
	{
		result_now_ = expandinglist[0].GetBack();
		result_now_.total_mass_ = expandinglist[0].op_index_.total_mass_;
		result_now_.time_cost_ = expandinglist[0].op_index_.time_cost;
		result_now_.total_observed_num_ = expandinglist[0].op_index_.observed_num_;
		result_now_.height_aver_ = expandinglist[0].op_index_.height_aver;
	}
	std::cout << "Removal total num, Time:  " << result_now_.total_observed_num_ << " " << result_now_.time_cost_ << std::endl;
	fout1 << "Removal total num, Time:  " << result_now_.total_observed_num_ << " " << result_now_.time_cost_ << std::endl;


	std::cout << expandinglist.size();


	for (int id_sat = 0; id_sat < TreeNum; id_sat++)
	{
		sort(result_now_.solution_[id_sat].node_info_.begin(), result_now_.solution_[id_sat].node_info_.end(), sort_by_out);
	}
		

	for (int id_sat = 0; id_sat < TreeNum; id_sat++)
	{
		for (int i = 2; i < result_now_.solution_[id_sat].node_info_.size(); i++)
		{
			double dt = result_now_.solution_[id_sat].node_info_[i].time_acc_ - result_now_.solution_[id_sat].node_info_[i - 1].time_acc_;

			if (fabs(dt) < 1e-5)
			{
				memcpy(result_now_.solution_[id_sat].node_info_[i].rv_acc_, result_now_.solution_[id_sat].node_info_[i - 1].rv_acc_, 6 * sizeof(double));
				result_now_.solution_[id_sat].node_info_[i].m_ = result_now_.solution_[id_sat].node_info_[i - 1].m_;
			}
		}
	}
	
	//按照格式输出最终结果
	
	fout1 << std::endl << std::endl << "清理个数：" << result_now_.total_observed_num_ << "当且节点数量" << expandinglist.size() << "当前指标" << result_now_.time_cost_ << "质量" << result_now_.total_mass_ << std::endl;
	for (int id_sat = 0; id_sat < TreeNum; id_sat++)
	{
		for (int i = 0; i < result_now_.solution_[id_sat].node_info_.size(); i++)
		{
			fout1 << std::fixed << i + 1 << " " << result_now_.solution_[id_sat].node_info_.back().branch_ << " "
				<< result_now_.solution_[id_sat].node_info_[i].action_ << " " << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].time_acc_ << " ";
			for (int j = 0; j < 6; j++)
				fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].rv_acc_[j] << " ";
			for (int j = 0; j < 3; j++)
				fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].dv_[j] << " ";
			fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].m_ << " " << result_now_.solution_[id_sat].node_info_[i].point_id_ << std::endl;
		}
	}



	std::cout << std::endl;
	fout1 << std::endl;

}


