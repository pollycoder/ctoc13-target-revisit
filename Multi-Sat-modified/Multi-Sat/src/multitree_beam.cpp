#include "multitree_beam.h"

#include "OrbitMath.h"
#include "single_impluse.h"


bool sort_by_out(const OutputResult& a, const OutputResult& b)
{
	bool result = false;
	if (fabs(a.time_acc_ - b.time_acc_) < 1.0e-5)
	{
		if (a.action_ == 1 && b.action_ == 2) return false;
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
* 功  能   : 排序函数：根据起始节点编号进行total_impulse指标从小到大排序
* 输 入    : 1.节点编号 2.节点编号
****************************************************************************/
inline bool MultiTree::SortTNC(const TNC& a, const TNC& b)
{
	if (a.op_index_.total_impulse_ - b.op_index_.total_impulse_ < -5.0e-4) {
		return true;
	}
	else if (a.op_index_.time_cost < b.op_index_.time_cost) {
		return true;
	}
	else { return false; }
		
	
}

/****************************************************************************
* 函数名   : EqualTNC(const TNC& a, const TNC& b)
* 功  能   : 判断两个tnc序列是否相等，速度增量极其接近
****************************************************************************/
inline bool MultiTree::EqualTNC(const TNC& a, const TNC& b)
{
	if ((fabs(a.op_index_.total_impulse_ - b.op_index_.total_impulse_) < 1.0e-10 )
			 && (a.op_index_.observed_num_ == b.op_index_.observed_num_))
	{
		//得到两个序列
		return true;
	}
	
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
		//if ( !EqualTNC(expandinglist[i],expandinglist[i+1])) //相等
		//{
		//	new_expandinglist.push_back(std::move(expandinglist[i]));
		//}
		new_expandinglist.push_back(std::move(expandinglist[i]));
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
		temp_output.rv_acc_;
		temp_output.action_ = 0;
		temp_output.time_acc_ = 0.0;
		double sats_rv0[6]; int flag;
		coe2rv(flag, sats_rv0, sats_coe0[i], mu_km_s);
		memcpy(temp_output.rv_acc_, sats_rv0, 6 * sizeof(double));

		multi_tree_[i].root_->problem_.node_info_.push_back(temp_output);
	}

	// 每棵树构建一个虚拟根节点
	/*for (int i = 0; i < TreeNum; i++)
	{
		OutputResult temp_output;
		temp_output.action_ = 0;
		temp_output.time_acc_ = 0.0;
		for (int i = 0; i < 6; i++) temp_output.rv_acc_[i] = 0;

		multi_tree_[i].root_->problem_.node_info_.push_back(temp_output);
	}*/

	// 虚拟节点下面扩展实际的初始轨道根数
	// Mesh：4 x 2 x 5 x 7 x 7 x 7
	// SMA：7200km-7350km
	// ECC：0-0.001
	// INC：120°-160°
	// RAAN：0-360°
	// AP：0-360°
	// TA：0-360°
	/*double sma_step = 50.0;
	double ecc_step = 0.001;
	double inc_step = 10.0 * D2R;
	double angle_step = 60.0 * D2R;

	double sma_start = 7200.0, sma_end = 7350.0;
	double ecc_start = 0.0, ecc_end = 0.001;
	double inc_start = 120.0 * D2R, inc_end = 160.0 * D2R;
	double angle_start = 0.0, angle_end = D2PI;*/

	//虚拟节点的TNC
	//std::vector<Node*> old_expand_nodes;
	//for (int i = 0; i < TreeNum; i++) {
	//	old_expand_nodes.push_back(multi_tree_[i].root_);
	//}

	//TNC temp(old_expand_nodes);
	//expandinglist.push_back(std::move(temp));

	//std::vector<Node_problem> temp;						// 备选的初始轨道根数

	//for (double sma = sma_start; sma <= sma_end; sma += sma_step) {
	//	for (double inc = inc_start; inc < inc_end; inc += inc_step) {
	//		for (double raan = angle_start; raan < angle_end; raan += angle_step) {	
	//			for (double ta = angle_start; ta < angle_end; ta += angle_step) {
	//				// 每个循环是要增加一个TNC，TNC要从实际根节点开始放入
	//				OutputResult temp_output;
	//				temp_output.rv_acc_;
	//				temp_output.action_ = 0;
	//				temp_output.time_acc_ = 0.0;
	//				double sats_rv0[6]; int flag;
	//				double coe_0[6] = { sma, 0.0, inc, raan, 0.0, ta };
	//				coe2rv(flag, sats_rv0, coe_0, mu_km_s);
	//				memcpy(temp_output.rv_acc_, sats_rv0, 6 * sizeof(double));
	//				Node_problem node_prob;
	//				node_prob.node_info_.push_back(temp_output);
	//				temp.push_back(node_prob);
	//			}		
	//		}
	//	}
	//}

	//每棵树都把这些备选节点都扩展一遍
	//const int visited[21]{ 0 };
	////const Node_problem initial_coe = temp;
	//
	//std::vector<std::vector<Node*>> expandnodes_alltree;			//8棵树的待选节点
	//for (int i = 0; i < TreeNum; i++) {
	//	expandnodes_alltree[i] = multi_tree_[i].ExpandNode(multi_tree_[i].root_, visited, temp);
	//}

	std::vector<Node*> old_expand_nodes;
	for (int i = 0; i < TreeNum; i++) {
		old_expand_nodes.push_back(multi_tree_[i].root_);
	}

	TNC tnctemp(old_expand_nodes);
	expandinglist.push_back(std::move(tnctemp));
	
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

	temp.total_impulse_ = op_index_.total_impulse_;
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
	double total_impulse = 0.0;
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
			else if (solution_temp.solution_[i].node_info_[j].action_ == 1) {
				total_impulse += V_Norm2(solution_temp.solution_[i].node_info_[j].dv_, 3);
			}
		}

	}

	op_index_.time_cost = time_last;
	op_index_.height_aver = height_total / total_number;
	op_index_.total_impulse_ = total_impulse;
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
		

		// 我们在problems_里面挑选要扩展的节点
		// 扩展应该是都要扩展的
		if (visited[problem_[i].node_info_.back().point_id_] > 0 || V_Norm2(problem_[i].node_info_[0].dv_, 3) > dv_max_) continue; //已经访问过 大于dv_max

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
			Node* temp_node = new Node(node, problem_[i]);
			expandnodes.push_back(temp_node);           //放入新节点
		
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

// TODO: 去掉AccessPointObjects的调用，改成逐一打靶
inline  void children_nodes(Node* node, const int* visited, std::vector<Node_problem>& child_node_problems)
{
	for (int j = 0; j < TargetNum; j++)
	{
		for (int bra = 0; bra < 2; bra++) {
			for (int NR = 0; NR < 20; NR++) {
				if (visited[j] > 0)
				{
					continue;
				}
				double t0, h0, rv0[6], coe0[6], lambda0, phi0, tf, dv[3], a, e, hmin, hmax;
				int flag;

				memcpy(rv0, node->problem_.node_info_.back().rv_acc_, 6 * sizeof(double));
				t0 = node->problem_.node_info_.back().time_acc_;
				tf = t0 + 10800.0;						//从3h开始扰动，瞎给的，用的是张刚论文的方法给初值，tf除了ship随便取

				rv2coe(flag, coe0, rv0, mu_km_s);
				a = coe0[0]; e = coe0[1];
				hmin = a * (1 - e) - Re_km;
				hmax = a * (1 + e) - Re_km;
				if (hmin < 201.0 || hmax > 999.0) continue;

				double rvf[6];

				// 打靶到下一个目标，适配visit_gap
				obs_shooting(flag, dv, tf, rvf, t0, rv0, j, NR, bra);

				// 打靶失败则不扩展
				if (!flag) continue;
				

				double rv1[6], rv2[6];
				memcpy(rv1, rv0, 6 * sizeof(double));
				for (int i = 0; i < 3; i++) rv1[i + 3] += dv[i];

				// 超过两天则只积分到172800s
				if (tf > 2.0 * 86400.0) {
					tf = 2.0 * 86400.0;
					propagate_j2(rv1, rvf, t0, tf);
				}


				//TODO: 将结果输出
				Node_problem temp;
				OutputResult temp_out, temp_out2, temp_out3; //机动的信息，机动后的信息

				temp_out.action_ = 1;
				for (int i = 0; i < 3; i++)  rv0[3 + i] += dv[i];
				memcpy(temp_out.rv_acc_, rv0, 6 * sizeof(double));
				memcpy(temp_out.dv_, dv, 3 * sizeof(double));
				temp_out.time_acc_ = t0;
				temp_out.point_id_ = 0;
				temp.node_info_.push_back(temp_out);

				temp_out2.action_ = 2;
				temp_out2.time_acc_ = tf;
				memcpy(temp_out2.rv_acc_, rvf, 6 * sizeof(double));
				for (int i = 0; i < 3; i++) temp_out2.dv_[i] = 0.0;
				if (fabs(tf - 86400.0 * 2.0) < 1.0e-5) {
					temp_out2.action_ = 3;
					temp_out2.point_id_ = 0;
				}
				else {
					temp_out2.point_id_ = j + 1;
				}

			
				//检测t0到tf时段内有没有其他能看到的目标
				/*double tf_60m = tf - 60.0;
				double tf_60a = tf + 60.0;
				std::vector<std::vector<double>> results;
				AccessPointObjects(rv0, tf_60m, tf_60a, 60.0, 21, results);
				for (int i = 0; i < TargetNum; i++) {
					if (i == temp_out2.point_id_ - 1) continue;
					for (int k = 0; k < results[i].size(); k++) {
						temp_out3.action_ = 2;
						temp_out3.time_acc_ = results[i][k];
						memcpy(temp_out3.rv_acc_, rvf, 6 * sizeof(double));
						for (int m = 0; m < 3; m++) temp_out2.dv_[m] = 0.0;
						temp_out3.point_id_ = i + 1;
						temp.node_info_.push_back(temp_out3);
					}
				}*/

				//最后一个放tf的子节点
				temp.node_info_.push_back(temp_out2);
				/*if (j == TargetNum - 1) {
					std::cout << "成功扩展海上目标，观测时间：" << temp_out2.time_acc_ << "s，脉冲：" << V_Norm2(temp_out.dv_, 3) << std::endl;
				}*/

				child_node_problems.push_back(temp);
			}
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

	int visited[TargetNum]{};                     //该tnc的已观测序列，初始化为0，完成重访任务后为1
	std::vector<std::vector<double>> visible_timelist(21, std::vector<double>(0));	//该tnc的目标点访问时刻表
	//std::cout << visible_timelist.size() << std::endl;
	

	//TODO：扩展规则会改，此处暂时不更新visited
	double tf_max = 0.0;
	double tf_min = 1.0e10;
	for (int j = 0; j < TreeNum; j++)             //将访问过的所有节点按目标位置+1
	{
		std::vector<Node*> solution_one_node;
		tnc.tnc_[j]->return_node_sequence(solution_one_node);
		Solution_one temp;
		tnc.tnc_[j]->getback_problem(solution_one_node, temp);
		if(temp.node_info_[temp.node_info_.size() - 1].time_acc_ < tf_min) tf_min = temp.node_info_[temp.node_info_.size() - 1].time_acc_;
		if(temp.node_info_[temp.node_info_.size() - 1].time_acc_ > tf_max) tf_max = temp.node_info_[temp.node_info_.size() - 1].time_acc_;

		for (int k = 0; k < temp.node_info_.size(); k++)
		{
			if (temp.node_info_[k].action_ == 2)// && fmod(temp.node_info_[k].time_acc_, 60.0) < 1e-4)
				//visited[temp.node_info_[k].point_id_] ++;
				visible_timelist[temp.node_info_[k].point_id_ - 1].push_back(temp.node_info_[k].time_acc_);
		}

		if (temp.node_info_[temp.node_info_.size() - 1].time_acc_ < tf_min) {
			tf_min = temp.node_info_[temp.node_info_.size() - 1].time_acc_;
		}

		if (temp.node_info_[temp.node_info_.size() - 1].time_acc_ > tf_max) {
			tf_max = temp.node_info_[temp.node_info_.size() - 1].time_acc_;
		}
	}

	// Push back之后对可见性时刻表排序
	for (int j = 0; j < TargetNum; j++) {
		std::sort(visible_timelist[j].begin(), visible_timelist[j].end());
	}

	bool ifsync = false;
	if (tf_max - tf_min < 14400.0 && tf_min > 0.0) {
		ifsync = true;
	}

	// 如果有目标点的最大重访时间已经不可能满足要求，立即停止扩展
	// 这一段应该先给8颗星的可能性都做一下验证，确认8颗星已扩展的时间相差不到1.5h为止
	std::vector<double> max_revisit_gap;
	max_revisit_interval_beforeEnd(max_revisit_gap, visible_timelist);
	bool ifpossible = true;
	int idx = 0;
	for (auto iter = max_revisit_gap.begin(); iter != max_revisit_gap.end(); iter++) {
		idx++;
		if (idx != TargetNum && ifsync) {
			if (*iter > 21600.0 && ifsync) {
				//std::cout << "目标" << idx << "的最大重访时间已经达到" << *iter << "sec" << std::endl;
				ifpossible = false;
			}
			else if (*iter == 0.0) {
				//std::cout << "目标" << idx << "还未被看到" << std::endl;
			}
		}
		else {
			if (*iter > 10800.0 && ifsync) {
				//std::cout << "目标" << idx << "的最大重访时间已经达到" << *iter << "sec" << std::endl;
				if (ifsync) ifpossible = false;
			}
			else if (*iter == 0.0) {
				//std::cout << "目标" << idx << "还未被看到" << std::endl;
			}
		}
	}
	//std::cout << std::endl;
	if (!ifpossible) {
		return;
	}

	// 扩展完之后更新该时刻表的最大重访时间，如果重访时间满足要求则visited置1
	std::vector<double> max_revisit_list;
	max_reseetime(visible_timelist, max_revisit_list);
	for (int i = 0; i < TargetNum; i++) {
		if (i != 20) {
			if (max_revisit_list[i] < 6.0) {
				visited[i]++;
				//std::cout << "目标" << i + 1 << "观测已完成" << std::endl;
			}
		}
		else {
			if (max_revisit_list[i] < 3.0) {
				visited[i]++;
				//std::cout << "目标" << i + 1 << "观测已完成" << std::endl;
			}
		}
	}

	bool ifdone = true;
	for (int i = 0; i < TargetNum; i++) {
		if (visited[i] == 0) ifdone = false;
	}
	ifFinished_ = ifdone;
	if(ifdone) std::cout << "观测任务已经全部完成" << std::endl;

	double time[TreeNum]; double time_min = 1.e20; int counter = 0; //扩展时间最小的那颗树
	for (int i = 0; i < TreeNum; i++)
	{
		time[i] = tnc.tnc_[i]->problem_.node_info_.back().time_acc_;

		if (time[i] < time_min - 1.0e-3) {
			time_min = time[i]; counter = i;
		}
	}

	if (fabs(time_min - 86400.0 * 2.0) < 1.0e-3) {
		ifFinished_ = true;
		return;
	}

	if (tnc.op_index_.total_impulse_ > 1.0) {
		return;// 脉冲超过，停止扩展该TNC
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
	bool ifPossible = false;
	
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
			if (!newTNCs.empty()) {
				ifPossible = true;
			}
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
	if (!ifPossible) {
		std::cout << "观测已经不可能完成，退出。" << std::endl;
		ifFinished_ = true;
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
	
	while (expandinglist.size() > 0 && !ifFinished_) //非空
	{
		Expansion(expandinglist);

		beam_discard(expandinglist);                    //排序并取前W_个

		//delete_redundant_nodes();                       //删除多余的子节点

		RecordBestResult(expandinglist,fout0);			//记录最好信息

		if (ifFinished_)
		{
			//std::cout << "观测任务已经全部完成" << std::endl;
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
		result_now_.total_impulse_ = expandinglist[0].op_index_.total_impulse_;
		result_now_.time_cost_ = expandinglist[0].op_index_.time_cost;
		result_now_.total_observed_num_ = expandinglist[0].op_index_.observed_num_;
		result_now_.height_aver_ = expandinglist[0].op_index_.height_aver;
	}
	std::cout << "Removal total num, Time:  " << result_now_.total_observed_num_ << " " << result_now_.time_cost_ << std::endl;
	//fout1 << "Removal total num, Time:  " << result_now_.total_observed_num_ << " " << result_now_.time_cost_ << std::endl;


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
			}
		}
	}
	
	//按照格式输出最终结果
	
	fout1 << std::endl << std::endl << "清理个数：" << result_now_.total_observed_num_ << "当且节点数量" << expandinglist.size() << "当前时间" << result_now_.time_cost_ << std::endl;
	for (int id_sat = 0; id_sat < TreeNum; id_sat++)
	{
		for (int i = 0; i < result_now_.solution_[id_sat].node_info_.size(); i++)
		{
			//fout1 << std::fixed << i + 1 << " " << result_now_.solution_[id_sat].node_info_[i].action_ << " " << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].time_acc_ << " ";
			if (result_now_.solution_[id_sat].node_info_[i].action_ == 0) {
				double rv0[6], coe[6]; int flag;
				memcpy(rv0, result_now_.solution_[id_sat].node_info_[i].rv_acc_, 6 * sizeof(double));
				rv2coe(flag, coe, rv0, mu_km_s);
				fout1 << std::fixed << std::setprecision(16) << coe[0] << " " << coe[1] << " " << coe[2] << " " << coe[3] << " " << coe[4] << " " << coe[5] << std::endl;
			}
			else if (result_now_.solution_[id_sat].node_info_[i].action_ == 1) {
				fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].time_acc_ << " ";
				for (int j = 0; j < 3; j++) {
					fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].dv_[j] << " ";
				}
				fout1 << std::endl;
			}
			else {
				fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].time_acc_ << " ";
				fout1 << std::fixed << std::setprecision(16) << result_now_.solution_[id_sat].node_info_[i].point_id_ << std::endl;
			}
		}
	}



	std::cout << std::endl;
	fout1 << std::endl;

}


