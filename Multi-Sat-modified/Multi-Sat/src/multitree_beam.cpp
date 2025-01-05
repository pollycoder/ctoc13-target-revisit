#include "multitree_beam.h"

#include <unordered_map>

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

//bool sort_Nodes(const Node_problem& a, const Node_problem& b) {
//	bool result = false;
//	// 
//}


/****************************************************************************
* ������   : SortNode(int a,int b)
* ��  ��   : ��������������ʼ�ڵ��Ž���total_impulseָ���С��������
* �� ��    : 1.�ڵ��� 2.�ڵ���
****************************************************************************/
inline bool MultiTree::SortTNC(const TNC& a, const TNC& b)
{
		//double a_ave_gap = a.op_index_.time_cost;
		//double b_ave_gap = b.op_index_.time_cost;
		//if (a_ave_gap < b_ave_gap) {
		//	// ͬ��ʱ�俴�ĵ������Ϊ��
		//	return true;
		//}
		//else if (a_ave_gap == b_ave_gap) {
		//	if (a.op_index_.time_min < b.op_index_.time_min) {
		//		return true;
		//	}
		//	else {
		//		return false;
		//	}
		//}
		//else {
		//	return false;
		//}

	return a.op_index_.total_impulse_ < b.op_index_.total_impulse_;
	//return a.op_index_.time_cost < b.op_index_.time_cost;
	
}

/****************************************************************************
* ������   : EqualTNC(const TNC& a, const TNC& b)
* ��  ��   : �ж�����tnc�����Ƿ���ȣ��ٶ���������ӽ�
****************************************************************************/
inline bool MultiTree::EqualTNC(const TNC& a, const TNC& b)
{
	return  false;
}

/****************************************************************************
* ������   : void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
* ��  ��   : ȥ�غ���(֮ǰ����������õ�ǰW����ȥ�ص�)
****************************************************************************/
void MultiTree::unique_remove(std::vector<TNC>& expandinglist)
{
	std::vector<TNC> new_expandinglist;
	new_expandinglist.reserve(W_);

	for (int i = 0; i<expandinglist.size(); i++) //ÿ��ɾ�����С�����
	{
		//if ( !EqualTNC(expandinglist[i],expandinglist[i+1])) //���
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
* ������   : beam_discard(std::vector<TNC>& expandinglist)
* ��  ��   : ����������ȡǰW_��
****************************************************************************/
void MultiTree::beam_discard(std::vector<TNC>& expandinglist)
{
	//����h�Ӵ�С����
	sort(expandinglist.begin(), expandinglist.end(), MultiTree::SortTNC);

	unique_remove(expandinglist);

}

/****************************************************************************
* ������   : Initialize(std::vector<TNC>& expandinglist)
* ��  ��   : ��ʼ������չ�ڵ�
****************************************************************************/
void MultiTree::Initialize(std::vector<TNC>& expandinglist)
{
	layer_ = 0;

	//����ÿ�����ĳ�ʼ�������
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

	

	std::vector<Node*> old_expand_nodes;
	for (int i = 0; i < TreeNum; i++) {
		old_expand_nodes.push_back(multi_tree_[i].root_);
	}

	TNC tnctemp(old_expand_nodes);
	expandinglist.push_back(std::move(tnctemp));
	
}


/****************************************************************************
* ������   : Tree::getback(Node* node)
* ��  ��   : ���ݺ���������һ���ڵ�ֱ�����ڵ����Ϣ
****************************************************************************/
inline void Node::getback_problem(std::vector<Node*>& solution_one_node, Solution_one & temp_solution_one)
{

	//���һ�����ʼ�Ľڵ�
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
	for (int i = solution_one_node.size() - 1; i > -1; i--) //����ѭ�����Ӹ��ڵ㿪ʼ
	{
		int j = solution_one_node.size() - 1 - i; //��ʾ�ԳƵ�λ�ã������ڵ�ʱ��ʾΪ0
		Node_problem & now_node_problem = solution_one_node[i]->problem_; //��ǰ�ڵ������ṹ������
	}
}

/****************************************************************************
* ������   : Remove()
* ��  ��   : ɾ��һ���ڵ㼰�丸�ڵ㣬ǰ���Ǹýڵ�û���ӽڵ�
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
* ������   : GetBack()
* ��  ��   : ����TNC�����нڵ��ȫ����Ϣ
****************************************************************************/
Solution TNC::GetBack()
{
	Solution temp;
	for (int i =0; i< TreeNum; i++)
	{
		std::vector<Node*> solution_one_node;
	    tnc_[i]->return_node_sequence(solution_one_node);             
		tnc_[i]->getback_problem(solution_one_node, temp.solution_[i]);//�õ���Ƭ����
	}

	temp.total_impulse_ = op_index_.total_impulse_;
	temp.total_observed_num_ = op_index_.observed_num_;

	return temp;
}

/****************************************************************************
* ������   : Calculate_op_index()
* ��  ��   : ����ָ��ṹ��
****************************************************************************/
void TNC::Calculate_op_index()
{
	//����ָ��
	double total_impulse = 0.0;
	int  total_number = 0;
	Solution solution_temp = GetBack();
	double height_total = 0.0;
	int total_ship_num = 0;
	double time_min = 1.0e10;
	
	double time_last = -1.0e10;
	for (int i= 0; i< TreeNum; i++)
	{
		if(solution_temp.solution_[i].node_info_.size() == 0) continue;
		
		double time_temp = solution_temp.solution_[i].node_info_.back().time_acc_; //����ʱ��
		if (time_temp > time_last)  time_last = time_temp;
		if (time_temp < time_min) time_min = time_temp;

		for(int j = 0; j < solution_temp.solution_[i].node_info_.size(); j++)      //�۲����
		{
			if (solution_temp.solution_[i].node_info_[j].action_ == 2)
			{
				total_number++;
				height_total += V_Norm2(solution_temp.solution_[i].node_info_[j].rv_acc_, 3) - Re_km;
				if (solution_temp.solution_[i].node_info_[j].point_id_ == 21) total_ship_num++;
				if (solution_temp.solution_[i].node_info_[j].point_id_ == 4) total_ship_num++;
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
	op_index_.total_ship_num = total_ship_num;
	op_index_.time_min = time_min;
}

/****************************************************************************
* ������   : vector<Node*> Tree::ExpandNode(Node* node, int* visited_debris)
* ��  ��   : ��չһ���ڵ�
* �� ��    : �ýڵ�node���Լ��Ѿ��۲������visited���۲����Ϊ1������Ϊ0
* ��    �� : vector<Node*> һ����չ���б�
****************************************************************************/
std::vector<Node*> Tree::ExpandNode(Node* node, const int* visited, const std::vector<Node_problem> & problem_) //����������չ�ӽڵ�,���������Ѿ��ظ��Ľڵ�Ͳ���Ҫ�ظ���
{

	std::vector<Node*> expandnodes;                 //��չ�ӽڵ�

	for (int i = 0; i < problem_.size(); i++)
	{
		// ������problems_������ѡҪ��չ�Ľڵ�
		// ��չӦ���Ƕ�Ҫ��չ��
		if (visited[problem_[i].node_info_.back().point_id_] > 0 || V_Norm2(problem_[i].node_info_[0].dv_, 3) > dv_max_) continue; //�Ѿ����ʹ� ����dv_max

		int id = problem_[i].node_info_.back().point_id_;
		int t_60_int = static_cast<int>(round(problem_[i].node_info_.back().time_acc_/60.0));
		auto ifchild = std::find_if(node->child_.begin(), node->child_.end(),
			[id, t_60_int](Node* a) {return a->key_ == id*1e4 + t_60_int; }); //�жϸõ��Ƿ������ӽڵ���

		if (ifchild == node->child_.end()) //�������ӽڵ�
		{
			Node* temp_node = new Node(node, problem_[i]);
			expandnodes.push_back(temp_node);           //�����½ڵ�
		}
		else //���ڼ����������ӽڵ�
		{
			Node* temp_node = new Node(node, problem_[i]);
			expandnodes.push_back(temp_node);           //�����½ڵ�
		}
	}

	return expandnodes;
}


/****************************************************************************
* ������   : Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
* ��  ��   : ��չһ��TNC
* �� ��    : const TNC& tnc ԭTNC
* ��    �� : std::vector<TNC>& newTNCs  ����չ��tncs
****************************************************************************/

//rv0 ��dvǰ
void calculate_node_problem_info(int target_id, double t0, double rv0[6], double tf, double dv[3], double rvf[6], Node_problem & temp)
{
	OutputResult temp_out, temp_out2; //��������Ϣ�����������Ϣ

	temp_out.action_ = 1;
	for (int i = 0; i < 3; i++)  rv0[3 + i] += dv[i];
	//�������dv�������rv��Ϣ
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
		temp_out2.point_id_ = target_id + 1;
	}

	//���һ����tf���ӽڵ�
	temp.node_info_.push_back(temp_out2);
}

// TODO: ȥ��AccessPointObjects�ĵ��ã��ĳ���һ���
inline  void children_nodes(Node* node, const int* visited, std::vector<Node_problem>& child_node_problems)
{
	//��չ���2�ֿ��ܣ�1. ֱ�ӹ۲� 2.������۲⣨����м��пɼ������Ϊ��һ���ɼ��㣩
	//ÿһ�ֹ۲ⶼ�ǣ�����1+�۲�2
	int last_id = 100;
	if (node->problem_.node_info_.back().action_== 2)
		last_id = node->problem_.node_info_.back().point_id_ - 1; 
	int last2_id = 100;
	if (node->parent_ != nullptr) {
		last2_id = node->parent_->problem_.node_info_.back().point_id_ - 1;
	}

	double t0, rv_sat[6];
	t0 = node->problem_.node_info_.back().time_acc_;
	memcpy(rv_sat, node->problem_.node_info_.back().rv_acc_, 6 * sizeof(double));

	double rv_sat_after60[6];
	propagate_j2(rv_sat, rv_sat_after60, t0, t0 + 60.0);

	const double dt = 3600.0 * 3;
	std::vector<std::vector<double>> accesstime_low; // ������ɼ��Խ���б�
	AccessPointObjects(rv_sat_after60, t0 + 60.0, std::min(t0+ dt,86400.0 * 2.0), 60.0, 21, accesstime_low, 19.99 * D2R);

	double angle_out = 45.0 * D2R;
	std::vector<std::vector<double>> accesstime_large; // ������ɼ��Խ���б�
	AccessPointObjects(rv_sat_after60, t0 + 60.0, std::min(t0 + dt, 86400.0 * 2.0), 60.0, 21, accesstime_large, angle_out);

	//3Сʱ�ڱ����ɼ��Ĳ�����
	//3Сʱ����angle_out�ɼ��ģ��Ż���
	std::vector<std::pair<int, double>> id_time_low;
	std::vector<std::pair<int, double>> id_time_large;
	for (int target_id = 0; target_id < TargetNum; target_id++)
	{
		if (accesstime_low[target_id].size() > 0) {
			id_time_low.push_back(std::make_pair(target_id, accesstime_low[target_id][0]));
		}
		if (accesstime_large[target_id].size() > 0 && accesstime_low[target_id].size() == 0) {  //�ٳ���ԭ���ɼ���
			id_time_large.push_back(std::make_pair(target_id, accesstime_large[target_id][0]));
		}
	}


	// ��չ1���ɼ��ԣ�ֻ��չ��һ��
	if (!id_time_low.empty())
	{
		// ����id_time_low��ʱ����С�ĵ�
		auto min_time_point = id_time_low[0];
		for (const auto& point : id_time_low) {
			//if (point.first != last_id && point.first != last2_id) //��ֹ��֮ǰ���������ظ���������һ��
			{
				if (point.second < min_time_point.second)
				{
					min_time_point = point;
				}
			}
		}
		int target_id = min_time_point.first;
		double tf = min_time_point.second;
		Node_problem node_problem;
		double dv[3] = { 0.0, 0.0,0.0 }; //Ϊ��������������һ��
		double rvf[6]; 
		propagate_j2(rv_sat, rvf, t0, tf);
		calculate_node_problem_info(target_id, t0, rv_sat, tf, dv, rvf, node_problem);
		child_node_problems.push_back(node_problem);
	}
	else
	{
		std::vector<std::vector<double>> accesstime_long; // ������ɼ��Խ���б�
		AccessPointObjects(rv_sat_after60, t0 + 60.0,86400.0 * 2.0, 60.0, 21, accesstime_long, 19.99 * D2R);

		// ����id_time_low��ʱ����С�ĵ�
		std::pair<int, double> min_time_point = std::make_pair(100, 1.0e10);
		for (int id_x = 0; id_x < accesstime_long.size(); id_x++)
		{
			if (accesstime_long[id_x].size() > 0)
			{
				if (accesstime_long[id_x][0] < min_time_point.second)
				{
					min_time_point = std::make_pair(id_x, accesstime_long[id_x][0]);
				}
			}
		}
		double dv[3] = { 0.0, 0.0,0.0 }; //Ϊ��������������һ��
		double rvf[6];
		propagate_j2(rv_sat_after60, rvf, t0 + 60.0, min_time_point.second);
		Node_problem node_problem;
		calculate_node_problem_info(min_time_point.first, t0, rv_sat, min_time_point.second, dv, rvf, node_problem);
		child_node_problems.push_back(node_problem);
	}

	// ��չ2��������ɼ�
	for (int index_id = 0; index_id < id_time_large.size(); index_id++)
	{
		int target_id = id_time_large[index_id].first;
		double t_f_possible = id_time_large[index_id].second;
		if (visited[target_id] > 0)
		{
			continue;
		}

		double target_r0[3];
		get_target_R(target_id, t_f_possible, target_r0);

		double rv0[6];
		memcpy(rv0, node->problem_.node_info_.back().rv_acc_, 6 * sizeof(double));
		t0 = node->problem_.node_info_.back().time_acc_;

		//���Ҫ�ؼ��
		int flag;
		double coe0[6], a, e, hmin, hmax;
		rv2coe(flag, coe0, rv0, mu_km_s);
		a = coe0[0]; e = coe0[1];
		hmin = a * (1 - e) - Re_km;
		hmax = a * (1 + e) - Re_km;
		if (hmin < 201.0 || hmax > 999.0) {
			continue;
		}

		//�������������
		auto dv_infos = obs_shooting_zzmodified(t0, rv0, target_id);

		//ע��˴�dv_infos�ж��ֿ��ܣ������Ҫ���治ͬʱ�����С�ٶ�������
		// ʹ�ù�ϣ��(key=ʱ��tf(ȡ����), value=��ǰʱ���¶�������С�� dv_shooding_info)
		std::unordered_map<int, dv_shooding_info> unique_dv_infos;
		unique_dv_infos.reserve(dv_infos.size());
		for (auto& info : dv_infos){
			auto it = unique_dv_infos.find(static_cast<int>(round(info.tf/60.0)));
			if (it == unique_dv_infos.end()){ // ��û�г��ֹ����ʱ�䣬��ֱ�Ӳ���
				unique_dv_infos[static_cast<int>(round(info.tf / 60.0))] = info;
			}
			else { // ������ͬʱ�䣬��˭���ٶ�������������С
				if (V_Norm2(info.dv, 3) < V_Norm2(it->second.dv, 3)) {
					it->second = info;
				}
			}
		}
		// �����ȥ�غ�Ľ����ŵ�һ���� vector
		dv_infos.clear();
		for (auto& kv : unique_dv_infos)
		{
			dv_infos.push_back(kv.second);
		}

		for(int i = 0 ; i< dv_infos.size(); i++)
		{
			double tf = dv_infos[i].tf;
			double dv[3];
			memcpy(dv, dv_infos[i].dv, 3 * sizeof(double));
			double rvf[6]; memcpy(rvf, rv0, 6 * sizeof(double));

			// ���ʱ�����
			double dt_max = 3.0 * 3600.0;  //��ֵ
			if (tf - t0 > dt_max) {continue;}
			double dv_max = 1.0; //1km/s
			if(V_Norm2(dv, 3) > dv_max) { continue; }

			// ����������ֻ���ֵ�172800s
			if (tf > 2.0 * 86400.0) {
				tf = 2.0 * 86400.0;
				double rv1[6];
				memcpy(rv1, rv0, 6 * sizeof(double));
				for (int i = 0; i < 3; i++) rv1[i + 3] += dv[i];
				propagate_j2(rv1, rvf, t0, tf);
			}

			//����м��Ƿ��пɼ������
			double rv1[6];
			memcpy(rv1, rv0, 6 * sizeof(double));
			for (int i = 0; i < 3; i++) rv1[i + 3] += dv[i];

			propagate_j2(rv1, rv1, t0, t0 + 60.0);

			std::vector<std::vector<double>> accesstime_between; // ������ɼ��Խ���б�
			AccessPointObjects(rv1, t0 + 60.0, tf-60.0, 60.0, 21, accesstime_between, 19.99 * D2R);

			// ����id_time_low��ʱ����С�ĵ�
			std::pair<int, double> min_time_point = std::make_pair(100, 1.0e10);
			for(int id_x = 0; id_x < accesstime_between.size(); id_x++)
			{
				if (accesstime_between[id_x].size() > 0)
				{
					if (accesstime_between[id_x][0] < min_time_point.second)
					{
						min_time_point = std::make_pair(id_x, accesstime_between[id_x][0]);
					}
				}
			}
			if(min_time_point.first < 100)
			{
				propagate_j2(rv1, rvf, t0 + 60.0, min_time_point.second);
				Node_problem node_problem;
				calculate_node_problem_info(min_time_point.first, t0, rv0, min_time_point.second, dv, rvf, node_problem);
				child_node_problems.push_back(node_problem);
			}
			else
			{
				Node_problem node_problem;
				calculate_node_problem_info(target_id, t0, rv0, tf, dv, rvf, node_problem);
				child_node_problems.push_back(node_problem);
			}
		}
	}



	/*if (child_node_problems.empty()) {
		std::cout << "�޽ڵ����չ��" << std::endl;
	}*/
	//std::sort(child_node_problems.begin(), child_node_problems.end(), sort_Nodes);
}

/****************************************************************************
* ������   : Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
* ��  ��   : ��չһ��TNC
* �� ��    : const TNC& tnc ԭTNC
* ��    �� : std::vector<TNC>& newTNCs  ����չ��tncs 
****************************************************************************/
void MultiTree::Expansion_one_TNC(const TNC& tnc, std::vector<TNC>& newTNCs)
{
	//���ǵ������нڵ㶼������չ������չ������tnc���п�����ͬ��
	newTNCs.clear();
	std::vector<TNC> delete_tnc;

	int visited[TargetNum]{};                     //��tnc���ѹ۲����У���ʼ��Ϊ0������ط������Ϊ1
	std::vector<std::vector<double>> visible_timelist(21, std::vector<double>(0));	//��tnc��Ŀ������ʱ�̱�
	//std::cout << visible_timelist.size() << std::endl;

	//TODO����չ�����ģ��˴���ʱ������visited
	double tf_max = 0.0;
	double tf_min = 1.0e10;
	for (int j = 0; j < TreeNum; j++)             //�����ʹ������нڵ㰴Ŀ��λ��+1
	{
		std::vector<Node*> solution_one_node;
		tnc.tnc_[j]->return_node_sequence(solution_one_node);
		Solution_one temp;
		tnc.tnc_[j]->getback_problem(solution_one_node, temp);
		if(temp.node_info_.back().time_acc_ < tf_min) tf_min = temp.node_info_.back().time_acc_;
		if(temp.node_info_.back().time_acc_ > tf_max) tf_max = temp.node_info_.back().time_acc_;

		for (int k = 0; k < temp.node_info_.size(); k++) {
			if (temp.node_info_[k].action_ == 2) {
				const int id = temp.node_info_[k].point_id_ - 1;
				double rv_standard[6], t_standard, dv[3], t_period[2];
				t_standard = temp.node_info_[k].time_acc_;
				
				memcpy(rv_standard, temp.node_info_[k].rv_acc_, 6 * sizeof(double));
				if (k < temp.node_info_.size() - 1 && temp.node_info_[k + 1].action_ == 1) {
					memcpy(dv, temp.node_info_[k + 1].dv_, 3 * sizeof(double));
				}
				else {
					for (int w = 0; w < 3; w++) dv[w] = 0.0;
				}

				AccessPeriodCertainObject(rv_standard, t_standard, dv, id, t_period);
				if (t_period[0] > 172800.0) {
					// ATK�޷���⣬����
					continue;
				}
				else {
					visible_timelist[id].push_back(t_period[0]);
					visible_timelist[id].push_back(t_period[1]);
				}

			}
		}

		if (temp.node_info_[temp.node_info_.size() - 1].time_acc_ < tf_min) {
			tf_min = temp.node_info_[temp.node_info_.size() - 1].time_acc_;
		}

		if (temp.node_info_[temp.node_info_.size() - 1].time_acc_ > tf_max) {
			tf_max = temp.node_info_[temp.node_info_.size() - 1].time_acc_;
		}
	}

	// Push back֮��Կɼ���ʱ�̱�����
	for (int j = 0; j < TargetNum; j++) {
		std::sort(visible_timelist[j].begin(), visible_timelist[j].end());
	}
	

	// �����Ŀ��������ط�ʱ���Ѿ�����������Ҫ������ֹͣ��չ
	// ��һ��Ӧ���ȸ�8���ǵĿ����Զ���һ����֤��ȷ��8��������չ��ʱ������1.5hΪֹ
	// ����ǰʱ����̵�����ʱ��Ϊֹ����������طã���ô�����Ѿ���������ɣ���TNCֹͣ��չ
	std::vector<double> max_revisit_gap;
	max_reseetime(visible_timelist, 0.0, tf_min, max_revisit_gap);			// ע���ط�ʱ�䵥λ�����h
	bool ifpossible = true;
	int idx = 0;
	for (auto iter = max_revisit_gap.begin(); iter != max_revisit_gap.end(); iter++) {
		idx++;
		if (idx != TargetNum) {
			if (*iter > 6.0) {
				//std::cout << "Ŀ��" << idx << "������ط�ʱ���Ѿ��ﵽ" << *iter << "h" << std::endl;
				ifpossible = false;
			}
			
		}
		else {
			if (*iter > 3.0) {
				//std::cout << "Ŀ��" << idx << "������ط�ʱ���Ѿ��ﵽ" << *iter << "h" << std::endl;
				ifpossible = false;
			}
			
		}
	}
	//std::cout << std::endl;
	if (!ifpossible) {
		return;
	}

	// ��չ��֮����¸�ʱ�̱������ط�ʱ�䣬����ط�ʱ������Ҫ����visited��1
	std::vector<double> max_revisit_list;
	max_reseetime(visible_timelist, max_revisit_list);
	for (int i = 0; i < TargetNum; i++) {
		if (i != 20) {
			if (max_revisit_list[i] < 6.0) {
				visited[i]++;
				//std::cout << "Ŀ��" << i + 1 << "�۲������" << std::endl;
			}
		}
		else {
			if (max_revisit_list[i] < 3.0) {
				visited[i]++;
				//std::cout << "Ŀ��" << i + 1 << "�۲������" << std::endl;
			}
		}
	}

	bool ifdone = true;
	for (int i = 0; i < TargetNum; i++) {
		if (visited[i] == 0) ifdone = false;
	}
	ifFinished_ = ifdone;
	if(ifdone) std::cout << "�۲������Ѿ�ȫ�����" << std::endl;

	double time[TreeNum]; double time_min = 1.e20; int counter = 0; //��չʱ����С���ǿ���
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

	if (tnc.op_index_.total_impulse_ > TreeNum * dv_max - 1.0e-3) {
		return;// ���峬����ֹͣ��չ��TNC
	}

	for (int i = 0; i < TreeNum; i++) //������չ
	{
		if (i != counter) continue;

		std::vector<Node_problem> child_node_problems;
		children_nodes(tnc.tnc_[i], visited, child_node_problems);     //�ų�����Щ�ܴ��ٶ���������Ƭ TODO: ��Ҫ��


		/*if (child_node_problems.empty()) {
			std::cout << "û�б�ѡ�ӽڵ�" << std::endl;
		}*/

		omp_set_lock(&tnc.tnc_[i]->lock);
		std::vector<Node*> leafnode = multi_tree_[i].ExpandNode(tnc.tnc_[i], visited, child_node_problems); //�ýڵ���չ���ӽڵ�
		omp_unset_lock(&tnc.tnc_[i]->lock);

		/*if (leafnode.empty()) {
			std::cout << "û�п���չ���ӽڵ�" << std::endl;
		}*/
		
		for (int j = 0; j < leafnode.size(); j++) //�ýڵ���չ���ӽڵ���Ϣ
		{
			Node* temp[TreeNum];
			for (int i = 0; i < TreeNum; i++) temp[i] = tnc.tnc_[i];
			temp[i] = leafnode[j];
			TNC temp_tnc(temp);
			newTNCs.push_back(std::move(temp_tnc));
		}
		
	}

	//��չ�󣬻����TNC��ָ��
	for (int i = 0; i < newTNCs.size(); i++)
	{
		newTNCs[i].Calculate_op_index();
	}

	//��ɸѡ
	sort(newTNCs.begin(), newTNCs.end(), MultiTree::SortTNC);
	//��ֵ��expand�б�
	if (newTNCs.size() > b_)
	{
		newTNCs.erase(newTNCs.begin() + b_, newTNCs.end());
	}
}

/****************************************************************************
* ������   : Expansion(std::vector<TNC>& expandinglist) ��չ�׶κ���
* ��  ��   : ������չ��Ϣ�õ��µ���չ�б�
****************************************************************************/
void MultiTree::Expansion(std::vector<TNC>& expandinglist)
{
	layer_++;                                             //��ǰ����
	std::vector< std::vector<TNC>> new_childnode_private;               //�ֲ�˽�б���
	bool ifPossible = false;
	
#pragma omp parallel
	{
		auto nthreads = omp_get_num_threads();
		auto id = omp_get_thread_num();

#pragma omp single
		{
			new_childnode_private.resize(nthreads);
		}

		//������Ҫ��չ�Ľڵ�
		
#pragma omp for schedule(dynamic) 
		for (int i = 0; i < expandinglist.size(); i++)
		{
			std::vector<TNC> newTNCs;
			Expansion_one_TNC(expandinglist[i], newTNCs);
			if (!newTNCs.empty()) {
				ifPossible = true;
				move(newTNCs.begin(), newTNCs.end(), back_inserter(new_childnode_private[id]));
			}
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
		std::cout << "�۲��Ѿ���������ɣ��˳���" << std::endl;
		ifFinished_ = true;
	}
}


/****************************************************************************
* ������   : Run()
* ��  ��   : ���������������ܺ���
****************************************************************************/
void MultiTree::Run()
{
	for (int i = 0; i < TreeNum; i++)  multi_tree_[i].dv_max_ = dv_max_; 

	std::ofstream fout0("../output_result/result.txt");
	
	std::vector<TNC> expandinglist;
	Initialize( expandinglist);                         //��ʼ���״���չ��
	
	while (expandinglist.size() > 0 && !ifFinished_) //�ǿ�
	{
		Expansion(expandinglist);

		beam_discard(expandinglist);                    //����ȡǰW_��

		//delete_redundant_nodes();                       //ɾ��������ӽڵ�

		RecordBestResult(expandinglist,fout0);			//��¼�����Ϣ

		if (ifFinished_)
		{
			//std::cout << "�۲������Ѿ�ȫ�����" << std::endl;
			std::vector<TNC> x;
			RecordBestResult(x, fout0);			         //��¼�����Ϣ
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
	std::cout << "Removal total num, Time, Impulse:  " << result_now_.total_observed_num_ << " " << result_now_.time_cost_ << " " << result_now_.total_impulse_ << std::endl;
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
	
	//���ո�ʽ������ս��
	
	fout1 << std::endl << std::endl << "���������" << result_now_.total_observed_num_ << " ��ǰ�ڵ�������" << expandinglist.size() << " ��ǰʱ�䣺" << result_now_.time_cost_ << " ��ǰ����������" << result_now_.total_impulse_ << std::endl;
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


