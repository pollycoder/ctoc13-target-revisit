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
* ������   : SortNode(int a,int b)
* ��  ��   : ��������������ʼ�ڵ��Ž���hָ��Ӵ�С����
* �� ��    : 1.�ڵ��� 2.�ڵ���
****************************************************************************/
inline bool MultiTree::SortTNC(const TNC& a, const TNC& b)
{
	double h_a, h_b;
	h_a = 1.0 / pow(a.op_index_.time_cost, 1);
	h_b = 1.0 / pow(b.op_index_.time_cost, 1);
	return  h_a > h_b;
}

/****************************************************************************
* ������   : EqualTNC(const TNC& a, const TNC& b)
* ��  ��   : �ж�����tnc�����Ƿ���ȣ��ٶ���������ӽ�
****************************************************************************/
inline bool MultiTree::EqualTNC(const TNC& a, const TNC& b)
{
	//if (fabs(a.op_index_.total_mass_ - b.op_index_.total_mass_) < 1.0e-10)
	//{
	//	//�õ���������
	//	return true;
	//}
	
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

	for (int i = 0; i<expandinglist.size()-1; i++) //ÿ��ɾ�����С�����
	{
		if ( !EqualTNC(expandinglist[i],expandinglist[i+1])) //���
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
		temp_output.branch_ = 0;
		temp_output.m_ = 100000.0;
		temp_output.rv_acc_;
		temp_output.action_ = 0;
		double sats_rv0[6]; int flag;
		coe2rv(flag, sats_rv0, sats_coe0[i], mu_km_s);
		memcpy(temp_output.rv_acc_, sats_rv0, 6 * sizeof(double));

		multi_tree_[i].root_->problem_.node_info_.push_back(temp_output);
	}

	//���tnc���
	std::vector<Node*> old_expand_nodes;
	for (int i=0; i< TreeNum; i++)
	{
		old_expand_nodes.push_back(multi_tree_[i].root_);
	}

	TNC temp(old_expand_nodes);
	expandinglist.push_back(std::move(temp));
	

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

	temp.total_mass_ = op_index_.total_mass_;
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
	double total_mass = 0.0;
	int  total_number = 0;
	Solution solution_temp = GetBack();
	double height_total = 0.0;
	
	double time_last = -1.0e10;
	for (int i= 0; i< TreeNum; i++)
	{
		if(solution_temp.solution_[i].node_info_.size() == 0) continue;
		
		double time_temp = solution_temp.solution_[i].node_info_.back().time_acc_; //����ʱ��
		if (time_temp > time_last)  time_last = time_temp;

		for(int j = 0; j < solution_temp.solution_[i].node_info_.size(); j++)      //�۲����
		{
			if (solution_temp.solution_[i].node_info_[j].action_ == 2)
			{
				total_number++;

				height_total += V_Norm2(solution_temp.solution_[i].node_info_[j].rv_acc_, 3) - Re_km;
			}
		}

		total_mass += solution_temp.solution_[i].node_info_.back().m_;             //��������
	}

	op_index_.time_cost = time_last;
	op_index_.height_aver = height_total / total_number;
	op_index_.total_mass_ = total_mass;
	op_index_.observed_num_ = total_number;
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
		if (visited[problem_[i].node_info_.back().point_id_] > 0 || V_Norm2(problem_[i].node_info_[0].dv_,3) > dv_max_) continue; //�Ѿ����ʹ� ����dv_max

		int id = problem_[i].node_info_.back().point_id_;
		
		auto ifchild = std::find_if(node->child_.begin(), node->child_.end(),
			[id](Node* a) {return a->key_ == id; }); //�жϸõ��Ƿ������ӽڵ���

		if (ifchild == node->child_.end()) //�������ӽڵ�
		{
			Node* temp_node = new Node(node, problem_[i]);
			expandnodes.push_back(temp_node);           //�����½ڵ�
		}
		else //���ڼ����������ӽڵ�
		{
			expandnodes.push_back(*ifchild);
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
inline  void children_nodes(Node* node, const int* visited, std::vector<Node_problem>& child_node_problems)
{
	for (int j = 0; j < GapNum; j++)
	{
		// �����ͽ�����Ҫ�㣬branchû���ڴ����ڵ�ʱָ����ֻ���������޸ģ�����child_node_problems�ڵĽڵ���ȷ����
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

			// ��γ������ָ��������visit_gap
			double next_target_geodetic[2];
			int id = static_cast<int>(visit_gap[j][0]);
			get_target_geogetic(id, visit_gap[j][1], next_target_geodetic);
			lambda0 = next_target_geodetic[1];
			phi0 = next_target_geodetic[0];

			double norm_r = V_Norm2(rv0, 3);      //����߶�Լ��
			if (norm_r - Re_km < 200.0) continue;

			double mass_end = 1000.0;
			double rvf[6];

			// ��е���һ��Ŀ�꣬����visit_gap
			shooting_target2target(t0, rv0, visit_gap[j][1], lambda0, phi0, tf, dv, rvf, branch);

			//TODO: ��������
			Node_problem temp;		// temp_D��������������ֶ�Ҫpush back
			OutputResult temp_out, temp_out2; //��������Ϣ�����������Ϣ

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

	int visited[GapNum]{};                     //��tnc���ѹ۲����У���ʼ��Ϊ0���۲��Ϊ1


	for (int j = 0; j < TreeNum; j++)             //�����ʹ������нڵ㰴Ŀ��λ��+1
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

	double time[TreeNum]; double time_min = 1.e20; int counter = 0; //��չʱ����С���ǿ���
	for (int i = 0; i < TreeNum; i++)
	{
		time[i] = tnc.tnc_[i]->problem_.node_info_.back().time_acc_;

		if (time[i] < time_min - 1.0e-3) {
			time_min = time[i]; counter = i;
		}
	}

	for (int i = 0; i < TreeNum; i++) //������չ
	{
		if (i != counter) continue;

		std::vector<Node_problem> child_node_problems;
		children_nodes(tnc.tnc_[i], visited, child_node_problems);     //�ų�����Щ�ܴ��ٶ���������Ƭ TODO: ��Ҫ��

		omp_set_lock(&tnc.tnc_[i]->lock);
		std::vector<Node*> leafnode = multi_tree_[i].ExpandNode(tnc.tnc_[i], visited, child_node_problems); //�ýڵ���չ���ӽڵ�
		omp_unset_lock(&tnc.tnc_[i]->lock);
		
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
* ������   : Run()
* ��  ��   : ���������������ܺ���
****************************************************************************/
void MultiTree::Run()
{
	for (int i = 0; i < TreeNum; i++)  multi_tree_[i].dv_max_ = dv_max_; 

	std::ofstream fout0("../output_result/result.txt");
	
	std::vector<TNC> expandinglist;
	Initialize( expandinglist);                         //��ʼ���״���չ��
	
	while (expandinglist.size() > 0 && layer_ < GapNum) //�ǿ�
	{
		Expansion(expandinglist);

		beam_discard(expandinglist);                    //����ȡǰW_��

		delete_redundant_nodes();                       //ɾ��������ӽڵ�

		RecordBestResult(expandinglist,fout0);			//��¼�����Ϣ

		if (layer_ == GapNum)
		{
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
	
	//���ո�ʽ������ս��
	
	fout1 << std::endl << std::endl << "���������" << result_now_.total_observed_num_ << "���ҽڵ�����" << expandinglist.size() << "��ǰָ��" << result_now_.time_cost_ << "����" << result_now_.total_mass_ << std::endl;
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


