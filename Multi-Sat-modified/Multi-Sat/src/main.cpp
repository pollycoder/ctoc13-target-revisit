/****************************************************************************
* Copyright (C), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
* ����: ���� zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* �ļ���: main.cpp
* ���ݼ�����GTOC9���������������򣬺��޸�����CTOC13����
*
* �ļ���ʷ��
* �汾��     ����         ����       ˵��
* 01       2021-05-12    ����      �����ļ�
* 02	   2024-09-28	 ��ܲ	   ���Թ۲��к���
* 03	   2024-09-30	 ��		   ����������
* 04	   2024-10-15	 ��		   ��ʼ��ֵ̽��
* 05	   2024-12-03	 ��		   �Զ�ȫ���Ż�����
****************************************************************************/
#include <stdlib.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"
#include "J2propagation.h"
#include "single_sat.h"
#include <iomanip>
#include "main.h"

const std::string space = " ";
std::vector<std::vector<std::vector<double>>> fixed_imp = { {}, {}, {}, {}, {}, {}, {}, {} };
double t_anchor = 45000.0;
double t_gap = 21600.0;

void output_result(std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, const double t)
{
	std::ofstream fout0("../output_result/result.txt", std::ios::app);
	fout0 << "Breakpoint " << ": " << t << "s" << std::endl;
	for (const auto& sat : sat_info_list) {
		// ��д�����
		for (const auto& coe : std::get<0>(sat)) {
			fout0 << coe << " ";
		}
		fout0 << std::endl;

		// ��д������
		if (std::get<1>(sat).empty()) {
			fout0 << std::setprecision(14) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << std::endl;
		}
		else {
			for (const auto& row : std::get<1>(sat)) {
				fout0 << std::setprecision(14) << row[0] << " " << row[1] << " " << row[2] << " " << row[3] << std::endl;
			}
		}
	}
	fout0 << std::endl;
}


// �����������������
void multitree_search() {
	MultiTree multi_tree(1000, 4, 100, dv_max);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "�ܺ�ʱ��" << duration_second << "��" << std::endl;
	//0.01s��չһ���ڵ�
	//W=10000����չһ��130s
	//����176��Լ22880s��Ԥ��7h����
	//������չԼ300-400�㣬һ���ڱ������
}


void multi_sat_opt() {
	std::ofstream fout0("../output_result/result.txt");
	fout0 << "" << std::endl;
	int var_num = std::accumulate(imp_num.begin(), imp_num.end(), 0) * 4;
	std::vector<double> X;
	double f;
	std::vector<std::vector<double>> AccessTable;
	std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>> sat_info_list;

	
	for (int i = 0; i < var_num; i++) {
		if (i % 4 == 0) {
			X.push_back(0.5);
		}
		else X.push_back(0.5);
		//X.push_back(0.5);
	}

	int idx = 0;
	while (true) {
		idx++;
		DE_parallel(obj_func, nullptr, X, f, var_num, 80, 3000, 100);
		nlopt_main(obj_func, nullptr, X, f, var_num, 0, 1);

		get_score_info(X, nullptr, f, sat_info_list, AccessTable);
		for (int i = 0; i < fixed_imp.size(); i++) {
			if (!std::get<1>(sat_info_list[i]).empty()) {
				fixed_imp[i].push_back(std::get<1>(sat_info_list[i]).back());
			}
		}
		t_anchor += t_gap;
		if (f <= 3.0) {
			std::cout << "Pipeline " << idx << " Succeeded." << std::endl;
		}
		else {
			std::cout << "Revisiting task failed." << std::endl;
			break;
		}


		std::cout << "f = " << f << std::endl;
		std::vector<double> revisit_gap;
		max_reseetime(AccessTable, revisit_gap);
		for (const auto& t : revisit_gap) {
			std::cout << t << std::endl;
		}

		output_result(sat_info_list, t_anchor - t_gap);
	}
	// 16�ˣ�2h�Ż�һ��
	// Ԥ��8h���
}



void max_revisit_verify() {
	std::vector<std::vector<double>> rv0_list;
	std::vector<std::vector<double>> result;
	for (int i = 0; i < TreeNum; i++) {
		std::vector<double> rv0;
		double rv0temp[6];
		int flag;
		coe2rv(flag, rv0temp, sats_coe0[i], mu_km_s);
		for (int j = 0; j < 6; j++) rv0.push_back(rv0temp[j]);
		rv0_list.push_back(rv0);
	}
	MultiSat_AccessPointObjects(rv0_list, 0.0, 172800.0, 60.0, 21, result);
	std::vector<double> max_revisit;
	max_reseetime(result, max_revisit);
	for (auto row = result.begin(); row != result.end(); row++) {
		for (auto it = row->begin(); it != row->end(); it++) {
			std::cout << *it << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
	std::cout << "Max revisit gap:" << std::endl;
	for (auto it = max_revisit.begin(); it != max_revisit.end(); it++) {
		std::cout << *it << std::endl;
	}
}
 
void shooting_test()
{
	double rv0[6], rvf[6], geo[2];
	int flag, NR;

	get_target_geogetic(3, 0.0, geo);

	for (int id = 0; id < TargetNum; id++) {
		std::cout << "Target " << id + 1 << std::endl;
		for (int i = 0; i < 4; i++) {
			coe2rv(flag, rv0, sats_coe0[i], mu_km_s);
			std::cout << "============================================" << std::endl;
			std::cout << "Satellite " << i << std::endl;
			std::cout << "--------------------------------------------" << std::endl;
			for (int N = 1; N < 21; N++) {
				for (int bra = 0; bra < 2; bra++) {
					double m0 = 1000.0, mf, tf, dv[3];
					NR = N;
					//single_imp(m0, 0.0, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, bra);
					obs_shooting(flag, dv, tf, rvf, 0.0, rv0, id, NR, bra);
					if (flag && V_Norm2(dv, 3) < 1.0 && V_Norm2(dv, 3) > 0.001) {
						std::cout << "NR = " << NR << ", branch = " << bra << std::endl;
						std::cout << "tf = " << tf << std::endl;
						std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
						std::cout << "imp = " << V_Norm2(dv, 3) << std::endl;
						std::cout << "--------------------------------------------" << std::endl;
					}
				}
			}
			std::cout << "============================================" << std::endl;
			std::cout << std::endl;
		}
	}
}

int main() {
	
	//multi_sat_opt();
	shooting_test();
	
	return 0;
}
