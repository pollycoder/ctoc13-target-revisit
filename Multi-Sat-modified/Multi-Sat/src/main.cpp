/****************************************************************************
* Copyright (C), 2020-2031 清华大学航天航空学院动力学与控制实验室
* 作者: 张众 zhong-zh19@mails.tsinghua.edu.cn
*            539977562@qq.com
* 文件名: main.cpp
* 内容简述：GTOC9并行算例的主程序，后被修改适配CTOC13丙题
*
* 文件历史：
* 版本号     日期         作者       说明
* 01       2021-05-12    张众      创建文件
* 02	   2024-09-28	 周懿	   测试观测打靶函数
* 03	   2024-09-30	 周		   测试树搜索
* 04	   2024-10-15	 周		   开始初值探索
* 05	   2024-12-03	 周		   自动全局优化脉冲
****************************************************************************/
#include <stdlib.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"
#include "J2propagation.h"
#include "single_sat.h"
#include <iomanip>
#include "main.h"
#include <chrono>
#include <unordered_map>

#include "visibility_4_targets.h"


const std::string space = " ";
std::vector<std::vector<std::vector<double>>> fixed_imp = { {}, {}, {}, {}, {}, {}, {}, {} };
double t_anchor = 45000.0;
double t_gap = 21600.0;

void output_result(std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list, const double t)
{
	std::ofstream fout0("../output_result/result.txt", std::ios::app);
	fout0 << "Breakpoint " << ": " << t << "s" << std::endl;
	for (const auto& sat : sat_info_list) {
		// 先写入根数
		for (const auto& coe : std::get<0>(sat)) {
			fout0 << coe << " ";
		}
		fout0 << std::endl;

		// 再写入脉冲
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


// 多棵树搜索（重启）
void multitree_search() {
	MultiTree multi_tree(1, 4, 100, dv_max);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "总耗时：" << duration_second << "秒" << std::endl;
	//0.01s扩展一个节点
	//W=10000，扩展一层130s
	//至少176层约22880s，预计7h到达
	//完整扩展约300-400层，一天内必须完成
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
	// 16核，2h优化一轮
	// 预计8h完成
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

	coe2rv(flag, rv0, sats_coe0[0], mu_km_s);
	std::vector<std::vector<double>> results_vector; // 输出：可见性结果列表

	AccessPointObjects(rv0, 0.0, 3600.0*3, 60.0, 21, results_vector, 45.0 * D2R);

	//AccessPointObjects_linearJ2(rv0, 0.0, 3600.0 * 3, 60.0, 21, results_vector, 19.5 * D2R);

	int id = 0;
	int N = 0;
	int bra = 1;
	//{  //从0开始
	//	get_target_geogetic(id, 0.0, geo);
	//	coe2rv(flag, rv0, sats_coe0[0], mu_km_s);




	//	double m0 = 1000.0, mf, tf, dv[3];
	//	NR = N;
	//	{
	//		single_imp(m0, 0.0, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, bra);
	//		single_imp_zzmodified(m0, 0.0, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, bra);
	//		//积分过去看情况
	//		double rvf_J2[6];
	//		double rv0_afterdv[6];
	//		memcpy(rv0_afterdv, rv0, 6 * sizeof(double));
	//		V_Add(rv0_afterdv + 3, rv0_afterdv + 3, dv, 3);
	//		propagate_j2(rv0_afterdv, rvf_J2, 0.0, tf, 1.0e-10, 1.0e-10);

	//		V_Multi(rvf_J2, rvf_J2, Re_km / V_Norm2(rvf_J2, 3), 3);
	//		double v_norm2 = V_Norm2(dv, 3);
	//		double rvf_target[6];
	//		get_target_R(id, tf, rvf_target);
	//		double result = 0;
	//	}
	//}


	for (int id = 0; id < TargetNum; id++) {
		std::cout << "Target " << id << std::endl;

		get_target_geogetic(id, 0.0, geo);
		for (int i = 0; i < 1; i++) {
			coe2rv(flag, rv0, sats_coe0[i], mu_km_s);

			//20241223 ZZ test:
			double dt = 0000.0;
			rv02rvf(flag, rv0, rv0, dt, mu_km_s);


			std::cout << "============================================" << std::endl;
			std::cout << "Satellite " << i << std::endl;
			std::cout << "--"
				"------------------------------------------" << std::endl;
			for (int bra = 1; bra < 2; bra++) {  //从0开始
				int N = 0;  { //ZZ TEST
				//for (int N = 0; N < 21; N++) {
					double m0 = 1000.0, mf, tf, dv[3];
					NR = N;
					//{
					//	single_imp(m0, dt, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, bra);
					//	single_imp_zzmodified(m0, dt, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, bra);
					//	//积分过去看情况
					//	double rvf_J2[6];
					//	double rv0_afterdv[6];
					//	memcpy(rv0_afterdv, rv0, 6 * sizeof(double));
					//	V_Add(rv0_afterdv + 3, rv0_afterdv + 3, dv, 3);
					//	propagate_j2(rv0_afterdv, rvf_J2, 0.0, tf, 1.0e-10, 1.0e-10);

					//	V_Multi(rvf_J2, rvf_J2, Re_km / V_Norm2(rvf_J2, 3), 3);
					//	double v_norm2 = V_Norm2(dv, 3);
					//	double rvf_target[6];
					//	get_target_R(id, tf, rvf_target);
					//	double result = 0;
					//}
					//{
					//	dv[0] = 20.0;
					//	single_imp_zzmodified(m0, 0.0, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, bra);
					//	//积分过去看情况
					//	double rvf_J2[6];
					//	double rv0_afterdv[6];
					//	memcpy(rv0_afterdv, rv0, 6 * sizeof(double));
					//	V_Add(rv0_afterdv + 3, rv0_afterdv + 3, dv, 3);
					//	propagate_j2(rv0_afterdv, rvf_J2, 0.0, tf, 1.0e-10, 1.0e-10);

					//	V_Multi(rvf_J2, rvf_J2, Re_km / V_Norm2(rvf_J2, 3), 3);

					//	double rvf_target[6];
					//	get_target_R(id, tf, rvf_target);
					//	double result = 0;
					//}


					//obs_shooting(flag, dv, tf, rvf, dt, rv0, id, NR, bra);

					auto dv_infos = obs_shooting_zzmodified(dt, rv0, id);

					//注意此处dv_infos有多种可能，因此需要保存不同时间的最小速度增量的
					// 使用哈希表(key=时间tf(取分钟), value=当前时间下二范数最小的 dv_shooding_info)
					std::unordered_map<int, dv_shooding_info> unique_dv_infos;
					unique_dv_infos.reserve(dv_infos.size());
					for (auto& info : dv_infos) {
						auto it = unique_dv_infos.find(static_cast<int>(round(info.tf / 60.0)));
						if (it == unique_dv_infos.end()) { // 还没有出现过这个时间，就直接插入
							unique_dv_infos[static_cast<int>(round(info.tf / 60.0))] = info;
						}
						else { // 已有相同时间，看谁的速度增量二范数更小
							if (V_Norm2(info.dv, 3) < V_Norm2(it->second.dv, 3)) {
								it->second = info;
							}
						}
					}
					// 整理出去重后的结果存放到一个新 vector
					dv_infos.clear();
					for (auto& kv : unique_dv_infos)
					{
						dv_infos.push_back(kv.second);
					}


					//if (flag && V_Norm2(dv, 3) < 100.0 && tf < 15000.0) {
					//	std::cout << "NR = " << NR << ", branch = " << bra << std::endl;
					//	std::cout << "tf = " << tf << std::endl;
					//	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
					//	std::cout << "imp = " << V_Norm2(dv, 3) << std::endl;
					//	std::cout << "--------------------------------------------" << std::endl;
					//}

					for(int i = 0; i < dv_infos.size() ; i++)
					{
						tf = dv_infos[i].tf;
						memcpy(dv, dv_infos[i].dv, 3 * sizeof(double));
						if (V_Norm2(dv, 3) < 1.0 && tf < 15000.0) {
							std::cout << "NR = " << NR << ", branch = " << bra << std::endl;
							std::cout << "tf = " << tf << std::endl;
							std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
							std::cout << "imp = " << V_Norm2(dv, 3) << std::endl;
							std::cout << "--------------------------------------------" << std::endl;
						}
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
	//shooting_test();
	multitree_search();

	double rv0[6];
	int flag;
	coe2rv(flag, rv0, sats_coe0[0], mu_km_s);
	double rvf_2body[6];
	double rvf_J2[6];
	
	//for (double t = 1.0 * 3600.0; t < 48.0 * 3600.0; t += 1.0* 3600.0) {
	//	
	//	propagate_j2(rv0, rvf_J2, 0.0, t, 1.0e-10, 1.0e-10);
	//	rv02rvf(flag, rvf_2body, rv0, t, mu_km_s);
	//	//propagate_j2linear(rv0, rvf_2body, 0.0, t);
	//	auto start_t = std::chrono::steady_clock::now();
	//	//propagate_linearJ2(rv0, rvf_2body, 0.0,t);

	//	auto end_t = std::chrono::steady_clock::now();
	//	double duration_t = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 1000.0;
	//	//std::cout << "Time is " << duration_t << "s ; which is " << duration_t / 60.0 << " mins;  which is " << duration_t / 3600.0 << " hours" << std::endl;
	//	//计算两者的差值
	//	double diff[3];
	//	for (int i = 0; i < 3; i++) {
	//		diff[i] = rvf_J2[i] - rvf_2body[i];
	//	}
	//	double diff_norm = V_Norm2(diff, 3);
	//	//output difference
	//	std::cout << "t" << t/3600.0 << "   Difference between J2 and 2-body propagation: " << diff_norm << std::endl;
	//}


	return 0;
}
