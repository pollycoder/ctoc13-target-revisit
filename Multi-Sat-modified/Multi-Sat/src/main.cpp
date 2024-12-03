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

const std::string space = " ";
std::vector<std::vector<std::vector<double>>> fixed_imp = { {}, {}, {}, {}, {}, {}, {}, {} };

void output_result(std::vector<std::tuple<std::vector<double>, std::vector<std::vector<double>>>>& sat_info_list)
{
	std::ofstream fout0("../output_result/result.txt");
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
}


// 多棵树搜索（重启）
void multitree_search() {
	MultiTree multi_tree(1000, 4, 100, dv_max);

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


// 读取数据库，筛选并生成新的数据库
void read_db() {
	const std::string input_filename1 = "../output_result/single_sat_21_R.bin";
	const std::string output_filename1 = "../output_result/single_better_21_R.bin";
	int threshold = 10;

	auto beforeTime = std::chrono::steady_clock::now();
	if (filter_bin_file(input_filename1, output_filename1, threshold)) {
		auto afterTime = std::chrono::steady_clock::now();
		double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
		std::cout << "1号数据库处理完成，总耗时：" << duration_second << "秒" << std::endl;
	}
	else {
		std::cerr << "处理失败" << std::endl;
	}
	// 0.001s读取一条数据，预计4.16h跑完一个数据库

	const std::string input_filename2 = "../output_result/single_sat_21_P.bin";
	const std::string output_filename2 = "../output_result/single_better_21_P.bin";

	beforeTime = std::chrono::steady_clock::now();
	if (filter_bin_file(input_filename2, output_filename2, threshold)) {
		auto afterTime = std::chrono::steady_clock::now();
		double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
		std::cout << "2号数据库处理完成，总耗时：" << duration_second << "秒" << std::endl;
	}
	else {
		std::cerr << "处理失败" << std::endl;
	}
	// 0.001s读取一条数据，预计4.16h跑完一个数据库

	// 16核，总计8.33h完成

}


void multi_sat_opt() {
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
	DE_parallel(obj_func, nullptr, X, f, var_num, 80, 3000, 100);
	nlopt_main(obj_func, nullptr, X, f, var_num, 0, 10000);

	get_score_info(X, nullptr, f, sat_info_list, AccessTable);
	for (int i = 0; i < fixed_imp.size(); i++) {
		fixed_imp[i].insert(fixed_imp[i].end(), std::get<1>(sat_info_list[i]).begin(), std::get<1>(sat_info_list[i]).end());
	}
	for (const auto& tdv : fixed_imp) {
		for (const auto& el : tdv) {
			std::cout << el[0] << space << el[1] << space << el[2] << space << el[3] << std::endl;
		}
	}
	

	std::cout << "f = " << f << std::endl;
	std::vector<double> revisit_gap;
	max_reseetime(AccessTable, revisit_gap);
	for (const auto& t : revisit_gap) {
		std::cout << t << std::endl;
	}

	/*std::ofstream fout0("../output_result/result.txt");
	for (const auto& coe : coelist) {
		for (const auto& oe : coe) {
			fout0 << std::setprecision(14) << oe << space;
		}
		fout0 << std::endl;
	}*/
	
	
	output_result(sat_info_list);
	// 16核，0.025s运行一轮
	// 预计1.67h完成
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

	get_target_geogetic(18, 0.0, geo);

	for (int i = 0; i < 4; i++) {
		coe2rv(flag, rv0, sats_coe0[i], mu_km_s);
		std::cout << "============================================" << std::endl;
		std::cout << "Satellite " << i << std::endl;
		std::cout << "--------------------------------------------" << std::endl;
		for (int N = 1; N < 21; N++) {
			for (int bra = 0; bra < 2; bra++) {
				double m0 = 1000.0, mf, tf, dv[3];
				NR = N;
				single_imp(m0, 0.0, rv0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, bra);
				if (flag && tf < 21600.0 && V_Norm2(dv, 3) < 1.0) {
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

int main() {
	
	multi_sat_opt();
	
	return 0;
}
