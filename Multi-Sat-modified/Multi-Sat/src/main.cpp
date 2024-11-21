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
****************************************************************************/
#include <stdlib.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"
#include "J2propagation.h"
#include "single_sat.h"
#include <iomanip>



const std::string space = " ";


// 多棵树搜索（重启）
void multitree_search() {
	MultiTree multi_tree(1000, 4, 50, dv_max);

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


/***********************************************************************/
// single_imp的失败算例
// 目标：12
// 预计观测时间：1620s
// 问题：tG估计错误（84000s），与圈数严重不匹配时方程无解
// 目标：13
// 预计观测时间：1800s
// 问题：升降交正确，但圈数不匹配（NR=3，实际1圈就能观测），且脉冲很大（dv = 0.11047 -3.51521 -2.95107）
// 目标：14
// 预计观测时间：8400s（过点非常精确）
// 问题：升降交正确，但圈数不匹配（NR=3，实际2圈就能观测），且脉冲不接近0（dv = 0.00431052 -0.137163 -0.11515）
// 从13000s第二次观测目标12开始，估算基本正确
/***********************************************************************/
void fail_example_singleimp() {
	double m0 = 1000.0, mf, t0 = 0.0, tf;
	double geo[2];
	get_target_geogetic(11, 0.0, geo);
	double dv[3];
	int NR = 4, branch = 0, flag;
	double coe_00[6] = { sats_coe0[0][0], sats_coe0[0][1], sats_coe0[0][2], sats_coe0[0][3], sats_coe0[0][4], sats_coe0[0][5] };
	double RV0[6];
	coe2rv(flag, RV0, coe_00, mu_km_s);
	single_imp(m0, t0, RV0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, branch);

	std::cout << "flag = " << flag << std::endl;
	std::cout << "tf = " << tf << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
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

void single_sat_opt() {
	double para[7] = { sats_coe0[1][0], sats_coe0[1][1], sats_coe0[1][2], sats_coe0[1][3], sats_coe0[1][4], sats_coe0[1][5], 86400.0 };
	std::vector<double> X = { 0.5, 0.5, 0.5, 0.5 };
	double f;

	// DE global
	auto beforeTime = std::chrono::steady_clock::now();
	DE_parallel(obj_single_sat, para, X, f, X.size(), 10, 5000, 100);
	auto afterTime = std::chrono::steady_clock::now();
	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "DE总耗时：" << duration_second << "秒" << std::endl;

	std::vector<double> max_revisit;
	double t_imp, dv[3];
	get_revisit(X, para, max_revisit, t_imp, dv, f);

	std::ofstream fout0("../output_result/result.txt");
	std::cout << "t = " << t_imp << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;

	fout0 << sats_coe0[0][0] << " " << sats_coe0[0][1] << " " << sats_coe0[0][2] << " " << sats_coe0[0][3] << " " << sats_coe0[0][4] << " " << sats_coe0[0][5] << std::endl;
	fout0 << t_imp << " " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;
	int index = 0;
	for (auto iter = max_revisit.begin(); iter != max_revisit.end(); iter++) {
		index++;
		fout0 << "Target" << index << ": " << *iter << std::endl;
		std::cout << *iter << std::endl;
	}

	// 种群：50
	// 迭代次数：10000
	// 0.025s目标函数运行一次
	// 16核，预计3.33h完成
}

void multi_sat_opt() {
	int para[impNum] = {4, 5, 6, 7};
	std::vector<double> X;
	std::vector<int> imp_idx;
	for (int i = 0; i < impNum; i++) {
		for (int j = 0; j < 4; j++) X.push_back(0.5);
		imp_idx.push_back(para[i]);
	}

	double f;

	// DE global
	/*auto beforeTime = std::chrono::steady_clock::now();
	DE_parallel(obj_multi_sat_certain, para, X, f, X.size(), 40, 10000, 50);
	nlopt_main(obj_multi_sat_certain, para, X, f, X.size(), 5000);
	auto afterTime = std::chrono::steady_clock::now();
	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "总耗时：" << duration_second << "秒" << std::endl;*/

	std::vector<double> max_revisit;
	std::vector<double> t_imp;
	std::vector<std::vector<double>> dv;
	
	//DE_parallel(obj_func_imp_coe, para, X, f, X.size(), 5 * X.size(), 10000, 100);
	get_revisit_certain(X, para, max_revisit, t_imp, dv, f);
	std::cout << "f = " << f << std::endl;

	std::ofstream fout0("../output_result/result.txt");
	int idx = 0;
	for (int i = 0; i < TreeNum; i++) {
		fout0 << sats_coe0[i][0] << " " << sats_coe0[i][1] << " " << sats_coe0[i][2] << " " << sats_coe0[i][3] << " " << sats_coe0[i][4] << " " << sats_coe0[i][5] << std::endl;
		if (std::find(imp_idx.begin(), imp_idx.end(), i) != imp_idx.end()) {
			fout0 << t_imp[idx] << " " << dv[idx][0] << " " << dv[idx][1] << " " << dv[idx][2] << std::endl;
			idx++;
		}
	}

	int index = 0;
	for (auto iter = max_revisit.begin(); iter != max_revisit.end(); iter++) {
		index++;
		fout0 << "Target" << index << ": " << *iter << std::endl;
		std::cout << "Target" << index << ": " << *iter << std::endl;
	}

	// 种群：160
	// 迭代次数：10000
	// 0.016s目标函数运行一次
	// 16核，预计8h完成
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
 
int main() {

	// 核数：16核
	// 任务：
	//	1. 树搜索，
	// 开始时间：15:00
	// 第一次验收时间：21:00
	multitree_search();
	//max_revisit_verify();
	//multi_sat_opt();
	return 0;
}
