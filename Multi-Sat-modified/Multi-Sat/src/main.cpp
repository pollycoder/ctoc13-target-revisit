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


// 多棵树搜索（有问题，搁置）
void multitree_search() {
	MultiTree multi_tree(1000, 4, 50, 0.3);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "总耗时：" << duration_second << "秒" << std::endl;
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
	get_target_geogetic(0, 0.0, geo);
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

 
int main() {
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

	return 0;
}
