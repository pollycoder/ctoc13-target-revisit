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
#include <omp.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"
#include "J2propagation.h"
#include "single_sat.h"
#include <iomanip>
const std::string space = " ";
 
int main() {

	/***********************************************************************/
	// 搜索
	/*MultiTree multi_tree(10, 4, 50, 0.3);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "总耗时：" << duration_second << "秒" << std::endl;*/


	/***********************************************************************/
	// single_imp的失败算例
	/*double m0 = 1000.0, mf, t0 = 0.0, tf;
	double geo[2];
	get_target_geogetic(2, 0.0, geo);
	double dv[3];
	int NR = 1, branch = 0, flag;
	double coe_00[6] = { sats_coe0[1][0], sats_coe0[1][1], sats_coe0[1][2], sats_coe0[1][3], sats_coe0[1][4], sats_coe0[1][5] };
	double RV0[6];
	coe2rv(flag, RV0, coe_00, mu_km_s);
	single_imp(m0, t0, RV0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, branch);

	std::cout << "flag = " << flag << std::endl;
	std::cout << "tf = " << tf << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;*/


	/***********************************************************************/
	// 输出不同根数的观测次数
	std::ofstream fout0;
	fout0.open("../output_result/single_sat.bin", std::ios::out | std::ios::binary);

	auto beforeTime = std::chrono::steady_clock::now();

	const double h_start = 500.0; const double h_end = 1000.0;
	const double Omega_start = 0.0; const double Omega_end = D2PI;
	const double inc_start = 0.0; const double inc_end = DPI;
	const double h_mesh = (h_end - h_start) / 20;
	const double Omega_mesh = (Omega_end - Omega_start) / 5;
	const double inc_mesh = (inc_end - inc_start) / 2;

	// Start the loop: 20 x 50 x 2
	const int num_loops = 20 * 5 * 2;
	std::vector<std::string> buffer(num_loops);

	int index = 0;
	for (int i = 0; i < 20; i++) {
		double Omega = Omega_start + i * Omega_mesh;
		for (int j = 0; j < 5; j++) {
			double inc = inc_start + j * inc_mesh;
			for (int k = 0; k < 2; k++) {
				// Each loop: 0.015s
				double h = h_start + k * h_mesh;
				double coe0[6] = { Re_km + h, 0.0, inc, Omega, 0.0, 0.0 };
				int num = single_sat_score(coe0);
				buffer[index++] = std::to_string(coe0[0]) + " "
					+ std::to_string(coe0[1]) + " "
					+ std::to_string(coe0[2]) + " "
					+ std::to_string(coe0[3]) + " "
					+ std::to_string(coe0[4]) + " "
					+ std::to_string(coe0[5]) + " "
					+ std::to_string(num) + "\n";
			}
		}
	}

	// Write the entire buffer to file at once
	if (fout0.is_open()) {
		for (const auto& line : buffer) {
			fout0 << line;
		}
		fout0.close();
	}
	else {
		std::cerr << "Unable to open file";
	}

	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "总耗时：" << duration_second << "秒" << std::endl;
	// 16核，预计7.81小时完成

	return 0;
}
