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

 
int main() {
	/***********************************************************************/
	// 搜索
	/*MultiTree multi_tree(10000, 4, 50, 0.3);

	auto beforeTime = std::chrono::steady_clock::now();
	multi_tree.Run();
	auto afterTime = std::chrono::steady_clock::now();

	double duration_second = std::chrono::duration<double>(afterTime - beforeTime).count();
	std::cout << "总耗时：" << duration_second << "秒" << std::endl;*/
	/***********************************************************************/


	/***********************************************************************/
	// single_imp的失败算例
	/*double m0 = 1000.0, mf, t0 = 0.0, tf;
	double geo[2];
	get_target_geogetic(2, 0.0, geo);
	double dv[3];
	int NR = 2, branch = 1, flag;
	double coe_00[6] = { sats_coe0[0][0], sats_coe0[0][1], sats_coe0[0][2], sats_coe0[0][3], sats_coe0[0][4], sats_coe0[0][5] };
	double RV0[6];
	coe2rv(flag, RV0, coe_00, mu_km_s);
	single_imp(m0, t0, RV0, geo[1], geo[0], 1, flag, mf, tf, dv, NR, branch);

	std::cout << "flag = " << flag << std::endl;
	std::cout << "tf = " << tf << std::endl;
	std::cout << "dv = " << dv[0] << " " << dv[1] << " " << dv[2] << std::endl;*/


	/***********************************************************************/
	// 输出不同根数的总观测次数
	// create_db_single();
	/***********************************************************************/

	/***********************************************************************/
	// 输出不同根数的特定目标观测次数
	// 当前目标：ship
	// 顺逆：0
	int target_id = 21;
	int branch = 0;
	std::string bra_str;
	if (branch) { bra_str = "R"; }		// Retrograde，逆行
	else { bra_str = "P"; }				// Prograde，顺行
	std::string id = std::to_string(target_id);
	std::string db_filepath = "../output_result/single_sat_" + id + "_" + bra_str + ".bin";
	std::ofstream fout0(db_filepath, std::ios::out | std::ios::binary);
	create_db_target(fout0, target_id);
	// 16核，平均一条数据0.002s
	// 500 * 200 * 150数据，预计8.33h完成
	/***********************************************************************/

	return 0;
}
