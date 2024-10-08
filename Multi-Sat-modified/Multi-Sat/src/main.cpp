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
****************************************************************************/
#include <stdlib.h>
#include "multitree_beam.h"
#include "single_impluse.h"
#include "J2Lambert.h"

const std::string space = " ";
 
int main() {
	MultiTree multi_tree(100, 3, 10, 1000.0);
	multi_tree.Run();

	/*std::vector<double> X = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	double f = 1000.0;
	DE_parallel(obj_func_initial_coe, nullptr, X, f, 6, 1000, 10000, 100);

	double coe0[6];
	get_initial_coe(X, coe0, f);
	std::cout << coe0[0] << space << coe0[1] << space << coe0[2] << space << coe0[3] << space << coe0[4] << space << coe0[5] << std::endl;*/

	test_lambert();


	return 0;
}
