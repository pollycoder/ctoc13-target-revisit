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
const std::string space = " ";
 
int main() {
	MultiTree multi_tree(10000, 4, 50, 0.1);
	multi_tree.Run();

	//test_lambert();
	//test_single_impulse();

	return 0;
}
